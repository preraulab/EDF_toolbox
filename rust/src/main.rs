//! convert_edf -- standalone EDF resample / rewrite tool.
//!
//! Usage:
//!   convert_edf input.edf -F 100 -o output.edf.gz
//!   convert_edf -F 100 -o /path/to/out /path/to/edfs
//!   convert_edf -F 100 -o /path/to/out -R /path/to/study_root
//!   convert_edf -F 100 -o /path/to/out a.edf b.edf c.edf
//!
//! Positional arguments may be files or directories, mixed. Pass `-R` to
//! recurse into subdirectories. Output structure mirrors the source tree
//! by default; pass `--flatten` to drop all outputs at the top of `--out`.
//!
//! Pipeline per file:
//!   read EDF (auto-detect .edf / .edf.gz / .edf.zst)
//!   for each non-annotation signal:
//!     int16 -> f32 (physical units)
//!     resample at target_rate
//!     f32 -> int16 (recompute physical_min/max from data range, by default)
//!   write EDF (auto-pick .edf / .edf.gz / .edf.zst by output extension)

mod edf;
mod resample;

use anyhow::{anyhow, bail, Context, Result};
use clap::{Parser, ValueEnum};
use rayon::prelude::*;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::Instant;

#[derive(Copy, Clone, Debug, ValueEnum)]
enum Compress {
    None,
    Gzip,
    Zstd,
}

#[derive(Copy, Clone, Debug, ValueEnum)]
enum AutoScale {
    /// Keep source physical_min/max. Filter ringing may clip to int16 range.
    Preserve,
    /// Recompute physical_min/max from resampled-channel data range. (default)
    Recompute,
}

#[derive(Parser, Debug)]
#[command(version, about = "Resample and rewrite EDF files")]
struct Args {
    /// Input files or directories (mixable). Directories are scanned for
    /// .edf / .edf.gz / .edf.zst; pass -R to recurse.
    inputs: Vec<PathBuf>,

    /// Backward-compat alias for passing a single directory positionally.
    /// Prefer the positional form: `convert_edf -F 100 --out OUT /path/to/dir`.
    #[arg(long, hide = true)]
    batch: Option<PathBuf>,

    /// Recurse into subdirectories of any directory inputs.
    #[arg(short = 'R', long)]
    recursive: bool,

    /// In recursive mode, write all outputs flat in --out instead of
    /// mirroring the source directory tree. Ignored if not recursive.
    #[arg(long)]
    flatten: bool,

    /// Target sampling rate (Hz). Required unless --read-bench is set.
    /// Short form is -F (sampling frequency, Fs). -r is intentionally not
    /// used to avoid visual confusion with -R (recursive). --rate is
    /// accepted as a long alias.
    #[arg(short = 'F', long = "target-rate", visible_alias = "rate", default_value_t = 0.0)]
    target_rate: f64,

    /// Read-only benchmark: read every input EDF (decompressing as needed)
    /// and print the elapsed wall time. No output is written.
    #[arg(long, default_value_t = false)]
    read_bench: bool,

    /// Output file (single-file input only) or output directory (multi-file
    /// or directory inputs). When mirroring directory structure, this is the
    /// root under which the source tree is reproduced.
    #[arg(short = 'o', long)]
    out: Option<PathBuf>,

    /// Compression for output. Default: gzip. With the gzp parallel-gzip
    /// encoder, gzip-9 is pareto-optimal in this pipeline -- nearly the same
    /// output size as zstd-9 (within ~1%) at roughly half the wall time.
    /// Pick zstd if you specifically want faster decode (gzip adds ~80 ms/file
    /// vs zstd's ~20 ms/file -- small in absolute terms either way).
    #[arg(long, value_enum, default_value_t = Compress::Gzip)]
    compress: Compress,

    /// Zstd compression level (1-22). zstd-3 is fast and a reasonable
    /// alternative to gzip-9 for decode-speed-sensitive use; zstd-19+ is for
    /// archival (smaller, much slower to compress).
    #[arg(long, default_value_t = 3)]
    zstd_level: i32,

    /// Gzip compression level (1-9). Default 9: smallest output of the fast
    /// codec tier in the parallel-gzip path.
    #[arg(long, default_value_t = 9)]
    gzip_level: u32,

    /// Auto-rescale physical_min/max to fit resampled data range.
    #[arg(long, value_enum, default_value_t = AutoScale::Recompute)]
    auto_scale: AutoScale,

    /// Number of parallel workers in batch mode (0 = use rayon default).
    #[arg(long, default_value_t = 0)]
    jobs: usize,

    /// Verbose output.
    #[arg(short = 'v', long)]
    verbose: bool,
}

fn main() -> Result<()> {
    let args = Args::parse();

    // Merge positional inputs and the legacy --batch alias.
    let mut roots: Vec<PathBuf> = args.inputs.clone();
    if let Some(b) = &args.batch {
        roots.push(b.clone());
    }
    if roots.is_empty() {
        bail!("no inputs given (pass FILE..., DIR..., or --batch DIR)");
    }
    if !args.read_bench && args.target_rate <= 0.0 {
        bail!("--target-rate / -F is required (or pass --read-bench to skip conversion)");
    }

    if args.jobs > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.jobs)
            .build_global()
            .ok();
    }

    // Collect (top_root, file_path) so we can mirror the source tree.
    // top_root is the user-supplied path under which file_path was found;
    // for a single file argument, top_root is the file's parent directory.
    let pairs = collect_inputs(&roots, args.recursive)?;

    if args.read_bench {
        return run_read_bench_pairs(&pairs);
    }

    if pairs.is_empty() {
        bail!("no .edf / .edf.gz / .edf.zst files found");
    }

    // Single-file shortcut: exactly one input file, --out points at a file
    // (not a directory). Preserve the historical "give me one file out"
    // ergonomic.
    if pairs.len() == 1 {
        if let Some(out) = &args.out {
            let is_dir = out.is_dir() || roots.iter().any(|r| r.is_dir());
            if !is_dir {
                let (_root, input) = &pairs[0];
                let t0 = Instant::now();
                convert_one(input, out, &args)?;
                if args.verbose {
                    println!("{} -> {} in {:.2}s", input.display(), out.display(),
                             t0.elapsed().as_secs_f64());
                }
                return Ok(());
            }
        } else {
            // No --out and one file: write next to the input.
            let (_root, input) = &pairs[0];
            let out = default_output_path(input, args.target_rate, args.compress);
            let t0 = Instant::now();
            convert_one(input, &out, &args)?;
            if args.verbose {
                println!("{} -> {} in {:.2}s", input.display(), out.display(),
                         t0.elapsed().as_secs_f64());
            }
            return Ok(());
        }
    }

    // Multi-file / directory mode. --out must be a directory.
    let out_dir = args
        .out
        .clone()
        .ok_or_else(|| anyhow!("--out DIR is required when converting multiple files or a directory"))?;
    if out_dir.is_file() {
        bail!(
            "--out points at a regular file ({}) but {} inputs were given; \
             pass a directory in multi-file mode",
            out_dir.display(),
            pairs.len()
        );
    }
    std::fs::create_dir_all(&out_dir).with_context(|| format!("creating {}", out_dir.display()))?;

    println!("convert_edf: {} files, target {} Hz", pairs.len(), args.target_rate);

    let t0 = Instant::now();
    let outcomes: Vec<(PathBuf, Result<f64>)> = pairs
        .par_iter()
        .map(|(root, input)| {
            let out = mirror_out_path(root, input, &out_dir, args.target_rate, args.compress, args.flatten);
            // Ensure parent dir of `out` exists (mirrored subdirs).
            if let Some(parent) = out.parent() {
                if let Err(e) = std::fs::create_dir_all(parent) {
                    return (input.clone(),
                        Err(anyhow!("creating {}: {}", parent.display(), e)));
                }
            }
            let t1 = Instant::now();
            match convert_one(input, &out, &args) {
                Ok(()) => (input.clone(), Ok(t1.elapsed().as_secs_f64())),
                Err(e) => (input.clone(), Err(e)),
            }
        })
        .collect();
    let total = t0.elapsed().as_secs_f64();

    let mut ok = 0usize;
    for (path, res) in &outcomes {
        match res {
            Ok(secs) => {
                ok += 1;
                if args.verbose {
                    println!("[ok]   {} ({:.2}s)", path.display(), secs);
                }
            }
            Err(e) => println!("[fail] {}: {:#}", path.display(), e),
        }
    }
    println!(
        "done: {}/{} ok in {:.1}s wall ({:.2} files/s)",
        ok,
        outcomes.len(),
        total,
        outcomes.len() as f64 / total
    );
    if ok < outcomes.len() {
        std::process::exit(1);
    }
    Ok(())
}

/// Read every input EDF in parallel (rayon), decompressing as needed,
/// loading all int16 samples into memory and dropping them. Prints total
/// wall time and aggregate throughput. Used to measure how compression
/// affects read speed end-to-end.
fn run_read_bench_pairs(pairs: &[(PathBuf, PathBuf)]) -> Result<()> {
    let inputs: Vec<PathBuf> = pairs.iter().map(|(_, f)| f.clone()).collect();
    if inputs.is_empty() {
        bail!("no inputs to read-bench");
    }
    let total_input_bytes: u64 = inputs
        .iter()
        .map(|p| std::fs::metadata(p).map(|m| m.len()).unwrap_or(0))
        .sum();

    let t0 = Instant::now();
    let outcomes: Vec<Result<usize>> = inputs
        .par_iter()
        .map(|p| {
            let (hdr, data) = edf::read_edf_path(p)?;
            // Touch every sample so the optimizer can't elide the read.
            let mut total = 0usize;
            for ch in &data {
                total += ch.len();
            }
            // Use hdr to make sure parsing happened
            let _ = hdr.num_signals;
            Ok(total)
        })
        .collect();
    let wall = t0.elapsed().as_secs_f64();

    let mut ok = 0usize;
    for r in &outcomes {
        if r.is_ok() {
            ok += 1;
        }
    }
    let mb_per_s = (total_input_bytes as f64) / (1024.0 * 1024.0) / wall;
    println!(
        "read-bench: {}/{} files in {:.2}s ({} bytes, {:.1} MB/s on-disk input)",
        ok,
        inputs.len(),
        wall,
        total_input_bytes,
        mb_per_s
    );
    Ok(())
}

/// True if the filename has an EDF extension we know about.
fn is_edf_filename(name: &str) -> bool {
    name.ends_with(".edf") || name.ends_with(".edf.gz") || name.ends_with(".edf.zst")
}

/// Resolve user-supplied root paths into a list of (root, file) pairs.
/// `root` is the original positional argument under which `file` was
/// discovered; for a file argument the root is the file's parent directory
/// (so the relative-path-from-root is just the basename).
fn collect_inputs(roots: &[PathBuf], recursive: bool) -> Result<Vec<(PathBuf, PathBuf)>> {
    let mut out = Vec::new();
    for root in roots {
        if !root.exists() {
            bail!("input not found: {}", root.display());
        }
        if root.is_file() {
            let name = root.file_name().and_then(|n| n.to_str()).unwrap_or("");
            if !is_edf_filename(name) {
                bail!("not an EDF file: {}", root.display());
            }
            let parent = root.parent().unwrap_or_else(|| Path::new("."));
            out.push((parent.to_path_buf(), root.clone()));
        } else if root.is_dir() {
            walk_dir(root, root, recursive, &mut out)?;
        } else {
            bail!("input is not a regular file or directory: {}", root.display());
        }
    }
    out.sort();
    Ok(out)
}

/// Recursive (or single-level) directory walker. Skips hidden entries
/// (names starting with '.') so we don't pick up macOS .DS_Store, .git,
/// etc. by accident.
fn walk_dir(
    root: &Path,
    dir: &Path,
    recursive: bool,
    out: &mut Vec<(PathBuf, PathBuf)>,
) -> Result<()> {
    for entry in std::fs::read_dir(dir).with_context(|| format!("read_dir {}", dir.display()))? {
        let entry = entry?;
        let p = entry.path();
        let name = entry.file_name();
        let name_str = name.to_string_lossy();
        if name_str.starts_with('.') {
            continue;
        }
        let ft = entry.file_type()?;
        if ft.is_file() {
            if is_edf_filename(&name_str) {
                out.push((root.to_path_buf(), p));
            }
        } else if ft.is_dir() && recursive {
            walk_dir(root, &p, recursive, out)?;
        }
    }
    Ok(())
}

/// Compute the output path for a discovered file.
/// - flatten: drop the source-tree path; place the renamed file directly
///   under `out_dir` (named like default_output_path's basename).
/// - mirror (default): replicate the path of `input` relative to `root`
///   under `out_dir`, then rename the leaf file.
fn mirror_out_path(
    root: &Path,
    input: &Path,
    out_dir: &Path,
    rate: f64,
    c: Compress,
    flatten: bool,
) -> PathBuf {
    let renamed = default_output_path(input, rate, c)
        .file_name()
        .map(|n| n.to_owned())
        .unwrap_or_default();
    if flatten {
        return out_dir.join(renamed);
    }
    let rel = input.strip_prefix(root).unwrap_or(input);
    let parent_rel = rel.parent().unwrap_or_else(|| Path::new(""));
    out_dir.join(parent_rel).join(renamed)
}

fn default_output_path(input: &Path, rate: f64, c: Compress) -> PathBuf {
    let stem = input
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("out")
        .trim_end_matches(".zst")
        .trim_end_matches(".gz")
        .trim_end_matches(".edf");
    let ext = match c {
        Compress::None => "edf",
        Compress::Gzip => "edf.gz",
        Compress::Zstd => "edf.zst",
    };
    let parent = input.parent().unwrap_or_else(|| Path::new("."));
    parent.join(format!("{}_{}Hz.{}", stem, rate as i64, ext))
}

fn compress_from_path(p: &Path) -> Option<Compress> {
    let name = p.file_name().and_then(|n| n.to_str()).unwrap_or("");
    if name.ends_with(".edf.gz") || name.ends_with(".gz") {
        Some(Compress::Gzip)
    } else if name.ends_with(".edf.zst") || name.ends_with(".zst") {
        Some(Compress::Zstd)
    } else if name.ends_with(".edf") {
        Some(Compress::None)
    } else {
        None
    }
}

fn convert_one(input: &Path, output: &Path, args: &Args) -> Result<()> {
    if args.verbose {
        eprintln!("opening {}", input.display());
    }
    let (mut hdr, data_i16) = edf::read_edf_path(input)?;
    if args.verbose {
        eprintln!(
            "  header: {} signals, {} records x {}s",
            hdr.num_signals, hdr.num_data_records, hdr.data_record_duration
        );
        for (i, s) in hdr.signals.iter().enumerate() {
            eprintln!(
                "  ch[{}] '{}' spr={} digital=[{},{}] phys=[{},{}]",
                i,
                s.label.trim(),
                s.samples_per_record,
                s.digital_min,
                s.digital_max,
                s.physical_min,
                s.physical_max
            );
        }
    }
    if hdr.num_data_records < 0 {
        bail!("could not infer num_data_records from file (still -1 after read)");
    }

    let target_rate = args.target_rate;
    let record_duration = hdr.data_record_duration;

    let new_spr_f = target_rate * record_duration;
    if (new_spr_f - new_spr_f.round()).abs() > 1e-9 {
        bail!(
            "target_rate * data_record_duration = {} not integer",
            new_spr_f
        );
    }
    let new_spr = new_spr_f.round() as u32;

    // --- Group signals by (orig_rate, input length) ---
    // Same-group channels are resampled together with one FftFixedIn instance,
    // amortizing FFT plan setup + buffer alloc across all channels of that rate.
    use std::collections::BTreeMap;
    let mut groups: BTreeMap<(u64, usize), Vec<usize>> = BTreeMap::new();
    for (si, sig) in hdr.signals.iter().enumerate() {
        if sig.is_annotations() {
            continue;
        }
        let orig_rate = sig.sampling_frequency(record_duration);
        let key = (orig_rate.to_bits(), data_i16[si].len());
        groups.entry(key).or_default().push(si);
    }

    // f32 resampled outputs, indexed by signal index. Annotations stay None.
    let mut yphys_per_signal: Vec<Option<Vec<f32>>> = (0..hdr.signals.len()).map(|_| None).collect();

    for ((orig_rate_bits, _len), idxs) in &groups {
        let orig_rate = f64::from_bits(*orig_rate_bits);
        // Decode i16 -> f32 for each channel of the group
        let xphys: Vec<Vec<f32>> = idxs
            .iter()
            .map(|&si| {
                let s = &hdr.signals[si];
                let scale = s.scale() as f32;
                let offset = s.offset() as f32;
                data_i16[si]
                    .iter()
                    .map(|&v| (v as f32) * scale + offset)
                    .collect()
            })
            .collect();

        let yphys: Vec<Vec<f32>> = if (orig_rate - target_rate).abs() < 1e-9 {
            xphys
        } else {
            resample::resample_channels(&xphys, orig_rate, target_rate)?
        };

        for (k, &si) in idxs.iter().enumerate() {
            yphys_per_signal[si] = Some(yphys[k].clone());
        }
    }

    // --- Quantize back to int16 + update header ---
    let mut new_signals = Vec::with_capacity(hdr.signals.len());
    let mut new_data: Vec<Vec<i16>> = Vec::with_capacity(hdr.signals.len());
    let mut min_records = i64::MAX;
    for (si, sig) in hdr.signals.iter().enumerate() {
        if sig.is_annotations() {
            new_signals.push(sig.clone());
            new_data.push(data_i16[si].clone());
            continue;
        }
        let yphys = yphys_per_signal[si]
            .as_ref()
            .expect("non-annotation signal missing resampled data");

        let nrec = (yphys.len() / new_spr as usize) as i64;
        if nrec < 1 {
            bail!("signal '{}' has no full records after resample", sig.label.trim());
        }
        if nrec < min_records {
            min_records = nrec;
        }
        let len = nrec as usize * new_spr as usize;
        let yphys = &yphys[..len];

        let (pmin, pmax) = match args.auto_scale {
            AutoScale::Preserve => (sig.physical_min, sig.physical_max),
            AutoScale::Recompute => {
                let (mut lo, mut hi) = (f32::INFINITY, f32::NEG_INFINITY);
                for &v in yphys {
                    if v < lo { lo = v; }
                    if v > hi { hi = v; }
                }
                if !lo.is_finite() || !hi.is_finite() || (hi - lo) < 1e-30 {
                    (sig.physical_min, sig.physical_max)
                } else {
                    (lo as f64, hi as f64)
                }
            }
        };

        let dmin: i32 = -32768;
        let dmax: i32 = 32767;
        let new_scale = (pmax - pmin) / (dmax - dmin) as f64;
        let new_scale = if new_scale == 0.0 { 1.0 } else { new_scale };
        let new_offset = pmin - new_scale * dmin as f64;
        let new_data_sig: Vec<i16> = yphys
            .iter()
            .map(|&v| {
                let q = ((v as f64 - new_offset) / new_scale).round();
                q.clamp(dmin as f64, dmax as f64) as i16
            })
            .collect();

        let mut new_sig = sig.clone();
        new_sig.physical_min = pmin;
        new_sig.physical_max = pmax;
        new_sig.digital_min = dmin;
        new_sig.digital_max = dmax;
        new_sig.samples_per_record = new_spr;
        new_signals.push(new_sig);
        new_data.push(new_data_sig);
    }

    if min_records == i64::MAX {
        min_records = hdr.num_data_records;
    }
    hdr.signals = new_signals;
    hdr.num_data_records = min_records;

    write_output(output, &hdr, &new_data, args)
}

fn write_output(path: &Path, hdr: &edf::Header, data: &[Vec<i16>], args: &Args) -> Result<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent)?;
        }
    }
    // If the output path's extension implies a different compression than
    // --compress, the path wins (so `-o foo.edf.gz` does what you'd expect
    // even with the zstd default).
    let compress = compress_from_path(path).unwrap_or(args.compress);

    let f = std::fs::File::create(path).with_context(|| format!("creating {}", path.display()))?;
    match compress {
        Compress::None => {
            let mut w = std::io::BufWriter::new(f);
            edf::write_header(&mut w, hdr)?;
            edf::write_data(&mut w, hdr, data)?;
            w.flush()?;
        }
        Compress::Gzip => {
            // Parallel gzip via gzp -- output is standard gzip format, readable
            // by gunzip / zlib / flate2 / MATLAB read_EDF.
            use gzp::deflate::Gzip;
            use gzp::par::compress::ParCompressBuilder;
            use gzp::Compression as GzpCompression;
            use gzp::ZWriter;
            let lvl = GzpCompression::new(args.gzip_level);
            let mut w = ParCompressBuilder::<Gzip>::new()
                .num_threads(num_cpus().max(1))
                .map_err(|e| anyhow!("gzp num_threads: {e:?}"))?
                .compression_level(lvl)
                .from_writer(f);
            edf::write_header(&mut w, hdr)?;
            edf::write_data(&mut w, hdr, data)?;
            w.finish().map_err(|e| anyhow!("gzp finish: {e:?}"))?;
        }
        Compress::Zstd => {
            let mut enc = zstd::stream::write::Encoder::new(f, args.zstd_level)?;
            enc.multithread(num_cpus().max(1) as u32).ok();
            let mut w = enc.auto_finish();
            edf::write_header(&mut w, hdr)?;
            edf::write_data(&mut w, hdr, data)?;
            w.flush()?;
        }
    }
    Ok(())
}

fn num_cpus() -> usize {
    std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1)
}
