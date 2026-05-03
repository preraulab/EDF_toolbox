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
use clap::{ArgGroup, Parser, ValueEnum};
use rayon::prelude::*;
use std::io::{BufRead, IsTerminal, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;
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

/// Policy for what to do when an intended output path already exists.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum ClobberPolicy {
    /// Default: skip the file and warn. Safe for shared NFS data dirs.
    Warn,
    /// `-n`: skip silently (cp -n idiom).
    NoClobber,
    /// `-f`: overwrite without asking.
    Force,
    /// `-i`: prompt with [y]es / [n]o / [A]ll / [N]one. Falls back to Warn
    /// if stdin is not a tty.
    Prompt,
}

/// Sticky decision from an interactive prompt. Once the user picks `A`
/// (yes-to-all) or `N` (no-to-all), we apply that to every later collision
/// in the same run without prompting again.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
enum StickyPrompt {
    Undecided,
    YesAll,
    NoAll,
}

/// Policy for deleting source EDFs after a successful convert (or, in
/// `--cleanup-only` mode, after confirming an output already exists).
#[derive(Copy, Clone, Debug, PartialEq, Eq, ValueEnum)]
enum DeletePolicy {
    /// Default: keep every source. No deletion happens.
    Never,
    /// Verify the output decompresses cleanly, then prompt
    /// [y]es / [n]o / [A]ll / [N]one before deleting each source.
    Ask,
    /// Verify the output decompresses cleanly, then delete silently.
    Silent,
    /// Skip verification AND prompt; delete after a successful rename.
    /// Fastest, useful only for trusted pipelines.
    Force,
}

/// `clap`'s `--version` string. Combines the Cargo package version with
/// the git SHA and commit date captured by `build.rs`. After this change
/// `convert_edf --version` prints e.g.
///     convert_edf 0.1.0 (d869c4d44b 2026-05-03)
/// so users can compare the running binary against `git ls-remote` without
/// having a source tree.
const VERSION: &str = concat!(
    env!("CARGO_PKG_VERSION"),
    " (",
    env!("GIT_SHA"),
    " ",
    env!("GIT_BUILD_DATE"),
    ")"
);

#[derive(Parser, Debug)]
#[command(version = VERSION, about = "Resample and rewrite EDF files")]
#[command(group = ArgGroup::new("clobber").multiple(false).args(["interactive", "force", "no_clobber"]))]
struct Args {
    /// Input files or directories (mixable). Directories are scanned for
    /// .edf / .edf.gz / .edf.zst; pass -R to recurse.
    inputs: Vec<PathBuf>,

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

    /// Number of parallel workers in batch mode. 0 (default) = auto:
    /// `min(num_cpus, 8)`. The cap exists because each worker holds a full
    /// decoded PSG in RAM (often 1-3 GB as f32), so on a 64-core box with
    /// no cap you'd OOM long before saturating CPU. Bump explicitly only
    /// if you know the files are small or the machine has lots of RAM.
    /// Spare cores are absorbed into per-file codec threads automatically.
    #[arg(long, default_value_t = 0)]
    jobs: usize,

    /// Prompt before each existing-file overwrite ([y]es / [n]o / [A]ll /
    /// [N]one). Falls back to skip-and-warn if stdin isn't a terminal.
    /// Mutually exclusive with -f / -n.
    #[arg(short = 'i', long)]
    interactive: bool,

    /// Force overwrite of existing output files without prompting.
    /// Mutually exclusive with -i / -n.
    #[arg(short = 'f', long)]
    force: bool,

    /// Skip silently if the output already exists (cp -n idiom). The default
    /// (no flag) is the same skip behavior but with a one-line warning per
    /// collision; -n suppresses the warning. Mutually exclusive with -i / -f.
    #[arg(short = 'n', long)]
    no_clobber: bool,

    /// Delete each source EDF after the convert successfully writes its
    /// compressed sibling. Modes:
    ///   never   — keep every source (default)
    ///   ask     — verify the output decompresses, then prompt y/n/A/N
    ///   silent  — verify the output decompresses, then delete silently
    ///   force   — skip verification AND prompt; delete after rename
    /// `ask` falls back to skip-and-warn if stdin is not a terminal.
    #[arg(long, value_enum, default_value_t = DeletePolicy::Never)]
    delete_source: DeletePolicy,

    /// Skip the conversion step entirely. For each input, if its expected
    /// output (`<stem>_<rate>Hz.edf.{gz,zst}`) already exists, apply the
    /// `--delete-source` policy to the source. Useful for finishing
    /// disk-cleanup after a previous run completed without
    /// `--delete-source`. Requires `--delete-source` not be `never`.
    #[arg(long)]
    cleanup_only: bool,

    /// Threads to use within each per-file codec (gzip via gzp, zstd via
    /// libzstd). 0 (default) = auto-divide available cores across rayon
    /// workers so total active threads stay near `nproc` instead of
    /// `nproc²`. Set explicitly if you have a specific reason; otherwise
    /// leave at 0.
    #[arg(long, default_value_t = 0)]
    codec_threads: usize,

    /// Verbose output.
    #[arg(short = 'v', long)]
    verbose: bool,
}

impl Args {
    fn clobber_policy(&self) -> ClobberPolicy {
        if self.force {
            ClobberPolicy::Force
        } else if self.no_clobber {
            ClobberPolicy::NoClobber
        } else if self.interactive {
            ClobberPolicy::Prompt
        } else {
            ClobberPolicy::Warn
        }
    }
}

fn main() -> Result<()> {
    let mut args = Args::parse();

    let roots: Vec<PathBuf> = args.inputs.clone();
    if roots.is_empty() {
        bail!("no inputs given (pass one or more FILE / DIR arguments)");
    }
    if !args.read_bench && args.target_rate <= 0.0 {
        bail!("--target-rate / -F is required (or pass --read-bench to skip conversion)");
    }

    // Resolve `--jobs 0` to a memory-aware default. See Args::jobs docs.
    const DEFAULT_JOBS_CAP: usize = 8;
    if args.jobs == 0 {
        args.jobs = num_cpus().min(DEFAULT_JOBS_CAP).max(1);
    }
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.jobs)
        .build_global()
        .ok();

    let clobber = args.clobber_policy();
    // Sticky decision shared across rayon workers for `-i` interactive mode.
    let prompt_state = Mutex::new(StickyPrompt::Undecided);

    let delete_policy = args.delete_source;
    let delete_prompt_state = Mutex::new(StickyPrompt::Undecided);
    if args.cleanup_only && delete_policy == DeletePolicy::Never {
        bail!(
            "--cleanup-only requires --delete-source <ask|silent|force>; without a delete \
             policy, cleanup-only mode would walk the tree and do nothing"
        );
    }

    // Collect (top_root, file_path) so we can mirror the source tree.
    // top_root is the user-supplied path under which file_path was found;
    // for a single file argument, top_root is the file's parent directory.
    let pairs = collect_inputs(&roots, args.recursive, args.target_rate)?;

    if args.read_bench {
        return run_read_bench_pairs(&pairs);
    }

    if pairs.is_empty() {
        bail!("no .edf / .edf.gz / .edf.zst files found");
    }

    // Auto-divide cores between rayon workers (one per file in flight) and
    // codec threads (gzp / zstd-mt inside each worker), so total active
    // threads stay near `nproc` rather than `nproc²`. Without this, on a
    // 64-core box every rayon worker spawns 64 codec threads => ~4096
    // threads contending => severe thrashing and OOM kills on large PSGs.
    if args.codec_threads == 0 {
        let cpus = num_cpus();
        let effective_workers = args.jobs.min(pairs.len().max(1)).max(1);
        args.codec_threads = (cpus / effective_workers).max(1);
        if args.verbose {
            println!(
                "thread budget: {} cores / {} rayon workers => {} codec threads each",
                cpus, effective_workers, args.codec_threads
            );
        }
    }

    // Single-file shortcut: exactly one input file, --out points at a file
    // (not a directory). Preserve the historical "give me one file out"
    // ergonomic.
    if pairs.len() == 1 {
        if let Some(out) = &args.out {
            let is_dir = out.is_dir() || roots.iter().any(|r| r.is_dir());
            if !is_dir {
                let (_root, input) = &pairs[0];
                if args.cleanup_only {
                    if out.exists() {
                        check_delete(input, out, delete_policy, &delete_prompt_state, args.verbose)?;
                    } else if args.verbose {
                        eprintln!("[cleanup-skip] {}: no output at {}", input.display(), out.display());
                    }
                    return Ok(());
                }
                if !check_clobber(out, clobber, &prompt_state, args.verbose) {
                    return Ok(());
                }
                let t0 = Instant::now();
                convert_one(input, out, &args)?;
                if args.verbose {
                    println!("{} -> {} in {:.2}s", input.display(), out.display(),
                             t0.elapsed().as_secs_f64());
                }
                check_delete(input, out, delete_policy, &delete_prompt_state, args.verbose)?;
                return Ok(());
            }
        } else {
            // No --out and one file: write next to the input.
            let (_root, input) = &pairs[0];
            let out = default_output_path(input, args.target_rate, args.compress);
            if args.cleanup_only {
                if out.exists() {
                    check_delete(input, &out, delete_policy, &delete_prompt_state, args.verbose)?;
                } else if args.verbose {
                    eprintln!("[cleanup-skip] {}: no output at {}", input.display(), out.display());
                }
                return Ok(());
            }
            if !check_clobber(&out, clobber, &prompt_state, args.verbose) {
                return Ok(());
            }
            let t0 = Instant::now();
            convert_one(input, &out, &args)?;
            if args.verbose {
                println!("{} -> {} in {:.2}s", input.display(), out.display(),
                         t0.elapsed().as_secs_f64());
            }
            check_delete(input, &out, delete_policy, &delete_prompt_state, args.verbose)?;
            return Ok(());
        }
    }

    // Multi-file / directory mode. --out is optional:
    //   - given:   write outputs under that directory (mirrored or flat)
    //   - omitted: write each output next to its input (in-situ)
    let out_dir = args.out.clone();
    if let Some(d) = &out_dir {
        if d.is_file() {
            bail!(
                "--out points at a regular file ({}) but {} inputs were given; \
                 pass a directory in multi-file mode (or omit --out to write in-situ)",
                d.display(),
                pairs.len()
            );
        }
        std::fs::create_dir_all(d).with_context(|| format!("creating {}", d.display()))?;
    }

    let mode_msg = match &out_dir {
        Some(d) if args.flatten => format!("--out {} (flat)", d.display()),
        Some(d)                 => format!("--out {} (mirrored)", d.display()),
        None                    => "in-situ (writing next to each input)".to_string(),
    };
    if args.cleanup_only {
        println!(
            "convert_edf: {} inputs, target {} Hz, {} -- cleanup-only ({})",
            pairs.len(), args.target_rate, mode_msg,
            match delete_policy {
                DeletePolicy::Ask => "verify+prompt",
                DeletePolicy::Silent => "verify+silent",
                DeletePolicy::Force => "no-verify+silent",
                DeletePolicy::Never => unreachable!(),  // bailed earlier
            }
        );
    } else {
        println!("convert_edf: {} files, target {} Hz, {}",
                 pairs.len(), args.target_rate, mode_msg);
    }

    let n = pairs.len();
    let t0 = Instant::now();
    // Live progress counter -- atomically incremented as each file
    // finishes. With many rayon workers, file completions arrive out of
    // order; this counter is the user-visible "N of total done" tally.
    let done = std::sync::atomic::AtomicUsize::new(0);
    // Sources actually deleted across the run (Inline + cleanup-only).
    let deleted = std::sync::atomic::AtomicUsize::new(0);
    // Approximate ETA from the wall time so far. Cheap and good enough.
    // Outcome value: Ok(Some(secs)) = converted (or, in cleanup-only, the
    // source was deleted), Ok(None) = skipped, Err = failed.
    let outcomes: Vec<(PathBuf, Result<Option<f64>>)> = pairs
        .par_iter()
        .map(|(root, input)| {
            let out = match &out_dir {
                Some(d) => mirror_out_path(root, input, d, args.target_rate, args.compress, args.flatten),
                None    => default_output_path(input, args.target_rate, args.compress),
            };

            // Cleanup-only branch: don't convert, just check whether the
            // expected output is already there and apply delete policy.
            if args.cleanup_only {
                let i = done.fetch_add(1, std::sync::atomic::Ordering::Relaxed) + 1;
                if !out.exists() {
                    if args.verbose {
                        println!("[{}/{}] cleanup-skip {}: no output at {}",
                                 i, n, input.display(), out.display());
                        std::io::stdout().flush().ok();
                    }
                    return (input.clone(), Ok(None));
                }
                match check_delete(input, &out, delete_policy, &delete_prompt_state, false) {
                    Ok(true) => {
                        deleted.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                        if args.verbose {
                            println!("[{}/{}] cleanup-deleted {}", i, n, input.display());
                            std::io::stdout().flush().ok();
                        }
                        (input.clone(), Ok(Some(0.0)))
                    }
                    Ok(false) => {
                        if args.verbose {
                            println!("[{}/{}] cleanup-kept {}", i, n, input.display());
                            std::io::stdout().flush().ok();
                        }
                        (input.clone(), Ok(None))
                    }
                    Err(e) => {
                        if args.verbose {
                            println!("[{}/{}] cleanup-fail {}: {:#}", i, n, input.display(), e);
                            std::io::stdout().flush().ok();
                        }
                        (input.clone(), Err(e))
                    }
                }
            } else {
                // Normal convert path.
                // Ensure parent dir of `out` exists (mirrored subdirs).
                if let Some(parent) = out.parent() {
                    if let Err(e) = std::fs::create_dir_all(parent) {
                        let i = done.fetch_add(1, std::sync::atomic::Ordering::Relaxed) + 1;
                        if args.verbose {
                            println!("[{}/{}] fail {}: creating {}: {}",
                                     i, n, input.display(), parent.display(), e);
                        }
                        return (input.clone(),
                            Err(anyhow!("creating {}: {}", parent.display(), e)));
                    }
                }
                // Honor the configured clobber policy. For `-i` (Prompt) the
                // mutex inside check_clobber serializes prompts across workers.
                if !check_clobber(&out, clobber, &prompt_state, args.verbose) {
                    let i = done.fetch_add(1, std::sync::atomic::Ordering::Relaxed) + 1;
                    if args.verbose {
                        println!("[{}/{}] skip {} (output exists)", i, n, out.display());
                        std::io::stdout().flush().ok();
                    }
                    return (input.clone(), Ok(None));
                }
                let t1 = Instant::now();
                let result = match convert_one_quiet(input, &out, &args) {
                    Ok(()) => Ok(Some(t1.elapsed().as_secs_f64())),
                    Err(e) => Err(e),
                };
                // On a successful convert, optionally delete the source.
                if let Ok(Some(_)) = &result {
                    if let Ok(true) = check_delete(input, &out, delete_policy, &delete_prompt_state, false) {
                        deleted.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    }
                }
                let i = done.fetch_add(1, std::sync::atomic::Ordering::Relaxed) + 1;
                if args.verbose {
                    let elapsed = t0.elapsed().as_secs_f64();
                    // ETA: extrapolate from completed/elapsed. Avoid div-by-zero.
                    let eta = if i > 0 && i < n {
                        let per = elapsed / i as f64;
                        per * (n - i) as f64
                    } else { 0.0 };
                    match &result {
                        Ok(Some(secs)) => println!(
                            "[{}/{}] ok   {} ({:.2}s) | elapsed {:.0}s, ETA {:.0}s",
                            i, n, input.display(), secs, elapsed, eta),
                        Ok(None) => {} // skip already announced above
                        Err(e) => println!(
                            "[{}/{}] fail {}: {:#}", i, n, input.display(), e),
                    }
                    std::io::stdout().flush().ok();
                }
                (input.clone(), result)
            }
        })
        .collect();
    let total = t0.elapsed().as_secs_f64();
    let total_deleted = deleted.load(std::sync::atomic::Ordering::Relaxed);

    let mut ok = 0usize;
    let mut skipped = 0usize;
    for (path, res) in &outcomes {
        match res {
            Ok(Some(_)) => ok += 1,
            Ok(None) => skipped += 1,
            // Always surface failures, even without -v (so silent runs
            // still tell you what didn't make it).
            Err(e) if !args.verbose => println!("[fail] {}: {:#}", path.display(), e),
            Err(_) => {}  // already printed live under -v
        }
    }
    let failed = n - ok - skipped;
    let mut parts: Vec<String> = vec![format!("{}/{} ok", ok, n)];
    if skipped > 0 {
        parts.push(format!("{} skipped", skipped));
    }
    if failed > 0 {
        parts.push(format!("{} failed", failed));
    }
    if total_deleted > 0 || delete_policy != DeletePolicy::Never {
        parts.push(format!("{} sources deleted", total_deleted));
    }
    println!(
        "done: {} in {:.1}s wall ({:.2} files/s)",
        parts.join(", "), total, n as f64 / total
    );
    if failed > 0 {
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

/// True if `name` looks like one of *our* output files for the given target
/// rate -- e.g. "foo_100Hz.edf.gz" when rate=100. Used to make in-situ
/// recursive runs idempotent: a second pass over a tree that already
/// contains "_100Hz.edf*" outputs from a previous run skips them instead of
/// producing "_100Hz_100Hz.edf*" garbage. Skip is silent.
///
/// Only filters during directory walks. Explicit positional file arguments
/// are always honored even if they match this pattern -- if the user names
/// the file by hand, they meant it.
fn looks_like_our_output_at_rate(name: &str, rate: f64) -> bool {
    let bare = name
        .strip_suffix(".zst")
        .or_else(|| name.strip_suffix(".gz"))
        .unwrap_or(name);
    if let Some(stem) = bare.strip_suffix(".edf") {
        let needle = format!("_{}Hz", rate as i64);
        stem.ends_with(&needle)
    } else {
        false
    }
}

/// Resolve user-supplied root paths into a list of (root, file) pairs.
/// `root` is the original positional argument under which `file` was
/// discovered; for a file argument the root is the file's parent directory
/// (so the relative-path-from-root is just the basename).
///
/// `target_rate` is used to filter our own past outputs out of directory
/// walks (idempotent re-runs). Pass 0.0 to disable the filter (e.g. for
/// --read-bench, which has no rate to compare against).
fn collect_inputs(
    roots: &[PathBuf],
    recursive: bool,
    target_rate: f64,
) -> Result<Vec<(PathBuf, PathBuf)>> {
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
            walk_dir(root, root, recursive, target_rate, &mut out)?;
        } else {
            bail!("input is not a regular file or directory: {}", root.display());
        }
    }
    out.sort();
    Ok(out)
}

/// Recursive (or single-level) directory walker. Skips hidden entries
/// (names starting with '.') so we don't pick up macOS .DS_Store, .git,
/// etc. by accident, and skips files matching our own output pattern at
/// `target_rate` to keep in-situ re-runs idempotent.
fn walk_dir(
    root: &Path,
    dir: &Path,
    recursive: bool,
    target_rate: f64,
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
            if is_edf_filename(&name_str)
                && !(target_rate > 0.0 && looks_like_our_output_at_rate(&name_str, target_rate))
            {
                out.push((root.to_path_buf(), p));
            }
        } else if ft.is_dir() && recursive {
            walk_dir(root, &p, recursive, target_rate, out)?;
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

/// Decide whether to (over)write `out` given the configured clobber policy.
/// Returns true if the caller should proceed with conversion + write.
///
/// `prompt_state` is shared across rayon workers so that the sticky
/// All / None decisions made in one prompt stick globally for the rest of
/// the run.
fn check_clobber(
    out: &Path,
    policy: ClobberPolicy,
    prompt_state: &Mutex<StickyPrompt>,
    verbose: bool,
) -> bool {
    if !out.exists() {
        return true;
    }
    match policy {
        ClobberPolicy::Force => true,
        ClobberPolicy::NoClobber => false,
        ClobberPolicy::Warn => {
            if verbose {
                eprintln!(
                    "[skip] {} already exists (use -f to overwrite, -i to prompt)",
                    out.display()
                );
            } else {
                eprintln!("[skip] {} already exists", out.display());
            }
            false
        }
        ClobberPolicy::Prompt => {
            // Hold the lock for the entire prompt so concurrent workers
            // don't interleave reads from stdin or writes of the prompt.
            let mut sticky = prompt_state.lock().unwrap();
            match *sticky {
                StickyPrompt::YesAll => return true,
                StickyPrompt::NoAll => return false,
                StickyPrompt::Undecided => {}
            }
            if !std::io::stdin().is_terminal() {
                eprintln!(
                    "[skip] {} already exists (-i requires a tty; non-interactive run)",
                    out.display()
                );
                return false;
            }
            print!(
                "overwrite '{}'? [y]es / [n]o / [A]ll / [N]one: ",
                out.display()
            );
            std::io::stdout().flush().ok();
            let mut buf = String::new();
            if std::io::stdin().lock().read_line(&mut buf).is_err() {
                return false;
            }
            match buf.trim() {
                "y" | "Y" | "yes" => true,
                "A" | "all" => {
                    *sticky = StickyPrompt::YesAll;
                    true
                }
                "N" | "none" => {
                    *sticky = StickyPrompt::NoAll;
                    false
                }
                _ => false,
            }
        }
    }
}

/// Decompress `path` end-to-end into `io::sink()` to confirm the codec
/// trailer (CRC + length for gzip; frame end + checksum for zstd) is
/// well-formed. A truncated / corrupt file -- e.g. the half-baked output
/// of an earlier kill-9 -- raises an error here, which then *blocks* the
/// matching source from being deleted.
///
/// Cheap because the read is decompress-only; nothing is written to disk.
/// Plain `.edf` (no compression) needs no integrity check beyond
/// existence -- nothing to verify.
fn verify_output(path: &Path) -> Result<()> {
    use std::io;
    let codec = compress_from_path(path).unwrap_or(Compress::None);
    let f = std::fs::File::open(path)
        .with_context(|| format!("opening {} for verify", path.display()))?;
    match codec {
        Compress::None => Ok(()),
        Compress::Gzip => {
            let mut dec = flate2::read::MultiGzDecoder::new(std::io::BufReader::new(f));
            io::copy(&mut dec, &mut io::sink())
                .with_context(|| format!("verify gunzip {}", path.display()))?;
            Ok(())
        }
        Compress::Zstd => {
            let mut dec = zstd::stream::read::Decoder::new(f)
                .with_context(|| format!("creating zstd decoder for {}", path.display()))?;
            io::copy(&mut dec, &mut io::sink())
                .with_context(|| format!("verify zstd {}", path.display()))?;
            Ok(())
        }
    }
}

/// Decide whether to delete `input` given the configured `DeletePolicy`
/// and the just-written or pre-existing `output`. Returns Ok(true) if the
/// source was deleted, Ok(false) if it was kept (either by policy, by
/// failed verify, or by user prompt).
///
/// Mirrors `check_clobber`: the `prompt_state` mutex serializes prompts
/// across rayon workers and tracks sticky All / None decisions.
///
/// Safety belt: refuses to delete if `input` and `output` canonicalize
/// to the same path (would happen only if a future caller passes
/// `-o <input>` -- but the check is one syscall so it's free insurance).
fn check_delete(
    input: &Path,
    output: &Path,
    policy: DeletePolicy,
    prompt_state: &Mutex<StickyPrompt>,
    verbose: bool,
) -> Result<bool> {
    if policy == DeletePolicy::Never {
        return Ok(false);
    }
    // Same-path safety belt.
    if let (Ok(a), Ok(b)) = (
        std::fs::canonicalize(input),
        std::fs::canonicalize(output),
    ) {
        if a == b {
            eprintln!(
                "[skip-delete] refusing to delete {}: input and output canonicalize to the same path",
                input.display()
            );
            return Ok(false);
        }
    }
    // Verification step (skipped only for Force).
    if policy != DeletePolicy::Force {
        if let Err(e) = verify_output(output) {
            eprintln!(
                "[skip-delete] {}: output {} failed integrity check: {:#}",
                input.display(),
                output.display(),
                e
            );
            return Ok(false);
        }
    }
    // Prompt branch (Ask only).
    if policy == DeletePolicy::Ask {
        let mut sticky = prompt_state.lock().unwrap();
        let go = match *sticky {
            StickyPrompt::YesAll => true,
            StickyPrompt::NoAll => false,
            StickyPrompt::Undecided => {
                if !std::io::stdin().is_terminal() {
                    eprintln!(
                        "[skip-delete] {} (--delete-source ask requires a tty; non-interactive run)",
                        input.display()
                    );
                    return Ok(false);
                }
                print!(
                    "delete source '{}'? [y]es / [n]o / [A]ll / [N]one: ",
                    input.display()
                );
                std::io::stdout().flush().ok();
                let mut buf = String::new();
                if std::io::stdin().lock().read_line(&mut buf).is_err() {
                    return Ok(false);
                }
                match buf.trim() {
                    "y" | "Y" | "yes" => true,
                    "A" | "all" => {
                        *sticky = StickyPrompt::YesAll;
                        true
                    }
                    "N" | "none" => {
                        *sticky = StickyPrompt::NoAll;
                        false
                    }
                    _ => false,
                }
            }
        };
        if !go {
            return Ok(false);
        }
    }
    // Do the delete.
    std::fs::remove_file(input)
        .with_context(|| format!("deleting source {}", input.display()))?;
    if verbose {
        eprintln!("[deleted] {}", input.display());
    }
    Ok(true)
}

fn convert_one(input: &Path, output: &Path, args: &Args) -> Result<()> {
    convert_one_inner(input, output, args, args.verbose)
}

/// Internal entry point that lets the batch mode suppress the per-file
/// header / channel dump (which is useful when debugging one file but
/// pure noise across hundreds). Live `[N/total]` progress lines come
/// from main()'s loop, not from here.
fn convert_one_quiet(input: &Path, output: &Path, args: &Args) -> Result<()> {
    convert_one_inner(input, output, args, false)
}

fn convert_one_inner(input: &Path, output: &Path, args: &Args, verbose: bool) -> Result<()> {
    if verbose {
        eprintln!("opening {}", input.display());
    }
    let (mut hdr, data_i16) = edf::read_edf_path(input)?;
    if verbose {
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

    // Atomic write: stream into <name>.partial, then rename to the final
    // path on success. A Ctrl-C or crash mid-write leaves only the .partial
    // sibling -- the final filename either does not exist or still holds
    // the previous good copy. Rename is atomic on a single POSIX filesystem.
    let tmp_name = format!(
        "{}.partial",
        path.file_name().and_then(|s| s.to_str()).unwrap_or(".out")
    );
    let tmp = path.with_file_name(tmp_name);

    let write_to_tmp = || -> Result<()> {
        let f = std::fs::File::create(&tmp)
            .with_context(|| format!("creating {}", tmp.display()))?;
        match compress {
            Compress::None => {
                let mut w = std::io::BufWriter::new(f);
                edf::write_header(&mut w, hdr)?;
                edf::write_data(&mut w, hdr, data)?;
                w.flush()?;
            }
            Compress::Gzip => {
                // Parallel gzip via gzp -- output is standard gzip format,
                // readable by gunzip / zlib / flate2 / MATLAB read_EDF.
                use gzp::deflate::Gzip;
                use gzp::par::compress::ParCompressBuilder;
                use gzp::Compression as GzpCompression;
                use gzp::ZWriter;
                let lvl = GzpCompression::new(args.gzip_level);
                let mut w = ParCompressBuilder::<Gzip>::new()
                    .num_threads(args.codec_threads.max(1))
                    .map_err(|e| anyhow!("gzp num_threads: {e:?}"))?
                    .compression_level(lvl)
                    .from_writer(f);
                edf::write_header(&mut w, hdr)?;
                edf::write_data(&mut w, hdr, data)?;
                w.finish().map_err(|e| anyhow!("gzp finish: {e:?}"))?;
            }
            Compress::Zstd => {
                let mut enc = zstd::stream::write::Encoder::new(f, args.zstd_level)?;
                enc.multithread(args.codec_threads.max(1) as u32).ok();
                let mut w = enc.auto_finish();
                edf::write_header(&mut w, hdr)?;
                edf::write_data(&mut w, hdr, data)?;
                w.flush()?;
            }
        }
        Ok(())
    };

    match write_to_tmp() {
        Ok(()) => {
            std::fs::rename(&tmp, path)
                .with_context(|| format!("renaming {} -> {}", tmp.display(), path.display()))?;
            Ok(())
        }
        Err(e) => {
            // Best-effort cleanup of the half-written .partial.
            let _ = std::fs::remove_file(&tmp);
            Err(e)
        }
    }
}

fn num_cpus() -> usize {
    std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1)
}
