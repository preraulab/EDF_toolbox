//! EDF / EDF+ format reader and writer.
//!
//! Spec summary:
//!   - 256-byte main header
//!   - n_signals * 256 bytes of signal headers (label, dims, scale, etc.)
//!   - n_records data blocks, each block = sum(samples_per_record[s]) * i16 LE
//!
//! Annotation signals are labeled "EDF Annotations" and carry TAL (time-stamped
//! annotation list) bytes packed two-per-i16. We pass them through verbatim
//! when writing.

use anyhow::{anyhow, bail, Context, Result};
use std::io::{Read, Write};

const MAIN_HEADER_LEN: usize = 256;
const SIGNAL_HEADER_LEN: usize = 256;
const ANNOTATIONS_LABEL: &str = "EDF Annotations";

#[derive(Debug, Clone)]
pub struct Header {
    pub version: String,
    pub patient_id: String,
    pub recording_id: String,
    pub start_date: String, // dd.mm.yy
    pub start_time: String, // hh.mm.ss
    pub bytes_in_header: u32,
    pub reserved: String, // EDF+C / EDF+D / blank
    pub num_data_records: i64,
    pub data_record_duration: f64, // seconds
    pub num_signals: u32,
    pub signals: Vec<SignalHeader>,
}

#[derive(Debug, Clone)]
pub struct SignalHeader {
    pub label: String,
    pub transducer: String,
    pub physical_dim: String,
    pub physical_min: f64,
    pub physical_max: f64,
    pub digital_min: i32,
    pub digital_max: i32,
    pub prefiltering: String,
    pub samples_per_record: u32,
    pub reserved: String,
}

impl SignalHeader {
    pub fn is_annotations(&self) -> bool {
        self.label.trim() == ANNOTATIONS_LABEL
    }
    pub fn sampling_frequency(&self, record_duration: f64) -> f64 {
        self.samples_per_record as f64 / record_duration
    }
    /// Map int16 sample -> physical units (mV / uV / etc).
    pub fn scale(&self) -> f64 {
        let dr = (self.digital_max - self.digital_min) as f64;
        if dr == 0.0 {
            1.0
        } else {
            (self.physical_max - self.physical_min) / dr
        }
    }
    pub fn offset(&self) -> f64 {
        self.physical_min - self.scale() * self.digital_min as f64
    }
}

fn read_str<R: Read>(r: &mut R, n: usize) -> Result<String> {
    let mut buf = vec![0u8; n];
    r.read_exact(&mut buf)?;
    let s = std::str::from_utf8(&buf)
        .map_err(|_| anyhow!("non-utf8 in EDF header"))?
        .to_string();
    Ok(s)
}

fn parse_field<T: std::str::FromStr>(s: &str, field: &str) -> Result<T>
where
    T::Err: std::fmt::Display,
{
    s.trim()
        .parse::<T>()
        .map_err(|e| anyhow!("parsing {field} from '{s}': {e}"))
}

pub fn read_header<R: Read>(r: &mut R) -> Result<Header> {
    let version = read_str(r, 8)?;
    let patient_id = read_str(r, 80)?;
    let recording_id = read_str(r, 80)?;
    let start_date = read_str(r, 8)?;
    let start_time = read_str(r, 8)?;
    let bytes_in_header: u32 = parse_field(&read_str(r, 8)?, "bytes_in_header")?;
    let reserved = read_str(r, 44)?;
    let num_data_records: i64 = parse_field(&read_str(r, 8)?, "num_data_records")?;
    let data_record_duration: f64 = parse_field(&read_str(r, 8)?, "data_record_duration")?;
    let num_signals: u32 = parse_field(&read_str(r, 4)?, "num_signals")?;
    if std::env::var("EDF_DEBUG").is_ok() {
        eprintln!(
            "DEBUG main: bih={} reserved='{}' nrec={} dur={} ns={}",
            bytes_in_header, reserved.trim(), num_data_records, data_record_duration, num_signals
        );
    }

    let ns = num_signals as usize;
    let mut labels = Vec::with_capacity(ns);
    for _ in 0..ns {
        labels.push(read_str(r, 16)?);
    }
    let mut transducers = Vec::with_capacity(ns);
    for _ in 0..ns {
        transducers.push(read_str(r, 80)?);
    }
    let mut phys_dims = Vec::with_capacity(ns);
    for _ in 0..ns {
        phys_dims.push(read_str(r, 8)?);
    }
    let mut phys_mins = Vec::with_capacity(ns);
    for _ in 0..ns {
        phys_mins.push(parse_field::<f64>(&read_str(r, 8)?, "physical_min")?);
    }
    let mut phys_maxs = Vec::with_capacity(ns);
    for _ in 0..ns {
        phys_maxs.push(parse_field::<f64>(&read_str(r, 8)?, "physical_max")?);
    }
    // Some EDFs encode digital_min/max as floats ("-300.000"); accept either.
    let mut dig_mins = Vec::with_capacity(ns);
    for _ in 0..ns {
        let s = read_str(r, 8)?;
        let v: f64 = parse_field(&s, "digital_min")?;
        dig_mins.push(v.round() as i32);
    }
    let mut dig_maxs = Vec::with_capacity(ns);
    for _ in 0..ns {
        let s = read_str(r, 8)?;
        let v: f64 = parse_field(&s, "digital_max")?;
        dig_maxs.push(v.round() as i32);
    }
    let mut prefilters = Vec::with_capacity(ns);
    for _ in 0..ns {
        prefilters.push(read_str(r, 80)?);
    }
    let mut spr = Vec::with_capacity(ns);
    for _ in 0..ns {
        spr.push(parse_field::<u32>(&read_str(r, 8)?, "samples_per_record")?);
    }
    if std::env::var("EDF_DEBUG").is_ok() {
        eprintln!("DEBUG: spr per signal = {:?}", spr);
        eprintln!("DEBUG: labels = {:?}", labels.iter().map(|s| s.trim().to_string()).collect::<Vec<_>>());
    }
    let mut reserveds = Vec::with_capacity(ns);
    for _ in 0..ns {
        reserveds.push(read_str(r, 32)?);
    }

    let signals = (0..ns)
        .map(|i| SignalHeader {
            label: labels[i].clone(),
            transducer: transducers[i].clone(),
            physical_dim: phys_dims[i].clone(),
            physical_min: phys_mins[i],
            physical_max: phys_maxs[i],
            digital_min: dig_mins[i],
            digital_max: dig_maxs[i],
            prefiltering: prefilters[i].clone(),
            samples_per_record: spr[i],
            reserved: reserveds[i].clone(),
        })
        .collect();

    Ok(Header {
        version,
        patient_id,
        recording_id,
        start_date,
        start_time,
        bytes_in_header,
        reserved,
        num_data_records,
        data_record_duration,
        num_signals,
        signals,
    })
}

/// Read all signals into Vec<Vec<i16>> -- one inner vec per signal, length =
/// num_data_records * samples_per_record[s].
///
/// If `hdr.num_data_records < 0` (EDF+ continuous "unknown record count"
/// sentinel), this reads until EOF and updates `hdr.num_data_records` to the
/// actual record count.
pub fn read_data<R: Read>(r: &mut R, hdr: &mut Header) -> Result<Vec<Vec<i16>>> {
    let ns = hdr.num_signals as usize;
    let record_size: usize = hdr
        .signals
        .iter()
        .map(|s| s.samples_per_record as usize * 2)
        .sum();
    if record_size == 0 {
        bail!("zero-length data record (no signals?)");
    }
    let mut rec_buf = vec![0u8; record_size];

    let known_nrec = hdr.num_data_records >= 0;
    let cap_per_signal: usize = if known_nrec {
        hdr.num_data_records as usize
    } else {
        // unknown: just hint a starting capacity, will grow as we go
        1024
    };
    let mut out: Vec<Vec<i16>> = (0..ns)
        .map(|s| Vec::with_capacity(cap_per_signal * hdr.signals[s].samples_per_record as usize))
        .collect();

    let mut nrec_read: i64 = 0;
    loop {
        // Try to read exactly one record; on partial/EOF, stop.
        let mut got = 0usize;
        while got < record_size {
            match r.read(&mut rec_buf[got..]) {
                Ok(0) => break,
                Ok(n) => got += n,
                Err(e) => {
                    if matches!(e.kind(), std::io::ErrorKind::Interrupted) {
                        continue;
                    }
                    return Err(e).context("reading EDF data record");
                }
            }
        }
        if got == 0 {
            break;
        }
        if got < record_size {
            bail!(
                "truncated final record: read {} of {} bytes (record {})",
                got,
                record_size,
                nrec_read
            );
        }
        let mut off = 0usize;
        for s in 0..ns {
            let n = hdr.signals[s].samples_per_record as usize;
            let bytes = &rec_buf[off..off + n * 2];
            out[s].extend(bytes.chunks_exact(2).map(|c| i16::from_le_bytes([c[0], c[1]])));
            off += n * 2;
        }
        nrec_read += 1;
        if known_nrec && nrec_read >= hdr.num_data_records {
            break;
        }
    }
    if !known_nrec {
        hdr.num_data_records = nrec_read;
    }
    Ok(out)
}

/// Convenience: read an entire EDF file path, returning (header, data).
pub fn read_edf_path(path: &std::path::Path) -> Result<(Header, Vec<Vec<i16>>)> {
    let f = std::fs::File::open(path).with_context(|| format!("opening {}", path.display()))?;
    // Choose decoder by extension
    let mut reader: Box<dyn Read> = match path.extension().and_then(|e| e.to_str()) {
        Some("gz") => Box::new(flate2::read::MultiGzDecoder::new(f)),
        Some("zst") => Box::new(zstd::stream::read::Decoder::new(f)?),
        _ => Box::new(f),
    };
    let mut hdr = read_header(&mut reader)?;
    let data = read_data(&mut reader, &mut hdr)?;
    Ok((hdr, data))
}

// ---------- WRITER ----------

fn fmt_field(s: &str, width: usize) -> Vec<u8> {
    let mut v = s.as_bytes().to_vec();
    if v.len() > width {
        v.truncate(width);
    } else {
        v.resize(width, b' ');
    }
    v
}

fn fmt_int(v: i64, width: usize) -> Vec<u8> {
    fmt_field(&format!("{}", v), width)
}

fn fmt_dbl(v: f64, width: usize) -> Vec<u8> {
    // EDF wants the most precision that fits in the field.
    let mut s = format!("{:.6}", v);
    if s.len() > width {
        // try fewer decimals
        for d in (0..6).rev() {
            s = format!("{:.*}", d, v);
            if s.len() <= width {
                break;
            }
        }
    }
    if s.len() > width {
        s.truncate(width);
    }
    fmt_field(&s, width)
}

pub fn write_header<W: Write>(w: &mut W, hdr: &Header) -> Result<()> {
    let ns = hdr.num_signals as usize;
    let bytes_in_header = MAIN_HEADER_LEN + ns * SIGNAL_HEADER_LEN;
    w.write_all(&fmt_field(&hdr.version, 8))?;
    w.write_all(&fmt_field(&hdr.patient_id, 80))?;
    w.write_all(&fmt_field(&hdr.recording_id, 80))?;
    w.write_all(&fmt_field(&hdr.start_date, 8))?;
    w.write_all(&fmt_field(&hdr.start_time, 8))?;
    w.write_all(&fmt_int(bytes_in_header as i64, 8))?;
    w.write_all(&fmt_field(&hdr.reserved, 44))?;
    w.write_all(&fmt_int(hdr.num_data_records, 8))?;
    w.write_all(&fmt_dbl(hdr.data_record_duration, 8))?;
    w.write_all(&fmt_int(ns as i64, 4))?;

    for s in &hdr.signals {
        w.write_all(&fmt_field(&s.label, 16))?;
    }
    for s in &hdr.signals {
        w.write_all(&fmt_field(&s.transducer, 80))?;
    }
    for s in &hdr.signals {
        w.write_all(&fmt_field(&s.physical_dim, 8))?;
    }
    for s in &hdr.signals {
        w.write_all(&fmt_dbl(s.physical_min, 8))?;
    }
    for s in &hdr.signals {
        w.write_all(&fmt_dbl(s.physical_max, 8))?;
    }
    for s in &hdr.signals {
        w.write_all(&fmt_int(s.digital_min as i64, 8))?;
    }
    for s in &hdr.signals {
        w.write_all(&fmt_int(s.digital_max as i64, 8))?;
    }
    for s in &hdr.signals {
        w.write_all(&fmt_field(&s.prefiltering, 80))?;
    }
    for s in &hdr.signals {
        w.write_all(&fmt_int(s.samples_per_record as i64, 8))?;
    }
    for s in &hdr.signals {
        w.write_all(&fmt_field(&s.reserved, 32))?;
    }
    Ok(())
}

pub fn write_data<W: Write>(w: &mut W, hdr: &Header, data: &[Vec<i16>]) -> Result<()> {
    if data.len() != hdr.num_signals as usize {
        bail!(
            "write_data: header has {} signals but got {} data vectors",
            hdr.num_signals,
            data.len()
        );
    }
    for r in 0..hdr.num_data_records as usize {
        for (s, sig) in hdr.signals.iter().enumerate() {
            let n = sig.samples_per_record as usize;
            let off = r * n;
            let end = off + n;
            if end > data[s].len() {
                bail!(
                    "write_data: signal {} length {} short of record {} (need {} samples)",
                    s,
                    data[s].len(),
                    r,
                    end
                );
            }
            for &v in &data[s][off..end] {
                w.write_all(&v.to_le_bytes())?;
            }
        }
    }
    Ok(())
}

