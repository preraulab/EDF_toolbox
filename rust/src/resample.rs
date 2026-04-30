//! Resampling using rubato (pure Rust, polyphase FIR with sinc interpolation).
//!
//! This is the convert pipeline:
//!   int16 samples -> f32 (apply scale+offset to physical units)
//!   -> rubato resampler (sinc, ratio P/Q)
//!   -> f32 output samples
//!   -> caller's responsibility to quantize back to int16

use anyhow::{Context, Result};
use rubato::{
    Resampler, SincFixedIn, SincInterpolationParameters, SincInterpolationType, WindowFunction,
};

/// Resample a single channel from `orig_rate` to `target_rate`.
/// Input is in physical units (f32). Output is also f32.
pub fn resample_channel(x: &[f32], orig_rate: f64, target_rate: f64) -> Result<Vec<f32>> {
    if (orig_rate - target_rate).abs() < 1e-9 {
        return Ok(x.to_vec());
    }

    let ratio = target_rate / orig_rate;

    // Kaiser-windowed sinc params, tuned to roughly match MATLAB resample's
    // default (10 zero-crossings per side, beta=5).
    let params = SincInterpolationParameters {
        sinc_len: 256,
        f_cutoff: 0.95,
        interpolation: SincInterpolationType::Linear,
        oversampling_factor: 256,
        window: WindowFunction::BlackmanHarris2,
    };

    let chunk_size = 1 << 14; // 16384 samples per chunk
    let mut resampler = SincFixedIn::<f32>::new(ratio, 2.0, params, chunk_size, 1)
        .context("creating SincFixedIn resampler")?;

    let mut output = Vec::<f32>::with_capacity((x.len() as f64 * ratio) as usize + 1024);
    let mut pos = 0usize;
    let n = x.len();

    // Feed full chunks
    while pos + chunk_size <= n {
        let inbuf = vec![x[pos..pos + chunk_size].to_vec()];
        let out = resampler.process(&inbuf, None).context("resampler chunk")?;
        output.extend_from_slice(&out[0]);
        pos += chunk_size;
    }
    // Pad the trailing block to chunk_size (rubato's SincFixedIn requires
    // exact-sized blocks). Trim the corresponding output to discard the
    // padding's filtered contribution.
    if pos < n {
        let remain = n - pos;
        let mut padded = Vec::with_capacity(chunk_size);
        padded.extend_from_slice(&x[pos..]);
        padded.resize(chunk_size, 0.0);
        let inbuf = vec![padded];
        let out = resampler.process(&inbuf, None).context("resampler tail")?;
        let keep = (remain as f64 * ratio).round() as usize;
        let take = keep.min(out[0].len());
        output.extend_from_slice(&out[0][..take]);
    }

    Ok(output)
}
