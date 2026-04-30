//! Resampling using rubato's FftFixedIn (FFT-based fixed-ratio resampler).
//!
//! For fixed integer ratios this is significantly faster than the
//! direct-form sinc resampler: an FFT-domain anti-alias multiply replaces
//! the per-sample convolution.
//!
//! [`resample_channels`] runs all channels through one resampler instance
//! in lockstep, which is meaningfully faster than instantiating a fresh
//! resampler per channel (each instance allocates FFT plans + buffers).

use anyhow::{Context, Result};
use rubato::{FftFixedIn, Resampler};

/// Resample several channels (each is a `Vec<f32>`, all the same length)
/// from `orig_rate` to `target_rate`. Returns one resampled `Vec<f32>` per
/// input channel.
pub fn resample_channels(
    channels: &[Vec<f32>],
    orig_rate: f64,
    target_rate: f64,
) -> Result<Vec<Vec<f32>>> {
    if channels.is_empty() {
        return Ok(Vec::new());
    }
    if (orig_rate - target_rate).abs() < 1e-9 {
        return Ok(channels.to_vec());
    }
    let in_rate = orig_rate.round() as usize;
    let out_rate = target_rate.round() as usize;
    if in_rate == 0 || out_rate == 0 {
        anyhow::bail!("non-integer-Hz resample not supported");
    }
    let n_in = channels[0].len();
    for c in channels {
        if c.len() != n_in {
            anyhow::bail!("resample_channels: all channels must have equal length");
        }
    }
    let nbr_channels = channels.len();

    let chunk_size_in = 1024;
    let sub_chunks = 2;
    let mut resampler = FftFixedIn::<f32>::new(
        in_rate,
        out_rate,
        chunk_size_in,
        sub_chunks,
        nbr_channels,
    )
    .context("creating FftFixedIn resampler")?;

    let needed = resampler.input_frames_next();
    let est_out = ((n_in as f64) * target_rate / orig_rate) as usize + needed;
    let mut output: Vec<Vec<f32>> =
        (0..nbr_channels).map(|_| Vec::with_capacity(est_out)).collect();

    // Per-call input buffer: outer = channels, inner = chunk samples
    let mut inbuf: Vec<Vec<f32>> = (0..nbr_channels).map(|_| vec![0.0f32; needed]).collect();

    let mut pos = 0usize;
    while pos < n_in {
        let real_in = (n_in - pos).min(needed);
        for (ci, ch) in channels.iter().enumerate() {
            inbuf[ci][..real_in].copy_from_slice(&ch[pos..pos + real_in]);
            if real_in < needed {
                inbuf[ci][real_in..].fill(0.0);
            }
        }
        let out = resampler.process(&inbuf, None).context("resampler block")?;
        if real_in == needed {
            for (ci, oc) in output.iter_mut().enumerate() {
                oc.extend_from_slice(&out[ci]);
            }
        } else {
            // Tail: output corresponds to real_in real samples plus padding.
            let real_out = (real_in as f64 * target_rate / orig_rate).round() as usize;
            for (ci, oc) in output.iter_mut().enumerate() {
                let take = real_out.min(out[ci].len());
                oc.extend_from_slice(&out[ci][..take]);
            }
            break;
        }
        pos += needed;
    }

    Ok(output)
}
