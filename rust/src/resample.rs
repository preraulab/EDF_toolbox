//! Resampling using rubato's FftFixedIn (FFT-based fixed-ratio resampler).
//!
//! For fixed integer ratios this is significantly faster than the
//! direct-form sinc resampler: an FFT-domain anti-alias multiply replaces
//! the per-sample convolution.
//!
//! [`resample_channels`] runs all channels through one resampler instance
//! in lockstep, which is meaningfully faster than instantiating a fresh
//! resampler per channel (each instance allocates FFT plans + buffers).
//!
//! Zero-phase alignment: `FftFixedIn` is a linear-phase FIR with a
//! constant group delay of `resampler.output_delay()` output samples.
//! We compensate by (a) discarding that many leading output samples and
//! (b) feeding zero-padded chunks past end-of-input so the trailing real
//! samples are flushed out of the filter pipeline. The final output has
//! length `round(n_in * target_rate / orig_rate)` and is time-aligned
//! with the input — matching MATLAB's `resample()` to <1% RMS.

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

    // Linear-phase group delay (in output samples) of the FFT FIR. We strip
    // this many leading samples from the assembled output to make the
    // resampler zero-phase relative to the input.
    let leading_delay = resampler.output_delay();
    // Final output length to return (matches MATLAB resample's length convention).
    let expected_out = ((n_in as f64) * (out_rate as f64) / (in_rate as f64)).round() as usize;
    // We need to collect leading_delay + expected_out output samples total.
    let target_collect = leading_delay + expected_out;

    let needed = resampler.input_frames_next();
    let mut output: Vec<Vec<f32>> = (0..nbr_channels)
        .map(|_| Vec::with_capacity(target_collect + needed))
        .collect();

    // Per-call input buffer: outer = channels, inner = chunk samples.
    let mut inbuf: Vec<Vec<f32>> = (0..nbr_channels).map(|_| vec![0.0f32; needed]).collect();

    let mut pos = 0usize;
    // Keep processing chunks until every channel has at least target_collect
    // output samples. After we run out of real input we feed zero-padded
    // chunks to flush the FIR pipeline (this is what recovers the trailing
    // real samples that would otherwise be lost to group delay).
    loop {
        let have = output[0].len();
        if have >= target_collect {
            break;
        }
        let need_in_this_chunk = resampler.input_frames_next();
        if need_in_this_chunk != needed {
            // FftFixedIn has constant input_frames_next, but be defensive.
            for ch_buf in inbuf.iter_mut() {
                ch_buf.resize(need_in_this_chunk, 0.0);
            }
        }
        let real_in = if pos < n_in {
            (n_in - pos).min(need_in_this_chunk)
        } else {
            0
        };
        for (ci, ch) in channels.iter().enumerate() {
            if real_in > 0 {
                inbuf[ci][..real_in].copy_from_slice(&ch[pos..pos + real_in]);
            }
            if real_in < need_in_this_chunk {
                inbuf[ci][real_in..].fill(0.0);
            }
        }
        let out = resampler.process(&inbuf, None).context("resampler block")?;
        for (ci, oc) in output.iter_mut().enumerate() {
            oc.extend_from_slice(&out[ci]);
        }
        pos += real_in;
        // Safety: if we are past EOF and the resampler produces nothing,
        // bail rather than spin forever.
        if real_in == 0 && out[0].is_empty() {
            break;
        }
    }

    // Trim leading group delay and truncate to expected length.
    for oc in output.iter_mut() {
        if oc.len() > leading_delay {
            oc.drain(..leading_delay);
        } else {
            oc.clear();
        }
        if oc.len() > expected_out {
            oc.truncate(expected_out);
        }
    }

    Ok(output)
}
