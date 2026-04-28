# EDF_toolbox

High-performance EDF / EDF+ reader **and writer** for MATLAB, with compiled MEX backends and pure-MATLAB fallbacks. Reads and writes plain `.edf`, gzip-compressed `.edf.gz`, and zstd-compressed `.edf.zst` (on the fly, no temp file). Bundled with `convert_EDF` (read → resample → write) and `batch_convert_EDF` (parallel multi-file pipeline), plus a `bin/convert_edf` shell CLI.

Handles the things that break off-the-shelf EDF readers in practice: malformed record counts in clinical files, EDF+ TAL annotations, A-B rereferencing at load time, de-identification, compressed archives, anti-aliased resampling, and large polysomnography recordings where a naive loop is slow.

Part of the Prerau Lab [`preraulab_utilities`](https://github.com/preraulab/preraulab_utilities) meta-repository. Can also be used standalone.

## Quick start

```matlab
% Read (.edf / .edf.gz / .edf.zst all work)
[header, signal_header, signal_cell, annotations] = read_EDF('sleep.edf');

% Write — output extension picks the format
write_EDF('sleep_copy.edf.gz',  header, signal_header, signal_cell, annotations);
write_EDF('sleep_copy.edf.zst', header, signal_header, signal_cell, annotations);

% Read → resample to 128 Hz → write (zstd-compressed, anti-aliased)
convert_EDF('sleep.edf', 128);
% writes sleep_128Hz.edf.zst next to the input

% Batch a directory of EDFs, parallel
result = batch_convert_EDF('/path/to/edfs', 128, 'OutputDir', '/path/out');
```

From the shell:

```sh
bin/convert_edf -r 128 -o /path/out /path/to/edfs       # zstd by default
bin/convert_edf -r 128 --compress-mode gzip --gzip-level 1 /path/to/edfs
```

- `header` — struct with main-header metadata (patient ID, record duration, start date/time, etc.)
- `signal_header` — struct array, one per channel (label, transducer, physical range, sample rate)
- `signal_cell` — cell array of channel vectors in physical units (scaled from digital)
- `annotations` — struct array of EDF+ events (onset, text)

## Features

| Feature | What it does |
|---|---|
| **Full header parsing** | 256-byte main header + 16-field-per-signal header per EDF spec |
| **MEX acceleration** | Compiled C reader; pure MATLAB fallback if MEX isn't available or is disabled |
| **Compressed input/output** | Reads and writes `.edf.gz` (zlib) and `.edf.zst` (libzstd) directly — no temp file, streaming on the fly. zstd is ~5× faster to compress than gzip-6 with smaller output and is the recommended format. |
| **Channel subsetting** | Load only specific channels by name |
| **A-B rereferencing** | `'EEG C3-A2'` — compute on-the-fly during load |
| **Epoch subsetting** | Load only a `[start_epoch end_epoch]` slice |
| **EDF+ annotations** | Parse TAL (Time-Annotation List) format |
| **Header repair** | Correct invalid `num_data_records` and save a `_fixed` copy (plain `.edf` only) |
| **De-identification** | Strip PHI fields and save a `_deidentified` copy (plain `.edf` only) |
| **Digital → physical scaling** | `phys_val = phys_min + (dig_val - dig_min) × (phys_max - phys_min) / (dig_max - dig_min)` — done correctly per-channel |
| **Streaming I/O** | Reads and writes record-by-record; peak memory is ~one record + the output arrays, even for multi-GB files |
| **Write side** | `write_EDF` mirrors `read_EDF`; same struct shapes round-trip through `read → write → read` to within ~1 digital LSB |
| **Resample pipeline** | `convert_EDF` reads, resamples every channel to a target rate (anti-aliased Kaiser FIR via SPT `resample`), and writes — defaults to `.edf.zst` output |
| **Parallel batch** | `batch_convert_EDF` runs `convert_EDF` over a list/dir/cell of files via `parfor`; per-file errors don't abort the batch |
| **Shell CLI** | `bin/convert_edf` — POSIX shell wrapper, single `matlab -batch` invocation per call |

## Usage

### Basic

```matlab
% Load all channels
[header, sig_hdr, signals, annot] = read_EDF('psg.edf');

% Load specific channels (names are matched case-insensitively)
[~, ~, signals] = read_EDF('psg.edf', 'Channels', {'EEG C3-A2', 'EEG O2-A1', 'EOG Left-Right'});
```

### A-B rereferencing

```matlab
% Rereference on load — EDF file contains both EEG C3 and EEG A2 separately
[~, ~, s] = read_EDF('psg.edf', 'Channels', {'EEG C3-A2'});
```

The reader parses `'C3-A2'`, loads both channels, subtracts, and returns a single rereferenced signal.

### Load a specific time window

```matlab
% Load epochs 30 to 60 (0-indexed, in units of the EDF record duration)
[~, ~, signals] = read_EDF('psg.edf', 'Epochs', [30 60]);
```

### Read a compressed EDF

Pass an `.edf.gz` or `.edf.zst` path. The reader streams through zlib / libzstd — no temp file is written.

```matlab
[header, sig_hdr, signals, annot] = read_EDF('psg.edf.gz');
[header, sig_hdr, signals, annot] = read_EDF('psg.edf.zst');
```

`RepairHeader` and `deidentify` cannot modify a compressed archive in place; they're disabled (with a warning) for both `.gz` and `.zst` inputs. Decompress to `.edf` first if you need either.

### Repair a broken file

Clinical EDFs sometimes have invalid `num_data_records` — the field says one number but the file actually contains a different amount of data. Passing `'RepairHeader', true` recomputes the correct value and writes it to a `<filename>_fixed.edf` file:

```matlab
[~, ~, ~] = read_EDF('broken.edf', 'RepairHeader', true);
% Writes broken_fixed.edf
```

### De-identify

```matlab
[~, ~, ~] = read_EDF('patient.edf', 'deidentify', true);
% Writes patient_deidentified.edf with PHI fields blanked
```

### Force pure-MATLAB path (disable MEX)

```matlab
[~, ~, s] = read_EDF('f.edf', 'forceMATLAB', true);
```

### Write an EDF

```matlab
% Round-trip
[h, sh, sc, ann] = read_EDF('in.edf');
write_EDF('out.edf', h, sh, sc, ann);

% Compressed output — choose the format by the file extension
write_EDF('out.edf.gz',  h, sh, sc, ann);                       % gzip, level 6
write_EDF('out.edf.gz',  h, sh, sc, ann, 'GzipLevel', 1);       % gzip, faster
write_EDF('out.edf.zst', h, sh, sc, ann);                       % zstd, level 3
write_EDF('out.edf.zst', h, sh, sc, ann, 'ZstdLevel', 9);       % zstd, smaller

% AutoScale: 'preserve' (default for write_EDF) keeps physical_min/max and
%            clips data to fit;
%            'recompute' sets physical_min/max from the data (no clipping).
write_EDF('out.edf', h, sh, sc, [], 'AutoScale', 'recompute');
```

The default `'preserve'` mode is required for lossless `read → write → read` round-trip; round-tripped signals match the originals to within one digital LSB per channel. Use `'recompute'` whenever the data may exceed the existing `physical_min/max` (e.g. after resampling — see `convert_EDF` below).

**Compression formats.** zstd is the recommended default: roughly 5× faster to compress than gzip level 6 with similar or better ratios, and decompression is also faster. Use `.edf.gz` only when you need to hand the file to a tool that does not understand `.zst`.

### Resample one EDF (read → resample → write)

```matlab
% Default: writes sleep_128Hz.edf.zst next to the input
convert_EDF('sleep.edf', 128);

% Choose a different compression mode
convert_EDF('sleep.edf', 128, 'CompressMode', 'gzip');           % .edf.gz output
convert_EDF('sleep.edf', 128, 'CompressMode', 'gzip', 'GzipLevel', 1);
convert_EDF('sleep.edf', 128, 'CompressMode', 'zstd', 'ZstdLevel', 9);

% No compression
convert_EDF('sleep.edf', 128, 'OutputName', '/tmp/out.edf', 'CompressMode', 'none');
```

Resampling uses SPT's `resample` (Kaiser-windowed sinc FIR, anti-aliased). The annotation channel is preserved; only signal channels are resampled. `target_rate * data_record_duration` must be an integer (the writer needs an integer `samples_in_record`); the function errors with the closest valid alternative if not.

`convert_EDF` defaults `'AutoScale'` to `'recompute'` (vs `write_EDF`'s `'preserve'` default). The anti-aliasing filter can briefly produce samples outside the source file's stored `physical_min/max`; `recompute` widens the range to the actual resampled data so those samples are not clipped on encode.

### Batch many EDFs

```matlab
% Cell array of explicit paths
result = batch_convert_EDF({'a.edf', 'b.edf'}, 128, ...
    'OutputDir', '/path/out', 'Parallel', true);

% Whole directory
result = batch_convert_EDF('/path/in', 128, 'Pattern', '*.edf');

% Text file with one path per line (.txt or .list)
result = batch_convert_EDF('list.txt', 128);
```

`result` is a table with columns `input`, `output`, `status` (`'ok'`/`'failed'`), `elapsed_s`, `error_message`. Per-file failures don't abort the batch.

**Compress a folder of EDFs** — the most common batch use: resample every `.edf` (and `.edf.gz` / `.edf.zst`) in a directory to a uniform target rate and write zstd-compressed copies to an output folder. Runs in parallel across files when a parpool is available.

```matlab
% Default: zstd compression, AutoScale='recompute', parfor across files
result = batch_convert_EDF('/data/edfs_in', 100, ...
    'OutputDir', '/data/edfs_out');

% Pick gzip-1 instead of zstd (e.g. for compatibility with tools that only read .gz)
result = batch_convert_EDF('/data/edfs_in', 100, ...
    'OutputDir',    '/data/edfs_out', ...
    'CompressMode', 'gzip', ...
    'GzipLevel',    1);

% Inspect what happened
nok = sum(strcmp(result.status, 'ok'));
fprintf('%d / %d files converted, total elapsed %.1f s\n', ...
        nok, height(result), sum(result.elapsed_s));
disp(result(strcmp(result.status, 'failed'), {'input', 'error_message'}));
```

### Shell CLI

`bin/convert_edf` is a POSIX shell wrapper for `batch_convert_EDF`. It auto-detects `matlab` on `$PATH` (or in standard install locations on macOS/Linux), invokes a single `matlab -batch`, and exits non-zero if any file failed.

```sh
bin/convert_edf --help

# Single file (zstd by default)
bin/convert_edf -r 128 sleep.edf

# A whole directory, custom output dir
bin/convert_edf -r 128 -o /path/out /path/in

# Switch to gzip output, level 1 (fastest gzip)
bin/convert_edf -r 128 --compress-mode gzip --gzip-level 1 /path/in

# Disable compression and parfor; run quietly
bin/convert_edf -r 128 --compress-mode none --no-parallel /path/file.edf
```

The CLI shells out to MATLAB once per invocation — startup cost (~10 s) is amortized across the whole batch, which is why a single `convert_edf` call across many files is much faster than scripting per-file invocations.

### Inspect the header in a GUI

```matlab
[header, sig_hdr] = read_EDF('f.edf');
[ht, st] = header_gui(header, sig_hdr);
% Opens a UIFIGURE with two tables — main header + per-signal headers
```

## Name-value pairs

| Name | Type | Default | Description |
|---|---|---|---|
| `Channels` | cell array of char | `{}` (all) | Which channels to load; supports A-B rereferencing syntax |
| `Epochs` | 1x2 double | `[]` (all) | `[start_epoch end_epoch]`, 0-indexed |
| `Verbose` | logical | `false` | Print progress messages |
| `RepairHeader` | logical | `false` | Fix invalid `num_data_records`, write `_fixed.edf` |
| `forceMATLAB` | logical | `false` | Disable MEX, use pure MATLAB reader |
| `debug` | logical | `false` | Verbose MEX diagnostics |
| `deidentify` | logical | `false` | Blank PHI fields, write `_deidentified.edf` |

## Files

| File | Role |
|---|---|
| `read_EDF.m` | Read entry point — dispatches to MEX or pure-MATLAB, handles rereferencing, annotations, gz decompression |
| `read_EDF_mex.c` | C source for the read MEX accelerator |
| `read_EDF_mex.mexmaca64` | Pre-built read MEX for Apple Silicon |
| `write_EDF.m` | Write entry point — mirror of read_EDF (gz-aware, MEX + MATLAB fallback) |
| `write_EDF_mex.c` | C source for the write MEX accelerator |
| `write_EDF_mex.mexmaca64` | Pre-built write MEX for Apple Silicon |
| `convert_EDF.m` | Read → resample → write helper (one file) |
| `batch_convert_EDF.m` | Multi-file driver with `parfor` and per-file error handling |
| `bin/convert_edf` | POSIX shell wrapper for `batch_convert_EDF` |
| `compile_edf_mex.m` | Shared auto-compile helper (vendored zlib + system libzstd) |
| `zlib/` | Vendored zlib 1.3.2 source (BSD-style license). Compiled into both MEX files so there's no system zlib dependency. |
| `header_gui.m` | Optional UI for inspecting header + signal-header tables |

If a pre-built MEX isn't available for your platform (Linux, Windows), `read_EDF.m` and `write_EDF.m` will auto-compile on first call. The build pulls in the vendored zlib (no system dep) and links against the system **libzstd** for `.edf.zst` support — the auto-compile script searches `/opt/homebrew` and `/usr/local` for it. To rebuild manually:

```matlab
% From the EDF_toolbox directory
compile_edf_mex(pwd, 'read_EDF_mex.c');
compile_edf_mex(pwd, 'write_EDF_mex.c');
```

## Install

```matlab
addpath('/path/to/EDF_toolbox');
```

When used as part of `preraulab_utilities`, the top-level path setup handles this automatically.

## Dependencies

- **MATLAB R2020a or later**
- **No required toolboxes** for the core `read_EDF.m`
- **libzstd** (only when (re)building the MEX — pre-built `.mexmaca64` files are committed). On macOS: `brew install zstd`. On Debian/Ubuntu: `apt install libzstd-dev`. On RHEL: `dnf install libzstd-devel`.
- `header_gui.m` requires `uifigure` (available in R2016a+ but most polished in R2020a+)
- A C compiler configured for MEX (`mex -setup C`) if you need to rebuild the MEX on a new platform

## Implementation notes

### Why MEX?

EDF files store samples as `int16` and require per-channel digital-to-physical conversion. The naive MATLAB loop is memory-bound and slow on large files. The MEX version:
- Streams records one at a time through a small reusable buffer (peak raw-bytes memory ≈ one record, regardless of file size)
- Decodes int16 → double in tight C loops
- Parses EDF+ annotations in the same pass as the signal data (no second seek-heavy pass)
- Reads/writes `.edf.gz` and `.edf.zst` directly via zlib / libzstd, on the fly with no temp file

Typical speedup for a full-night PSG file: 5-20× over the pure-MATLAB path, depending on channel count.

### gzip vs zstd

zstd at level 3 is the default for `convert_EDF` because on real EDF data it is roughly 5× faster than gzip level 6 with similar or slightly better compression ratios. gzip level 1 is roughly 2× faster than gzip level 6 with a small (~few percent) size penalty, so it is a reasonable fast-path if you must keep `.edf.gz` for compatibility with downstream tools. All three modes produce bit-identical decoded signals.

### Pure-MATLAB fallback

Kicks in when the MEX binary isn't available for your platform and auto-compile fails, or when `'forceMATLAB', true` is passed. Produces bit-identical output; just slower. For `.edf.gz` and `.edf.zst` inputs the MATLAB path decompresses to a temp file first (cleaned up automatically when the call returns) — for `.zst` it shells out to the `zstd` CLI, which must be on `$PATH` for the fallback. Rewriting the MATLAB reader to consume from an in-memory buffer wasn't worth the complexity for the slow path.

## License

BSD 3-Clause. See [`LICENSE`](LICENSE).
