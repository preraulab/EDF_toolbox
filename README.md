# read_EDF

High-performance EDF / EDF+ reader for MATLAB, with a compiled MEX backend and a pure-MATLAB fallback. Reads plain `.edf` and gzip-compressed `.edf.gz` (on the fly, no temp file).

Handles the things that break off-the-shelf EDF readers in practice: malformed record counts in clinical files, EDF+ TAL annotations, A-B rereferencing at load time, de-identification, compressed archives, and large polysomnography recordings where a naive loop is slow.

Part of the Prerau Lab [`preraulab_utilities`](https://github.com/preraulab/preraulab_utilities) meta-repository. Can also be used standalone.

## Quick start

```matlab
[header, signal_header, signal_cell, annotations] = read_EDF('sleep.edf');
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
| **Compressed input** | Reads `.edf.gz` directly via vendored zlib — no temp file, streaming decompression |
| **Channel subsetting** | Load only specific channels by name |
| **A-B rereferencing** | `'EEG C3-A2'` — compute on-the-fly during load |
| **Epoch subsetting** | Load only a `[start_epoch end_epoch]` slice |
| **EDF+ annotations** | Parse TAL (Time-Annotation List) format |
| **Header repair** | Correct invalid `num_data_records` and save a `_fixed` copy (plain `.edf` only) |
| **De-identification** | Strip PHI fields and save a `_deidentified` copy (plain `.edf` only) |
| **Digital → physical scaling** | `phys_val = phys_min + (dig_val - dig_min) × (phys_max - phys_min) / (dig_max - dig_min)` — done correctly per-channel |
| **Streaming I/O** | Reads record-by-record; peak memory is ~one record + the output arrays, even for multi-GB files |

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

### Read a gzip-compressed EDF

Pass an `.edf.gz` path. The reader streams through zlib — no temp file is written.

```matlab
[header, sig_hdr, signals, annot] = read_EDF('psg.edf.gz');
```

`RepairHeader` and `deidentify` cannot modify a gzip archive in place; they're disabled (with a warning) for `.gz` inputs. Decompress to `.edf` first if you need either.

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
| `read_EDF.m` | Main entry point — dispatches to MEX or pure-MATLAB, handles rereferencing, annotations, gz decompression |
| `read_EDF_mex.c` | C source for the MEX accelerator (uses bundled zlib for `.edf.gz`) |
| `read_EDF_mex.mexmaca64` | Pre-built MEX binary for Apple Silicon |
| `zlib/` | Vendored zlib 1.3.2 source (BSD-style license). Compiled into the MEX so there's no system zlib dependency on any platform. |
| `header_gui.m` | Optional UI for inspecting header + signal-header tables |

If a pre-built MEX isn't available for your platform (Linux, Windows), `read_EDF.m` will auto-compile on first call. The build pulls in the vendored zlib, so it works the same on macOS, Linux, and Windows with no system dependencies. To rebuild manually:

```matlab
% From the read_EDF directory
mex -O -largeArrayDims -Izlib read_EDF_mex.c zlib/*.c
```

## Install

```matlab
addpath('/path/to/read_EDF');
```

When used as part of `preraulab_utilities`, the top-level path setup handles this automatically.

## Dependencies

- **MATLAB R2020a or later**
- **No required toolboxes** for the core `read_EDF.m`
- `header_gui.m` requires `uifigure` (available in R2016a+ but most polished in R2020a+)
- A C compiler configured for MEX (`mex -setup C`) if you need to rebuild the MEX on a new platform

## Implementation notes

### Why MEX?

EDF files store samples as `int16` and require per-channel digital-to-physical conversion. The naive MATLAB loop is memory-bound and slow on large files. The MEX version:
- Streams records one at a time through a small reusable buffer (peak raw-bytes memory ≈ one record, regardless of file size)
- Decodes int16 → double in tight C loops
- Parses EDF+ annotations in the same pass as the signal data (no second seek-heavy pass)
- Reads `.edf.gz` directly via zlib, decompressing on the fly with no temp file

Typical speedup for a full-night PSG file: 5-20× over the pure-MATLAB path, depending on channel count. Reading `.edf.gz` adds the cost of zlib decompression (~150-300 MB/s of uncompressed throughput, single-core) but avoids the I/O cost of reading the larger plain file.

### Pure-MATLAB fallback

Kicks in when the MEX binary isn't available for your platform and auto-compile fails, or when `'forceMATLAB', true` is passed. Produces bit-identical output; just slower. For `.edf.gz` inputs the MATLAB path decompresses to a temp file first (cleaned up automatically when the call returns), since rewriting the MATLAB reader to consume from an in-memory buffer wasn't worth the complexity for the slow path.

## License

BSD 3-Clause. See [`LICENSE`](LICENSE).
