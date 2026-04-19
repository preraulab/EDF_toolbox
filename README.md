# read_EDF

High-performance EDF / EDF+ reader for MATLAB, with a compiled MEX backend and a pure-MATLAB fallback.

Handles the things that break off-the-shelf EDF readers in practice: malformed record counts in clinical files, EDF+ TAL annotations, A-B rereferencing at load time, de-identification, and large polysomnography recordings where a naive loop is slow.

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
| **Channel subsetting** | Load only specific channels by name |
| **A-B rereferencing** | `'EEG C3-A2'` — compute on-the-fly during load |
| **Epoch subsetting** | Load only a `[start_epoch end_epoch]` slice |
| **EDF+ annotations** | Parse TAL (Time-Annotation List) format |
| **Header repair** | Correct invalid `num_data_records` and save a `_fixed` copy |
| **De-identification** | Strip PHI fields and save a `_deidentified` copy |
| **Digital → physical scaling** | `phys_val = phys_min + (dig_val - dig_min) × (phys_max - phys_min) / (dig_max - dig_min)` — done correctly per-channel |

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
| `read_EDF.m` | Main entry point — dispatches to MEX or pure-MATLAB, handles rereferencing, annotations |
| `read_EDF_mex.c` | C source for the MEX accelerator |
| `read_EDF_mex.mexmaca64` | Pre-built MEX binary for Apple Silicon |
| `read_EDF_mex.mexa64` | Pre-built MEX binary for Linux x86_64 |
| `header_gui.m` | Optional UI for inspecting header + signal-header tables |

No pre-built Windows MEX is shipped; MATLAB will transparently fall back to `read_EDF.m` (pure MATLAB), or you can build your own:

```matlab
mex read_EDF_mex.c
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
- Reads file regions in one `fread` per record
- Does the int16 → double conversion in tight C loops
- Handles annotation channels separately (they're not numeric data)
- Uses less peak memory because it doesn't build MATLAB arrays until the end

Typical speedup for a full-night PSG file: 5-20× over the pure-MATLAB path, depending on channel count.

### Pure-MATLAB fallback

Kicks in when the MEX binary isn't available for your platform (Windows unless you build it) or when `'forceMATLAB', true` is passed. Produces bit-identical output; just slower.

## License

BSD 3-Clause. See [`LICENSE`](LICENSE).
