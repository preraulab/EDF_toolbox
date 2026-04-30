#!/usr/bin/env bash
# Compare the standalone Rust binary against the MATLAB pipeline on the
# same input directory. Run from the EDF_toolbox directory.
#
#   ./compare_rust_vs_matlab.sh [INDIR]
# Defaults: INDIR=/scratch/test_raw_edf

set -uo pipefail

INDIR="${1:-/scratch/test_raw_edf}"
RUST_OUT="/scratch/edf_rust_out"
MATLAB_OUT="/scratch/edf_matlab_out"

echo "Pulling latest..."
git pull --ff-only

# --- Toolchain check ------------------------------------------------------
if ! command -v cargo >/dev/null 2>&1; then
  echo
  echo "ERROR: cargo not found. Install Rust first:"
  echo "  curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y"
  echo "  source \"\$HOME/.cargo/env\""
  exit 1
fi
echo "rustc: $(rustc --version)"
echo "cargo: $(cargo --version)"

# --- Build ----------------------------------------------------------------
echo
echo "=== Building convert_edf (release) ==="
( cd rust && cargo build --release ) || { echo "build failed"; exit 1; }
BIN="$(pwd)/rust/target/release/convert_edf"
[ -x "$BIN" ] || { echo "binary not found at $BIN"; exit 1; }
echo "binary: $BIN"

# --- Run ------------------------------------------------------------------
mkdir -p "$RUST_OUT" "$MATLAB_OUT"

echo
echo "=== Rust: 1 thread (zstd-3) ==="
rm -f "$RUST_OUT"/*.edf*
( time "$BIN" --batch "$INDIR" -r 100 --out "$RUST_OUT" --jobs 1 --compress zstd --zstd-level 3 ) 2>&1 | tail -5

echo
echo "=== Rust: 4 threads (zstd-3) ==="
rm -f "$RUST_OUT"/*.edf*
( time "$BIN" --batch "$INDIR" -r 100 --out "$RUST_OUT" --jobs 4 --compress zstd --zstd-level 3 ) 2>&1 | tail -5

echo
echo "=== Rust: all-cores (zstd-3) ==="
rm -f "$RUST_OUT"/*.edf*
( time "$BIN" --batch "$INDIR" -r 100 --out "$RUST_OUT" --jobs 0 --compress zstd --zstd-level 3 ) 2>&1 | tail -5

echo
echo "=== MATLAB: serial ==="
rm -rf "$MATLAB_OUT" && mkdir -p "$MATLAB_OUT"
matlab -batch "addpath(genpath(pwd)); delete(gcp('nocreate')); t=tic; r=batch_convert_EDF('$INDIR',100,'Parallel',false,'OutputDir','$MATLAB_OUT'); fprintf('matlab serial: %.1fs ok=%d/%d\n',toc(t),sum(strcmp(r.status,'ok')),height(r))" 2>&1 | grep -E "matlab|Error" | tail -5

echo
echo "=== MATLAB: W=4 WT=2 ==="
rm -rf "$MATLAB_OUT" && mkdir -p "$MATLAB_OUT"
matlab -batch "addpath(genpath(pwd)); delete(gcp('nocreate')); parpool('Processes',4); t=tic; r=batch_convert_EDF('$INDIR',100,'Parallel',true,'WorkerThreads',2,'OutputDir','$MATLAB_OUT'); fprintf('matlab w4 wt2: %.1fs ok=%d/%d\n',toc(t),sum(strcmp(r.status,'ok')),height(r)); delete(gcp('nocreate'))" 2>&1 | grep -E "matlab|Error" | tail -5

echo
echo "DONE"
