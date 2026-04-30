#!/usr/bin/env bash
# Sweep compression codec/level for the Rust convert_edf pipeline (all cores).
# Reports wall time and total output bytes for each setting.
#
# Run from the EDF_toolbox directory:
#   ./compare_compression.sh [INDIR]
# Defaults: INDIR=/scratch/test_raw_edf

set -uo pipefail

INDIR="${1:-/scratch/test_raw_edf}"
OUT_DIR="/scratch/edf_rust_cmp_out"

echo "Pulling latest..."
git pull --ff-only

if ! command -v cargo >/dev/null 2>&1; then
  echo "ERROR: cargo not found"; exit 1
fi

echo
echo "=== Building convert_edf (release) ==="
( cd rust && cargo build --release ) || { echo "build failed"; exit 1; }
BIN="$(pwd)/rust/target/release/convert_edf"
[ -x "$BIN" ] || { echo "binary not found at $BIN"; exit 1; }

# Total input size for context
in_bytes=$(du -sb "$INDIR" 2>/dev/null | awk '{print $1}')
in_human=$(du -sh "$INDIR" 2>/dev/null | awk '{print $1}')
echo "input dir: $INDIR  ($in_human, $in_bytes bytes)"
echo

results_file=$(mktemp)
printf "%-15s %10s %12s %10s %10s\n" "codec" "wall_s" "out_bytes" "out_size" "ratio_in" | tee "$results_file"

run_one() {
    local label="$1"; shift
    rm -rf "$OUT_DIR" && mkdir -p "$OUT_DIR"
    # Time the pipeline; use seconds-with-decimal via `date +%s.%N`.
    local t0 t1 wall out_bytes out_human ratio
    t0=$(date +%s.%N)
    "$BIN" --batch "$INDIR" -r 100 --out "$OUT_DIR" --jobs 0 "$@" >/dev/null 2>&1
    local rc=$?
    t1=$(date +%s.%N)
    wall=$(awk -v a="$t1" -v b="$t0" 'BEGIN{printf "%.2f", a-b}')
    if [ "$rc" -ne 0 ]; then
        printf "%-15s %10s %12s %10s %10s\n" "$label" "$wall" "FAIL" "-" "-" | tee -a "$results_file"
        return
    fi
    out_bytes=$(du -sb "$OUT_DIR" | awk '{print $1}')
    out_human=$(du -sh "$OUT_DIR" | awk '{print $1}')
    ratio=$(awk -v out="$out_bytes" -v in_b="$in_bytes" 'BEGIN{ if(in_b>0) printf "%.3f", out/in_b; else print "-"}')
    printf "%-15s %10s %12s %10s %10s\n" "$label" "$wall" "$out_bytes" "$out_human" "$ratio" | tee -a "$results_file"
}

echo
echo "=== Sweeping codecs (all cores) ==="
run_one "none"           --compress none
run_one "gzip-1"         --compress gzip --gzip-level 1
run_one "gzip-6"         --compress gzip --gzip-level 6
run_one "gzip-9"         --compress gzip --gzip-level 9
run_one "zstd-1"         --compress zstd --zstd-level 1
run_one "zstd-3"         --compress zstd --zstd-level 3
run_one "zstd-9"         --compress zstd --zstd-level 9
run_one "zstd-19"        --compress zstd --zstd-level 19
run_one "zstd-22"        --compress zstd --zstd-level 22

echo
echo "Summary:"
cat "$results_file"
rm -f "$results_file"

echo
echo "DONE"
