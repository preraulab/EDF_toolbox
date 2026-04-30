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
REF_DIR="/scratch/edf_rust_cmp_ref"

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

# --- Build the lossless reference (uncompressed) once ---
# Every codec's decompressed output must be byte-identical to this. If any
# codec produces different bytes the compression itself is corrupting data
# (which would be a bug in the codec lib, not in our pipeline) -- so this
# also catches non-determinism in the resample/quantize path.
echo
echo "=== Generating uncompressed reference ==="
rm -rf "$REF_DIR" && mkdir -p "$REF_DIR"
"$BIN" --batch "$INDIR" -r 100 --out "$REF_DIR" --jobs 0 --compress none >/dev/null 2>&1 || {
    echo "reference run failed"; exit 1;
}
ref_count=$(ls "$REF_DIR"/*.edf 2>/dev/null | wc -l)
echo "  reference: $ref_count uncompressed EDFs in $REF_DIR"

# --- Decoder helpers ---
# zstd CLI is needed to compare zst outputs. Most dev boxes have it.
HAVE_ZSTD=0
if command -v zstd >/dev/null 2>&1; then HAVE_ZSTD=1; fi
HAVE_ZCAT=0
if command -v zcat >/dev/null 2>&1; then HAVE_ZCAT=1; fi

# Compare every output file to its corresponding REF_DIR file. Returns
# "ok" or "<n>/<m> mismatch".
check_accuracy() {
    local mismatch=0 total=0
    for out in "$OUT_DIR"/*; do
        [ -f "$out" ] || continue
        local base=$(basename "$out")
        # Strip compression extension to find ref name
        local ref_name="$base"
        case "$base" in
            *.edf.gz)  ref_name="${base%.gz}" ;;
            *.edf.zst) ref_name="${base%.zst}" ;;
        esac
        local ref="$REF_DIR/$ref_name"
        [ -f "$ref" ] || { echo "(no ref for $base)"; continue; }
        total=$((total + 1))
        case "$out" in
            *.gz)
                if [ "$HAVE_ZCAT" -eq 1 ]; then
                    if ! zcat "$out" | cmp -s - "$ref"; then mismatch=$((mismatch+1)); fi
                else
                    echo "(no zcat; cannot verify $base)"; continue
                fi ;;
            *.zst)
                if [ "$HAVE_ZSTD" -eq 1 ]; then
                    if ! zstd -d -q -c "$out" | cmp -s - "$ref"; then mismatch=$((mismatch+1)); fi
                else
                    echo "(no zstd CLI; cannot verify $base)"; continue
                fi ;;
            *)
                if ! cmp -s "$out" "$ref"; then mismatch=$((mismatch+1)); fi ;;
        esac
    done
    if [ "$total" -eq 0 ]; then
        echo "no-files"
    elif [ "$mismatch" -eq 0 ]; then
        echo "ok($total)"
    else
        echo "FAIL($mismatch/$total)"
    fi
}

results_file=$(mktemp)
printf "%-15s %10s %12s %10s %10s %12s\n" "codec" "wall_s" "out_bytes" "out_size" "ratio_in" "accuracy" | tee "$results_file"

run_one() {
    local label="$1"; shift
    rm -rf "$OUT_DIR" && mkdir -p "$OUT_DIR"
    local t0 t1 wall out_bytes out_human ratio accuracy
    t0=$(date +%s.%N)
    "$BIN" --batch "$INDIR" -r 100 --out "$OUT_DIR" --jobs 0 "$@" >/dev/null 2>&1
    local rc=$?
    t1=$(date +%s.%N)
    wall=$(awk -v a="$t1" -v b="$t0" 'BEGIN{printf "%.2f", a-b}')
    if [ "$rc" -ne 0 ]; then
        printf "%-15s %10s %12s %10s %10s %12s\n" "$label" "$wall" "FAIL" "-" "-" "-" | tee -a "$results_file"
        return
    fi
    out_bytes=$(du -sb "$OUT_DIR" | awk '{print $1}')
    out_human=$(du -sh "$OUT_DIR" | awk '{print $1}')
    ratio=$(awk -v out="$out_bytes" -v in_b="$in_bytes" 'BEGIN{ if(in_b>0) printf "%.3f", out/in_b; else print "-"}')
    accuracy=$(check_accuracy)
    printf "%-15s %10s %12s %10s %10s %12s\n" "$label" "$wall" "$out_bytes" "$out_human" "$ratio" "$accuracy" | tee -a "$results_file"
}

echo
echo "=== Sweeping codecs (all cores; accuracy = decompressed bytes match uncompressed reference) ==="
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
