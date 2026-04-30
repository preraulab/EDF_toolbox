#!/usr/bin/env bash
# Diagnostic sweep for batch_convert_EDF parallel performance.
#
# Probes the system, runs disk + iostat/vmstat monitors, and runs a MATLAB
# sweep across (Workers, WorkerThreads, StageDir) configurations.
#
# The script lives in EDF_toolbox/, so run it from here:
#   ./diag_edf.sh [INDIR]
# Defaults: INDIR=/scratch/edf
# Output dir: /scratch/edf_diag_<timestamp>

set -uo pipefail

TOOLBOX="$(pwd)"
INDIR="${1:-/scratch/edf}"
OUT="/scratch/edf_diag_$(date +%Y%m%d_%H%M%S)"

mkdir -p "$OUT"
echo "Toolbox    : $TOOLBOX"
echo "Input dir  : $INDIR"
echo "Output dir : $OUT"

# ============================================================================
# Phase A: system probes
# ============================================================================
{
  echo "==== date ===="; date
  echo; echo "==== hostname ===="; hostname
  echo; echo "==== uname -a ===="; uname -a
  echo; echo "==== lscpu ===="; lscpu | head -25
  echo; echo "==== free -h ===="; free -h
  echo; echo "==== nproc ===="; nproc
  echo; echo "==== mount (auto/nfs) ===="; mount | grep -E 'auto|nfs' || echo "(none)"
  echo; echo "==== df -hT (relevant paths) ===="
  for p in "$INDIR" /scratch /dev/shm /tmp; do
    df -hT "$p" 2>/dev/null | tail -n +2
  done
  echo; echo "==== readlink resolves ===="
  for p in "$INDIR" /scratch /dev/shm "$TOOLBOX"; do
    printf "  %-40s -> %s\n" "$p" "$(readlink -f "$p" 2>/dev/null || echo MISSING)"
  done
  echo; echo "==== input listing ===="
  ls -lah "$INDIR" | head -40
  echo; echo "==== input total size ===="; du -sh "$INDIR"
  echo; echo "==== /dev/shm capacity ===="; df -h /dev/shm
} | tee "$OUT/system_info.txt"

# ============================================================================
# Phase B: raw disk read bandwidth (two reads, second should be cache-warm)
# ============================================================================
FIRST=$(ls "$INDIR"/*.edf 2>/dev/null | head -1)
if [ -n "$FIRST" ]; then
  echo
  echo "==== dd read 1 ===="
  ( time dd if="$FIRST" of=/dev/null bs=8M ) 2>&1 | tee "$OUT/dd_run1.txt"
  echo "==== dd read 2 (warm cache) ===="
  ( time dd if="$FIRST" of=/dev/null bs=8M ) 2>&1 | tee "$OUT/dd_run2.txt"
fi

# ============================================================================
# Phase C: iostat + vmstat monitors during MATLAB run
# ============================================================================
if command -v iostat >/dev/null 2>&1; then
  iostat -x 5 > "$OUT/iostat.log" 2>&1 &
  IOPID=$!
else
  IOPID=
fi
vmstat 5 > "$OUT/vmstat.log" 2>&1 &
VMPID=$!
cleanup() { [ -n "$IOPID" ] && kill "$IOPID" 2>/dev/null; kill "$VMPID" 2>/dev/null; true; }
trap cleanup EXIT

# ============================================================================
# Phase D: write the diag_sweep.m and run it under MATLAB
# ============================================================================
cat > "$OUT/diag_sweep.m" <<'MATLAB_EOF'
function diag_sweep(outdir, toolbox, indir)
addpath(genpath(toolbox));
fprintf('=== diag_sweep starting at %s ===\n', char(datetime));

files = dir(fullfile(indir, '*.edf'));
files = arrayfun(@(x) fullfile(x.folder, x.name), files, 'uni', false);
fprintf('outdir=%s\n  toolbox=%s\n  indir=%s\n  %d files\n\n', ...
    outdir, toolbox, indir, numel(files));

% --- Pre-compile MEX on head node ---
% read_EDF / write_EDF auto-compile MEX on first call when the platform
% .mex file is missing. Without this step, every parfor worker would
% race to write the same .mexa64 file -> corruption + N x compile cost.
fprintf('--- pre-compiling MEX (head node only) ---\n');
read_mex  = fullfile(toolbox, ['read_EDF_mex.'  mexext]);
write_mex = fullfile(toolbox, ['write_EDF_mex.' mexext]);
if ~isfile(read_mex)
    compile_edf_mex(toolbox, 'read_EDF_mex.c');
end
if ~isfile(write_mex)
    compile_edf_mex(toolbox, 'write_EDF_mex.c');
end
fprintf('  read_mex  : %s (%d bytes)\n', read_mex,  file_size_safe(read_mex));
fprintf('  write_mex : %s (%d bytes)\n', write_mex, file_size_safe(write_mex));

% --- BLAS pinning sanity check ---
fprintf('--- BLAS pinning check ---\n');
delete(gcp('nocreate'));
parpool('Processes', 4);
fut = parfevalOnAll(gcp, @maxNumCompThreads, 1); wait(fut);
fprintf('default worker BLAS threads: %s\n', mat2str(fetchOutputs(fut)));
fut = parfevalOnAll(gcp, @maxNumCompThreads, 0, 1); wait(fut);
fut = parfevalOnAll(gcp, @maxNumCompThreads, 1); wait(fut);
fprintf('after WT=1, worker BLAS threads: %s\n', mat2str(fetchOutputs(fut)));
delete(gcp('nocreate'));

% --- Read-only contention test (isolates I/O from compute) ---
fprintf('\n--- read-only contention ---\n');
t = tic;
for k = 1:numel(files); [~,~] = read_EDF(files{k}); end
t_ser = toc(t);
fprintf('serial read-only: %.1fs (%.2fs/file)\n', t_ser, t_ser/numel(files));
read_rows = {1, t_ser, 1.0};
for W = [2 4 8]
    delete(gcp('nocreate'));
    parpool('Processes', W);
    fut = parfevalOnAll(gcp, @maxNumCompThreads, 0, 1); wait(fut);
    t = tic;
    parfor k = 1:numel(files); [~,~] = read_EDF(files{k}); end
    tw = toc(t);
    fprintf('parallel read-only W=%d: %.1fs (speedup=%.2fx)\n', W, tw, t_ser/tw);
    read_rows(end+1, :) = {W, tw, t_ser/tw}; %#ok<AGROW>
end
writetable(cell2table(read_rows, 'VariableNames', {'workers','wall_s','speedup'}), ...
    fullfile(outdir, 'read_only.csv'));
delete(gcp('nocreate'));

% --- Full convert sweep ---
fprintf('\n--- full convert sweep ---\n');
configs = {
    'serial',     1, 1, '';
    'W2_WT1',     2, 1, '';
    'W4_WT1',     4, 1, '';
    'W8_WT1',     8, 1, '';
    'W2_WT2',     2, 2, '';
    'W2_WT4',     2, 4, '';
    'W4_WT2',     4, 2, '';
    'W2_WT1_shm', 2, 1, '/dev/shm';
    'W4_WT1_shm', 4, 1, '/dev/shm';
};
sweep_rows = {};
for i = 1:size(configs,1)
    name = configs{i,1}; W = configs{i,2}; WT = configs{i,3}; stage = configs{i,4};
    delete(gcp('nocreate'));
    if W > 1, parpool('Processes', W); end
    out_dir = fullfile(outdir, ['out_' name]);
    if isfolder(out_dir), rmdir(out_dir, 's'); end
    mkdir(out_dir);
    args = {'Parallel', W>1, 'WorkerThreads', WT, ...
            'OutputDir', out_dir, 'CompressMode', 'zstd', ...
            'AutoScale', 'recompute'};
    if ~isempty(stage)
        args = [args, {'StageLocal', true, 'StageDir', stage}];
    end
    t = tic;
    try
        r = batch_convert_EDF(indir, 100, args{:});
        wall = toc(t);
        nok = sum(strcmp(r.status,'ok')); ntot = height(r);
    catch ME
        wall = toc(t); nok = -1; ntot = -1;
        fprintf('  ERROR in %s: %s\n', name, ME.message);
    end
    fprintf('[%s] W=%d WT=%d stage=%s wall=%.1fs ok=%d/%d\n', ...
        name, W, WT, stage, wall, nok, ntot);
    sweep_rows(end+1, :) = {name, W, WT, stage, wall, nok, ntot}; %#ok<AGROW>
    if isfolder(out_dir), rmdir(out_dir, 's'); end
end
delete(gcp('nocreate'));

T = cell2table(sweep_rows, 'VariableNames', ...
    {'config','workers','worker_threads','stage','wall_s','n_ok','n_total'});
writetable(T, fullfile(outdir, 'sweep_results.csv'));
disp(T);
fprintf('\n=== diag_sweep done at %s ===\n', char(datetime));
end

function n = file_size_safe(p)
if isfile(p), d = dir(p); n = d.bytes; else, n = -1; end
end
MATLAB_EOF

matlab -batch "addpath('$OUT'); diag_sweep('$OUT','$TOOLBOX','$INDIR')" 2>&1 | tee "$OUT/matlab.log"

cleanup
echo
echo "DONE. Results in $OUT/"
ls -la "$OUT/"
