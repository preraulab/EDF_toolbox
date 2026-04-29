function results = bench_convert_modes(in_dir, out_dir, target_rate, varargin)
%BENCH_CONVERT_MODES  Time single-stage convert_EDF in 3 modes on the same corpus.
%
%   results = bench_convert_modes(in_dir, out_dir, target_rate)
%   results = bench_convert_modes(..., 'Workers', N, 'StageDir', D, 'Runs', R)
%
%   Modes tested (single-stage = read+resample+compress in one MEX pass):
%     1. serial                    : Parallel=false, no staging
%     2. parallel                  : Parallel=true,  no staging
%     3. parallel + StageLocal     : Parallel=true,  StageLocal=true
%     4. serial   + StageLocal     : Parallel=false, StageLocal=true (control)
%
%   For NFS-backed inputs the 'parallel' mode is the one that typically
%   suffers from server contention. Mode 3 should claw the speedup back.
%
%   Inputs:
%       in_dir      : directory of input .edf / .edf.gz / .edf.zst files
%       out_dir     : where outputs go. Cleared between modes.
%       target_rate : target sample rate (Hz)
%
%   Name-value pairs:
%       'Workers'      : parpool size (default: feature('numcores'))
%       'WorkerThreads': BLAS threads per worker (default 1). See
%                        batch_convert_EDF docs for why 1 is the right
%                        default.
%       'StageDir'     : local scratch dir (default: '' -> tempdir)
%       'Runs'         : runs per mode, takes min (default: 2)
%       'CompressMode' : 'zstd' (default) | 'gzip' | 'none'
%       'AutoScale'    : 'recompute' (default) | 'preserve'
%
%   Output:
%       results : table of (mode, run, wall_seconds, n_ok, n_failed)
%   Also writes <out_dir>/bench_modes_raw.csv and bench_modes_summary.csv.

p = inputParser;
addRequired(p, 'in_dir');
addRequired(p, 'out_dir');
addRequired(p, 'target_rate', @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'Workers',       feature('numcores'), @(x) isnumeric(x) && x>=1);
addParameter(p, 'WorkerThreads', 1, @(x) isnumeric(x) && x>=1);
addParameter(p, 'StageDir',      '', @ischar);
addParameter(p, 'Runs',          2, @(x) isnumeric(x) && x>=1);
addParameter(p, 'CompressMode',  'zstd', @(s) any(strcmpi(s,{'gzip','zstd','none'})));
addParameter(p, 'AutoScale',     'recompute', @(s) any(strcmpi(s,{'preserve','recompute'})));
parse(p, in_dir, out_dir, target_rate, varargin{:});

in_dir         = char(p.Results.in_dir);
out_dir        = char(p.Results.out_dir);
target_rate    = double(p.Results.target_rate);
n_workers      = double(p.Results.Workers);
worker_threads = double(p.Results.WorkerThreads);
stage_dir      = p.Results.StageDir;
n_runs         = double(p.Results.Runs);
compress_mode  = p.Results.CompressMode;
autoscale      = p.Results.AutoScale;

% --- Resolve corpus -----------------------------------------------------
d = [dir(fullfile(in_dir,'*.edf')); ...
     dir(fullfile(in_dir,'*.edf.gz')); ...
     dir(fullfile(in_dir,'*.edf.zst'))];
files = arrayfun(@(x) fullfile(x.folder, x.name), d, 'UniformOutput', false);
if isempty(files)
    error('bench_convert_modes:NoFiles', 'No .edf* files in %s', in_dir);
end
n_files = numel(files);

if ~isfolder(out_dir), mkdir(out_dir); end
total_in = sum(arrayfun(@(x) x.bytes, d));
fprintf('Corpus    : %d files, %.1f MB total\n', n_files, total_in/1e6);
fprintf('Target    : %g Hz, CompressMode=%s, AutoScale=%s\n', target_rate, compress_mode, autoscale);
fprintf('Workers   : %d  WorkerThreads: %d  StageDir: %s\n', ...
    n_workers, worker_threads, ifempty(stage_dir, '<tempdir>'));
fprintf('Runs/mode : %d (best taken)\n\n', n_runs);

% --- Set up parpool of the requested size up front ----------------------
need_pool = true;
pool = gcp('nocreate');
if ~isempty(pool) && pool.NumWorkers == n_workers
    need_pool = false;
end
if need_pool
    if ~isempty(pool), delete(pool); end
    fprintf('Starting parpool with %d workers...\n', n_workers);
    parpool(n_workers);
end

% --- Modes to test ------------------------------------------------------
modes = {
    'serial',                 false, false;
    'parallel',               true,  false;
    'parallel+stage_local',   true,  true;
    'serial+stage_local',     false, true;
};

raw_rows = {};
fprintf('\n%-22s | %-3s | %10s | %5s | %6s\n', 'mode','run','wall_s','n_ok','n_fail');
fprintf('%s\n', repmat('-', 1, 60));

for m = 1:size(modes, 1)
    mode_name   = modes{m, 1};
    do_parallel = modes{m, 2};
    stage_local = modes{m, 3};

    for run = 1:n_runs
        clear_outputs(out_dir);
        % Sync filesystem to flush previous-run dirty pages
        try, system('sync'); catch, end

        t0 = tic;
        opts = {'OutputDir', out_dir, ...
                'CompressMode', compress_mode, ...
                'AutoScale', autoscale, ...
                'Parallel', do_parallel, ...
                'WorkerThreads', worker_threads, ...
                'StageLocal', stage_local};
        if stage_local && ~isempty(stage_dir)
            opts = [opts, {'StageDir', stage_dir}];
        end
        r = batch_convert_EDF(in_dir, target_rate, opts{:});
        wall_s = toc(t0);

        n_ok   = sum(strcmp(r.status, 'ok'));
        n_fail = height(r) - n_ok;

        fprintf('%-22s | %3d | %10.2f | %5d | %6d\n', ...
            mode_name, run, wall_s, n_ok, n_fail);

        raw_rows(end+1, :) = {mode_name, run, wall_s, n_ok, n_fail, ...
            do_parallel, stage_local}; %#ok<AGROW>
    end
    fprintf('%s\n', repmat('-', 1, 60));
end

raw = cell2table(raw_rows, 'VariableNames', ...
    {'mode','run','wall_s','n_ok','n_fail','parallel','stage_local'});

% --- Best per mode ------------------------------------------------------
mode_names = unique(raw.mode, 'stable');
sum_rows = {};
for k = 1:numel(mode_names)
    rows = raw(strcmp(raw.mode, mode_names{k}), :);
    [best_s, idx] = min(rows.wall_s);
    sum_rows(end+1, :) = {mode_names{k}, best_s, rows.n_ok(idx), ...
        rows.n_fail(idx), (total_in/1e6)/best_s}; %#ok<AGROW>
end
summary = cell2table(sum_rows, 'VariableNames', ...
    {'mode','best_wall_s','n_ok','n_fail','MB_per_s_in'});

% Speedup vs serial baseline
serial_row = summary(strcmp(summary.mode, 'serial'), :);
if ~isempty(serial_row)
    summary.speedup_vs_serial = serial_row.best_wall_s ./ summary.best_wall_s;
end

fprintf('\n=========  SUMMARY (best of %d, sorted)  =========\n', n_runs);
disp(sortrows(summary, 'best_wall_s'));

% --- Write CSVs ---------------------------------------------------------
raw_csv     = fullfile(out_dir, 'bench_modes_raw.csv');
summary_csv = fullfile(out_dir, 'bench_modes_summary.csv');
writetable(raw, raw_csv);
writetable(summary, summary_csv);
fprintf('\nSaved: %s\n', raw_csv);
fprintf('Saved: %s\n', summary_csv);

results = raw;
end


% =========================================================================
function clear_outputs(out_dir)
% Remove any *.edf, *.edf.gz, *.edf.zst already in out_dir so each run
% starts from a clean slate. Keeps the dir itself.
if ~isfolder(out_dir), return; end
for pat = {'*.edf','*.edf.gz','*.edf.zst'}
    d = dir(fullfile(out_dir, pat{1}));
    for k = 1:numel(d)
        try, delete(fullfile(d(k).folder, d(k).name)); catch, end
    end
end
end


% =========================================================================
function s = ifempty(s, dflt)
if isempty(s), s = dflt; end
end
