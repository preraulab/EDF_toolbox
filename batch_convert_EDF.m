function result = batch_convert_EDF(input_spec, target_rate, varargin)
%BATCH_CONVERT_EDF  Resample-and-rewrite a batch of EDF files in parallel.
%
%   result = batch_convert_EDF(input_spec, target_rate, ...)
%
%   Inputs:
%       input_spec  : one of
%                       - cell array of file paths
%                       - directory string (uses 'Pattern' to glob)
%                       - text file path with one filename per line
%       target_rate : target sample rate (Hz) for all signals
%
%   Name-value pairs:
%       'OutputDir'    : '' (default: alongside each input file)
%       'CompressMode' : 'zstd' (default) | 'gzip' | 'none'
%       'Compress'     : (deprecated) true -> 'zstd', false -> 'none'.
%                        'CompressMode' takes precedence if both are passed.
%       'GzipLevel'    : integer 1..9 (default 6)
%       'ZstdLevel'    : integer 1..22 (default 9)
%       'Parallel'     : true (default if a parpool exists). When true and no
%                        pool exists, attempts to start one; falls back to
%                        serial on failure.
%       'WorkerThreads': integer >= 1 (default 1). BLAS / computational
%                        threads per parpool worker. Default 1 prevents the
%                        N_workers x N_cores oversubscription that
%                        otherwise tanks resample's FIR filter performance
%                        under parfor on multi-core boxes. Increase only if
%                        your workload is not BLAS-heavy and you have spare
%                        cores. Ignored when 'Parallel' is false.
%       'StageLocal'   : false (default). When true, each file is first
%                        copied to a local scratch directory, converted
%                        there, then the output is moved to its final
%                        location. Use this when inputs/outputs live on
%                        NFS / shared storage where concurrent reads from
%                        many parfor workers would contend at the file
%                        server. Disk space required per worker = one
%                        input file + one output file at a time.
%       'StageDir'     : '' (default: system temp dir from `tempname`).
%                        Only used when 'StageLocal' is true.
%       'Pattern'      : '*.edf' (default; matches .edf, .edf.gz, .edf.zst)
%       'Verbose'      : false (default)
%       'AutoScale'    : 'recompute' (default) | 'preserve'. See convert_EDF
%                        for the rationale on why 'recompute' is the default.
%
%   Output:
%       result : table with columns
%                  input          : char (full path)
%                  output         : char (full path) or '' on failure
%                  status         : 'ok' or 'failed'
%                  elapsed_s      : double (seconds, NaN on early failure)
%                  error_message  : char ('' on success)
%
%   Per-file errors are caught and recorded; the batch always completes
%   for the remaining files.
%
%   See also: convert_EDF, read_EDF, write_EDF.

p = inputParser;
addRequired(p,  'input_spec');
addRequired(p,  'target_rate', @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'OutputDir',    '', @ischar);
addParameter(p, 'CompressMode', '', @(s) ischar(s) && (isempty(s) || any(strcmpi(s, {'gzip','zstd','none'}))));
addParameter(p, 'Compress',     true, @islogical);
addParameter(p, 'GzipLevel',    6, @(x) isnumeric(x) && isscalar(x) && x >= 1 && x <= 9);
addParameter(p, 'ZstdLevel',    9, @(x) isnumeric(x) && isscalar(x) && x >= 1 && x <= 22);
addParameter(p, 'Parallel',     [], @(x) isempty(x) || islogical(x));
addParameter(p, 'WorkerThreads',1, @(x) isnumeric(x) && isscalar(x) && x >= 1);
addParameter(p, 'StageLocal',   false, @islogical);
addParameter(p, 'StageDir',     '', @ischar);
addParameter(p, 'Pattern',      '*.edf', @ischar);
addParameter(p, 'Verbose',      false, @islogical);
addParameter(p, 'AutoScale',    'recompute', @(s) any(strcmpi(s, {'preserve','recompute'})));
parse(p, input_spec, target_rate, varargin{:});

target_rate  = double(p.Results.target_rate);
out_dir      = p.Results.OutputDir;
compress_mode= lower(p.Results.CompressMode);
compress     = p.Results.Compress;
gzip_level   = double(p.Results.GzipLevel);
zstd_level   = double(p.Results.ZstdLevel);
do_parallel  = p.Results.Parallel;
worker_threads = round(double(p.Results.WorkerThreads));
stage_local  = p.Results.StageLocal;
stage_dir    = p.Results.StageDir;
pattern      = p.Results.Pattern;
verbose      = p.Results.Verbose;
autoscale    = p.Results.AutoScale;

if isempty(compress_mode)
    if compress, compress_mode = 'zstd'; else, compress_mode = 'none'; end
end

% --- Resolve input list ---------------------------------------------------
files = resolve_inputs(input_spec, pattern);
if isempty(files)
    warning('batch_convert_EDF:Empty', 'No input files resolved.');
    result = empty_result_table();
    return
end

if verbose
    fprintf('batch_convert_EDF: %d files @ %g Hz\n', numel(files), target_rate);
end

% --- Decide on parallel ---------------------------------------------------
if isempty(do_parallel)
    do_parallel = ~isempty(gcp('nocreate'));
end
if do_parallel
    pool = gcp('nocreate');
    if isempty(pool)
        try
            pool = parpool;
            if verbose, fprintf('  started parpool with %d workers\n', pool.NumWorkers); end
        catch ME
            if verbose, fprintf('  parpool failed (%s); running serial\n', ME.message); end
            do_parallel = false;
        end
    end
    % Pin each worker to worker_threads BLAS threads. Without this,
    % every Process-profile worker tries to use all cores for BLAS, so
    % N_workers x N_cores threads compete for N_cores -- resample tanks
    % and parfor runs slower than serial. Default worker_threads=1.
    %
    % parfevalOnAll is asynchronous, so we MUST wait() on the future
    % before launching the parfor below -- otherwise the parfor may
    % start while workers still have their default (all-cores) BLAS
    % thread count, producing the very oversubscription this guards
    % against. Skipping the wait was the root cause of bimodal "fast
    % some runs, very slow other runs" timings under parfor.
    if do_parallel && ~isempty(pool)
        try
            fut = parfevalOnAll(pool, @maxNumCompThreads, 0, worker_threads);
            wait(fut);
        catch ME
            if verbose
                fprintf('  could not set worker BLAS threads (%s); continuing\n', ME.message);
            end
        end
    end
end

% --- Run -----------------------------------------------------------------
n = numel(files);
out_paths      = cell(n, 1);
statuses       = cell(n, 1);
elapsed_s      = nan(n, 1);
error_messages = cell(n, 1);
[out_paths{:}]      = deal('');
[statuses{:}]       = deal('failed');
[error_messages{:}] = deal('');

if do_parallel
    parfor k = 1:n
        [out_paths{k}, statuses{k}, elapsed_s(k), error_messages{k}] = ...
            run_one(files{k}, target_rate, out_dir, compress_mode, ...
                    gzip_level, zstd_level, autoscale, ...
                    stage_local, stage_dir, false);
    end
else
    for k = 1:n
        if verbose, fprintf('  [%d/%d] %s\n', k, n, files{k}); end
        [out_paths{k}, statuses{k}, elapsed_s(k), error_messages{k}] = ...
            run_one(files{k}, target_rate, out_dir, compress_mode, ...
                    gzip_level, zstd_level, autoscale, ...
                    stage_local, stage_dir, verbose);
    end
end

result = table(files(:), out_paths, statuses, elapsed_s, error_messages, ...
    'VariableNames', {'input', 'output', 'status', 'elapsed_s', 'error_message'});

if verbose
    nok = sum(strcmp(statuses, 'ok'));
    fprintf('batch_convert_EDF: %d/%d ok, %d failed\n', nok, n, n - nok);
end
end


%% =========================================================================
function [out_path, status, elapsed, errmsg] = run_one(in_fname, target_rate, ...
    out_dir, compress_mode, gzip_level, zstd_level, autoscale, ...
    stage_local, stage_dir, verbose)

t0 = tic;
out_path = '';
status = 'failed';
errmsg = '';

try
    [dir_, base, ext] = fileparts(in_fname);
    if strcmpi(ext, '.gz') || strcmpi(ext, '.zst')
        [~, base, ~] = fileparts(base);
    end
    switch compress_mode
        case 'gzip', out_ext = '.edf.gz';
        case 'zstd', out_ext = '.edf.zst';
        otherwise,   out_ext = '.edf';
    end
    if isempty(out_dir), out_dir = dir_; end
    final_out = fullfile(out_dir, sprintf('%s_%dHz%s', base, round(target_rate), out_ext));

    if ~isfolder(out_dir), mkdir(out_dir); end

    if stage_local
        % Copy input to local scratch, convert there, move output back.
        % Avoids per-worker NFS read contention.
        if isempty(stage_dir)
            scratch = tempname();          % uses system $TMPDIR / tempdir()
        else
            scratch = tempname(stage_dir); % unique subdir under user-supplied root
        end
        mkdir(scratch);
        cleanup = onCleanup(@() rmdir(scratch, 's')); %#ok<NASGU>

        [~, in_base, in_ext] = fileparts(in_fname);
        local_in = fullfile(scratch, [in_base in_ext]);
        copyfile(in_fname, local_in);

        local_out = fullfile(scratch, ...
            sprintf('%s_%dHz%s', base, round(target_rate), out_ext));

        if verbose
            fprintf('    staged to %s\n', scratch);
        end

        convert_EDF(local_in, target_rate, ...
            'OutputName',   local_out, ...
            'CompressMode', compress_mode, ...
            'GzipLevel',    gzip_level, ...
            'ZstdLevel',    zstd_level, ...
            'Verbose',      verbose, ...
            'AutoScale',    autoscale);

        movefile(local_out, final_out, 'f');
    else
        convert_EDF(in_fname, target_rate, ...
            'OutputName',   final_out, ...
            'CompressMode', compress_mode, ...
            'GzipLevel',    gzip_level, ...
            'ZstdLevel',    zstd_level, ...
            'Verbose',      verbose, ...
            'AutoScale',    autoscale);
    end

    out_path = final_out;
    status   = 'ok';
catch ME
    errmsg = ME.message;
end
elapsed = toc(t0);
end


%% =========================================================================
function files = resolve_inputs(input_spec, pattern)
files = {};

if iscell(input_spec)
    files = cellfun(@char, input_spec, 'UniformOutput', false);
    return
end

if ~ischar(input_spec) && ~isstring(input_spec)
    error('batch_convert_EDF:InputType', ...
        'input_spec must be a cell array, directory string, or list-file path.');
end
input_spec = char(input_spec);

if isfolder(input_spec)
    d1 = dir(fullfile(input_spec, pattern));
    d2 = dir(fullfile(input_spec, [pattern '.gz']));
    d  = [d1; d2];
    files = arrayfun(@(x) fullfile(x.folder, x.name), d, 'UniformOutput', false);
    return
end

if isfile(input_spec)
    [~, ~, ext] = fileparts(input_spec);
    if strcmpi(ext, '.txt') || strcmpi(ext, '.list')
        % list file
        fid = fopen(input_spec, 'r');
        if fid < 0
            error('batch_convert_EDF:ListOpen', 'Cannot open list file: %s', input_spec);
        end
        c = textscan(fid, '%s', 'Delimiter', '\n', 'CommentStyle', '#');
        fclose(fid);
        files = c{1};
        files = files(~cellfun(@isempty, strtrim(files)));
        return
    end
    % Single file path passed as a string
    files = {input_spec};
    return
end

error('batch_convert_EDF:NotFound', ...
    'input_spec not found: %s', input_spec);
end


%% =========================================================================
function t = empty_result_table()
t = table('Size', [0 5], ...
    'VariableTypes', {'string', 'string', 'string', 'double', 'string'}, ...
    'VariableNames', {'input', 'output', 'status', 'elapsed_s', 'error_message'});
end
