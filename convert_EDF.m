function out_path = convert_EDF(in_fname, target_rate, varargin)
%CONVERT_EDF  Read an EDF, resample all signals to a uniform rate, write it.
%
%   out_path = convert_EDF(in_fname, target_rate)
%   out_path = convert_EDF(in_fname, target_rate, 'OutputName', name, ...)
%
%   Reads an EDF / .edf.gz / .edf.zst file with read_EDF, resamples every
%   non-annotation channel to target_rate Hz using SPT's anti-aliased
%   resample, and writes a new EDF (zstd-compressed by default — change
%   with 'CompressMode').
%
%   Inputs:
%       in_fname    : path to a .edf, .edf.gz, or .edf.zst file
%       target_rate : target sample rate (Hz) for all signals
%
%   Name-value pairs:
%       'OutputName'   : explicit output path. Default:
%                        <dir>/<basename>_<target_rate>Hz<ext>
%                        where <ext> depends on 'CompressMode'.
%       'CompressMode' : 'zstd' (default) | 'gzip' | 'none'. Sets the
%                        default output extension when 'OutputName' is not
%                        provided (.edf.zst / .edf.gz / .edf).
%       'Compress'     : (deprecated) true -> 'zstd', false -> 'none'.
%                        'CompressMode' takes precedence if both are passed.
%       'GzipLevel'    : integer 1..9 (default 6). Used for .gz outputs.
%       'ZstdLevel'    : integer 1..22 (default 3). Used for .zst outputs.
%       'Verbose'      : logical (default false)
%       'AutoScale'    : 'recompute' (default) | 'preserve' (passed to
%                        write_EDF). 'recompute' resets each channel's
%                        physical_min/max to the actual data range so the
%                        anti-aliasing filter's ringing isn't clipped at
%                        the source file's narrower stored range.
%       'MaxResampleBlockMB' : 4096 (default). Cap on the size of any
%                        single resample() matrix call. Channels of the
%                        same source rate and length are grouped and
%                        resampled together (one matrix call instead of
%                        N per-channel calls); this cap chunks the group
%                        if the stacked matrix would exceed the budget.
%                        Lower it on memory-constrained machines.
%
%   Output:
%       out_path : full path of the written file
%
%   Validation:
%       target_rate * data_record_duration must be a positive integer for
%       every signal (so that samples_in_record is integer). If not,
%       errors with the closest valid alternative.
%
%   Behavior:
%       Annotation channels are passed through untouched (their byte
%       layout per record is preserved). Each non-annotation channel is
%       resampled with SPT's resample (Kaiser-windowed sinc FIR; this is
%       the anti-aliasing path -- naive downsample is NOT used).
%
%   Why AutoScale defaults to 'recompute':
%       The anti-aliasing filter can briefly produce samples outside the
%       source file's declared physical_min/max (filter ringing at sharp
%       transients). With 'preserve' those samples get clipped to the
%       int16 range on write. 'recompute' widens the stored range to the
%       actual resampled data range, eliminating the clipping.
%
%   See also: read_EDF, write_EDF, batch_convert_EDF.

p = inputParser;
addRequired(p,  'in_fname');
addRequired(p,  'target_rate', @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'OutputName',   '', @ischar);
addParameter(p, 'CompressMode', '', @(s) ischar(s) && (isempty(s) || any(strcmpi(s, {'gzip','zstd','none'}))));
addParameter(p, 'Compress',     true, @islogical);
addParameter(p, 'GzipLevel',    6, @(x) isnumeric(x) && isscalar(x) && x >= 1 && x <= 9);
addParameter(p, 'ZstdLevel',    3, @(x) isnumeric(x) && isscalar(x) && x >= 1 && x <= 22);
addParameter(p, 'Verbose',      false, @islogical);
addParameter(p, 'AutoScale',    'recompute', @(s) any(strcmpi(s, {'preserve','recompute'})));
addParameter(p, 'MaxResampleBlockMB', 4096, @(x) isnumeric(x) && isscalar(x) && x > 0);
parse(p, in_fname, target_rate, varargin{:});

in_fname     = char(p.Results.in_fname);
target_rate  = double(p.Results.target_rate);
out_name     = p.Results.OutputName;
compress_mode= lower(p.Results.CompressMode);
compress     = p.Results.Compress;
gzip_level   = double(p.Results.GzipLevel);
zstd_level   = double(p.Results.ZstdLevel);
verbose      = p.Results.Verbose;
autoscale    = p.Results.AutoScale;
max_block_bytes = double(p.Results.MaxResampleBlockMB) * 1024 * 1024;

% Resolve compression mode. CompressMode wins; otherwise legacy Compress flag.
if isempty(compress_mode)
    if compress, compress_mode = 'zstd'; else, compress_mode = 'none'; end
end

if ~isfile(in_fname)
    error('convert_EDF:FileNotFound', 'EDF file not found: %s', in_fname);
end

% --- Default output name --------------------------------------------------
if isempty(out_name)
    [dir_, base, ext] = fileparts(in_fname);
    if strcmpi(ext, '.gz') || strcmpi(ext, '.zst')
        [~, base, ~] = fileparts(base);  % strip the .edf
    end
    switch compress_mode
        case 'gzip', out_ext = '.edf.gz';
        case 'zstd', out_ext = '.edf.zst';
        otherwise,   out_ext = '.edf';
    end
    out_name = fullfile(dir_, sprintf('%s_%dHz%s', base, round(target_rate), out_ext));
end

if verbose
    fprintf('convert_EDF: %s -> %s @ %g Hz\n', in_fname, out_name, target_rate);
end

% --- Read --------------------------------------------------------------
[header, signal_header, signal_cell, annotations] = read_EDF(in_fname, ...
    'Verbose', verbose);

record_duration = header.data_record_duration;
n_signals = numel(signal_header);

% Validate target_rate * record_duration is a positive integer
new_spr = target_rate * record_duration;
if abs(new_spr - round(new_spr)) > 1e-9
    suggested = round(new_spr) / record_duration;
    error('convert_EDF:NonIntegerSPR', ...
        ['target_rate * data_record_duration = %g, which is not an integer ' ...
         'samples_in_record. Closest valid target_rate is %g Hz.'], ...
        new_spr, suggested);
end
new_spr = round(new_spr);

% --- Locate annotation channel ---------------------------------------------
annot_idx = 0;
for i = 1:n_signals
    if strcmp(strtrim(signal_header(i).signal_labels), 'EDF Annotations')
        annot_idx = i;
        break;
    end
end

% --- Resample each non-annotation channel ----------------------------------
% Filter design is cached by (P,Q) across calls so we don't re-run firls /
% kaiser for every channel and every file. The cache is `persistent`, so
% under parfor each worker builds up its own cache as it processes files.
persistent resample_filter_cache
if isempty(resample_filter_cache)
    resample_filter_cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
end

new_signal_cell = signal_cell;
new_signal_header = signal_header;

% --- Group channels by (orig_rate, length) ---------------------------------
% Channels of the same source rate AND same length can be stacked into a
% single matrix and resampled with one resample() call (column-wise),
% avoiding per-channel function-call overhead and giving MATLAB's IPP /
% BLAS path a larger contiguous workload to dispatch on.
group_map = containers.Map('KeyType','char','ValueType','any');
for i = 1:n_signals
    if i == annot_idx, continue; end
    orig_rate = signal_header(i).sampling_frequency;
    if isempty(orig_rate) || ~isfinite(orig_rate) || orig_rate <= 0
        orig_rate = signal_header(i).samples_in_record / record_duration;
    end
    if abs(orig_rate - target_rate) < 1e-9, continue; end  % already at target
    L = numel(signal_cell{i});
    gkey = sprintf('%.10g_%d', orig_rate, L);
    if isKey(group_map, gkey)
        g = group_map(gkey);
        g.idx(end+1) = i;
        group_map(gkey) = g;
    else
        g.orig_rate = orig_rate;
        g.L         = L;
        g.idx       = i;
        group_map(gkey) = g;
    end
end

% --- Resample each group, chunked by RAM budget ----------------------------
% RAM headroom: a resample() call on an L x C double matrix peaks at
% roughly 3-4x the input matrix bytes (input + output + internal scratch).
safety_factor = 4;

gkeys = group_map.keys;
for gk = 1:numel(gkeys)
    g = group_map(gkeys{gk});
    [P, Q] = rat(target_rate / g.orig_rate, 1e-6);
    cache_key = sprintf('%d_%d', P, Q);
    if isKey(resample_filter_cache, cache_key)
        b = resample_filter_cache(cache_key);
    else
        [~, b] = resample(zeros(20,1), P, Q);
        resample_filter_cache(cache_key) = b;
    end

    bytes_per_ch = g.L * 8;
    chunk_ch = max(1, floor(max_block_bytes / (bytes_per_ch * safety_factor)));
    if verbose
        fprintf('  group: %g Hz -> %g Hz (ratio %d/%d), %d ch x %d samp, chunks of %d ch\n', ...
                g.orig_rate, target_rate, P, Q, numel(g.idx), g.L, chunk_ch);
    end

    for cs = 1:chunk_ch:numel(g.idx)
        ce = min(cs + chunk_ch - 1, numel(g.idx));
        chunk = g.idx(cs:ce);

        % Stack channels as columns of an L x nchunk matrix
        X = zeros(g.L, numel(chunk));
        for c = 1:numel(chunk)
            X(:, c) = signal_cell{chunk(c)}(:);
        end

        Y = resample(X, P, Q, b);

        for c = 1:numel(chunk)
            i = chunk(c);
            y = Y(:, c)';
            n_records = floor(numel(y) / new_spr);
            if n_records < 1
                error('convert_EDF:TooShort', ...
                    'After resampling, signal %d (%s) has no complete records.', ...
                    i, signal_header(i).signal_labels);
            end
            if numel(y) > n_records * new_spr
                y = y(1 : n_records * new_spr);
            end
            new_signal_cell{i}                       = y;
            new_signal_header(i).samples_in_record   = new_spr;
            new_signal_header(i).sampling_frequency  = target_rate;
        end
        clear X Y
    end
end

% --- Update header.num_data_records ----------------------------------------
% Use minimum across non-annotation channels (write_EDF will also do this,
% but be explicit).
min_records = inf;
for i = 1:n_signals
    if i == annot_idx, continue; end
    r = floor(numel(new_signal_cell{i}) / new_signal_header(i).samples_in_record);
    if r < min_records, min_records = r; end
end
header.num_data_records = min_records;

% --- Write -----------------------------------------------------------------
write_EDF(out_name, header, new_signal_header, new_signal_cell, annotations, ...
          'Verbose', verbose, 'AutoScale', autoscale, ...
          'GzipLevel', gzip_level, 'ZstdLevel', zstd_level);

out_path = out_name;
end
