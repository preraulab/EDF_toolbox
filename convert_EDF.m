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
new_signal_cell = signal_cell;
new_signal_header = signal_header;

for i = 1:n_signals
    if i == annot_idx, continue; end

    orig_rate = signal_header(i).sampling_frequency;
    if isempty(orig_rate) || ~isfinite(orig_rate) || orig_rate <= 0
        orig_rate = signal_header(i).samples_in_record / record_duration;
    end

    if abs(orig_rate - target_rate) < 1e-9
        % already at target rate
        continue
    end

    [P, Q] = rat(target_rate / orig_rate, 1e-6);
    if verbose
        fprintf('  resample %s: %g Hz -> %g Hz (ratio %d/%d)\n', ...
                strtrim(signal_header(i).signal_labels), orig_rate, target_rate, P, Q);
    end

    x = signal_cell{i};
    y = resample(x(:), P, Q);
    y = y(:)';   % keep row-vector convention used by read_EDF outputs

    % Trim to integer multiple of new_spr
    n_records = floor(numel(y) / new_spr);
    if n_records < 1
        error('convert_EDF:TooShort', ...
            'After resampling, signal %d (%s) has no complete records.', ...
            i, signal_header(i).signal_labels);
    end
    if numel(y) > n_records * new_spr
        if verbose
            fprintf('    trimming %d trailing samples to fit record boundary\n', ...
                    numel(y) - n_records * new_spr);
        end
        y = y(1 : n_records * new_spr);
    end

    new_signal_cell{i}                       = y;
    new_signal_header(i).samples_in_record   = new_spr;
    new_signal_header(i).sampling_frequency  = target_rate;
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
