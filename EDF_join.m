function varargout = EDF_join(input_spec, varargin)
%EDF_JOIN  Join split EDF / EDF+ segments into one continuous recording.
%
%   out_path = EDF_join(filelist)
%   out_path = EDF_join(filelist, 'OutputName', name, ...)
%   [out_path, header, signal_header, signal_cell, annotations] = EDF_join(...)
%
%   Reads a set of EDF files that are pieces of a single recording (a
%   clinical acquisition split into multiple files, with or without gaps
%   between them), verifies that they belong together, orders them by
%   recording start date/time, and writes ONE continuous EDF spanning the
%   full time range. Any gap between the end of one segment and the start
%   of the next is filled with a constant "blank" value (0 in physical
%   units by default) so that every sample from the earliest start to the
%   latest end is present and correctly time-aligned.
%
%   Inputs:
%       input_spec  : one of
%                       - cell array of file paths (.edf/.edf.gz/.edf.zst)
%                       - directory string (globbed with 'Pattern')
%                       - text file path with one filename per line
%
%   Name-value pairs:
%       'OutputName'   : explicit output path. Default:
%                        <dir>/<first_basename>_joined<ext>, where <ext>
%                        follows 'CompressMode'.
%       'CompressMode' : 'zstd' (default) | 'gzip' | 'none'. Sets the
%                        default output extension (.edf.zst / .edf.gz /
%                        .edf) when 'OutputName' is not given.
%       'GzipLevel'    : integer 1..9  (default 6). Used for .gz output.
%       'ZstdLevel'    : integer 1..22 (default 9). Used for .zst output.
%       'AutoScale'    : 'preserve' (default) | 'recompute'. Passed to
%                        write_EDF. 'preserve' keeps each channel's stored
%                        physical_min/max (lossless when all segments share
%                        the same scaling); 'recompute' resets the range
%                        from the joined data. If the segments do NOT share
%                        identical per-channel physical/digital ranges, or
%                        the fill value falls outside a channel's stored
%                        range where a gap exists, 'recompute' is the safe
%                        choice (a warning points this out).
%       'FillValue'    : constant physical-unit value written into gaps
%                        (default 0). Applies to every non-annotation
%                        channel.
%       'StartTimes'   : override the per-file start datetime instead of
%                        parsing it from the headers. Cell array or
%                        datetime array with one entry per input file, in
%                        the SAME order as the resolved file list (see
%                        'Verbose' to print that order). Entries may be
%                        datetime, datenum, or a parseable date string.
%                        Use this when header dates are missing or wrong.
%       'Tolerance'    : alignment/overlap tolerance in seconds
%                        (default 1e-3). A segment must start an integer
%                        number of records after the global start to sit
%                        on the continuous record grid; deviations larger
%                        than this are an error. Segment time ranges that
%                        overlap by more than this are also an error.
%       'Pattern'      : '*.edf' (default) glob when input_spec is a
%                        directory (also matches .edf.gz / .edf.zst).
%       'Verbose'      : logical (default false). Prints the resolved,
%                        sorted segment list with start times, gaps, and
%                        the joined total duration.
%       'forceMATLAB'  : logical (default false). Forwarded to read_EDF /
%                        write_EDF to bypass the MEX backends.
%
%   Outputs:
%       out_path      : full path of the written joined file
%       header        : joined file-level header struct (as written)
%       signal_header : joined per-signal header struct array
%       signal_cell   : joined per-channel physical-unit vectors
%       annotations   : merged EDF+ annotations (onsets shifted to the
%                       joined timeline)
%
%   -------------------------------------------------------------------------
%   What is checked (all hard errors unless noted)
%       • Ordering       : segments are sorted by start datetime. Start
%                           datetime is taken from the EDF+ 'Startdate
%                           DD-MMM-YYYY' token in local_rec_id when present
%                           (4-digit year, unambiguous), otherwise from the
%                           main-header recording_startdate 'dd.mm.yy'
%                           (EDF 1985 century rule) plus recording_starttime.
%                           A missing/invalid start time is an error unless
%                           supplied via 'StartTimes'.
%       • Record grid    : data_record_duration must be identical across
%                           segments, and each segment must begin an integer
%                           number of records after the earliest start (so
%                           samples land on one shared record grid). A
%                           non-integer offset larger than 'Tolerance' is an
%                           error — the segments cannot be tiled exactly.
%       • No overlap     : segment [start, end) time ranges must not
%                           overlap (touching end-to-start is fine, that is
%                           a seamless join). Any overlap beyond 'Tolerance'
%                           is an error, since overlapping samples would be
%                           ambiguous.
%       • Identical      : the ordered list of non-annotation channel
%         montage          labels must match across all segments, and each
%                           channel's samples_in_record (hence sampling
%                           rate) must match. Differences are an error.
%                           Differing physical/digital ranges are a WARNING
%                           (use 'recompute' to avoid re-quantization loss).
%
%   Gaps between segments are filled with 'FillValue' (physical units) on
%   every non-annotation channel. EDF+ annotation channels are regenerated
%   by write_EDF for the full continuous span; annotations from every
%   segment are merged with their onsets shifted onto the joined timeline.
%
%   -------------------------------------------------------------------------
%   Examples:
%
%       % Join three explicitly-listed segments (zstd output next to file 1)
%       EDF_join({'night_part1.edf', 'night_part2.edf', 'night_part3.edf'});
%
%       % Join every EDF in a directory, plain .edf output at a set path
%       EDF_join('/data/split_night', ...
%           'CompressMode', 'none', 'OutputName', '/data/night_full.edf', ...
%           'Verbose', true);
%
%       % Headers have unreliable dates — supply start times explicitly
%       EDF_join({'a.edf','b.edf'}, ...
%           'StartTimes', {datetime(2024,3,1,22,0,0), datetime(2024,3,1,23,5,0)});
%
%       % Also capture the joined structures without re-reading
%       [out, hdr, shdr, sc, ann] = EDF_join(files);
%
%   Memory note: every segment is read fully into memory and the joined
%   recording is assembled in memory before writing, so peak use is roughly
%   the input total plus the output. For very large cohorts, join in groups.
%
%   See also: read_EDF, write_EDF, convert_EDF, batch_convert_EDF.

%% ---------------- INPUT PARSING ----------------
if nargin < 1
    error('EDF_join:InvalidInput', 'A list of EDF files is required.');
end

p = inputParser;
p.FunctionName    = 'EDF_join';
p.KeepUnmatched   = false;
p.PartialMatching = false;
p.CaseSensitive   = false;
addRequired(p,  'input_spec');
addParameter(p, 'OutputName',   '', @ischar);
addParameter(p, 'CompressMode', 'zstd', @(s) ischar(s) && any(strcmpi(s, {'gzip','zstd','none'})));
addParameter(p, 'GzipLevel',    6, @(x) isnumeric(x) && isscalar(x) && x >= 1 && x <= 9);
addParameter(p, 'ZstdLevel',    9, @(x) isnumeric(x) && isscalar(x) && x >= 1 && x <= 22);
addParameter(p, 'AutoScale',    'preserve', @(s) any(strcmpi(s, {'preserve','recompute'})));
addParameter(p, 'FillValue',    0, @(x) isnumeric(x) && isscalar(x) && isfinite(x));
addParameter(p, 'StartTimes',   [], @(x) isempty(x) || iscell(x) || isdatetime(x) || isnumeric(x));
addParameter(p, 'Tolerance',    1e-3, @(x) isnumeric(x) && isscalar(x) && x >= 0);
addParameter(p, 'Pattern',      '*.edf', @ischar);
addParameter(p, 'Verbose',      false, @islogical);
addParameter(p, 'forceMATLAB',  false, @islogical);
parse(p, input_spec, varargin{:});

out_name      = p.Results.OutputName;
compress_mode = lower(p.Results.CompressMode);
gzip_level    = double(p.Results.GzipLevel);
zstd_level    = double(p.Results.ZstdLevel);
autoscale     = lower(p.Results.AutoScale);
fill_value    = double(p.Results.FillValue);
start_override = p.Results.StartTimes;
tol_sec       = double(p.Results.Tolerance);
pattern       = p.Results.Pattern;
verbose       = p.Results.Verbose;
force_matlab  = p.Results.forceMATLAB;

%% ---------------- RESOLVE FILE LIST ----------------
files = resolve_input_list(input_spec, pattern);
if isempty(files)
    error('EDF_join:Empty', 'No input files resolved from input_spec.');
end
n_files = numel(files);

for k = 1:n_files
    if ~isfile(files{k})
        error('EDF_join:FileNotFound', 'EDF file not found: %s', files{k});
    end
end

% Normalize any StartTimes override to a datetime vector aligned to `files`.
start_override = normalize_start_override(start_override, n_files);

if verbose
    fprintf('EDF_join: %d input file(s)\n', n_files);
end

%% ---------------- READ EVERY SEGMENT ----------------
% Each segment is read in full (header + per-signal physical data +
% annotations). read_EDF handles .edf / .edf.gz / .edf.zst and MEX
% acceleration transparently.
seg(n_files) = struct('header', [], 'signal_header', [], 'signal_cell', [], ...
                      'annotations', [], 'start_dt', NaT, 'num_records', 0);
for k = 1:n_files
    if verbose
        fprintf('  reading [%d/%d] %s\n', k, n_files, files{k});
    end
    [hdr, shdr, sc, ann] = read_EDF(files{k}, 'Verbose', false, ...
                                    'forceMATLAB', force_matlab);
    seg(k).header        = hdr;
    seg(k).signal_header = shdr;
    seg(k).signal_cell   = sc;
    seg(k).annotations   = ann;
    seg(k).num_records   = hdr.num_data_records;

    if isempty(start_override)
        seg(k).start_dt = parse_edf_datetime(hdr, files{k});
    else
        seg(k).start_dt = start_override(k);
    end
end

%% ---------------- VALIDATE START TIMES / SORT ----------------
start_dts = [seg.start_dt];
if any(isnat(start_dts))
    bad = find(isnat(start_dts));
    error('EDF_join:NoStartTime', ...
        ['Could not determine a start date/time for %d file(s): %s. ' ...
         'Provide them explicitly via the ''StartTimes'' name-value pair.'], ...
        numel(bad), strjoin(files(bad), ', '));
end

[~, order] = sort(start_dts);
seg   = seg(order);
files = files(order);
start_dts = start_dts(order);
T0 = start_dts(1);

%% ---------------- VALIDATE RECORD DURATION ----------------
rec_durs = arrayfun(@(s) s.header.data_record_duration, seg);
D = rec_durs(1);
if ~isfinite(D) || D <= 0
    error('EDF_join:BadRecordDuration', ...
        'First segment has a non-positive data_record_duration (%g).', D);
end
if any(abs(rec_durs - D) > 1e-9)
    error('EDF_join:RecordDurationMismatch', ...
        ['data_record_duration differs across segments (found: %s s). ' ...
         'All segments must share one record duration to tile onto a common grid.'], ...
        strjoin(compose('%g', rec_durs), ', '));
end

%% ---------------- MONTAGE / SAMPLING-RATE CHECK ----------------
% Compare the ordered non-annotation channel layout of every segment to
% the first (template) segment. Annotation channels are handled separately.
tmpl_shdr      = seg(1).signal_header;
tmpl_annot_idx = find_annotation_indices(tmpl_shdr);
tmpl_data_idx  = setdiff(1:numel(tmpl_shdr), tmpl_annot_idx, 'stable');
tmpl_labels    = arrayfun(@(s) strtrim(s.signal_labels), tmpl_shdr(tmpl_data_idx), ...
                          'UniformOutput', false);
tmpl_spr       = arrayfun(@(s) s.samples_in_record, tmpl_shdr(tmpl_data_idx));

if isempty(tmpl_data_idx)
    error('EDF_join:NoChannels', 'The first segment has no non-annotation signals.');
end

% Per-segment map: data channel c (template order) -> index into that
% segment's signal arrays. Built from the montage check below.
data_map = zeros(n_files, numel(tmpl_data_idx));
ranges_match = true;

for k = 1:n_files
    shdr      = seg(k).signal_header;
    annot_idx = find_annotation_indices(shdr);
    data_idx  = setdiff(1:numel(shdr), annot_idx, 'stable');
    labels    = arrayfun(@(s) strtrim(s.signal_labels), shdr(data_idx), ...
                         'UniformOutput', false);
    spr       = arrayfun(@(s) s.samples_in_record, shdr(data_idx));

    if numel(labels) ~= numel(tmpl_labels) || ...
            ~all(strcmpi(labels, tmpl_labels))
        error('EDF_join:MontageMismatch', ...
            ['Segment %s has a different channel montage than %s.\n' ...
             '  expected: %s\n  found:    %s'], ...
            files{k}, files{1}, strjoin(tmpl_labels, ', '), strjoin(labels, ', '));
    end
    if ~isequal(spr(:), tmpl_spr(:))
        error('EDF_join:RateMismatch', ...
            ['Segment %s has different samples_in_record (sampling rate) ' ...
             'than %s for one or more channels.'], files{k}, files{1});
    end

    data_map(k, :) = data_idx;

    % Physical/digital range agreement is not required, but a mismatch
    % under AutoScale=preserve means each segment's samples get re-quantized
    % into the template's stored range on write, which can clip or lose
    % precision. Flag it so the caller can switch to 'recompute'.
    for c = 1:numel(tmpl_data_idx)
        a = tmpl_shdr(tmpl_data_idx(c));
        b = shdr(data_idx(c));
        if a.physical_min ~= b.physical_min || a.physical_max ~= b.physical_max || ...
                a.digital_min ~= b.digital_min || a.digital_max ~= b.digital_max
            ranges_match = false;
        end
    end
end

if ~ranges_match && strcmp(autoscale, 'preserve')
    warning('EDF_join:RangeMismatch', ...
        ['Segments have differing per-channel physical/digital ranges. With ' ...
         'AutoScale=''preserve'' their samples will be re-quantized into the ' ...
         'first segment''s range on write (possible clipping/precision loss). ' ...
         'Pass ''AutoScale'',''recompute'' to size the range from the joined data.']);
end

%% ---------------- PLACE SEGMENTS ON THE RECORD GRID ----------------
% Offset of each segment, measured in whole records from the global start.
offset_sec  = seconds(start_dts - T0);
offset_recs = offset_sec / D;
rec_start   = round(offset_recs);

for k = 1:n_files
    resid_sec = abs(offset_sec(k) - rec_start(k) * D);
    if resid_sec > tol_sec
        error('EDF_join:GridMisalignment', ...
            ['Segment %s starts %.6g s after the earliest segment, which is ' ...
             'not an integer number of %g s records (residual %.6g s > tolerance ' ...
             '%.6g s). The segments do not lie on a shared record grid.'], ...
            files{k}, offset_sec(k), D, resid_sec, tol_sec);
    end
end

seg_records = arrayfun(@(s) s.num_records, seg);
rec_end     = rec_start(:) + seg_records(:);   % exclusive end, in records

% Overlap check on the sorted segments (start times are non-decreasing).
tol_recs = tol_sec / D;
for k = 2:n_files
    prev_end = max(rec_end(1:k-1));
    if rec_start(k) < prev_end - tol_recs
        overlap_sec = (prev_end - rec_start(k)) * D;
        error('EDF_join:Overlap', ...
            ['Segment %s overlaps an earlier segment by %.6g s. Overlapping ' ...
             'samples are ambiguous; EDF_join will not merge them.'], ...
            files{k}, overlap_sec);
    end
end

total_records = max(rec_end);

%% ---------------- ASSEMBLE CONTINUOUS CHANNELS ----------------
% Build one physical-unit vector per non-annotation channel: allocate the
% full continuous length, prefill with FillValue, then drop each segment's
% samples at its record-aligned offset. Untouched regions remain the fill.
n_data = numel(tmpl_data_idx);
joined_cell = cell(1, n_data);
gap_present = false;

for c = 1:n_data
    spr_c  = tmpl_spr(c);
    n_tot  = total_records * spr_c;
    outvec = fill_value * ones(1, n_tot);

    for k = 1:n_files
        sidx = data_map(k, c);
        x    = seg(k).signal_cell{sidx};
        off  = rec_start(k) * spr_c;
        L    = min(numel(x), seg_records(k) * spr_c);
        outvec(off + 1 : off + L) = x(1:L);
    end

    joined_cell{c} = outvec;
end

% Was any record left as fill (a real gap, not just the joins)?
covered = false(1, total_records);
for k = 1:n_files
    covered(rec_start(k) + 1 : rec_end(k)) = true;
end
gap_present = any(~covered);

% Under 'preserve', warn if the fill value can't be represented in a
% channel's stored physical range wherever a gap exists (it would clip).
if gap_present && strcmp(autoscale, 'preserve')
    for c = 1:n_data
        a = tmpl_shdr(tmpl_data_idx(c));
        if fill_value < min(a.physical_min, a.physical_max) || ...
                fill_value > max(a.physical_min, a.physical_max)
            warning('EDF_join:FillOutOfRange', ...
                ['FillValue %g is outside the stored physical range [%g, %g] of ' ...
                 'channel ''%s''; gap samples will clip to the nearest bound under ' ...
                 'AutoScale=''preserve''. Use ''recompute'' or a different ''FillValue''.'], ...
                fill_value, a.physical_min, a.physical_max, tmpl_labels{c});
        end
    end
end

%% ---------------- MERGE ANNOTATIONS ----------------
merged_annotations = merge_annotations(seg, rec_start, D);

%% ---------------- BUILD OUTPUT HEADER / SIGNAL HEADERS ----------------
% Template = first (earliest) segment. Keep its full channel layout,
% including any annotation channel, so channel positions are preserved.
out_header = seg(1).header;
out_header = update_start_fields(out_header, T0);
out_header.data_record_duration = D;
out_header.num_data_records     = total_records;

out_shdr = tmpl_shdr;
out_cell = cell(1, numel(out_shdr));

% Non-annotation channels get their assembled continuous vectors.
for c = 1:n_data
    out_cell{tmpl_data_idx(c)} = joined_cell{c};
end

% Annotation channel(s): write_EDF regenerates the TAL content per record
% from the record timestamps and the merged annotations, so the payload
% here is an unused placeholder of the right length.
for j = 1:numel(tmpl_annot_idx)
    ai = tmpl_annot_idx(j);
    spr_a = out_shdr(ai).samples_in_record;
    if isempty(spr_a) || ~isfinite(spr_a) || spr_a <= 0
        spr_a = 60;
        out_shdr(ai).samples_in_record = spr_a;
    end
    out_cell{ai} = zeros(1, spr_a * total_records);
end

% No annotation channel in the template, but annotations exist -> add a
% standard EDF+ annotation channel so they are not silently dropped.
if isempty(tmpl_annot_idx) && ~isempty(merged_annotations)
    [out_shdr, out_cell] = append_annotation_channel(out_shdr, out_cell, total_records);
end

out_header.num_signals      = numel(out_shdr);
out_header.num_header_bytes = 256 * (1 + numel(out_shdr));

%% ---------------- RESOLVE OUTPUT NAME ----------------
if isempty(out_name)
    [dir_, base, ext] = fileparts(files{1});
    if strcmpi(ext, '.gz') || strcmpi(ext, '.zst')
        [~, base, ~] = fileparts(base);   % strip the .edf
    end
    switch compress_mode
        case 'gzip', out_ext = '.edf.gz';
        case 'zstd', out_ext = '.edf.zst';
        otherwise,   out_ext = '.edf';
    end
    out_name = fullfile(dir_, sprintf('%s_joined%s', base, out_ext));
end

%% ---------------- VERBOSE SUMMARY ----------------
if verbose
    fprintf('EDF_join: %d segments on a %g s record grid\n', n_files, D);
    for k = 1:n_files
        gap_before = '';
        if k > 1
            g = (rec_start(k) - max(rec_end(1:k-1))) * D;
            if g > 0, gap_before = sprintf('   [gap before: %g s]', g); end
        end
        fprintf('  %s  start=%s  records=%d  (offset %g s)%s\n', ...
            files{k}, char(start_dts(k)), seg_records(k), ...
            rec_start(k) * D, gap_before);
    end
    total_sec = total_records * D;
    fprintf('  joined: %d records, %g s (%s), %d signals -> %s\n', ...
        total_records, total_sec, char(duration(0, 0, total_sec)), ...
        numel(out_shdr), out_name);
    if ~gap_present
        fprintf('  (segments are contiguous; no gap fill needed)\n');
    end
end

%% ---------------- WRITE JOINED FILE ----------------
write_EDF(out_name, out_header, out_shdr, out_cell, merged_annotations, ...
          'Verbose', verbose, 'AutoScale', autoscale, ...
          'GzipLevel', gzip_level, 'ZstdLevel', zstd_level, ...
          'forceMATLAB', force_matlab);

%% ---------------- OUTPUT ----------------
varargout{1} = out_name;
if nargout > 1, varargout{2} = out_header; end
if nargout > 2, varargout{3} = out_shdr; end
if nargout > 3, varargout{4} = out_cell; end
if nargout > 4, varargout{5} = merged_annotations; end

end


%% =========================================================================
%  INPUT LIST RESOLUTION
%  -----------------------------------------------------------------------
%  Accept a cell array of paths, a directory (globbed with Pattern), or a
%  text/list file with one path per line. Mirrors batch_convert_EDF's
%  resolver so the two share input conventions.
% =========================================================================
function files = resolve_input_list(input_spec, pattern)
files = {};

if iscell(input_spec)
    files = cellfun(@char, input_spec, 'UniformOutput', false);
    files = files(~cellfun(@isempty, strtrim(files)));
    return
end

if isstring(input_spec) && ~isscalar(input_spec)
    files = cellstr(input_spec(:));
    files = files(~cellfun(@isempty, strtrim(files)));
    return
end

if ~ischar(input_spec) && ~isstring(input_spec)
    error('EDF_join:InputType', ...
        'input_spec must be a cell array, directory string, or list-file path.');
end
input_spec = char(input_spec);

if isfolder(input_spec)
    d1 = dir(fullfile(input_spec, pattern));
    d2 = dir(fullfile(input_spec, [pattern '.gz']));
    d3 = dir(fullfile(input_spec, [pattern '.zst']));
    d  = [d1; d2; d3];
    files = arrayfun(@(x) fullfile(x.folder, x.name), d, 'UniformOutput', false);
    return
end

if isfile(input_spec)
    [~, ~, ext] = fileparts(input_spec);
    if strcmpi(ext, '.txt') || strcmpi(ext, '.list')
        fid = fopen(input_spec, 'r');
        if fid < 0
            error('EDF_join:ListOpen', 'Cannot open list file: %s', input_spec);
        end
        c = textscan(fid, '%s', 'Delimiter', '\n', 'CommentStyle', '#');
        fclose(fid);
        files = c{1};
        files = files(~cellfun(@isempty, strtrim(files)));
        return
    end
    files = {input_spec};   % a single EDF path
    return
end

error('EDF_join:NotFound', 'input_spec not found: %s', input_spec);
end


%% =========================================================================
%  START-TIME OVERRIDE NORMALIZATION
% =========================================================================
function dt = normalize_start_override(start_override, n_files)
if isempty(start_override)
    dt = datetime.empty;
    return
end

if isdatetime(start_override)
    dt = start_override(:).';
elseif iscell(start_override)
    dt = NaT(1, numel(start_override));
    for k = 1:numel(start_override)
        v = start_override{k};
        if isdatetime(v)
            dt(k) = v;
        elseif isnumeric(v) && isscalar(v)
            dt(k) = datetime(v, 'ConvertFrom', 'datenum');
        elseif ischar(v) || isstring(v)
            dt(k) = datetime(char(v));
        else
            error('EDF_join:BadStartTimes', ...
                'StartTimes entry %d is not a datetime, datenum, or date string.', k);
        end
    end
elseif isnumeric(start_override)
    dt = datetime(start_override(:).', 'ConvertFrom', 'datenum');
else
    error('EDF_join:BadStartTimes', 'Unrecognized StartTimes type.');
end

if numel(dt) ~= n_files
    error('EDF_join:StartTimesCount', ...
        'StartTimes has %d entries but %d files were resolved.', numel(dt), n_files);
end
end


%% =========================================================================
%  EDF START DATE/TIME PARSING
%  -----------------------------------------------------------------------
%  Prefer the EDF+ 'Startdate DD-MMM-YYYY' token in local_rec_id (4-digit
%  year, unambiguous). Fall back to the main-header recording_startdate
%  'dd.mm.yy' with the EDF century rule (yy 85..99 -> 1985..1999,
%  00..84 -> 2000..2084). Time comes from recording_starttime 'hh.mm.ss'.
%  Returns NaT when no usable date can be recovered.
% =========================================================================
function dt = parse_edf_datetime(header, fname)
dt = NaT;

[H, M, S, time_ok] = parse_hms(getfield_default(header, 'recording_starttime', ''));

ymd = parse_startdate_plus(getfield_default(header, 'local_rec_id', ''));
if isempty(ymd)
    ymd = parse_startdate_basic(getfield_default(header, 'recording_startdate', ''));
end
if isempty(ymd)
    return
end

if ~time_ok
    % Date recovered but time unreadable: assume midnight and warn, since a
    % wrong time-of-day would misplace the segment on the grid.
    H = 0; M = 0; S = 0;
    warning('EDF_join:NoStartTimeField', ...
        'Unreadable recording_starttime in %s; assuming 00.00.00.', fname);
end

dt = datetime(ymd(1), ymd(2), ymd(3), H, M, S);
end

function [H, M, S, ok] = parse_hms(str)
H = 0; M = 0; S = 0; ok = false;
str = strtrim(char(str));
if isempty(str), return; end
parts = regexp(str, '[.:]', 'split');
if numel(parts) < 3, return; end
h = str2double(parts{1}); m = str2double(parts{2}); s = str2double(parts{3});
if any(isnan([h m s])), return; end
if h < 0 || h > 23 || m < 0 || m > 59 || s < 0 || s >= 61, return; end
H = h; M = m; S = s; ok = true;
end

function ymd = parse_startdate_basic(str)
ymd = [];
str = strtrim(char(str));
if isempty(str), return; end
parts = regexp(str, '[.:/-]', 'split');
if numel(parts) < 3, return; end
d = str2double(parts{1}); mo = str2double(parts{2}); yy = str2double(parts{3});
if any(isnan([d mo yy])), return; end
if yy < 100
    if yy >= 85, yy = 1900 + yy; else, yy = 2000 + yy; end   % EDF century rule
end
if mo < 1 || mo > 12 || d < 1 || d > 31, return; end
ymd = [yy mo d];
end

function ymd = parse_startdate_plus(local_rec_id)
ymd = [];
str = char(local_rec_id);
tok = regexpi(str, 'Startdate\s+(\d{1,2})-([A-Za-z]{3})-(\d{2,4})', 'tokens', 'once');
if isempty(tok), return; end
d  = str2double(tok{1});
mo = month_from_abbrev(tok{2});
yy = str2double(tok{3});
if isnan(d) || isnan(mo) || isnan(yy), return; end
if yy < 100
    if yy >= 85, yy = 1900 + yy; else, yy = 2000 + yy; end
end
if d < 1 || d > 31, return; end
ymd = [yy mo d];
end

function mo = month_from_abbrev(abbr)
names = {'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'};
mo = find(strcmpi(abbr, names), 1);
if isempty(mo), mo = NaN; end
end


%% =========================================================================
%  HEADER START-FIELD UPDATE
%  -----------------------------------------------------------------------
%  Stamp the joined recording's start into the main header: the plain-EDF
%  date/time fields, and the EDF+ 'Startdate DD-MMM-YYYY' token inside
%  local_rec_id when that field carries one.
% =========================================================================
function header = update_start_fields(header, T0)
header.recording_startdate = char(datetime(T0, 'Format', 'dd.MM.yy'));
header.recording_starttime = char(datetime(T0, 'Format', 'HH.mm.ss'));

if isfield(header, 'local_rec_id') && ~isempty(header.local_rec_id)
    upper_mon = upper(char(datetime(T0, 'Format', 'dd-MMM-yyyy')));
    new_id = regexprep(char(header.local_rec_id), ...
        'Startdate\s+\d{1,2}-[A-Za-z]{3}-\d{2,4}', ...
        ['Startdate ' upper_mon], 'ignorecase');
    header.local_rec_id = new_id;
end
end


%% =========================================================================
%  ANNOTATION HELPERS
% =========================================================================
function idx = find_annotation_indices(signal_header)
%FIND_ANNOTATION_INDICES  Indices of 'EDF Annotations' channels.
labels = arrayfun(@(s) strtrim(s.signal_labels), signal_header, ...
                  'UniformOutput', false);
idx = find(strcmp(labels, 'EDF Annotations'));
end

function merged = merge_annotations(seg, rec_start, D)
%MERGE_ANNOTATIONS  Combine per-segment annotations onto the joined timeline.
%   Each segment's onsets are relative to its own start; shift them by the
%   segment's record offset (rec_start * D seconds) and concatenate.
merged = struct('onset', {}, 'text', {});
for k = 1:numel(seg)
    ann = seg(k).annotations;
    if isempty(ann), continue; end
    shift = rec_start(k) * D;
    for a = 1:numel(ann)
        j = numel(merged) + 1;
        merged(j).onset = ann(a).onset + shift;
        merged(j).text  = ann(a).text;
    end
end
if ~isempty(merged)
    [~, ord] = sort([merged.onset]);
    merged = merged(ord);
end
end

function [shdr, scell] = append_annotation_channel(shdr, scell, total_records)
%APPEND_ANNOTATION_CHANNEL  Add a standard EDF+ annotation channel.
%   Used only when the template segment lacks one but annotations were
%   collected from other segments. write_EDF fills the TAL payload.
spr_a = 60;   % 120 bytes/record; write_EDF drops overflow with a warning
a = struct( ...
    'signal_labels',      'EDF Annotations', ...
    'transducer_type',    '', ...
    'physical_dimension', '', ...
    'physical_min',       -1, ...
    'physical_max',        1, ...
    'digital_min',        -32768, ...
    'digital_max',         32767, ...
    'prefiltering',       '', ...
    'samples_in_record',  spr_a, ...
    'reserve_2',          '');

% Match whatever field set the existing signal_header carries so the
% struct array stays homogeneous.
fn = fieldnames(shdr);
new_entry = shdr(1);
for i = 1:numel(fn)
    if isfield(a, fn{i})
        new_entry.(fn{i}) = a.(fn{i});
    else
        new_entry.(fn{i}) = [];
    end
end
if isfield(new_entry, 'sampling_frequency')
    new_entry.sampling_frequency = [];
end

shdr(end+1)  = new_entry;
scell{end+1} = zeros(1, spr_a * total_records);
end


%% =========================================================================
%  MISC
% =========================================================================
function v = getfield_default(s, fn, dflt)
if isfield(s, fn) && ~isempty(s.(fn))
    v = s.(fn);
else
    v = dflt;
end
end
