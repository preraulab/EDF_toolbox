function write_EDF(out_fname, header, signal_header, signal_cell, varargin)
%WRITE_EDF  Write an EDF / EDF+ file (with optional gzip output)
%
%   WRITE_EDF is the structural mirror of READ_EDF. It accepts the same
%   header / signal_header / signal_cell / annotations shapes that
%   READ_EDF returns, and produces a valid EDF or EDF+ file. Output is
%   written via a compiled MEX (with vendored zlib for .edf.gz) when
%   available, or a pure MATLAB fallback otherwise.
%
%   Usage:
%       write_EDF(out_fname, header, signal_header, signal_cell)
%       write_EDF(out_fname, header, signal_header, signal_cell, annotations)
%       write_EDF(..., 'AutoScale', 'preserve' | 'recompute', ...)
%
%   Inputs:
%       out_fname     : output file path. Ends in '.gz' -> gzipped output.
%       header        : 1x1 struct (matches read_EDF's first output)
%       signal_header : 1xN struct array (matches read_EDF's second output)
%       signal_cell   : 1xN cell array of physical-unit channel vectors
%       annotations   : Nx1 struct array of {onset, text} (optional, [] OK)
%
%   Name-value pairs:
%       'Verbose'     : logical (default false)
%       'forceMATLAB' : logical (default false)
%       'AutoScale'   : 'preserve' (default) -> keep existing physical_min/max
%                       'recompute'         -> set from data (lossless if data
%                                              fits within int16 dynamic range)
%       'debug'       : logical (default false)
%
%   Output naming:
%       If out_fname ends in '.gz' (case-insensitive), the writer streams
%       through zlib and produces a gzipped EDF on the fly (no temp file).
%
%   AutoScale notes:
%       'preserve' is the default because it is required for lossless
%       round-trip: read_EDF -> write_EDF reproduces the original file
%       up to ~1 digital LSB per sample (limited by float rounding in
%       physical-units storage).
%       'recompute' may shift the physical range if the input data
%       exceeds the original physical_min/max -- prevents clipping at
%       the cost of changing the file's stored scaling.
%
%   See also: read_EDF, convert_EDF, batch_convert_EDF.

%% ---------------- INPUT PARSING ----------------
if nargin < 4
    error('write_EDF:InvalidInput', 'write_EDF requires at least 4 args.');
end

% Annotations: optional positional or first name-value pair
annotations = [];
if ~isempty(varargin)
    if isstruct(varargin{1}) || isempty(varargin{1})
        annotations = varargin{1};
        varargin = varargin(2:end);
    end
end

p = inputParser;
addParameter(p, 'Verbose',     false, @islogical);
addParameter(p, 'forceMATLAB', false, @islogical);
addParameter(p, 'AutoScale',   'preserve', @(s) any(strcmpi(s, {'preserve','recompute'})));
addParameter(p, 'debug',       false, @islogical);
parse(p, varargin{:});

verbose      = p.Results.Verbose;
force_matlab = p.Results.forceMATLAB;
autoscale    = lower(p.Results.AutoScale);
debug        = p.Results.debug;

out_fname = char(out_fname);

if ~isstruct(header)
    error('write_EDF:Input', 'header must be a struct.');
end
if ~isstruct(signal_header)
    error('write_EDF:Input', 'signal_header must be a struct array.');
end
if ~iscell(signal_cell)
    error('write_EDF:Input', 'signal_cell must be a cell array.');
end
if numel(signal_cell) ~= numel(signal_header)
    error('write_EDF:Mismatch', ...
        'signal_cell has %d elements but signal_header has %d.', ...
        numel(signal_cell), numel(signal_header));
end

% Strip total_data_seconds / total_data_hms (read_EDF adds these for
% convenience; the MEX ignores unknown fields, but be defensive).
fields_to_drop = {'total_data_seconds', 'total_data_hms'};
for k = 1:numel(fields_to_drop)
    if isfield(header, fields_to_drop{k})
        header = rmfield(header, fields_to_drop{k});
    end
end

%% ---------------- MEX HANDLING ----------------
script_dir = fileparts(mfilename('fullpath'));
mex_file   = fullfile(script_dir, ['write_EDF_mex.' mexext]);
mex_exists = isfile(mex_file);

if ~force_matlab
    if ~mex_exists
        compile_edf_mex(script_dir, 'write_EDF_mex.c');
    end
    try
        autoscale_mode = double(strcmp(autoscale, 'recompute'));
        write_EDF_mex(out_fname, header, signal_header, signal_cell, ...
                      annotations, autoscale_mode, double(verbose), double(debug));
        return
    catch ME
        if verbose
            fprintf('write_EDF MEX failed (%s). Falling back to MATLAB writer.\n', ME.message);
        end
    end
else
    if verbose
        fprintf('forceMATLAB = true, using MATLAB writer.\n');
    end
end

%% ---------------- MATLAB FALLBACK ----------------
write_EDF_matlab(out_fname, header, signal_header, signal_cell, annotations, ...
                 autoscale, verbose);

end


%% =========================================================================
%  PURE MATLAB EDF WRITER
% =========================================================================
function write_EDF_matlab(out_fname, header, signal_header, signal_cell, ...
                          annotations, autoscale, verbose)

% gz output: write plain to a temp file, then gzip and delete temp.
gz_cleanup = []; %#ok<NASGU>
final_out = out_fname;
if endsWith(out_fname, '.gz', 'IgnoreCase', true)
    tmpdir = tempname;
    mkdir(tmpdir);
    gz_cleanup = onCleanup(@() rmdir(tmpdir, 's'));
    [~, base, ~] = fileparts(out_fname);
    [~, base2, ~] = fileparts(base);  % strip .edf
    if isempty(base2), base2 = base; end
    out_fname = fullfile(tmpdir, [base2 '.edf']);
end

n_signals = numel(signal_header);

% Locate annotation channel
annot_idx = 0;
for i = 1:n_signals
    if strcmp(strtrim(signal_header(i).signal_labels), 'EDF Annotations')
        annot_idx = i;
        break;
    end
end

% Determine num_data_records (truncate to common minimum across non-annot
% channels).
num_records = inf;
for i = 1:n_signals
    if i == annot_idx, continue; end
    spr = signal_header(i).samples_in_record;
    if spr <= 0
        error('write_EDF:Input', ...
            'signal %d (%s) has samples_in_record <= 0', i, signal_header(i).signal_labels);
    end
    r = floor(numel(signal_cell{i}) / spr);
    if r < num_records, num_records = r; end
end
if ~isfinite(num_records) || num_records == 0
    error('write_EDF:NoData', 'No complete data records to write.');
end

% Annotation channel default if not configured
if annot_idx > 0 && signal_header(annot_idx).samples_in_record <= 0
    signal_header(annot_idx).samples_in_record = 30;  % 60 bytes/record
end

% Optional autoscale
if strcmp(autoscale, 'recompute')
    for i = 1:n_signals
        if i == annot_idx, continue; end
        x = signal_cell{i}(1 : signal_header(i).samples_in_record * num_records);
        mn = min(x); mx = max(x);
        if mn == mx, mn = mn - 1; mx = mx + 1; end
        signal_header(i).physical_min = mn;
        signal_header(i).physical_max = mx;
    end
end

% Build & write header
data_record_duration = header.data_record_duration;
num_header_bytes = 256 * (1 + n_signals);

fid = fopen(out_fname, 'wb', 'ieee-le');
if fid < 0
    error('write_EDF:File', 'Cannot open %s for writing.', out_fname);
end
clean_fid = onCleanup(@() fclose_safe(fid)); %#ok<NASGU>

% Main 256-byte header
mh = repmat(' ', 1, 256);
mh = pad_field_m(mh,   1,   8, getf(header, 'edf_ver',              '0'));
mh = pad_field_m(mh,   9,  80, getf(header, 'patient_id',            'X X X X'));
mh = pad_field_m(mh,  89,  80, getf(header, 'local_rec_id',          'Startdate X X X X'));
mh = pad_field_m(mh, 169,   8, getf(header, 'recording_startdate',   '01.01.01'));
mh = pad_field_m(mh, 177,   8, getf(header, 'recording_starttime',   '00.00.00'));
mh = pad_field_m(mh, 185,   8, fmt_int(num_header_bytes));
mh = pad_field_m(mh, 237,   8, fmt_int(num_records));
mh = pad_field_m(mh, 245,   8, fmt_dbl(data_record_duration));
mh = pad_field_m(mh, 253,   4, fmt_int(n_signals));
fwrite(fid, mh, 'char');

% Per-signal header (transposed)
write_field = @(width, fn) deal_field(fid, signal_header, n_signals, width, fn);
write_field(16, 'signal_labels');
write_field(80, 'transducer_type');
write_field( 8, 'physical_dimension');
write_doubles(fid, signal_header, n_signals, 8, 'physical_min');
write_doubles(fid, signal_header, n_signals, 8, 'physical_max');
write_ints   (fid, signal_header, n_signals, 8, 'digital_min');
write_ints   (fid, signal_header, n_signals, 8, 'digital_max');
write_field(80, 'prefiltering');
write_ints   (fid, signal_header, n_signals, 8, 'samples_in_record');
fwrite(fid, repmat(' ', 1, 32 * n_signals), 'char');

% Build annotation per-record buffers
annot_bytes_per_rec = 0;
arecs = {};
if annot_idx > 0
    annot_bytes_per_rec = signal_header(annot_idx).samples_in_record * 2;
    arecs = cell(num_records, 1);
    for r = 1:num_records
        buf = zeros(1, annot_bytes_per_rec, 'uint8');
        ts = sprintf('%+g', (r-1) * data_record_duration);
        b = uint8([ts 20 20 0]);
        n = numel(b);
        if n <= annot_bytes_per_rec
            buf(1:n) = b;
        end
        arecs{r}.buf  = buf;
        arecs{r}.used = n;
    end
    if isstruct(annotations) && ~isempty(annotations)
        dropped = 0;
        for a = 1:numel(annotations)
            onset = annotations(a).onset;
            target = floor(onset / data_record_duration) + 1;
            target = max(1, min(num_records, target));
            txt = annotations(a).text;
            if ~iscell(txt), txt = {txt}; end
            for k = 1:numel(txt)
                placed = false;
                for r = target:num_records
                    [arecs{r}, ok] = append_tal_m(arecs{r}, onset, txt{k});
                    if ok, placed = true; break; end
                end
                if ~placed, dropped = dropped + 1; break; end
            end
        end
        if dropped > 0 && verbose
            fprintf('write_EDF: %d annotation(s) dropped (no record had room).\n', dropped);
        end
    end
end

% Per-signal byte offsets within a record
samples_per_record = arrayfun(@(s) s.samples_in_record, signal_header);
sample_offsets = [0 cumsum(samples_per_record(1:end-1))];

% Stream records
for r = 1:num_records
    rec_bytes = zeros(1, sum(samples_per_record) * 2, 'uint8');
    for i = 1:n_signals
        spr = samples_per_record(i);
        byte_off = sample_offsets(i) * 2;
        if i == annot_idx && annot_idx > 0
            n = arecs{r}.used;
            n = min(n, annot_bytes_per_rec);
            rec_bytes(byte_off + 1 : byte_off + n) = arecs{r}.buf(1:n);
            continue;
        end
        x = signal_cell{i}((r-1)*spr + 1 : r*spr);
        scale = (signal_header(i).physical_max - signal_header(i).physical_min) / ...
                (signal_header(i).digital_max  - signal_header(i).digital_min);
        if scale == 0, scale = 1; end
        dig = round((x - signal_header(i).physical_min) / scale + signal_header(i).digital_min);
        dig = max(min(dig, signal_header(i).digital_max), signal_header(i).digital_min);
        dig = max(min(dig, 32767), -32768);
        v = typecast(int16(dig), 'uint8');
        rec_bytes(byte_off + 1 : byte_off + spr*2) = v;
    end
    fwrite(fid, rec_bytes, 'uint8');
end

clear clean_fid;  % triggers fclose

% Compress if requested
if ~strcmp(out_fname, final_out)
    if verbose, fprintf('Gzipping output...\n'); end
    gzip(out_fname, fileparts(final_out));
    % gzip writes <out_fname>.gz; rename if needed
    written = [out_fname '.gz'];
    if ~strcmp(written, final_out)
        movefile(written, final_out, 'f');
    end
end

if verbose
    fprintf('write_EDF (MATLAB): %d signals, %d records -> %s\n', ...
            n_signals, num_records, final_out);
end

end


%% =========================================================================
%  HELPERS
% =========================================================================

function s = pad_field_m(s, start_pos, width, val)
val = char(val);
n = numel(val);
if n > width, n = width; end
s(start_pos : start_pos + n - 1) = val(1:n);
% trailing space-pad already present (mh was init'd to spaces)
end

function v = getf(s, fn, dflt)
if isfield(s, fn) && ~isempty(s.(fn))
    v = char(s.(fn));
else
    v = dflt;
end
end

function s = fmt_int(v)
s = sprintf('%d', round(double(v)));
end

function s = fmt_dbl(v)
v = double(v);
if v == floor(v) && abs(v) < 1e15
    s = sprintf('%d', int64(v));
else
    s = strip_trailing_zeros(sprintf('%.6f', v));
end
end

function s = strip_trailing_zeros(s)
if any(s == '.')
    while numel(s) > 1 && s(end) == '0', s(end) = []; end
    if s(end) == '.', s(end) = []; end
end
end

function deal_field(fid, sh, n, width, fname)
buf = repmat(' ', 1, width * n);
for i = 1:n
    val = sh(i).(fname);
    if isempty(val), continue; end
    val = char(val);
    k = min(numel(val), width);
    buf((i-1)*width + 1 : (i-1)*width + k) = val(1:k);
end
fwrite(fid, buf, 'char');
end

function write_doubles(fid, sh, n, width, fname)
buf = repmat(' ', 1, width * n);
for i = 1:n
    val = fmt_dbl(sh(i).(fname));
    k = min(numel(val), width);
    buf((i-1)*width + 1 : (i-1)*width + k) = val(1:k);
end
fwrite(fid, buf, 'char');
end

function write_ints(fid, sh, n, width, fname)
buf = repmat(' ', 1, width * n);
for i = 1:n
    val = fmt_int(sh(i).(fname));
    k = min(numel(val), width);
    buf((i-1)*width + 1 : (i-1)*width + k) = val(1:k);
end
fwrite(fid, buf, 'char');
end

function [rec, ok] = append_tal_m(rec, onset, text)
ts = sprintf('%+g', onset);
b = uint8([ts 20 uint8(text) 20 0]);
need = numel(b);
if rec.used + need > numel(rec.buf)
    ok = false;
    return
end
rec.buf(rec.used + 1 : rec.used + need) = b;
rec.used = rec.used + need;
ok = true;
end

function fclose_safe(fid)
if fid > 0
    try, fclose(fid); catch, end
end
end
