function varargout = read_EDF(edf_fname, varargin)
%READ_EDF  Load EDF / EDF+ / EDF.gz / EDF.zst file with full metadata, annotations, and MEX acceleration
%
%   READ_EDF reads European Data Format (EDF / EDF+) files using a compiled
%   MEX reader when available, and a pure MATLAB fallback otherwise. Plain
%   '.edf', gzip-compressed '.edf.gz', and zstd-compressed '.edf.zst' files
%   are all accepted; compressed files are streamed directly through zlib
%   or libzstd with no temp file (MEX path). The function provides full
%   access to header metadata, per-signal scaling (digital-to-physical
%   conversion), and EDF+ annotations.
%
%   Usage:
%       [header, signal_header, signal_cell, annotations] = read_EDF(filename, 'Channels', {'EEG Fpz-Cz'})
%
%   Inputs:
%       edf_fname      : string - path to a .edf or .edf.gz file
%       'Channels'     : cell array - subset of channels to read (default: all)
%                        Use 'ChA-ChB' to load ChA minus ChB as a rereferenced channel.
%       'Epochs'       : 1x2 vector [start_epoch end_epoch] (0-indexed, default: all)
%       'Verbose'      : logical - print progress and status info (default: false)
%       'RepairHeader' : logical - correct invalid record counts and save with _fixed suffix
%                        (default: false; ignored with a warning for .gz inputs)
%       'forceMATLAB'  : logical - disable MEX usage (default: false). For .gz inputs
%                        the MATLAB path decompresses to a temp file (auto-cleaned up).
%       'debug'        : logical - debug mode for MEX (default: false)
%       'deidentify'   : logical - overwrite PHI fields and save with _deidentified suffix
%                        (default: false; ignored with a warning for .gz inputs)
%
%   Outputs:
%       header        : structure containing EDF file-level metadata
%       signal_header : structure array of per-signal headers
%       signal_cell   : cell array containing each signal vector (in physical units)
%       annotations   : structure array of EDF+ annotations with onset and text
%
%   Examples:
%       % Plain EDF
%       [hdr, shdr, sc, ann] = read_EDF('sleep.edf', 'Channels', {'EEG C3-A2'});
%
%       % Gzip-compressed EDF (streamed via zlib, no temp file)
%       [hdr, shdr, sc, ann] = read_EDF('sleep.edf.gz');
%
%   -------------------------------------------------------------------------
%   EDF File Specification Summary:
%       • Each EDF file begins with a fixed-length 256-byte main header
%       • Followed by per-signal headers (16 fields x N signals)
%       • Digital samples stored as int16 are scaled to physical units:
%
%             phys_val = phys_min + (dig_val - dig_min) * (phys_max - phys_min) / (dig_max - dig_min)
%
%       • EDF+ annotation channels (labelled 'EDF Annotations') contain onset
%         times and event texts in TAL (Time-Annotation List) format.
%
%   -------------------------------------------------------------------------
%   Compressed input notes:
%       • Files ending in '.gz' (case-insensitive) are read directly through
%         the bundled zlib (vendored in the 'zlib/' subdirectory).
%       • Peak memory is one record's worth of raw bytes plus the per-signal
%         output arrays — files do not need to fit decompressed on disk.
%       • RepairHeader and deidentify both modify the file in place, which is
%         not meaningful on a gzip archive; both are silently disabled (with a
%         warning) for .gz inputs. Decompress to .edf first if you need either.

%% ---------------- INPUT PARSING ----------------
if nargin < 1
    error('read_EDF:InvalidInput', 'EDF filename required.');
end

edf_fname = char(edf_fname);
if ~isfile(edf_fname)
    error('read_EDF:FileNotFound', 'EDF file not found: %s', edf_fname);
end

p = inputParser;
addParameter(p, 'Channels', {}, @iscell);
addParameter(p, 'Epochs', [], @isnumeric);
addParameter(p, 'Verbose', false, @islogical);
addParameter(p, 'RepairHeader', false, @islogical);
addParameter(p, 'forceMATLAB', false, @islogical);
addParameter(p, 'debug', false, @islogical);
addParameter(p, 'deidentify', false, @islogical);
parse(p, varargin{:});

channels        = p.Results.Channels;
epochs          = p.Results.Epochs;
verbose         = p.Results.Verbose;
repair_header   = p.Results.RepairHeader;
force_matlab    = p.Results.forceMATLAB;
debug           = p.Results.debug;
deidentify      = p.Results.deidentify;

%% ---------------- COMPRESSED INPUT GATING ----------------
% RepairHeader and deidentify both require modifying the file in place,
% which is not meaningful on a compressed archive. Force them off and warn.
is_gz  = endsWith(edf_fname, '.gz',  'IgnoreCase', true);
is_zst = endsWith(edf_fname, '.zst', 'IgnoreCase', true);
is_compressed = is_gz || is_zst;
if is_compressed
    suffix = '.gz'; if is_zst, suffix = '.zst'; end
    if repair_header
        warning('read_EDF:RepairOnCompressed', ...
            'RepairHeader is not supported for %s inputs; in-memory header will still be corrected.', suffix);
        repair_header = false;
    end
    if deidentify
        warning('read_EDF:DeidentifyOnCompressed', ...
            'deidentify is not supported for %s inputs; decompress to .edf first if you need a deidentified copy.', suffix);
        deidentify = false;
    end
end

%% ---------------- MEX HANDLING ----------------
script_dir = fileparts(mfilename('fullpath'));
mex_file = fullfile(script_dir, ['read_EDF_mex.' mexext]);
mex_exists = isfile(mex_file);

if ~force_matlab
    if ~mex_exists
        compile_edf_mex(script_dir, 'read_EDF_mex.c');
    end
    try
        [varargout{1:nargout}] = read_EDF_mex(edf_fname, channels, epochs, verbose, repair_header, debug);

        if nargout>0
            %Add the total data in seconds
            total_seconds = varargout{1}.num_data_records * varargout{1}.data_record_duration;
            varargout{1}.total_data_seconds = total_seconds;
            varargout{1}.total_data_hms = char(duration(0,0,total_seconds));
        end

        %Remove the whitespace around the labels
        if nargout>1
            new_signal_labels = cellfun(@strip,{varargout{2}.signal_labels},'UniformOutput', false);
            [varargout{2}.signal_labels] = deal(new_signal_labels{:});
        end

        % Deidentify post-MEX if requested (MEX does not handle this)
        if deidentify
            deidentify_edf(edf_fname, repair_header, verbose);
        end
        return
    catch ME
        if verbose
            fprintf('MEX failed (%s). Falling back to MATLAB reader.\n', ME.message);
        end
    end

else
    if verbose
        fprintf('forceMATLAB = true, using MATLAB reader.\n');
    end
end

%% ---------------- MATLAB FALLBACK ----------------
[varargout{1:nargout}] = read_EDF_matlab(edf_fname, channels, epochs, verbose, repair_header, deidentify);

end


%% =========================================================================
%  PURE MATLAB EDF READER
% =========================================================================
function varargout = read_EDF_matlab(edf_fname, channels, epochs, verbose, repair_header, deidentify)

% If input is compressed, decompress to a temp .edf and run the rest
% against it. Temp directory is removed when this function returns.
comp_cleanup = []; %#ok<NASGU>
if endsWith(edf_fname, '.gz', 'IgnoreCase', true)
    if verbose
        fprintf('Decompressing %s to temp file for MATLAB reader...\n', edf_fname);
    end
    tmpdir = tempname;
    mkdir(tmpdir);
    comp_cleanup = onCleanup(@() rmdir(tmpdir, 's'));
    gunzipped = gunzip(edf_fname, tmpdir);
    edf_fname = gunzipped{1};
elseif endsWith(edf_fname, '.zst', 'IgnoreCase', true)
    if verbose
        fprintf('Decompressing %s to temp file for MATLAB reader...\n', edf_fname);
    end
    tmpdir = tempname;
    mkdir(tmpdir);
    comp_cleanup = onCleanup(@() rmdir(tmpdir, 's'));
    [~, base, ~] = fileparts(edf_fname);   % strip .zst -> base.edf
    out_path = fullfile(tmpdir, base);
    cmd = sprintf('zstd -q -d -f -o %s %s', shell_quote(out_path), shell_quote(edf_fname));
    [rc, msg] = system(cmd);
    if rc ~= 0
        error('read_EDF:ZstdDecompress', ...
            'zstd command failed for %s: %s', edf_fname, strtrim(msg));
    end
    edf_fname = out_path;
end

fid = fopen(edf_fname, 'r', 'ieee-le');
if fid < 0
    error('read_EDF:FileError', 'Cannot open file.');
end

%% ---------------- READ MAIN HEADER ----------------
A = fread(fid, 256, 'uint8=>char')';
field_sizes = [8,80,80,8,8,8,44,8,8,4];
fields = {'edf_ver','patient_id','local_rec_id',...
    'recording_startdate','recording_starttime',...
    'num_header_bytes','reserve_1','num_data_records',...
    'data_record_duration','num_signals'};

loc = [0; cumsum(field_sizes(:))];
header = struct();

for k = 1:numel(fields)
    str = strtrim(A(loc(k)+1:loc(k+1)));
    if any(strcmp(fields{k},{'num_header_bytes','num_data_records','data_record_duration','num_signals'}))
        header.(fields{k}) = str2double(str);
    else
        header.(fields{k}) = str;
    end
end

%% ---------------- READ SIGNAL HEADERS ----------------
fseek(fid, 256, 'bof');

num_signals = header.num_signals;
sig_sizes = [16,80,8,8,8,8,8,80,8,32];
sig_fields = {'signal_labels','transducer_type','physical_dimension',...
    'physical_min','physical_max','digital_min','digital_max',...
    'prefiltering','samples_in_record','reserve_2'};

A = fread(fid, header.num_header_bytes-256, 'uint8=>char')';
loc = [0; cumsum(sig_sizes(:)*num_signals)];

signal_header = struct();
for f = 1:numel(sig_fields)
    block = A(loc(f)+1:loc(f+1));
    for s = 1:num_signals
        idx1 = (s-1)*sig_sizes(f)+1;
        idx2 = s*sig_sizes(f);
        val = strtrim(block(idx1:idx2));
        if any(strcmp(sig_fields{f},{'physical_min','physical_max','digital_min','digital_max','samples_in_record'}))
            val = str2double(val);
        end
        signal_header(s).(sig_fields{f}) = val;
    end
end

new_signal_labels = cellfun(@strip,{signal_header.signal_labels},'UniformOutput', false);
[signal_header.signal_labels] = deal(new_signal_labels{:});

%% ---------------- CHANNEL PLAN ----------------
% Inspect all requested channels. Identify plain channel requests and A-B
% rereference pairs. Determine the superset of raw signals to decode.
all_labels = {signal_header.signal_labels};

if isempty(channels)
    signal_indices = 1:num_signals;
    plain_chs   = {};
    reref_pairs = {};
else
    [plain_chs, reref_pairs, load_labels] = parse_channel_plan(channels, all_labels);
    signal_indices = find(ismember(lower(all_labels), lower(load_labels)));
    if isempty(signal_indices)
        error('read_EDF:NoChannels', 'No valid channels found.');
    end
end

%% ---------------- READ DATA ----------------
samples_per_record = [signal_header.samples_in_record];
record_size = sum(samples_per_record);

fseek(fid, header.num_header_bytes, 'bof');
raw = fread(fid, inf, 'int16=>double');
fclose(fid);

if isempty(raw)
    error('EDF contains no data samples.');
end

actual_records = floor(length(raw) / record_size);

if actual_records < 1
    error('File does not contain complete data records.');
end

if actual_records ~= header.num_data_records
    if verbose
        fprintf('Header record count mismatch. Using actual value: %d\n', actual_records);
    end

    if repair_header
        % Mirror MEX fix_num_records: write corrected count as left-justified
        % 8-char ASCII string at byte offset 236 (0-based), same as MEX does.
        fw = fopen(edf_fname, 'r+', 'ieee-le');
        if fw < 0
            warning('read_EDF:RepairFailed', ...
                'Could not open file for header repair: %s', edf_fname);
        else
            rec_str = sprintf('%-8d', actual_records);   % left-justified, space-padded to 8 chars
            fseek(fw, 236, 'bof');
            fwrite(fw, rec_str, 'char');
            fclose(fw);
            if verbose
                fprintf('Header repaired on disk: num_data_records written as %d.\n', actual_records);
            end
        end
    end

    header.num_data_records = actual_records;
end

%% ---------------- EPOCH VALIDATION ----------------
if isempty(epochs)
    start_epoch = 1;
    end_epoch = header.num_data_records;
else
    start_epoch = epochs(1) + 1;
    end_epoch   = epochs(2);
end

start_epoch = max(start_epoch,1);
end_epoch   = min(end_epoch, header.num_data_records);

if end_epoch < start_epoch
    error('Invalid epoch selection.');
end

num_epochs = end_epoch - start_epoch + 1;

%% ---------------- TOTAL DURATION ----------------
total_seconds = header.num_data_records * header.data_record_duration;
header.total_data_seconds = total_seconds;
header.total_data_hms = char(duration(0,0,total_seconds));

%% ---------------- SAMPLING FREQUENCY ----------------
for ii = 1:length(signal_header)
    signal_header(ii).sampling_frequency = ...
        signal_header(ii).samples_in_record / header.data_record_duration;
end

%% ---------------- DIGITAL TO PHYSICAL ----------------
% Decode only the channels in signal_indices (plain + reref constituents).
signal_cells = cell(1, length(signal_indices));

for i = 1:length(signal_indices)

    sidx = signal_indices(i);
    offset = sum(samples_per_record(1:sidx-1));
    samples_per_epoch = samples_per_record(sidx);
    total_samples = samples_per_epoch * num_epochs;

    sig = zeros(1, total_samples);

    dig_min = signal_header(sidx).digital_min;
    dig_max = signal_header(sidx).digital_max;
    phy_min = signal_header(sidx).physical_min;
    phy_max = signal_header(sidx).physical_max;

    scale = (phy_max - phy_min) / (dig_max - dig_min);

    for r = 1:num_epochs
        rec = start_epoch + r - 1;
        idx1 = (rec-1)*record_size + offset + 1;
        idx2 = idx1 + samples_per_epoch - 1;

        raw_vals = raw(idx1:idx2);
        sig((r-1)*samples_per_epoch+1 : r*samples_per_epoch) = ...
            phy_min + (raw_vals - dig_min)*scale;
    end

    signal_cells{i} = sig;

end

% Trim signal_header to only the loaded subset (signal_cells{i} <-> signal_header(i)).
signal_header = signal_header(signal_indices);

%% ---------------- REREFERENCING ----------------
% Compute A-B difference signals and remove any constituent channels that
% were not independently requested.
if ~isempty(channels)
    [signal_header, signal_cells] = apply_channel_plan(signal_header, signal_cells, plain_chs, reref_pairs);
end

annotations = extractAnnotations(edf_fname, header, signal_header);

%% ---------------- DEIDENTIFY ----------------
if deidentify
    deidentify_edf(edf_fname, repair_header, verbose);

    % Also scrub the returned header struct so it matches the written file
    header.patient_id        = 'X X X X';
    header.local_rec_id      = 'Startdate X X X X';
    header.recording_startdate = '01.01.01';
end

%% ---------------- OUTPUT ----------------
varargout{1} = header;
if nargout > 1, varargout{2} = signal_header; end
if nargout > 2, varargout{3} = signal_cells; end
if nargout > 3, varargout{4} = annotations; end

end


%% =========================================================================
%  CHANNEL PLAN PARSER
% =========================================================================
function [plain_chs, reref_pairs, load_labels] = parse_channel_plan(channels, all_labels)
%PARSE_CHANNEL_PLAN  Classify requested channels as plain or A-B reref pairs.
%
%   For each entry in channels:
%     - If it directly matches a label (case-insensitive) → plain channel.
%     - Otherwise try each '-' as a split point. The first split where both
%       sides are valid label names is treated as a rereference request.
%     - Anything else is passed through as a plain name (will error on load).
%
%   Returns:
%     plain_chs   : cell array of directly-requested channel name strings
%     reref_pairs : cell array of {chA, chB, 'chA-chB'} rows
%     load_labels : unique set of labels to actually decode (plain + constituents)

all_labels_lower = lower(cellfun(@strtrim, all_labels, 'UniformOutput', false));
plain_chs   = {};
reref_pairs = {};

for k = 1:numel(channels)
    ch = strtrim(channels{k});

    % Direct label match → plain channel
    if ismember(lower(ch), all_labels_lower)
        plain_chs{end+1} = ch;
        continue;
    end

    % Try each '-' as a split point, leftmost first
    dashes = strfind(ch, '-');
    found  = false;
    for d = dashes
        chA = strtrim(ch(1:d-1));
        chB = strtrim(ch(d+1:end));
        if ~isempty(chA) && ~isempty(chB) && ...
                ismember(lower(chA), all_labels_lower) && ...
                ismember(lower(chB), all_labels_lower)
            reref_pairs{end+1} = {chA, chB, ch};
            found = true;
            break;
        end
    end

    if ~found
        plain_chs{end+1} = ch;  % unknown; will produce an error on lookup
    end
end

% Superset of labels to decode = plain channels + both sides of every reref pair
load_labels = plain_chs;
for k = 1:numel(reref_pairs)
    load_labels{end+1} = reref_pairs{k}{1};
    load_labels{end+1} = reref_pairs{k}{2};
end
load_labels = unique(load_labels, 'stable');
end


%% =========================================================================
%  APPLY CHANNEL PLAN
% =========================================================================
function [sh_out, sc_out] = apply_channel_plan(signal_header, signal_cell, plain_chs, reref_pairs)
%APPLY_CHANNEL_PLAN  Filter to requested channels and compute A-B reref signals.
%
%   signal_header / signal_cell contain the loaded superset (plain + constituents).
%   Output keeps only the explicitly requested plain channels, then appends one
%   entry per reref pair. Constituent-only channels are dropped.

labels = {signal_header.signal_labels};

% Indices of plain channels in the loaded data (preserve request order)
plain_idx = zeros(1, numel(plain_chs));
for k = 1:numel(plain_chs)
    plain_idx(k) = find(strcmpi(labels, strtrim(plain_chs{k})), 1);
end
plain_idx = plain_idx(plain_idx > 0);

% Compute rereferenced signals
n = numel(reref_pairs);
reref_sc = cell(1, n);
if n > 0
    reref_sh(1:n) = signal_header(1);   % preallocate with matching struct fields
    for k = 1:n
        chA  = reref_pairs{k}{1};
        chB  = reref_pairs{k}{2};
        name = reref_pairs{k}{3};
        idxA = find(strcmpi(labels, chA), 1);
        idxB = find(strcmpi(labels, chB), 1);
        reref_sh(k)               = signal_header(idxA);
        reref_sh(k).signal_labels = name;
        reref_sc{k}               = signal_cell{idxA} - signal_cell{idxB};
    end
    sh_out = [signal_header(plain_idx), reref_sh];
else
    sh_out = signal_header(plain_idx);
end
sc_out = [signal_cell(plain_idx), reref_sc];
end


%% =========================================================================
%  DEIDENTIFY HELPER — copy file with PHI fields blanked in the header
% =========================================================================
function deidentify_edf(edf_fname, repair_header, verbose)
% Build output filename with appropriate suffix
[fdir, fname, fext] = fileparts(edf_fname);
if repair_header
    out_fname = fullfile(fdir, [fname '_fixed_deidentified' fext]);
else
    out_fname = fullfile(fdir, [fname '_deidentified' fext]);
end

% Copy the entire file first, then patch header fields in the copy
copyfile(edf_fname, out_fname);

fw = fopen(out_fname, 'r+', 'ieee-le');
if fw < 0
    warning('read_EDF:DeidentifyFailed', ...
        'Could not open output file for deidentification: %s', out_fname);
    return
end

% EDF main header field byte offsets and widths (0-based offsets):
%   patient_id          offset=8,   width=80
%   local_rec_id        offset=88,  width=80
%   recording_startdate offset=168, width=8
phi_fields = {
    'X X X X',          8,   80;   % patient_id
    'Startdate X X X X', 88,  80;   % local_rec_id
    '01.01.01',         168,  8;   % recording_startdate
    };

for k = 1:size(phi_fields, 1)
    val     = phi_fields{k,1};
    offset  = phi_fields{k,2};
    width   = phi_fields{k,3};

    % Left-justify value, space-pad to exact field width
    padded = sprintf('%-*s', width, val);
    padded = padded(1:width);   % truncate if somehow over width

    fseek(fw, offset, 'bof');
    fwrite(fw, padded, 'char');
end

fclose(fw);

if verbose
    fprintf('Deidentified file written: %s\n', out_fname);
end
end


%% =========================================================================
%  EDF+ ANNOTATION PARSER — exclude blank annotation texts
% =========================================================================
function annotations = extractAnnotations(edf_fname, header, signal_header)
annotations = struct('onset', {}, 'text', {});
annotIdx = find(strcmp(strtrim({signal_header.signal_labels}), 'EDF Annotations'), 1);
if isempty(annotIdx), return; end

fid = fopen(edf_fname, 'rb');
if fid < 0, return; end

annotBytePos = sum([signal_header(1:annotIdx-1).samples_in_record]) * 2;
annotBytesPerRecord = signal_header(annotIdx).samples_in_record * 2;
totalBytesPerRecord = sum([signal_header.samples_in_record]) * 2;

fseek(fid, header.num_header_bytes, 'bof');
allAnnotations = [];
while ~feof(fid)
    data = fread(fid, totalBytesPerRecord, 'uint8=>uint8');
    if length(data) < totalBytesPerRecord, break; end
    annotData = data(annotBytePos + 1 : annotBytePos + annotBytesPerRecord);
    tals = parseTALs(annotData);
    allAnnotations = [allAnnotations; tals];
end
fclose(fid);

% Keep only annotations with non-empty text
if ~isempty(allAnnotations)
    mask = arrayfun(@(a) any(~cellfun(@isempty, strtrim(a.texts))), allAnnotations);
    valid = allAnnotations(mask);
    if ~isempty(valid)
        annotations = struct('onset', num2cell([valid.onset])', ...
            'text', {valid.texts}');
    end
end
end

%% =========================================================================
%  PARSE EDF+ TAL STRINGS
% =========================================================================
function tals = parseTALs(data)
tals = [];
idx = 1;
while idx <= length(data)
    while idx <= length(data) && data(idx) == 0
        idx = idx + 1;
    end
    if idx > length(data), break; end

    if data(idx) ~= 43 && data(idx) ~= 45  % '+' or '-'
        idx = idx + 1; continue;
    end
    onsetStart = idx;
    while idx <= length(data) && data(idx) ~= 20
        idx = idx + 1;
    end
    if idx > length(data), break; end
    onset = str2double(char(data(onsetStart:idx-1))');
    idx = idx + 1;

    texts = {};
    while idx <= length(data) && data(idx) ~= 0
        textStart = idx;
        while idx <= length(data) && data(idx) ~= 20 && data(idx) ~= 0
            idx = idx + 1;
        end
        texts{end+1} = strtrim(char(data(textStart:idx-1))');
        if idx <= length(data) && data(idx) == 20
            idx = idx + 1;
        end
    end
    if any(~cellfun(@isempty, strtrim(texts)))
        tals = [tals; struct('onset', onset, 'texts', {texts})];
    end
    idx = idx + 1;
end
end

%% =========================================================================
%  SHELL QUOTING (POSIX, for system() calls)
% =========================================================================
function s = shell_quote(s)
s = ['''' strrep(s, '''', '''\''''') ''''];
end
