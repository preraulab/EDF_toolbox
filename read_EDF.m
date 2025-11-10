function varargout = read_EDF(edf_fname, varargin)
%READ_EDF  Load EDF or EDF+ file with full metadata, annotations, and MEX acceleration
%
%   READ_EDF reads European Data Format (EDF/EDF+) files using a compiled MEX
%   reader when available, and a pure MATLAB fallback otherwise. The function
%   provides full access to header metadata, per-signal scaling (digital-to-
%   physical conversion), and EDF+ annotations.
%
%   Usage:
%       [hdr, sh, sig, ann] = read_EDF(filename, 'Channels', {'EEG Fpz-Cz'})
%
%   Inputs:
%       edf_fname      : string - EDF or EDF+ file path
%       'Channels'     : cell array - subset of channels to read (default: all)
%       'Epochs'       : 1x2 vector [start_epoch end_epoch] (0-indexed, default: all)
%       'Verbose'      : logical - print progress and status info (default: false)
%       'RepairHeader' : logical - correct invalid record counts (default: false)
%       'forceMATLAB'  : logical - disable MEX usage (default: false)
%
%   Outputs:
%       hdr   : structure containing EDF file-level metadata
%       sh    : structure array of per-signal headers
%       sig   : cell array containing each signal vector (in physical units)
%       ann   : structure array of EDF+ annotations with onset and text
%
%   Example:
%       [hdr, sh, sig, ann] = read_EDF('sleep.edf', 'Channels', {'EEG C3-A2'});
%
%   -------------------------------------------------------------------------
%   EDF File Specification Summary:
%       • Each EDF file begins with a fixed-length 256-byte main header
%       • Followed by per-signal headers (16 fields × N signals)
%       • Digital samples stored as int16 are scaled to physical units:
%
%             phys_val = phys_min + (dig_val - dig_min) * (phys_max - phys_min) / (dig_max - dig_min)
%
%       • EDF+ annotation channels (labelled 'EDF Annotations') contain onset
%         times and event texts in TAL (Time-Annotation List) format.


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
parse(p, varargin{:});
channels = p.Results.Channels;
epochs = p.Results.Epochs;
verbose = p.Results.Verbose;
repair_header = p.Results.RepairHeader;
force_matlab = p.Results.forceMATLAB;

%% ---------------- MEX HANDLING ----------------
script_dir = fileparts(mfilename('fullpath'));
mex_file = fullfile(script_dir, ['read_EDF_mex.' mexext]);
mex_exists = isfile(mex_file);

if ~force_matlab
    if ~mex_exists
        try
            if verbose, fprintf('Compiling read_EDF_mex.c...\n'); end
            cd(script_dir);
            mex('read_EDF_mex.c');
            cd(script_dir);
            mex_exists = isfile(mex_file);
        catch
            if verbose, fprintf('MEX compilation failed. Using MATLAB fallback.\n'); end
            mex_exists = false;
        end
    end
    if mex_exists
        try
            [varargout{1:nargout}] = read_EDF_mex(edf_fname, channels, epochs, verbose, repair_header);
            return;
        catch ME
            if verbose
                fprintf('MEX failed (%s). Falling back to MATLAB reader.\n', ME.message);
            end
        end
    end
else
    if verbose, fprintf('forceMATLAB = true, skipping MEX reader\n'); end
end

%% ---------------- MATLAB FALLBACK ----------------
[varargout{1:nargout}] = read_EDF_matlab(edf_fname, channels, epochs, verbose, repair_header);
end


%% =========================================================================
%  PURE MATLAB EDF READER
% =========================================================================
function varargout = read_EDF_matlab(edf_fname, channels, epochs, verbose, repair_header)
fid = fopen(edf_fname, 'r', 'ieee-le');
if fid < 0
    error('read_EDF:FileError', 'Cannot open file: %s', edf_fname);
end

if verbose, fprintf('Reading EDF header...\n'); end

% ---------------- Read EDF Main Header ----------------
A = fread(fid, 256, 'uint8=>char')';
header_fields = {'edf_ver','patient_id','local_rec_id',...
    'recording_startdate','recording_starttime',...
    'num_header_bytes','reserve_1','num_data_records',...
    'data_record_duration','num_signals'};
field_sizes = [8,80,80,8,8,8,44,8,8,4];
loc = [0; cumsum(field_sizes(:))];
header = struct();
for k = 1:numel(header_fields)
    str = strtrim(A(loc(k)+1:loc(k+1)));
    switch header_fields{k}
        case {'num_header_bytes','num_data_records','data_record_duration','num_signals'}
            header.(header_fields{k}) = str2double(str);
        otherwise
            header.(header_fields{k}) = str;
    end
end

% ---------------- Add total data length fields ----------------
total_seconds = header.num_data_records * header.data_record_duration;
header.total_data_seconds = total_seconds;
dur = seconds(total_seconds);
dur.Format="hh:mm:ss.SSS";
header.total_data_hms = char(dur);


if verbose
    fprintf('Total data length: %.2f sec (%s)\n', header.total_data_seconds, header.total_data_hms);
end

% ---------------- Read Per-Signal Headers ----------------
fseek(fid, 256, 'bof');
num_signals = header.num_signals;
sig_header_size = header.num_header_bytes - 256;
A = fread(fid, sig_header_size, 'uint8=>char')';
sig_fields = {'signal_labels','transducer_type','physical_dimension',...
    'physical_min','physical_max','digital_min','digital_max',...
    'prefiltering','samples_in_record','reserve_2'};
sig_field_sizes = [16,80,8,8,8,8,8,80,8,32];
loc = [0; cumsum(sig_field_sizes(:)*num_signals)];
signal_header = struct();
for f = 1:numel(sig_fields)
    block = A(loc(f)+1:loc(f+1));
    for s = 1:num_signals
        start_idx = (s-1)*sig_field_sizes(f)+1;
        end_idx = s*sig_field_sizes(f);
        val = strtrim(block(start_idx:end_idx));
        if any(strcmp(sig_fields{f},{'physical_min','physical_max','digital_min','digital_max','samples_in_record'}))
            val = str2double(val);
        end
        signal_header(s).(sig_fields{f}) = val;
    end
end

% ---------------- Validate Channel List ----------------
labels = strtrim({signal_header.signal_labels});
if isempty(channels)
    signal_indices = 1:num_signals;
else
    signal_indices = find(ismember(labels, channels));
    invalid_channels = setdiff(channels, labels);
    if ~isempty(invalid_channels)
        warning('read_EDF:InvalidChannel', ...
            'Some channels not found: %s\nValid channels: %s', ...
            strjoin(invalid_channels, ', '), strjoin(labels, ', '));
    end
    if isempty(signal_indices)
        error('read_EDF:NoValidChannels', ...
            'No valid channels found.\nRequested: %s\nAvailable: %s', ...
            strjoin(channels, ', '), strjoin(labels, ', '));
    end
end

% ---------------- Prepare Data Read ----------------
samples_per_record = [signal_header.samples_in_record];
record_size = sum(samples_per_record);
fseek(fid, header.num_header_bytes, 'bof');
raw = fread(fid, inf, 'int16=>double', 'ieee-le');
fclose(fid);

num_records = header.num_data_records;
if isempty(epochs)
    start_epoch = 1;
    end_epoch = num_records;
else
    start_epoch = epochs(1) + 1;
    end_epoch = epochs(2);
end
num_epochs = end_epoch - start_epoch + 1;

% ---------------- EDF Digital-to-Physical Conversion ----------------
signal_cell = cell(1, length(signal_indices));
for i = 1:length(signal_indices)
    sidx = signal_indices(i);
    sig_offset = sum(samples_per_record(1:sidx-1));
    num_samples = samples_per_record(sidx) * num_epochs;
    sig = zeros(1, num_samples, 'double');

    dig_min = double(signal_header(sidx).digital_min);
    dig_max = double(signal_header(sidx).digital_max);
    phy_min = double(signal_header(sidx).physical_min);
    phy_max = double(signal_header(sidx).physical_max);

    if dig_max < dig_min, [dig_min, dig_max] = deal(dig_max, dig_min); end
    if phy_max < phy_min, [phy_min, phy_max] = deal(phy_max, phy_min); end
    if dig_max == dig_min
        scale = 0;
    else
        scale = (phy_max - phy_min) / (dig_max - dig_min);
    end

    if verbose, fprintf('Reading signal %s...\n', signal_header(sidx).signal_labels); end

    for r = 1:num_epochs
        record_idx = start_epoch + r - 1;
        rec_start = (record_idx - 1) * record_size + sig_offset + 1;
        rec_end = rec_start + samples_per_record(sidx) - 1;
        raw_vals = double(raw(rec_start:rec_end));
        sig_start = (r - 1) * samples_per_record(sidx) + 1;
        sig_end = r * samples_per_record(sidx);
        sig(sig_start:sig_end) = phy_min + (raw_vals - dig_min) * scale;
    end
    signal_cell{i} = sig;
end

signal_header_out = signal_header(signal_indices);
annotations = extractAnnotations(edf_fname, header, signal_header);

% ---------------- Output Assignment ----------------
varargout{1} = header;
if nargout > 1, varargout{2} = signal_header_out; end
if nargout > 2, varargout{3} = signal_cell; end
if nargout > 3, varargout{4} = annotations; end
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
