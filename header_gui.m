function [header_tbl, signal_tbl, uifig, uitable_header, uitable_signal] = header_gui(header, signal_header, varargin)
%HEADER_GUI  Display EDF header and signal header information in UI tables
%
%   Usage:
%       [header_tbl, signal_tbl] = header_gui(header, signal_header)
%       [header_tbl, signal_tbl, uifig, uitable_header, uitable_signal] = header_gui(header, signal_header)
%       [...] = header_gui(header, signal_header, 'Parent', parentHandle)
%       [...] = header_gui(header, signal_header, 'CreateGUI', false)
%
%   Input:
%       header: struct - EDF main header structure returned by read_EDF -- required
%       signal_header: struct array - EDF per-signal header structures -- required
%
%   Name-Value Pairs:
%       'Parent': graphics handle - UIFIGURE or UIPANEL to contain the tables
%                 (default: new UIFIGURE)
%       'CreateGUI': logical - create GUI and uitables (default: true)
%
%   Output:
%       header_tbl: table - prettified table version of header struct
%       signal_tbl: table - prettified table version of signal_header struct
%       uifig: UIFIGURE handle (empty if GUI not created or Parent supplied)
%       uitable_header: uitable handle displaying header fields
%       uitable_signal: uitable handle displaying signal header fields
%
%   See also:
%       read_EDF, uitable, uifigure, struct2table

%************************************************************
%                 DEFINE REQUIRED EDF FIELDS
%************************************************************
required_header_fields = { ...
    'edf_ver', 'patient_id', 'local_rec_id', ...
    'recording_startdate', 'recording_starttime', ...
    'num_header_bytes', 'num_data_records', ...
    'data_record_duration', 'num_signals', ...
    'total_data_seconds', 'total_data_hms'};

required_signal_fields = { ...
    'signal_labels', 'transducer_type', 'physical_dimension', ...
    'physical_min', 'physical_max', ...
    'digital_min', 'digital_max', ...
    'prefiltering', 'samples_in_record', ...
    'sampling_frequency'};

%************************************************************
%                 PARSE INPUTS AND OPTIONS
%************************************************************
p = inputParser;
p.FunctionName = 'header_gui';

addRequired(p, 'header', @(x) validateHeaderStruct(x, required_header_fields));
addRequired(p, 'signal_header', @(x) validateSignalHeaderStruct(x, required_signal_fields));

addParameter(p, 'Parent', [], @(x) isempty(x) || isgraphics(x));
addParameter(p, 'CreateGUI', true, @(x) islogical(x) && isscalar(x));

parse(p, header, signal_header, varargin{:});

parent    = p.Results.Parent;
createGUI = p.Results.CreateGUI;

%************************************************************
%              CONVERT STRUCTS TO TABLES
%************************************************************
header_tbl = struct2table(header, 'AsArray', true);
header_tbl.Properties.VariableNames = prettifyNames(header_tbl.Properties.VariableNames);

signal_tbl = struct2table(signal_header);
if ~isempty(signal_tbl)
    signal_tbl = addvars(signal_tbl, (1:height(signal_tbl))', ...
        'Before', 1, 'NewVariableNames', 'Channel');
    signal_tbl.Properties.VariableNames = prettifyNames(signal_tbl.Properties.VariableNames);
end

%************************************************************
%              EXIT EARLY IF GUI NOT REQUESTED
%************************************************************
if ~createGUI
    uifig = [];
    uitable_header = [];
    uitable_signal = [];
    return;
end

%************************************************************
%               HANDLE EMPTY SIGNAL HEADER
%************************************************************
if isempty(signal_header)
    if isempty(parent)
        uifig = uifigure('Name','EDF Header', 'Position',[100 100 800 200]);
        parent = uifig;
    else
        uifig = [];
    end
    uilabel(parent, ...
        'Text','No signal header data available.', ...
        'Position',[20 20 400 30]);
    uitable_header = [];
    uitable_signal = [];
    return;
end

%************************************************************
%              CREATE UIFIGURE IF NEEDED
%************************************************************
if isempty(parent)
    screenSize = get(0, 'ScreenSize');
    figWidth  = min(1500, screenSize(3));
    figHeight = min(700,  screenSize(4));
    x = (screenSize(3) - figWidth)/2;
    y = (screenSize(4) - figHeight)/2;
    uifig = uifigure('Name','EDF Header Information', ...
        'Position',[x y figWidth figHeight]);
    parent = uifig;
else
    uifig = [];
end

%************************************************************
%            POSITION TABLES IN PIXELS
%************************************************************
originalUnits = parent.Units;
parent.Units = 'pixels';
parentPos = parent.Position;
parentWidth  = parentPos(3);
parentHeight = parentPos(4);

padding = 20;
rowHeight = 22;
headerHeight = rowHeight + 35;
signalHeight = parentHeight - headerHeight - 3*padding;

uitable_header = uitable(parent, ...
    'Data', header_tbl, ...
    'Position', [padding, parentHeight - headerHeight - padding, ...
                 parentWidth - 2*padding, headerHeight], ...
    'ColumnWidth','auto');

uitable_signal = uitable(parent, ...
    'Data', signal_tbl, ...
    'Position', [padding, padding, ...
                 parentWidth - 2*padding, signalHeight], ...
    'ColumnWidth','auto');

%************************************************************
%            CONVERT TO NORMALIZED UNITS
%************************************************************
uitable_header.Units = 'normalized';
uitable_signal.Units = 'normalized';
parent.Units = originalUnits;

end

%************************************************************
%          VALIDATE MAIN EDF HEADER STRUCTURE
%************************************************************
function tf = validateHeaderStruct(h, required_fields)
tf = false;
if ~isstruct(h) || ~isscalar(h)
    error('header_gui:InvalidHeader', ...
        'header must be a scalar struct returned by read_EDF.');
end
missing = setdiff(required_fields, fieldnames(h));
if ~isempty(missing)
    error('header_gui:MissingHeaderFields', ...
        'Header missing required fields: %s', strjoin(missing, ', '));
end
tf = true;
end

%************************************************************
%        VALIDATE PER-SIGNAL HEADER STRUCTURE
%************************************************************
function tf = validateSignalHeaderStruct(sh, required_fields)
tf = false;
if ~isstruct(sh)
    error('header_gui:InvalidSignalHeader', ...
        'signal_header must be a struct array returned by read_EDF.');
end
if isempty(sh)
    tf = true;
    return;
end
fn = fieldnames(sh);
missing = setdiff(required_fields, fn);
if ~isempty(missing)
    error('header_gui:MissingSignalHeaderFields', ...
        'Signal header missing required fields: %s', strjoin(missing, ', '));
end
tf = true;
end

%************************************************************
%                PRETTIFY VARIABLE NAMES
%************************************************************
function names = prettifyNames(names)
for k = 1:numel(names)
    name = strrep(names{k}, '_', ' ');
    parts = strsplit(lower(name), ' ');
    for p = 1:numel(parts)
        if ~isempty(parts{p})
            parts{p}(1) = upper(parts{p}(1));
        end
    end
    names{k} = strjoin(parts, ' ');
end
end
