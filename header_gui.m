function [header_tbl, signal_tbl, uifig, uitable_header, uitable_signal] = header_gui(header, signal_header, varargin)
%HEADER_GUI  Display EDF header and signal header information in UI tables
%
%   Usage:
%       [header_tbl, signal_tbl] = header_gui(header, signal_header)
%       [header_tbl, signal_tbl, uifig, uitable_header, uitable_signal] = header_gui(header, signal_header)
%       [...] = header_gui(header, signal_header, 'Parent', parentHandle)
%
%   Input:
%       header: struct - EDF main header structure -- required
%       signal_header: struct array - EDF per-signal header structures -- required
%
%   Name-Value Pairs:
%       'Parent': graphics handle - UIFIGURE or UIPANEL to contain the tables
%                 (default: new UIFIGURE)
%
%   Output:
%       header_tbl: table - prettified table version of header struct
%       signal_tbl: table - prettified table version of signal_header struct
%       uifig: UIFIGURE handle (empty if Parent is supplied)
%       uitable_header: uitable handle displaying header fields
%       uitable_signal: uitable handle displaying signal header fields
%
%   Notes:
%       - If only header_tbl and/or signal_tbl are requested, no GUI is created.
%       - If additional outputs are requested, or if no outputs are requested,
%         the GUI and uitables are created.
%       - Tables automatically resize with the parent container.
%
%   See also:
%       uitable, uifigure, struct2table

%************************************************************
%                 PARSE NAME-VALUE PAIRS
%************************************************************
p = inputParser;
addParameter(p, 'Parent', [], @(x) isempty(x) || isgraphics(x));
parse(p, varargin{:});
parent = p.Results.Parent;

%************************************************************
%               VALIDATE INPUT STRUCTURES
%************************************************************
if ~isstruct(header) || isempty(header)
    error('header_gui:InvalidInput', 'Header must be a non-empty structure.');
end
if ~isstruct(signal_header)
    error('header_gui:InvalidInput', 'Signal header must be a structure array.');
end

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
%        EARLY RETURN IF ONLY TABLES ARE REQUESTED
%************************************************************
if nargout <= 2
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
