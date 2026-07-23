function compile_edf_mex(script_dir, src_name)
%COMPILE_EDF_MEX  Build a MEX file in this directory against zlib and zstd
%
%   Usage:
%       compile_edf_mex(script_dir, src_name)
%
%   Inputs:
%       script_dir : char - directory containing the .c source and the
%                    bundled 'zlib/' and 'zstd/' subdirs -- required
%       src_name   : char - .c filename (e.g. 'read_EDF_mex.c') -- required
%
%   Outputs:
%       none (side effects only: writes the compiled MEX into script_dir)
%
%   Notes:
%       Used by read_EDF.m and write_EDF.m to auto-compile their MEX
%       backends on first call.
%
%       Linkage:
%           zlib : bundled source in script_dir/zlib if present, else
%                  system -lz.
%           zstd : bundled single-file amalgamation in script_dir/zstd if
%                  present (works on all platforms, including Windows,
%                  with no system library). Falls back to system libzstd:
%                  searches /opt/homebrew (Apple Silicon Homebrew) and
%                  /usr/local (Intel Homebrew / Linux), else plain -lzstd.
%
%       Requires a configured C compiler (see `mex -setup C`). On Windows,
%       install the free "MATLAB Support for MinGW-w64 C/C++ Compiler"
%       add-on (Home tab > Add-Ons) if no compiler is found; the error
%       raised here includes these instructions.
%
%   See also: read_EDF, write_EDF
%
%   ∿∿∿  Prerau Laboratory · sleepEEG.org  ∿∿∿

zlib_dir = fullfile(script_dir, 'zlib');
zstd_dir = fullfile(script_dir, 'zstd');
mex_src  = fullfile(script_dir, src_name);

if ~isfile(mex_src)
    error('compile_edf_mex:NoSource', 'Source file not found: %s', mex_src);
end

% --- Check for a configured C compiler ----------------------------------
% mex.getCompilerConfigurations returns empty when no C compiler is set
% up. Failing here with instructions beats the opaque error mex() throws.
try
    cc = mex.getCompilerConfigurations('C', 'Selected');
catch
    cc = [];  % very old MATLAB without this API; let mex() try anyway
end
if isempty(cc)
    if ispc
        error('compile_edf_mex:NoCompiler', ...
            ['No C compiler is configured for MEX. On Windows, install the free\n' ...
             '"MATLAB Support for MinGW-w64 C/C++ Compiler" add-on:\n' ...
             '    MATLAB Home tab > Add-Ons > Get Add-Ons > search "MinGW"\n' ...
             'then run:  mex -setup C\n' ...
             '(read_EDF/write_EDF will fall back to the pure-MATLAB backend\n' ...
             'until a compiler is available, at reduced speed.)']);
    else
        error('compile_edf_mex:NoCompiler', ...
            ['No C compiler is configured for MEX. Install one (macOS: Xcode\n' ...
             'Command Line Tools via `xcode-select --install`; Linux: gcc)\n' ...
             'then run:  mex -setup C']);
    end
end

% --- Locate zstd --------------------------------------------------------
% Prefer the bundled single-file amalgamation (script_dir/zstd/zstd.c):
% compiles anywhere, no system libzstd needed. Otherwise try Homebrew
% prefixes, falling back to letting the linker find it via -lzstd.
zstd_bundled  = isfile(fullfile(zstd_dir, 'zstd.c'));
zstd_inc_args = {};
zstd_lib_args = {};
zstd_src_args = {};
if zstd_bundled
    zstd_inc_args = {['-I' zstd_dir]};
    zstd_src_args = {fullfile(zstd_dir, 'zstd.c')};
    % The amalgamation bakes in ZSTD_MULTITHREAD: Win32 threads on
    % Windows (no extra flags), pthreads elsewhere. -lpthread is
    % harmless where pthreads live in libc (macOS, glibc >= 2.34).
    if ~ispc
        zstd_lib_args = {'-lpthread'};
    end
else
    zstd_prefixes = {'/opt/homebrew', '/usr/local'};
    for k = 1:numel(zstd_prefixes)
        pfx = zstd_prefixes{k};
        if isfile(fullfile(pfx, 'include', 'zstd.h'))
            zstd_inc_args = {['-I' fullfile(pfx, 'include')]};
            zstd_lib_args = {['-L' fullfile(pfx, 'lib')]};
            break;
        end
    end
    zstd_lib_args = [zstd_lib_args, {'-lzstd'}];
    if ispc
        warning('compile_edf_mex:NoBundledZstd', ...
            ['Bundled zstd sources not found in %s and Windows has no\n' ...
             'system libzstd by default; the build will likely fail.\n' ...
             'Restore the toolbox''s zstd/ directory (zstd.c + zstd.h).'], ...
            zstd_dir);
    end
end

% --- Build args ---------------------------------------------------------
common_args = {'-O', '-largeArrayDims', '-outdir', script_dir};

if isfolder(zlib_dir)
    fprintf('Compiling %s with bundled zlib + %s zstd...\n', src_name, ...
        ternary(zstd_bundled, 'bundled', 'system'));
    zlib_files = dir(fullfile(zlib_dir, '*.c'));
    zlib_paths = arrayfun(@(f) fullfile(f.folder, f.name), zlib_files, ...
        'UniformOutput', false);
    args = [common_args, {['-I' zlib_dir]}, zstd_inc_args, ...
            {mex_src}, zlib_paths(:)', zstd_src_args, zstd_lib_args];
else
    fprintf('Compiling %s with system zlib + %s zstd...\n', src_name, ...
        ternary(zstd_bundled, 'bundled', 'system'));
    args = [common_args, zstd_inc_args, ...
            {mex_src}, zstd_src_args, {'-lz'}, zstd_lib_args];
end

mex(args{:});
end

function out = ternary(cond, a, b)
%TERNARY  Return a if cond is true, else b
%
%   Inputs:
%       cond : logical - condition to test -- required
%       a    : any - value returned when cond is true -- required
%       b    : any - value returned when cond is false -- required
%
%   Outputs:
%       out : any - a or b
if cond, out = a; else, out = b; end
end
