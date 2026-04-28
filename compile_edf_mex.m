function compile_edf_mex(script_dir, src_name)
%COMPILE_EDF_MEX  Build a MEX file in this directory against zlib and libzstd.
%
%   Used by read_EDF.m and write_EDF.m to auto-compile their MEX backends
%   on first call.
%
%   Inputs:
%       script_dir : directory containing the .c source and (optionally)
%                    the bundled 'zlib/' subdir
%       src_name   : .c filename (e.g. 'read_EDF_mex.c')
%
%   Linkage:
%       zlib  : bundled in script_dir/zlib if present, else system -lz
%       zstd  : system libzstd. Searches /opt/homebrew (Apple Silicon
%               Homebrew) and /usr/local (Intel Homebrew / Linux) for
%               headers and libraries; falls back to plain -lzstd.

zlib_dir = fullfile(script_dir, 'zlib');
mex_src  = fullfile(script_dir, src_name);

if ~isfile(mex_src)
    error('compile_edf_mex:NoSource', 'Source file not found: %s', mex_src);
end

% --- Locate libzstd ------------------------------------------------------
% Try Homebrew prefixes; fall back to letting the linker find it via -lzstd.
zstd_inc_args = {};
zstd_lib_args = {};
zstd_prefixes = {'/opt/homebrew', '/usr/local'};
for k = 1:numel(zstd_prefixes)
    pfx = zstd_prefixes{k};
    if isfile(fullfile(pfx, 'include', 'zstd.h'))
        zstd_inc_args = {['-I' fullfile(pfx, 'include')]};
        zstd_lib_args = {['-L' fullfile(pfx, 'lib')]};
        break;
    end
end

% --- Build args ---------------------------------------------------------
common_args = {'-O', '-largeArrayDims', '-outdir', script_dir};

if isfolder(zlib_dir)
    fprintf('Compiling %s with bundled zlib + system libzstd...\n', src_name);
    zlib_files = dir(fullfile(zlib_dir, '*.c'));
    zlib_paths = arrayfun(@(f) fullfile(f.folder, f.name), zlib_files, ...
        'UniformOutput', false);
    args = [common_args, {['-I' zlib_dir]}, zstd_inc_args, ...
            zstd_lib_args, {mex_src}, zlib_paths(:)', {'-lzstd'}];
else
    fprintf('Compiling %s with system zlib + libzstd...\n', src_name);
    args = [common_args, zstd_inc_args, zstd_lib_args, ...
            {mex_src}, {'-lz', '-lzstd'}];
end

mex(args{:});
end
