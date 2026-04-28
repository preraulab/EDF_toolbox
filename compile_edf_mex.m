function compile_edf_mex(script_dir, src_name)
%COMPILE_EDF_MEX  Build a MEX file in this directory against the bundled zlib.
%
%   Used by read_EDF.m and write_EDF.m to auto-compile their MEX backends
%   on first call. Falls back to system zlib (-lz) if the vendored
%   'zlib/' directory is absent.
%
%   Inputs:
%       script_dir : directory containing the .c source and (optionally)
%                    the bundled 'zlib/' subdir
%       src_name   : .c filename (e.g. 'read_EDF_mex.c')

zlib_dir = fullfile(script_dir, 'zlib');
mex_src  = fullfile(script_dir, src_name);

if ~isfile(mex_src)
    error('compile_edf_mex:NoSource', 'Source file not found: %s', mex_src);
end

if isfolder(zlib_dir)
    fprintf('Compiling %s with bundled zlib...\n', src_name);
    zlib_files = dir(fullfile(zlib_dir, '*.c'));
    zlib_paths = arrayfun(@(f) fullfile(f.folder, f.name), zlib_files, ...
        'UniformOutput', false);
    args = [{'-O', '-largeArrayDims', ['-I' zlib_dir], '-outdir', script_dir, mex_src}, ...
            zlib_paths(:)'];
    mex(args{:});
else
    fprintf('Compiling %s with system zlib (-lz)...\n', src_name);
    mex('-O', '-largeArrayDims', '-outdir', script_dir, mex_src, '-lz');
end
end
