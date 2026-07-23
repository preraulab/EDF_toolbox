# Bundled Zstandard (zstd) v1.5.7

Single-file amalgamation of the Zstandard compression library, vendored so
the MEX backends build with **no system libzstd required** (mirrors the
bundled `zlib/` directory). Dual-licensed BSD-3-Clause / GPLv2 — see
`LICENSE` (BSD text; both licenses permit this use).

Contents:

- `zstd.c` — full library amalgamation, generated from the v1.5.7 release
  with `build/single_file_libs/create_single_file_library.sh`
  (equivalently: `python combine.py -r ../../lib -x legacy/zstd_legacy.h
  -o zstd.c zstd-in.c`). Assembly is disabled in the amalgamation
  (`ZSTD_DISABLE_ASM`), so it compiles as plain C on any platform/compiler
  MATLAB's `mex` supports, including MSVC and MinGW on Windows.
  `ZSTD_MULTITHREAD` is baked in: uses native Win32 threads on Windows and
  pthreads elsewhere (compile_edf_mex.m adds `-lpthread` on unix).
- `zstd.h`, `zstd_errors.h` — public API headers (`lib/zstd.h` and
  `lib/zstd_errors.h` from the same release; `zstd.h` includes the latter).
- `LICENSE` — upstream BSD license.

To upgrade: download a newer release from
https://github.com/facebook/zstd/releases, run the script above, and
replace `zstd.c` / `zstd.h` here. Do not edit these files by hand.
