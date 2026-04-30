#!/usr/bin/env bash
# Quick comparison: serial-no-batch vs serial-batched vs parallel (W=4 WT=2).
#
# Run from the EDF_toolbox directory (script is in the toolbox root):
#   ./compare_modes.sh [INDIR]
# Defaults: INDIR=/scratch/test_raw_edf
# Pulls latest from origin first, then runs the three configs back-to-back.

set -uo pipefail

INDIR="${1:-/scratch/test_raw_edf}"
OUT_DIR="/scratch/edf_cmp_out"

echo "Pulling latest..."
git pull --ff-only

echo
echo "Toolbox  : $(pwd)"
echo "Input    : $INDIR"
echo "Out dir  : $OUT_DIR"
echo

matlab -batch "addpath(genpath(pwd)); \
in_dir='$INDIR'; out_dir='$OUT_DIR'; \
if isfolder(out_dir), rmdir(out_dir,'s'); end; \
fprintf('\n--- serial, double precision ---\n'); \
delete(gcp('nocreate')); mkdir(out_dir); \
t=tic; r=batch_convert_EDF(in_dir,100,'Parallel',false,'OutputDir',out_dir); \
fprintf('serial double: %.1fs ok=%d/%d\n',toc(t),sum(strcmp(r.status,'ok')),height(r)); \
rmdir(out_dir,'s'); \
fprintf('\n--- serial, single precision ---\n'); \
mkdir(out_dir); \
t=tic; r=batch_convert_EDF(in_dir,100,'Parallel',false,'OutputDir',out_dir,'Precision','single'); \
fprintf('serial single: %.1fs ok=%d/%d\n',toc(t),sum(strcmp(r.status,'ok')),height(r)); \
rmdir(out_dir,'s'); \
fprintf('\n--- W=4 WT=2, double precision ---\n'); \
parpool('Processes',4); mkdir(out_dir); \
t=tic; r=batch_convert_EDF(in_dir,100,'Parallel',true,'WorkerThreads',2,'OutputDir',out_dir); \
fprintf('W4 WT2 double: %.1fs ok=%d/%d\n',toc(t),sum(strcmp(r.status,'ok')),height(r)); \
rmdir(out_dir,'s'); \
fprintf('\n--- W=4 WT=2, single precision ---\n'); \
mkdir(out_dir); \
t=tic; r=batch_convert_EDF(in_dir,100,'Parallel',true,'WorkerThreads',2,'OutputDir',out_dir,'Precision','single'); \
fprintf('W4 WT2 single: %.1fs ok=%d/%d\n',toc(t),sum(strcmp(r.status,'ok')),height(r)); \
delete(gcp('nocreate'))"

echo
echo "DONE"
