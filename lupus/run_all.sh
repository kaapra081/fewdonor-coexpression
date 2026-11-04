#!/bin/bash
# PLAN7 runner inside the Modal sandbox. usage: bash run_all.sh <jobs file> <parallel>
set -u
mkdir -p out code; [ -d code/repo ] || tar xzf code.tgz -C code
export LUPUS_REPO=$PWD/code/repo LUPUS_DATA=$PWD LUPUS_OUT=$PWD/out OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
snap() { tar czf out/_snap.tgz -C out $(cd out && ls *.json *.npz 2>/dev/null) 2>/dev/null; }
trap 'snap' EXIT TERM
( while true; do sleep 300; snap; done ) &
cat "$1" | xargs -P "$2" -L 1 python code/lupus/lupus_run.py
snap
