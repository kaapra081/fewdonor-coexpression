#!/bin/bash
# usage: bash run_all.sh  (from recomb_v2/). Runs every job in jobs.txt, nproc at a time, single-threaded BLAS.
cd "$(dirname "$0")"; mkdir -p logs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
N=${NPAR:-$(nproc)}
cat jobs.txt | xargs -P "$N" -I{} bash -c '{}; echo "$? {}" >> logs/_done.txt'
echo ALLDONE >> logs/_done.txt
