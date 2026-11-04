#!/bin/bash
# Plan 8 runner inside the Modal sandbox: bash run8.sh
set -u
mkdir -p out code; tar xzf code8.tgz -C code
export REV5_OUT=$PWD/out LUPUS_DATA=$PWD LUPUS_OUT=$PWD/out/_lupus_unused OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
snap() { (cd out && tar czf _snap.tgz f_*.json 2>/dev/null); }
trap 'snap' EXIT TERM
( while true; do sleep 240; snap; done ) &
python -c "import numpy, scipy; print('numpy', numpy.__version__, 'scipy', scipy.__version__)"
cat jobs8.txt | xargs -P 15 -L 1 python code/repo/rev5/rev5_f.py > out/run8.log 2>&1
echo "xargs exit $?"; snap; ls out | grep -c "^f_"
