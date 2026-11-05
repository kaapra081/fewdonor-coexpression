#!/bin/bash
# Plan 9 runner inside the Modal sandbox: bash run9.sh
set -u
mkdir -p out data code; tar xzf code9.tgz -C code
export OUT=$PWD/out BDATA=$PWD/data REPO=$PWD/code/repo OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
snap() { (cd out && tar czf _snap.tgz *.json *.csv *.log 2>/dev/null); }
trap 'snap' EXIT TERM
( while true; do sleep 240; snap; done ) &
pip install -q "cellxgene-census==1.18.0" "numpy==1.26.4" > out/pip.log 2>&1; echo "pip exit $?"
python -c "import cellxgene_census, numpy, scipy, tiledbsoma; print('census', cellxgene_census.__version__, 'numpy', numpy.__version__, 'scipy', scipy.__version__, 'soma', tiledbsoma.__version__)"
date -u; python code/repo/breadth/breadth_job.py select > out/select.log 2>&1; echo "select exit $?"; tail -3 out/select.log; date -u
NS=$(python -c "import json; print(len(json.load(open('out/selection.json'))['settings']))"); echo "settings $NS"
seq 1 $NS | xargs -P 8 -I{} sh -c 'python code/repo/breadth/breadth_job.py extract {} > out/extract_{}.log 2>&1; echo "extract {} exit $?"'
date -u
for i in $(seq 1 $NS); do for k in 0 5 10 15; do [ -f data/S${i}_X.npz ] && echo "null $i $k $((k+5))"; done; done > jobs9.txt; wc -l jobs9.txt
cat jobs9.txt | xargs -P 16 -L 1 sh -c 'python code/repo/breadth/breadth_job.py "$0" "$1" "$2" "$3" >> out/null.log 2>&1'
date -u; snap; ls out | grep -c "^b_S"
