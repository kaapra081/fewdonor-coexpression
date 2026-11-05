#!/bin/bash
# Plan 10 runner: bash run10.sh <ORDER fwd|rev> <WORKERS> <CENSUS_PYTHON_SETUP modal|colab>
# Each worker takes one setting: extract it, then run its 20 null reassignments as 4 parallel chunks of 5.
set -u
ORDER=$1; WORKERS=$2; HOSTKIND=$3
mkdir -p out data code; tar xzf code10.tgz -C code
export PLAN=10 MIN_DONORS=100 PER_COLL=3 MAX_SET=24 CAP=0 SEED0=93000000
export OUT=$PWD/out BDATA=$PWD/data REPO=$PWD/code/repo OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
snap() { (cd out && tar czf _snap.tgz *.json *.csv *.log 2>/dev/null); }
trap 'snap' EXIT TERM
( while true; do sleep 180; snap; done ) &
python3 -m venv $PWD/cenv && ($PWD/cenv/bin/pip install -q "cellxgene-census==1.18.0" scipy > out/pip.log 2>&1 || $PWD/cenv/bin/pip install -q cellxgene-census scipy >> out/pip.log 2>&1); echo "pip exit $?"
export CPY=$PWD/cenv/bin/python APY=$(command -v python3)
$CPY -c "import cellxgene_census, numpy, tiledbsoma; print('census', cellxgene_census.__version__, 'numpy', numpy.__version__, 'soma', tiledbsoma.__version__)"
$APY -c "import numpy, scipy; print('analysis numpy', numpy.__version__, 'scipy', scipy.__version__)"
date -u; $CPY code/repo/breadth/breadth_job.py select > out/select.log 2>&1; echo "select exit $?"; date -u
NS=$($CPY -c "import json; print(len(json.load(open('out/selection.json'))['settings']))"); echo "settings $NS"
if [ "$ORDER" = rev ]; then seq $NS -1 1; else seq 1 $NS; fi > order.txt
do_setting() {
  i=$1
  if [ ! -f data/S${i}_X.npz ]; then $CPY code/repo/breadth/breadth_job.py extract $i > out/extract_$i.log 2>&1 || { echo "extract $i failed"; return 0; }; fi
  for k in 0 5 10 15; do $APY code/repo/breadth/breadth_job.py null $i $k $((k+5)) >> out/null_$i.log 2>&1 & done; wait
  echo "setting $i done $(date -u +%H:%M:%S)"
}
export -f do_setting
cat order.txt | xargs -P $WORKERS -I{} bash -c 'do_setting {}'
date -u; snap; ls out | grep -c "^b_S"
