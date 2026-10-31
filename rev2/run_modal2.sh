#!/bin/bash
set -u
W=$(pwd); mkdir -p out /vol/res/rev2_out /vol/res/onek_out /vol/res/logs
if [ -f kb.tgz ]; then cp -f kb.tgz /vol/kb.tgz; fi
tar xzf /vol/kb.tgz && cd recomb_v2
# seed already-finished results so resumable scripts skip them
cp -n /vol/res/rev2_out/* rev2/out/ 2>/dev/null; cp -rn /vol/res/rev2_out/. rev2/out/ 2>/dev/null; cp -rn /vol/res/onek_out/. onek/out/ 2>/dev/null
rm -rf logs; ln -s /vol/res/logs logs
( while true; do date >> /vol/res/logs/mem.log; free -m | sed -n 2p >> /vol/res/logs/mem.log; ps -eo rss,args --sort=-rss | head -6 >> /vol/res/logs/mem.log; sleep 30; done ) &
( while true; do sleep 60; cp -ru rev2/out/. /vol/res/rev2_out/ 2>/dev/null; cp -ru onek/out/. /vol/res/onek_out/ 2>/dev/null; done ) &
nproc; free -g | head -2
NPAR=${NPAR:-16} bash run_all.sh
cp -ru rev2/out/. /vol/res/rev2_out/; cp -ru onek/out/. /vol/res/onek_out/
sync; sleep 10
tar czf "$W/out/results.tgz" -C /vol/res .
tail -3 logs/_done.txt; grep -c "^0 " logs/_done.txt; grep -v "^0 " logs/_done.txt | head
