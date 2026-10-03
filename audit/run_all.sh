#!/bin/sh
# Re-run every audit script; outputs land in audit/outputs/. Run from audit/.
#   cd audit && sh run_all.sh
# a2 (~4 min) and a5 (~15 min) are the slow ones.
PY="${PY:-python}"   # e.g. PY="sh /mnt/c/Metabolomics_project/scripts/py" on the author machine
for s in a1_pxd001064 a2_pxd004242 a3_tool a4_leak a5_margin a6_robust a7_fig6b a8_ci a9_claims56; do
  echo "== $s"
  if [ "$s" = a6_robust ]; then
    $PY $s.py "$PWD/outputs/cli_out" > outputs/$s.txt 2>&1
  else
    $PY $s.py > outputs/$s.txt 2>&1
  fi
done
