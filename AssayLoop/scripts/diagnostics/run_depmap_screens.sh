#!/usr/bin/env bash
# Screen the DepMap features on both screen sets. Needs achilles/Expression/Model
# in --depmap-dir; run depmap_to_npz.py first or the CSV parse dominates.
set -u
D=${1:?usage: run_depmap_screens.sh <depmap-dir>}
for S in paper_validation paper_test; do
  echo "=== $S"
  uv run python -I scripts/diagnostics/screen_depmap.py --depmap-dir "$D" --screen-set "$S" 2>&1 \
    | grep -v "^INFO" | tail -13
done
