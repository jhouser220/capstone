#!/usr/bin/env bash
# Profile-block screens on both screen sets, with the deployable column.
set -u
D=${1:?usage: run_profile_screens.sh <depmap-dir>}
echo "######## paper_validation ########"
uv run python scripts/diagnostics/screen_profiles.py --depmap-dir "$D" \
  --screen-set paper_validation --base-sweep sweep-fd2017be \
  --base-dir output/runs/baseline/gemini --k 50 --per-screen 2>&1 | grep -v "^INFO"
echo "######## paper_test ########"
uv run python scripts/diagnostics/screen_profiles.py --depmap-dir "$D" \
  --screen-set paper_test --base-sweep sweep-00866f93 \
  --base-dir output/runs/baseline/gemini --k 50 --per-screen 2>&1 | grep -v "^INFO"
