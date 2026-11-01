#!/usr/bin/env bash
# Score the confirmation arms, excluding the curated screens used to pick the weight.
set -u
SET=${1:-test}
case $SET in test) EX=paper_test;; validation) EX=paper_validation;; *) EX="";; esac
for ARM in ess expr; do
  echo "######## $ARM profile, $SET fold ########"
  uv run python scripts/diagnostics/codep_gated.py --screen-set "$SET" \
    --off "sweep-cf-$SET-w0" --on "sweep-cf-$SET-$ARM-w05" \
    ${EX:+--exclude-set $EX} --quiet 2>&1 | grep -v "^INFO"
done
