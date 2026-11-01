#!/usr/bin/env bash
# ATLAS GEOMETRY (§12 item 3 / §4.44): PC removal and amplitude weighting on the
# dependency profiles, on the validation FOLD.
#
# The idea. Every one of §4's eleven feature candidates varied the DATA SOURCE. The
# geometry of the space the kernel measures similarity in was never varied: cosine on
# raw 808-dim profiles. Difference-of-means cancels shared factors in the DIRECTION,
# which is part of why it beats the centroid, but nothing cancels them in the SPACE,
# and those are two separate corrections. Removing leading components is standard in
# the co-essentiality literature.
#
# BUT MEASURE THE PREMISE FIRST, and it is weak here. The leading component of the
# per-gene-centred Achilles matrix carries only 2.2% of the variance; top 3 is 5.9%,
# top 20 is 19.9%, and the mean random-pair cosine barely moves (+0.0024 -> +0.0017)
# after dropping 3. There is no dominant shared factor to remove, almost certainly
# because the Chronos release already models screen quality, growth rate and library
# effects and ships a corrected gene-effect matrix -- the literature advice applies to
# less-processed data. So the DOSES here are set by the variance curve rather than by
# the usual 1/3/5: pc20 is the first dose that removes a fifth of the variance.
#
# Amplitude has little to work with either: median-scaled amplitude runs p5 0.79 to
# p95 1.74, a narrow distribution. The reason to try it anyway is §4.25, which measured
# the hard screens' hits at amplitude 0.138 against 0.137 for random genes -- a flat
# profile's cosine is noise whatever its value, and those genes still compete at the
# top-100 cut. NOTE this knob cuts AGAINST the hard screens by construction, so read
# the difficulty split, not only the overall figure.
#
# FIVE PRE-SPECIFIED ARMS, both axes nested at their off value (drop_pcs=0,
# amp_power=0 is bit-identical to the deployed configuration, test-pinned):
#   pc3 / pc10 / pc20   drop_pcs, at the deployed 0.5/0.5
#   amp05 / amp10       amp_power, at the deployed 0.5/0.5
#
# Baseline and reference already on disk, not re-run: sweep-dm-validation-w0 and
# sweep-dm-validation-dm05knn. No API calls.
set -u
D=${1:?usage: run_geometry_sweep.sh <depmap-dir> [screen-set]}
# DEPMAP ATLAS, 2026-10-28: the ABLATED matrix is now the default. achilles_noeval.npz
# keeps 768 of the 808 Achilles cell lines, dropping every line used by any of the four
# evaluation sets, so no screen's own cell line can be in the atlas by construction.
# Override with ATLAS=achilles.npz to go back to all 808 lines.
ATLAS=${ATLAS:-achilles_noeval.npz}
SET=${2:-validation}
TH=output/analysis/train_hits
CKPT=/home/jeff/.cache/huggingface/hub/models--Genentech--assayformer/snapshots/4927703e5903984ae0254a9fe8e93ec0b29c27e4
LOG=output/analysis/geometry_sweep_logs
mkdir -p "$LOG"
# device defaults to cpu: GPU0 is in a driver error state and `device=cuda` fails
# every screen while the runner still exits 0 (§8, §4.43). DEV=cuda to override.
PY=.venv/bin/assayloop
[ -x "$PY" ] || { echo "no $PY -- run 'uv sync --extra torch' first"; exit 1; }

run() { id=$1; shift
  if compgen -G "output/runs/${id}-[0-9]*" > /dev/null; then
    echo "SKIP $id ($(ls -d output/runs/${id}-[0-9]* | wc -l) run dirs present)"; return
  fi
  echo "=== $id  $(date +%H:%M:%S)"
  "$PY" run --model amortized_ranker --acq codep_seed \
    --screen-set "$SET" --n-steps 10 --batch-size 100 \
    --full-genome --min-screen-freq 2 \
    --model-param checkpoint=$CKPT --model-param ckpt_file=model_last.pt \
    --model-param device=${DEV:-cpu} \
    --acq-param profile_npz="$D/$ATLAS" --acq-param fuse=rrf \
    --acq-param train_hits=$TH --acq-param query=diffmeans \
    --acq-param weight=0.5 --acq-param knn_weight=0.5 "$@" \
    --sweep-id "$id" --tag "$id" > "$LOG/$id.log" 2>&1
  echo "  exit $? | $(ls -d output/runs/${id}-[0-9]* 2>/dev/null | wc -l) run dirs"
  grep -iE "error|traceback" "$LOG/$id.log" | head -2 | sed 's/^/    /'
}

P=geom-$SET
run "sweep-$P-pc3"   --acq-param drop_pcs=3
run "sweep-$P-pc10"  --acq-param drop_pcs=10
run "sweep-$P-pc20"  --acq-param drop_pcs=20
run "sweep-$P-amp05" --acq-param amp_power=0.5
run "sweep-$P-amp10" --acq-param amp_power=1.0
echo "ALLDONE $(date +%H:%M:%S)"
