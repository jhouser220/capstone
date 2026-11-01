#!/usr/bin/env bash
# Joint weight sweep for the two fusion scorers, on the 218-screen `validation`
# FOLD rather than the curated 20s.
#
# Why the fold and not paper_validation. On 20 screens the family-bootstrap CIs are
# ~0.6 wide (see the recorded arms at +0.497 / +0.562 / +0.566) while the differences
# here are plausibly +0.02 to +0.05 -- the curated sets are ~10x too coarse to
# resolve them, and tuning a 20th configuration on them is how +0.562 became +0.175.
#
# The grid is a pre-specified CROSS through the deployed (0.5, 0.5), not a full grid,
# to keep multiplicity honest:
#   knn_weight axis  0.25 / 0.75 / 1.0  at weight 0.5   -- never tested above 0.5,
#                                                          and 1.0 is where the kNN
#                                                          scorer gains the ability to
#                                                          seat its own top pick at the
#                                                          head of a batch (§4.13)
#   weight axis      0.25 / 0.75        at knn_weight 0.5 -- the codep optimum was
#                                                          located with NO kNN present
#
# Already on disk and NOT re-run (all 218 screens each): sweep-dm-validation-w0 (the
# nested baseline) and sweep-dm-validation-dm05knn (the deployed 0.5/0.5 reference).
# Score with scripts/diagnostics/score_weight_sweep.py, which excludes the 20 curated
# screens and reports PAIRED contrasts against the reference.
#
# No API calls: the shipped text-embedding cache covers all 218 validation-fold
# screens (unlike the test fold, §7).
set -u
D=${1:?usage: run_weight_sweep.sh <depmap-dir> [screen-set]}
# DEPMAP ATLAS, 2026-10-28: the ABLATED matrix is now the default. achilles_noeval.npz
# keeps 768 of the 808 Achilles cell lines, dropping every line used by any of the four
# evaluation sets, so no screen's own cell line can be in the atlas by construction.
# Override with ATLAS=achilles.npz to go back to all 808 lines.
ATLAS=${ATLAS:-achilles_noeval.npz}
SET=${2:-validation}
TH=output/analysis/train_hits
CKPT=/home/jeff/.cache/huggingface/hub/models--Genentech--assayformer/snapshots/4927703e5903984ae0254a9fe8e93ec0b29c27e4
LOG=output/analysis/weight_sweep_logs
mkdir -p "$LOG"
# The venv's console script, NOT `uv run`. `uv run` exits 120 with no output
# intermittently in this environment (§8), which inside a long unattended sweep
# looks exactly like a completed arm that produced nothing.
PY=.venv/bin/assayloop
[ -x "$PY" ] || { echo "no $PY -- run 'uv sync --extra torch' first"; exit 1; }

run() { id=$1; shift
  # match any index width: the fold has 218 screens, so `-[0-9][0-9]-` would only
  # catch 00-99 and an interrupted arm would be skipped on the retry (§8)
  if compgen -G "output/runs/${id}-[0-9]*" > /dev/null; then
    n=$(ls -d output/runs/${id}-[0-9]* 2>/dev/null | wc -l)
    echo "SKIP $id ($n run dirs already present)"; return
  fi
  echo "=== $id  $(date +%H:%M:%S)"
  # stderr goes to a LOG, never to /dev/null: 286 of 334 test-fold screens once
  # failed silently on the embedding cache and the sweep looked complete (§8)
  "$PY" run --model amortized_ranker --acq codep_seed \
    --screen-set "$SET" --n-steps 10 --batch-size 100 \
    --full-genome --min-screen-freq 2 \
    --model-param checkpoint=$CKPT --model-param ckpt_file=model_last.pt \
    --model-param device=cuda \
    --acq-param profile_npz="$D/$ATLAS" --acq-param fuse=rrf \
    --acq-param train_hits=$TH --acq-param query=diffmeans "$@" \
    --sweep-id "$id" --tag "$id" > "$LOG/$id.log" 2>&1
  echo "  exit $? | $(ls -d output/runs/${id}-[0-9]* 2>/dev/null | wc -l) run dirs | tail:"
  tail -3 "$LOG/$id.log" | sed 's/^/    /'
}

P=wt-$SET
# knn_weight axis, codep weight held at the deployed 0.5
run "sweep-$P-c05k025" --acq-param weight=0.5  --acq-param knn_weight=0.25
run "sweep-$P-c05k075" --acq-param weight=0.5  --acq-param knn_weight=0.75
run "sweep-$P-c05k10"  --acq-param weight=0.5  --acq-param knn_weight=1.0
# codep weight axis, knn_weight held at the deployed 0.5
run "sweep-$P-c025k05" --acq-param weight=0.25 --acq-param knn_weight=0.5
run "sweep-$P-c075k05" --acq-param weight=0.75 --acq-param knn_weight=0.5
echo "ALLDONE $(date +%H:%M:%S)"
