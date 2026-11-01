#!/usr/bin/env bash
# SECOND HELD-OUT CONFIRMATION of the deployed method, on the `test` fold (§10 item 1).
#
# The headline +0.175 EF per screen / +0.331 per family is confirmed on the 198 fresh
# `validation`-fold screens. This repeats that measurement on 314 fresh `test`-fold
# screens, which chose nothing: no weight, no parameter, no gate was selected on them.
# Two independent confirmations instead of one.
#
# PREREQUISITE: every test-fold description must be in the embedding cache.
#   uv run python scripts/diagnostics/embed_screen_set.py --screen-set test
# 48 of 334 ship pre-embedded; the rest cost well under a cent, one-time. Without them
# each uncached screen dies with `RuntimeError: Text not found in the embedding cache`.
#
# TWO ARMS, the same nested pair the validation confirmation used, so the contrast is
# like-for-like:
#   w0        weight=0 knn_weight=0 -- reproduces AssayFormer pick for pick
#   dm05knn   the deployed weight=0.5 knn_weight=0.5
#
# NEW sweep ids (`sweep-tf-test-*`), deliberately NOT the existing `sweep-cf-test-w0`:
# that one holds 48 run dirs from the earlier run that silently failed on 286 uncached
# screens (§8), and mixing 48 old dirs with 286 new ones would make provenance
# unverifiable while the launcher's guard would SKIP the arm entirely.
#
# Score with:
#   score_weight_sweep.py --screen-set test --exclude-set paper_test \
#     --baseline sweep-tf-test-w0 --reference sweep-tf-test-dm05knn
#   score_resolution.py  --screen-set test --exclude-set paper_test ... (per-round + ANCOVA)
# Quote the fold figure, and the genome-scale subset (--min-library 1001) separately:
# 98 of the validation fold's screens were metric-degenerate and the test fold has its
# own share (§4.29).
set -u
D=${1:?usage: run_testfold_confirm.sh <depmap-dir>}
# DEPMAP ATLAS, 2026-10-28: the ABLATED matrix is now the default. achilles_noeval.npz
# keeps 768 of the 808 Achilles cell lines, dropping every line used by any of the four
# evaluation sets, so no screen's own cell line can be in the atlas by construction.
# Override with ATLAS=achilles.npz to go back to all 808 lines.
ATLAS=${ATLAS:-achilles_noeval.npz}
TH=output/analysis/train_hits
CKPT=/home/jeff/.cache/huggingface/hub/models--Genentech--assayformer/snapshots/4927703e5903984ae0254a9fe8e93ec0b29c27e4
LOG=output/analysis/testfold_logs
mkdir -p "$LOG"
# device defaults to cpu: GPU0 is in a driver error state and device=cuda fails every
# screen while the runner still exits 0 (§4.43). DEV=cuda once the GPU is back.
PY=.venv/bin/assayloop
[ -x "$PY" ] || { echo "no $PY -- run 'uv sync --extra torch' first"; exit 1; }

# Fail fast rather than burning an hour producing a half-fold: the embedding cache must
# already cover the set. This is the §8 lesson -- 286 of 334 screens once failed
# silently and the sweep looked complete.
.venv/bin/python scripts/diagnostics/embed_screen_set.py --screen-set test --dry-run \
  | tee "$LOG/precheck.log"
grep -q "nothing to do" "$LOG/precheck.log" || {
  echo "ABORT: test-fold descriptions are not all cached. Run:"
  echo "  .venv/bin/python scripts/diagnostics/embed_screen_set.py --screen-set test"
  exit 1; }

run() { id=$1; shift
  if compgen -G "output/runs/${id}-[0-9]*" > /dev/null; then
    echo "SKIP $id ($(ls -d output/runs/${id}-[0-9]* | wc -l) run dirs present)"; return
  fi
  echo "=== $id  $(date +%H:%M:%S)"
  "$PY" run --model amortized_ranker --acq codep_seed \
    --screen-set test --n-steps 10 --batch-size 100 \
    --full-genome --min-screen-freq 2 \
    --model-param checkpoint=$CKPT --model-param ckpt_file=model_last.pt \
    --model-param device=${DEV:-cpu} \
    --acq-param profile_npz="$D/$ATLAS" --acq-param fuse=rrf \
    --acq-param train_hits=$TH --acq-param query=diffmeans "$@" \
    --sweep-id "$id" --tag "$id" > "$LOG/$id.log" 2>&1
  echo "  exit $? | $(ls -d output/runs/${id}-[0-9]* 2>/dev/null | wc -l) run dirs"
  grep -icE "text not found|error|traceback" "$LOG/$id.log" \
    | xargs -I{} echo "    {} error lines in the log -- inspect before scoring"
}

run "sweep-tf-test-w0"      --acq-param weight=0   --acq-param knn_weight=0
run "sweep-tf-test-dm05knn" --acq-param weight=0.5 --acq-param knn_weight=0.5
echo "ALLDONE $(date +%H:%M:%S)"
