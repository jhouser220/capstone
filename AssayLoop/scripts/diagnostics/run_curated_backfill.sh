#!/usr/bin/env bash
# BACKFILL this session's arms onto the two CURATED 20-screen sets, so the notebook can
# compare the same methods across paper_validation / paper_test / validation / test.
#
# Those two sets are the ones the PAPER reports, so they are the only place a direct
# comparison to its claimed EF 5.67 is possible -- but they are also the sets that chose
# this project's settings, so their numbers run ~3x the held-out figures (§1). Treat them
# as the paper comparison, never as evidence.
#
# A gap worth knowing: paper_validation had NO deployed codep+kNN campaign at all, only
# codep-alone weight arms (sweep-cdv-*). That is arm 1 below.
set -u
D=${1:-/home/jeff/depmap}
# DEPMAP ATLAS, 2026-10-28: the ABLATED matrix is now the default. achilles_noeval.npz
# keeps 768 of the 808 Achilles cell lines, dropping every line used by any of the four
# evaluation sets, so no screen's own cell line can be in the atlas by construction.
# Override with ATLAS=achilles.npz to go back to all 808 lines.
ATLAS=${ATLAS:-achilles_noeval.npz}
TH=output/analysis/train_hits
PPI=output/analysis/bioplex.npz
CKPT=/home/jeff/.cache/huggingface/hub/models--Genentech--assayformer/snapshots/4927703e5903984ae0254a9fe8e93ec0b29c27e4
LOG=output/analysis/backfill_logs
mkdir -p "$LOG"
PY=.venv/bin/assayloop

one() { id=$1; set_=$2; shift 2
  if compgen -G "output/runs/${id}-[0-9]*" > /dev/null; then
    echo "SKIP $id"; return; fi
  echo "=== $id ($set_)  $(date +%H:%M:%S)"
  "$PY" run --model amortized_ranker --acq codep_seed \
    --screen-set "$set_" --n-steps 10 --batch-size 100 \
    --full-genome --min-screen-freq 2 \
    --model-param checkpoint=$CKPT --model-param ckpt_file=model_last.pt \
    --model-param device=${DEV:-cpu} \
    --acq-param profile_npz="$D/$ATLAS" --acq-param fuse=rrf \
    --acq-param train_hits=$TH --acq-param query=diffmeans "$@" \
    --sweep-id "$id" --tag "$id" > "$LOG/$id.log" 2>&1
  echo "  exit $? | $(ls -d output/runs/${id}-[0-9]* 2>/dev/null | wc -l) dirs"
}

# 1. the missing deployed arm on paper_validation
one sweep-pv-dm05knn paper_validation --acq-param weight=0.5 --acq-param knn_weight=0.5

# 2. BioPlex, both curated sets
for S in paper_validation paper_test; do
  P=$([ "$S" = paper_validation ] && echo pv || echo pt)
  one "sweep-ppi-$P-w025"   "$S" --acq-param weight=0.5 --acq-param knn_weight=0.5 \
      --acq-param ppi_npz=$PPI --acq-param ppi_weight=0.25
  one "sweep-ppi-$P-w025dr" "$S" --acq-param weight=0.5 --acq-param knn_weight=0.5 \
      --acq-param ppi_npz=$PPI --acq-param ppi_weight=0.25 \
      --acq-param ppi_degree_residual=true
done

# 3. description variants, both curated sets, kernel OFF to isolate the text channel
for S in paper_validation paper_test; do
  P=$([ "$S" = paper_validation ] && echo pv || echo pt)
  for V in notes dir all; do
    if compgen -G "output/runs/sweep-desc-$P-$V-[0-9]*" > /dev/null; then
      echo "SKIP sweep-desc-$P-$V"; continue; fi
    echo "=== sweep-desc-$P-$V ($S, ASSAYLOOP_DESC_EXTRA=$V)"
    ASSAYLOOP_DESC_EXTRA=$V "$PY" run --model amortized_ranker --acq codep_seed \
      --screen-set "$S" --n-steps 10 --batch-size 100 \
      --full-genome --min-screen-freq 2 \
      --model-param checkpoint=$CKPT --model-param ckpt_file=model_last.pt \
      --model-param device=${DEV:-cpu} \
      --acq-param profile_npz="$D/$ATLAS" --acq-param fuse=rrf \
      --acq-param train_hits=$TH --acq-param query=diffmeans \
      --acq-param weight=0 --acq-param knn_weight=0 \
      --sweep-id "sweep-desc-$P-$V" --tag "sweep-desc-$P-$V" \
      > "$LOG/sweep-desc-$P-$V.log" 2>&1
    echo "  exit $? | $(ls -d output/runs/sweep-desc-$P-$V-[0-9]* 2>/dev/null | wc -l) dirs"
  done
done
echo "ALLDONE-PART1 $(date +%H:%M:%S)"
