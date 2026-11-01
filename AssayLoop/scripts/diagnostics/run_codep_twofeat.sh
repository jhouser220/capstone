#!/usr/bin/env bash
# Campaigns with BOTH scorers: AssayFormer + DepMap co-essentiality + screen-kNN.
# Offline the two features are near-perfectly additive (combined oracle margin
# +0.260 validation / +0.733 test, about the sum of the parts).
# No API calls. Baselines are the existing weight-0 arms.
set -u
D=${1:?usage: run_codep_twofeat.sh <depmap-dir> [screen-set]}
# DEPMAP ATLAS, 2026-10-28: the ABLATED matrix is now the default. achilles_noeval.npz
# keeps 768 of the 808 Achilles cell lines, dropping every line used by any of the four
# evaluation sets, so no screen's own cell line can be in the atlas by construction.
# Override with ATLAS=achilles.npz to go back to all 808 lines.
ATLAS=${ATLAS:-achilles_noeval.npz}
SET=${2:-paper_test}
TH=output/analysis/train_hits          # from build_train_hits.py
CKPT=/home/jeff/.cache/huggingface/hub/models--Genentech--assayformer/snapshots/4927703e5903984ae0254a9fe8e93ec0b29c27e4

run() { id=$1; shift
  if compgen -G "output/runs/${id}-[0-9][0-9]-*" > /dev/null; then echo "SKIP $id (exists)"; return; fi
  echo "=== $id"
  uv run assayloop run --model amortized_ranker --acq codep_seed \
    --screen-set "$SET" --n-steps 10 --batch-size 100 \
    --full-genome --min-screen-freq 2 \
    --model-param checkpoint=$CKPT --model-param ckpt_file=model_last.pt \
    --model-param device=cuda \
    --acq-param profile_npz="$D/$ATLAS" --acq-param fuse=rrf \
    --acq-param train_hits=$TH "$@" \
    --sweep-id "$id" --tag "$id" 2>&1 | tail -2
}

case $SET in paper_test) P=tf-pt;; paper_validation) P=tf-pv;; *) P="tf-$SET";; esac
run "sweep-$P-w0"          --acq-param weight=0.0  --acq-param knn_weight=0.0   # nested baseline
run "sweep-$P-knn05"       --acq-param weight=0.0  --acq-param knn_weight=0.5   # kNN alone
run "sweep-$P-both05"      --acq-param weight=0.5  --acq-param knn_weight=0.5   # both
run "sweep-$P-both-c5k25"  --acq-param weight=0.5  --acq-param knn_weight=0.25
echo ALLDONE
