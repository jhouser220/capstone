#!/usr/bin/env bash
# Campaigns with the difference-of-means query direction (the 2026-10-22 default).
# Offline it roughly doubles the hit-centroid's oracle margin on both screen sets.
# Baselines are weight-0 arms, which reproduce AssayFormer pick for pick.
set -u
D=${1:?usage: run_codep_diffmeans.sh <depmap-dir> [screen-set]}
# DEPMAP ATLAS, 2026-10-28: the ABLATED matrix is now the default. achilles_noeval.npz
# keeps 768 of the 808 Achilles cell lines, dropping every line used by any of the four
# evaluation sets, so no screen's own cell line can be in the atlas by construction.
# Override with ATLAS=achilles.npz to go back to all 808 lines.
ATLAS=${ATLAS:-achilles_noeval.npz}
SET=${2:-paper_test}
TH=output/analysis/train_hits
CKPT=/home/jeff/.cache/huggingface/hub/models--Genentech--assayformer/snapshots/4927703e5903984ae0254a9fe8e93ec0b29c27e4

run() { id=$1; shift
  if compgen -G "output/runs/${id}-[0-9][0-9]-*" > /dev/null; then echo "SKIP $id"; return; fi
  echo "=== $id"
  uv run assayloop run --model amortized_ranker --acq codep_seed \
    --screen-set "$SET" --n-steps 10 --batch-size 100 \
    --full-genome --min-screen-freq 2 \
    --model-param checkpoint=$CKPT --model-param ckpt_file=model_last.pt \
    --model-param device=cuda \
    --acq-param profile_npz="$D/$ATLAS" --acq-param fuse=rrf \
    --acq-param train_hits=$TH "$@" \
    --sweep-id "$id" --tag "$id" 2>&1 | tail -1
}

case $SET in paper_test) P=dm-pt;; paper_validation) P=dm-pv;; *) P="dm-$SET";; esac
run "sweep-$P-w0"      --acq-param weight=0.0 --acq-param knn_weight=0.0
run "sweep-$P-dm05"    --acq-param weight=0.5 --acq-param knn_weight=0.0  --acq-param query=diffmeans
run "sweep-$P-dm05knn" --acq-param weight=0.5 --acq-param knn_weight=0.5  --acq-param query=diffmeans
run "sweep-$P-dm10"    --acq-param weight=1.0 --acq-param knn_weight=0.0  --acq-param query=diffmeans
echo ALLDONE
