#!/usr/bin/env bash
# Real 10-round campaigns with AssayFormer fused with DepMap co-essentiality.
# No API calls. weight=0 reproduces plain greedy, so the sweep contains its own
# nested baseline.
set -u
NPZ=${1:?usage: run_codep_campaigns.sh <path to achilles.npz>}
CKPT=/home/jeff/.cache/huggingface/hub/models--Genentech--assayformer/snapshots/4927703e5903984ae0254a9fe8e93ec0b29c27e4

run() { id=$1; set_=$2; shift 2
  if [ -d "output/runs/${id}-00-"* ] 2>/dev/null; then echo "SKIP $id"; return; fi
  echo "=== $id ($set_)"
  uv run assayloop run --model amortized_ranker --acq codep_seed \
    --screen-set "$set_" --n-steps 10 --batch-size 100 \
    --full-genome --min-screen-freq 2 \
    --model-param checkpoint=$CKPT --model-param ckpt_file=model_last.pt \
    --model-param device=cuda \
    --acq-param profile_npz="$NPZ" "$@" \
    --sweep-id "$id" --tag "$id" >/dev/null 2>&1 && echo "  ok" || echo "  FAILED"
}

for S in paper_validation paper_test; do
  case $S in paper_validation) P=cdv;; *) P=cdt;; esac
  run "sweep-$P-w0"      "$S" --acq-param weight=0.0                          # nested baseline
  run "sweep-$P-w05"     "$S" --acq-param weight=0.5
  run "sweep-$P-w1"      "$S" --acq-param weight=1.0
  run "sweep-$P-w2"      "$S" --acq-param weight=2.0
  run "sweep-$P-w1-zsum" "$S" --acq-param weight=1.0 --acq-param fuse=zsum
done
echo ALLDONE
