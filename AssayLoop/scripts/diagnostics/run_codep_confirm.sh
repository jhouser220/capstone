#!/usr/bin/env bash
# Pre-specified confirmation of codep_seed, plus the expression-profile variant.
#
# Weight 0.5 + RRF was chosen on paper_validation/paper_test, so those 20-screen
# curated sets are NOT held out. This runs the whole `test` fold (334 screens) and
# scores only the 314 that are not in paper_test. No API calls; ~7 min per arm on
# the laptop GPU.
#
#   ./scripts/diagnostics/run_codep_confirm.sh <depmap-dir> [screen-set]
set -u
D=${1:?usage: run_codep_confirm.sh <depmap-dir> [screen-set]}
# DEPMAP ATLAS, 2026-10-28: the ABLATED matrix is now the default. achilles_noeval.npz
# keeps 768 of the 808 Achilles cell lines, dropping every line used by any of the four
# evaluation sets, so no screen's own cell line can be in the atlas by construction.
# Override with ATLAS=achilles.npz to go back to all 808 lines.
ATLAS=${ATLAS:-achilles_noeval.npz}
SET=${2:-test}
CKPT=/home/jeff/.cache/huggingface/hub/models--Genentech--assayformer/snapshots/4927703e5903984ae0254a9fe8e93ec0b29c27e4

run() { id=$1; shift
  if compgen -G "output/runs/${id}-0*" > /dev/null; then echo "SKIP $id (exists)"; return; fi
  echo "=== $id"
  uv run assayloop run --model amortized_ranker --acq codep_seed \
    --screen-set "$SET" --n-steps 10 --batch-size 100 \
    --full-genome --min-screen-freq 2 \
    --model-param checkpoint=$CKPT --model-param ckpt_file=model_last.pt \
    --model-param device=cuda \
    --acq-param fuse=rrf "$@" \
    --sweep-id "$id" --tag "$id" >/dev/null 2>&1 && echo "  ok" || echo "  FAILED"
}

# weight 0 is the nested baseline and is identical whichever profile is loaded
run "sweep-cf-$SET-w0"        --acq-param profile_npz="$D/$ATLAS"   --acq-param weight=0.0
run "sweep-cf-$SET-ess-w05"   --acq-param profile_npz="$D/$ATLAS"   --acq-param weight=0.5
run "sweep-cf-$SET-expr-w05"  --acq-param profile_npz="$D/Expression.npz" --acq-param weight=0.5
echo "ALLDONE"
