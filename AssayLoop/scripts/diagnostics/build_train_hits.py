"""Cache the raw training hit matrix: 1349 screens x genes, as a sparse .npz.

This is the information AssayFormer's gene factors are a 10-dimensional
compression OF. The ceiling analysis (HANDOFF §4.9) says only 3.5% of a screen's
hit structure survives that compression, so the raw matrix is the obvious place
to look for signal outside span(V).

Using training labels at inference is not leakage: they are training data, and the
screen being evaluated is never in this matrix (the folds are disjoint).

    uv run python scripts/diagnostics/build_train_hits.py <out.npz>
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix, save_npz

sys.path.insert(0, "src")
from assayloop.tasks import load_screens   # noqa: E402

out = Path(sys.argv[1] if len(sys.argv) > 1 else "output/analysis/train_hits.npz")
out.parent.mkdir(parents=True, exist_ok=True)

screens = load_screens(target_set="train")
print(f"loaded {len(screens)} training screens")

genes: dict[str, int] = {}
rows, cols = [], []
for i, s in enumerate(screens):
    for g, h in zip(s.genes, s.hits):
        if not h:
            continue
        j = genes.setdefault(g, len(genes))
        rows.append(i); cols.append(j)

M = csr_matrix((np.ones(len(rows), dtype=np.float32), (rows, cols)),
               shape=(len(screens), len(genes)))
save_npz(out.with_suffix(".mat.npz"), M)
np.savez_compressed(out, gene_names=np.array(list(genes), dtype=object),
                    screen_names=np.array([s.dataset_name for s in screens], dtype=object),
                    n_hits=np.asarray(M.sum(axis=1)).ravel())
print(f"{M.shape[0]} screens x {M.shape[1]} genes, {M.nnz:,} hit entries "
      f"({100 * M.nnz / (M.shape[0] * M.shape[1]):.2f}% dense)")
print(f"wrote {out.with_suffix('.mat.npz')} and {out}")
