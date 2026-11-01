"""Pickle the 348 held-out benchmark screens so the EXTERNAL driver can run them.

The point of this is to put the benchmark through the identical code path as the
six external datasets: run_external.py sets ``universe_genes`` to each screen's
own library, so no pick can land outside it, nothing is forgiven, and the
effective budget is exactly 1000. On the benchmark's own runs the candidate pool
is the fold-wide union instead, and 562 of every 1000 picks are forgiven, so the
effective budget there is ~437. This isolates that difference.

Filters match the comparison notebook: genome-scale only, curated screens
excluded from the folds, and screens with no hits dropped.
"""
from __future__ import annotations
import argparse, logging, pickle
from pathlib import Path
from assayloop.tasks import load_screens

BUDGET = 1000

def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()

    curated = {x.dataset_name for n in ("paper_validation", "paper_test")
               for x in load_screens(target_set=n)}
    keep = []
    for sset in ("validation", "test"):
        for s in load_screens(target_set=sset):
            if len(s.genes) <= BUDGET + 1:       continue
            if s.dataset_name in curated:        continue
            if s.total_hits <= 0:                continue
            keep.append(s)
    # dataset_name must stay unique across the two folds for the driver's run ids
    seen = set()
    uniq = []
    for s in keep:
        if s.dataset_name in seen:
            continue
        seen.add(s.dataset_name)
        uniq.append(s)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "wb") as fh:
        pickle.dump(uniq, fh)
    import numpy as np
    H = np.array([s.total_hits for s in uniq]); G = np.array([len(s.genes) for s in uniq])
    print(f"{len(uniq)} screens -> {a.out}")
    print(f"  genes  min {G.min()} median {int(np.median(G))} max {G.max()}")
    print(f"  hits   min {H.min()} median {int(np.median(H))} max {H.max()}")


if __name__ == "__main__":
    main()
