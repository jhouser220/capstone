"""Does AssayBench's confirmed DepMap gain live entirely in its SPARSE screens?

The external datasets say the AF + DepMap gain tracks HIT SPARSITY: pooled over 23
external directional screens, hit fraction predicts the delta (Spearman -0.705,
p=0.0004), screens with fewer than 500 hits gain +0.595 [+0.195, +1.088] (6/6
positive) and screens with 500 or more lose -0.140 [-0.313, -0.009] (3/17). If that
is a real property of the method rather than an artifact of four small datasets, it
has to show up on the benchmark's own held-out folds, where n is 198 and 314 instead
of 1 to 14 -- and AssayBench's curated test screens have a median of 220 hits, i.e.
they sit in the sparse regime that would explain the +0.485 headline.

This script tests that on the folds, reusing the comparison notebook's own
conventions rather than a reimplementation: ``adjusted_ef_value`` and
``ef_curve``'s prefix accounting, ``GENOME_SCALE_ONLY``, the curated-screen
exclusion that makes the folds 198/314 fresh screens, ``gene_universe(...,
min_screen_freq=2)``, and the FAMILY bootstrap (single-linkage at hit-set
Jaccard >= 0.3) that this project uses for every paired delta.

Arms, both on disk:
  AF alone     sweep-dm-validation-w0 / sweep-tf-test-w0   (codep at weight 0)
  AF + DepMap  sweep-dmne-cd-validation / sweep-dmne-cd-test

The AF-alone arms read the FULL atlas at weight 0, where it cannot contribute, so
they are atlas-neutral and the reduced-atlas-only rule does not bite -- which is
exactly how the notebook's registry classifies them (``atlas="none"``).
"""

from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts" / "diagnostics"))

from assayloop.metrics.hits_auc import adjusted_ef_value          # noqa: E402
from assayloop.tasks import gene_universe, load_screens           # noqa: E402
from run_paths import run_dir  # noqa: E402

RUNS = ROOT / "output" / "runs"
BUDGET = 1000
GENOME_SCALE_ONLY = True
EXCLUDE_CURATED_FROM_FOLDS = True

ARMS = {
    "validation": ("sweep-dm-validation-w0", "sweep-dmne-cd-validation"),
    "test":       ("sweep-tf-test-w0",       "sweep-dmne-cd-test"),
}


def ef_at_budget(path: Path, lib: set[str], uni: set[str],
                 total_hits: int, domain: int) -> float | None:
    """Adjusted EF at the full budget, by the notebook's prefix accounting."""
    rd = json.loads(path.read_text())
    steps = rd.get("steps", [])
    if not steps:
        return None
    n1 = n2 = hits = 0
    for st in steps:
        for g, h in zip(st.get("acquired_batch") or [], st.get("hits") or []):
            if g in lib:
                n1 += 1
                hits += bool(h)
            elif g in uni:
                n2 += 1
    return adjusted_ef_value(hits, n1, n2, domain, total_hits, BUDGET)


def collect(sset: str) -> list[dict]:
    base_sid, test_sid = ARMS[sset]
    screens = load_screens(target_set=sset)
    uni = set(gene_universe(screens, min_screen_freq=2))
    curated = {x.dataset_name for n in ("paper_validation", "paper_test")
               for x in load_screens(target_set=n)}
    from score_ctx import families
    fam = families(screens)

    out = []
    for i, s in enumerate(screens):
        if GENOME_SCALE_ONLY and len(s.genes) <= BUDGET + 1:
            continue
        if EXCLUDE_CURATED_FROM_FOLDS and s.dataset_name in curated:
            continue
        if s.total_hits <= 0:
            continue
        lib = set(s.genes)
        vals = {}
        for arm, sid in (("base", base_sid), ("test", test_sid)):
            rp = run_dir(RUNS, sid, i, s.dataset_name) / "result.json"
            if not rp.is_file():
                break
            v = ef_at_budget(rp, lib, uni, s.total_hits, len(lib))
            if v is None:
                break
            vals[arm] = v
        if len(vals) != 2:
            continue
        out.append(dict(sset=sset, screen=s.dataset_name,
                        family=fam.get(s.dataset_name, s.dataset_name),
                        H=s.total_hits, D=len(lib),
                        frac=s.total_hits / len(lib),
                        base=vals["base"], test=vals["test"],
                        d=vals["test"] - vals["base"]))
    return out


def boot_families(rows: list[dict], n: int = 20000, seed: int = 0) -> tuple:
    """Mean and 95% CI of the paired delta, resampling FAMILIES."""
    if not rows:
        return float("nan"), float("nan"), float("nan"), 0
    g = defaultdict(list)
    for r in rows:
        g[(r["sset"], r["family"])].append(r["d"])
    fm = np.array([np.mean(v) for v in g.values()])
    rng = np.random.default_rng(seed)
    draws = fm[rng.integers(0, len(fm), size=(n, len(fm)))].mean(axis=1)
    return fm.mean(), np.percentile(draws, 2.5), np.percentile(draws, 97.5), len(fm)


def line(label: str, rows: list[dict]) -> None:
    o, lo, hi, nf = boot_families(rows)
    pos = sum(1 for r in rows if r["d"] > 0)
    tag = "excludes 0" if (lo > 0 or hi < 0) else "crosses 0 "
    print(f"  {label:26s} n={len(rows):4d} fam={nf:4d}  {o:+7.3f}  "
          f"[{lo:+.3f}, {hi:+.3f}]  {tag}  up {pos:4d}/{len(rows)}")


def main() -> None:
    import logging
    logging.disable(logging.CRITICAL)

    rows = []
    for sset in ARMS:
        r = collect(sset)
        print(f"{sset}: {len(r)} screens with both arms on disk")
        rows += r
    if not rows:
        raise SystemExit("no paired runs found")

    H = np.array([r["H"] for r in rows])
    print(f"\nhits per screen: min {H.min()}  p25 {np.percentile(H,25):.0f}  "
          f"median {np.median(H):.0f}  p75 {np.percentile(H,75):.0f}  max {H.max()}")

    print("\n=== THE TEST: paired dEF by hit count, family bootstrap ===")
    line("ALL", rows)
    print()
    line("SPARSE  H < 500", [r for r in rows if r["H"] < 500])
    line("DENSE   H >= 500", [r for r in rows if r["H"] >= 500])

    print("\n--- finer bins ---")
    edges = [0, 100, 250, 500, 1000, 2000, 10**9]
    names = ["H < 100", "100-250", "250-500", "500-1000", "1000-2000", "H >= 2000"]
    for nm, lo, hi in zip(names, edges[:-1], edges[1:]):
        line(nm, [r for r in rows if lo <= r["H"] < hi])

    print("\n--- by fold, sparse vs dense (does it replicate independently?) ---")
    for sset in ARMS:
        sub = [r for r in rows if r["sset"] == sset]
        line(f"{sset}: SPARSE H<500", [r for r in sub if r["H"] < 500])
        line(f"{sset}: DENSE  H>=500", [r for r in sub if r["H"] >= 500])

    print("\n--- rank correlation of dEF against sparsity ---")
    from scipy.stats import pearsonr, spearmanr
    y = np.array([r["d"] for r in rows])
    for nm, key in (("hit count H", "H"), ("hit fraction", "frac"),
                    ("AF-alone EF", "base")):
        x = np.array([r[key] for r in rows], dtype=float)
        pr, pp = pearsonr(x, y)
        sr, sp = spearmanr(x, y)
        print(f"  {nm:14s} pearson {pr:+.3f} (p={pp:.2g})   "
              f"spearman {sr:+.3f} (p={sp:.2g})")

    out = ROOT / "output" / "runs" / "external" / "benchmark_stratified.json"
    out.write_text(json.dumps(rows, indent=1))
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
