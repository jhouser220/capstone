"""Summarise an external-dataset AF vs AF + DepMap run.

**On the unit of analysis.** This repo's paired deltas bootstrap over hit-set
FAMILIES -- single-linkage clusters at Jaccard >= 0.3 -- because AssayBench
holds many screens that are near-duplicates of each other. That convention does
not transfer here: these are all knockout *fitness* screens in related cancer
lines, they share their core-essential hits (median pairwise Jaccard 0.481, max
0.804), and the rule collapses 13 of the 14 into ONE family. A family bootstrap
would be n=2.

So the primary interval here clusters by CELL LINE (n=7), which is the real
experimental unit: each line contributes a near-euploid control plus one or two
engineered aneuploid clones, and those are the screens that resemble each other
most. The per-screen interval (n=14) is reported beside it and is too narrow,
because it treats a control and its own derived clone as independent. Quote the
cell-line figure; this is a 7-unit experiment wearing 14 screens.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

B = 10000


def cluster(name: str) -> str:
    """Cell line, which is the resampling unit."""
    return name.replace("SCHUKKEN_", "").split("_")[0]


def boot(deltas: dict[str, list[float]], rng: np.random.Generator) -> tuple:
    """Cluster bootstrap: resample GROUPS, then average their member deltas."""
    keys = list(deltas)
    flat = [d for k in keys for d in deltas[k]]
    obs = float(np.mean(flat))
    draws = np.empty(B)
    for b in range(B):
        pick = rng.integers(0, len(keys), len(keys))
        vals = [d for i in pick for d in deltas[keys[i]]]
        draws[b] = np.mean(vals)
    return obs, float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True, type=Path)
    ap.add_argument("--base", default="AF alone")
    ap.add_argument("--test", default="AF + DepMap")
    a = ap.parse_args()

    rows = json.loads(a.results.read_text())
    by = {(r["arm"], r["screen"]): r for r in rows}
    screens = sorted({r["screen"] for r in rows})
    rng = np.random.default_rng(0)

    print(f"{len(screens)} screens | {a.base}  vs  {a.test}\n")
    hdr = (f"{'screen':24s} {'H':>5} | {'hits':>5} {'rate':>6} {'EF':>6} | "
           f"{'hits':>5} {'rate':>6} {'EF':>6} | {'dHits':>6} {'dEF':>7}")
    print(hdr)
    print("-" * len(hdr))

    d_ef: dict[str, list[float]] = {}
    d_hits: dict[str, list[float]] = {}
    tot = {a.base: 0, a.test: 0}
    for s in screens:
        b, t = by.get((a.base, s)), by.get((a.test, s))
        if not b or not t:
            continue
        c = cluster(s)
        d_ef.setdefault(c, []).append(t["ef"] - b["ef"])
        d_hits.setdefault(c, []).append(t["hits"] - b["hits"])
        tot[a.base] += b["hits"]
        tot[a.test] += t["hits"]
        print(f"{s.replace('SCHUKKEN_',''):24s} {b['total_hits']:5d} | "
              f"{b['hits']:5d} {b['hit_rate']:6.3f} {b['ef']:6.3f} | "
              f"{t['hits']:5d} {t['hit_rate']:6.3f} {t['ef']:6.3f} | "
              f"{t['hits']-b['hits']:+6d} {t['ef']-b['ef']:+7.3f}")

    n = len([s for s in screens if (a.base, s) in by and (a.test, s) in by])
    print("-" * len(hdr))
    mean = lambda arm, k: np.mean([by[(arm, s)][k] for s in screens if (arm, s) in by])
    print(f"{'MEAN':24s} {'':5s} | "
          f"{tot[a.base]/n:5.1f} {mean(a.base,'hit_rate'):6.3f} {mean(a.base,'ef'):6.3f} | "
          f"{tot[a.test]/n:5.1f} {mean(a.test,'hit_rate'):6.3f} {mean(a.test,'ef'):6.3f} | "
          f"{(tot[a.test]-tot[a.base])/n:+6.1f} "
          f"{mean(a.test,'ef')-mean(a.base,'ef'):+7.3f}")
    print(f"{'TOTAL hits, 10 rounds':24s} {'':5s} | {tot[a.base]:5d} {'':6s} {'':6s} | "
          f"{tot[a.test]:5d} {'':6s} {'':6s} | {tot[a.test]-tot[a.base]:+6d}")

    print("\nPaired deltas, cluster bootstrap by CELL LINE "
          f"(n={len(d_ef)} lines, {n} screens, B={B}):")
    for label, d in (("EF", d_ef), ("total hits", d_hits)):
        o, lo, hi = boot(d, rng)
        star = "excludes 0" if lo > 0 or hi < 0 else "CROSSES 0"
        print(f"  {label:11s} {o:+8.3f}  [{lo:+.3f}, {hi:+.3f}]   {star}")

    print(f"\nFor contrast, per-SCREEN bootstrap (n={n}, too narrow -- "
          "treats a clone as independent of its own control):")
    for label, d in (("EF", d_ef), ("total hits", d_hits)):
        flat = {f"s{i}": [v] for i, v in enumerate(x for k in d for x in d[k])}
        o, lo, hi = boot(flat, rng)
        print(f"  {label:11s} {o:+8.3f}  [{lo:+.3f}, {hi:+.3f}]")

    wins = sum(1 for k in d_ef for v in d_ef[k] if v > 0)
    print(f"\nScreens where {a.test} wins on EF: {wins}/{n}")

    print("\nPer-round cumulative hits (mean over screens):")
    print(f"{'round':>6} {'picks':>7} " + " ".join(f"{x:>14s}" for x in (a.base, a.test))
          + f" {'delta':>8}")
    for r in range(1, 11):
        got = {}
        for arm in (a.base, a.test):
            vals = [pr["cum_hits"] for s in screens if (arm, s) in by
                    for pr in by[(arm, s)]["per_round"] if pr["round"] == r]
            got[arm] = np.mean(vals) if vals else float("nan")
        print(f"{r:6d} {r*100:7d} {got[a.base]:14.1f} {got[a.test]:14.1f} "
              f"{got[a.test]-got[a.base]:+8.1f}")

    print("\nPer-round cumulative EF (mean over screens):")
    for r in range(1, 11):
        got = {}
        for arm in (a.base, a.test):
            vals = [pr["ef"] for s in screens if (arm, s) in by
                    for pr in by[(arm, s)]["per_round"] if pr["round"] == r]
            got[arm] = np.mean(vals) if vals else float("nan")
        print(f"{r:6d} {r*100:7d} {got[a.base]:14.3f} {got[a.test]:14.3f} "
              f"{got[a.test]-got[a.base]:+8.3f}")


if __name__ == "__main__":
    main()
