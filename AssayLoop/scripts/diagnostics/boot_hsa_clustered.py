"""Paired bootstrap for HSA arms, at BOTH screen level and publication-cluster level.

paper_validation has 20 screens over only 10 publications -- one accounts for 8 -- so a
screen-level bootstrap treats near-duplicate screens as independent and understates the
interval. The cluster bootstrap resamples PUBLICATIONS, which is the honest unit, and is
the same principle as `codep_gated.py --exclude-set` bootstrapping over hit-set families.

Usage:
  python3 scripts/diagnostics/boot_hsa_clustered.py \
      --results output/runs/external/paper_validation_hsa.json \
      --manifest .venv/.../assayloop-validation.yaml \
      [--base "AF alone" ...]
"""
import argparse, json
from collections import defaultdict

import numpy as np
import yaml

ap = argparse.ArgumentParser()
ap.add_argument("--results", nargs="+", required=True)
ap.add_argument("--manifest", required=True)
ap.add_argument("--base", nargs="+", default=["AF alone", "AF + DepMap", "HSA only (ours)"])
ap.add_argument("--n", type=int, default=20000)
ap.add_argument("--min-hits", type=int, default=0,
                help="keep only screens with total_hits >= this. 50 is the curated sets' own "
                     "num_hits_min; the 218-screen fold applies no such floor, which is why "
                     "1-hit lottery screens dominate its EF (HANDOFF 24.18).")
ap.add_argument("--max-genes", type=int, default=0, help="curated criterion: 22000")
ap.add_argument("--max-hit-rate", type=float, default=0.0, help="curated criterion: 0.15")
ap.add_argument("--adjusted", action="store_true",
                help="score the project's CANONICAL domain-adjusted EF instead of the stored "
                     "`ef` field. The stored field is the runner's `n_hits_vs_random` = "
                     "hits/(1000*H/pool): full budget, POOL prevalence. The canonical metric "
                     "(full_genome_table._adj_ef_from_run, and what the paper's 4.83 is) is "
                     "hits/(eff*H/D) with D = the screen's own library and eff = 1000 minus "
                     "the off-library-but-in-universe picks, which are FORGIVEN. HSA arms "
                     "spend fewer slots off-library, so the two conventions disagree BY ARM "
                     "and paired contrasts move too -- see HANDOFF 24.29. Use this flag for "
                     "anything that will sit in a table beside a sweep-era number.")
ap.add_argument("--min-genes", type=int, default=0,
                help="keep only screens with n_genes >= this. 1002 is the repo's "
                     "GENOME_SCALE_ONLY rule (len(genes) > BUDGET+1); 18000 matches the "
                     "curated sets' own num_genes_min, i.e. genuinely genome-wide.")
a = ap.parse_args()

rows = [r for f in a.results for r in json.load(open(f))]
if a.adjusted:
    # eff = 1000 - (off-library, in-universe) = n_measured, because every pick comes from
    # the f2 universe (so there are no out-of-universe picks) and the budget is always
    # filled -- both asserted here rather than assumed.
    for r in rows:
        assert r["picks"] == 1000, f"{r['arm']}/{r['screen']}: {r['picks']} picks, not 1000"
        assert r["n_measured"] + r["n_unknown"] == r["picks"], r["screen"]
        r["ef"] = r["hits"] / (r["n_measured"] * r["total_hits"] / r["n_genes"])
    print("metric: CANONICAL domain-adjusted EF (library domain, off-library picks forgiven)")
else:
    print("metric: stored `ef` = n_hits_vs_random (pool domain, full 1000 budget). "
          "NOT comparable to sweep-era tables -- see --adjusted")
by = defaultdict(dict)
for r in rows:
    by[r["arm"]][r["screen"]] = r
pub = {s["dataset_name"]: str(s["source_id"])
       for s in yaml.safe_load(open(a.manifest))["screens"]}
S = sorted(by[a.base[0]])
def _keep(s):
    r0 = by[a.base[0]][s]
    rate = r0["total_hits"] / r0["n_genes"] if r0["n_genes"] else 1.0
    return ((not a.min_genes or r0["n_genes"] >= a.min_genes)
            and (not a.max_genes or r0["n_genes"] <= a.max_genes)
            and (not a.min_hits or r0["total_hits"] >= a.min_hits)
            and (not a.max_hit_rate or rate <= a.max_hit_rate))

if a.min_genes or a.max_genes or a.min_hits or a.max_hit_rate:
    n0 = len(S)
    S = [s for s in S if _keep(s)]
    crit = [f"n_genes>={a.min_genes}" if a.min_genes else "",
            f"n_genes<={a.max_genes}" if a.max_genes else "",
            f"total_hits>={a.min_hits}" if a.min_hits else "",
            f"hit_rate<={a.max_hit_rate}" if a.max_hit_rate else ""]
    print("filter " + ", ".join(c for c in crit if c) + f": kept {len(S)} of {n0} screens")
clusters = defaultdict(list)
for s in S:
    clusters[pub[s]].append(s)
ck = sorted(clusters)
print(f"{len(S)} screens over {len(ck)} publications; "
      f"largest cluster {max(len(v) for v in clusters.values())} screens")
rng = np.random.default_rng(0)

def ci(d):
    d = np.asarray(d, float)
    m = d[rng.integers(0, len(d), (a.n, len(d)))].mean(1)
    return d.mean(), np.percentile(m, 2.5), np.percentile(m, 97.5)

def ci_cluster(per_screen):
    """Resample publications; each contributes the MEAN of its screens."""
    cm = np.array([np.mean([per_screen[s] for s in clusters[k]]) for k in ck])
    m = cm[rng.integers(0, len(cm), (a.n, len(cm)))].mean(1)
    return cm.mean(), np.percentile(m, 2.5), np.percentile(m, 97.5)

arms = sorted(by)
print(f"\n{'arm':<26}{'meanEF':>8}{'medEF':>8}{'hits':>7}"
      f"{'  per-pub mean EF':>18}")
for x in arms:
    e = [by[x][s]["ef"] for s in S]
    pe = np.mean([np.mean([by[x][s]["ef"] for s in clusters[k]]) for k in ck])
    print(f"{x:<26}{np.mean(e):>8.3f}{np.median(e):>8.3f}"
          f"{sum(by[x][s]['hits'] for s in S):>7}{pe:>18.3f}")

for b in a.base:
    if b not in by:
        continue
    print(f"\n== vs {b} " + "=" * (84 - len(b)))
    print(f"{'contrast':<26}{'dEF':>8}{'screen-level 95% CI':>24}"
          f"{'publication-level 95% CI':>28}{'d_hits':>8}{'w/l':>8}")
    for x in arms:
        if x == b:
            continue
        d = {s: by[x][s]["ef"] - by[b][s]["ef"] for s in S}
        m, lo, hi = ci(list(d.values()))
        cm, clo, chi = ci_cluster(d)
        dh = sum(by[x][s]["hits"] - by[b][s]["hits"] for s in S)
        w = sum(1 for v in d.values() if v > 1e-9); l = sum(1 for v in d.values() if v < -1e-9)
        s1 = "*" if (lo > 0 or hi < 0) else " "
        s2 = "*" if (clo > 0 or chi < 0) else " "
        print(f"{x:<26}{m:>+8.3f}  [{lo:>+7.3f},{hi:>+7.3f}]{s1}"
              f"     {cm:>+7.3f} [{clo:>+7.3f},{chi:>+7.3f}]{s2}{dh:>+8}{f'{w}/{l}':>8}")
print("\n* = interval excludes zero. The PUBLICATION-level interval is the one to quote.")
