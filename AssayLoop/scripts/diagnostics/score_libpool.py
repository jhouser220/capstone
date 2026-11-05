# The point estimate MUST be the clustered mean that `boot` returns, not a screen
# mean: the interval is built on publication means, so printing a screen mean beside
# it puts the centre and the interval on different units. They agree on balanced
# sets and diverge badly on unbalanced ones -- on the 146-screen test fold the dm.50
# contrast reads -0.659 by screen and -0.225 by publication, and only the second is
# consistent with its own [-0.769, +0.238]. `dscreen` is still printed, labelled.
"""paper_test scored on each screen's OWN library — the run with no metric ambiguity.

Why this run exists. The f2 pool (21,147 genes) is not the native setting: AssayBench hands
each screen its own library unless `universe_genes` says otherwise, and the paper's f2 table
exists only "with full-genome candidate pool for fair LLM comparison"
(`full_genome_table.py`). That pool costs a forgiveness term — off-library-but-in-pool picks
are refunded — and the refund is per-arm, so it reordered `AF + DepMap` and `DepMap + HSA
dm.50` (SUMMARY section 7 rule 1). With the library as the pool no pick can land off-library,
so n2 = n3 = 0, nothing is refunded, nothing is charged unfairly, and

    stored `ef` == canonical adjusted EF == hits / (1000 * H/D)

are the same number. That identity is ASSERTED below per row, not assumed.

The price: absolute EF here is NOT comparable to an f2 number (different denominator, section 7
rule 4). The f2 column is printed only to show how much each arm moved, never to be quoted
beside it.
"""
import argparse
import json
import statistics as st
from collections import defaultdict

import numpy as np
import yaml

ap = argparse.ArgumentParser()
ap.add_argument("--results", default="output/runs/external/paper_test_hsa_libpool.json")
ap.add_argument("--f2", nargs="+",
                default=["output/runs/external/paper_test_hsa_realmean.json",
                         "output/runs/external/paper_test_hsa.json"],
                help="the f2-pool run, for the 'moved by' column only. ORDER MATTERS: the "
                     "FIRST file listed wins for an arm present in several, so the REAL-MEAN "
                     "file goes first. The own-library runs use `banks_real_*`, so comparing "
                     "them against the mean=0 file mixes bank variants -- that mistake cost "
                     "-0.10 EF on paper_test's `HSA only` row. `--f2-provenance` prints the "
                     "file each arm was taken from.")
ap.add_argument("--f2-provenance", action="store_true",
                help="print which f2 file supplied each arm, then continue")
ap.add_argument("--manifest",
                default=".venv/lib/python3.11/site-packages/assaybench/data/"
                        "screen_sets/assayloop-test.yaml")
ap.add_argument("--base", default="AF + DepMap")
ap.add_argument("--n", type=int, default=20000)
a = ap.parse_args()

pub = {s["dataset_name"]: str(s["source_id"])
       for s in yaml.safe_load(open(a.manifest))["screens"]}
BUDGET = 1000
SHORT = {"AF alone": "AF", "AF + DepMap": "AF + DepMap", "HSA only (ours)": "HSA",
         "DepMap + HSA dm.25": "DepMap + HSA .25", "DepMap + HSA dm.50": "DepMap + HSA .50"}


def load(files):
    """First file listed wins for an arm it shares with a later one."""
    by, src = defaultdict(dict), {}
    for f in files:
        for r in json.load(open(f)):
            if r["screen"] not in by[r["arm"]]:
                by[r["arm"]][r["screen"]] = r
                src.setdefault(r["arm"], f)
    load.source = src
    return by


L = load([a.results])
F = load(a.f2)
F_SRC = dict(load.source)
if a.f2_provenance:
    print("f2 rows taken from:")
    for k in sorted(F_SRC):
        print(f"  {k:22s} <- {F_SRC[k]}")
    print()
S = sorted(L[a.base])
ARMS = [m for m in SHORT if m in L]

# the whole point of this run: one metric, verified row by row
for m in ARMS:
    for s in S:
        r = L[m][s]
        assert r["pool_mode"] == "library", f"{s}: {r['pool_mode']} pool, not library"
        assert r["n_unknown"] == 0 and r["n_measured"] == r["picks"] == BUDGET, r
        own = r["hits"] / (BUDGET * r["total_hits"] / r["n_genes"])
        assert abs(own - r["ef"]) < 2e-3, f"{m}/{s}: stored {r['ef']:.4f} vs own-library {own:.4f}"
print(f"{len(ARMS)} arms x {len(S)} screens: every row has 1,000 measured picks, "
      f"0 off-library, and stored ef == hits/(1000*H/D). ONE metric.\n")


def boot(vals):
    cl = defaultdict(list)
    for s, v in vals.items():
        cl[pub[s]].append(v)
    k = list(cl)
    mm = np.array([np.mean(cl[x]) for x in k])
    rng = np.random.default_rng(0)
    draws = [mm[rng.integers(0, len(k), len(k))].mean() for _ in range(a.n)]
    return float(mm.mean()), *np.percentile(draws, [2.5, 97.5])


print(f"{len(S)} screens over {len(set(pub[s] for s in S))} publications, "
      f"pool = own library (mean {st.mean(L[a.base][s]['pool'] for s in S):,.0f} genes, "
      f"bank-axis restricted)\n")
print(f"{'arm':20s} {'EF':>8s} {'hits':>7s} {'hit rate':>9s} | "
      f"{'f2 EF':>7s} {'f2 hits':>8s} {'moved':>7s}")
for m in sorted(ARMS, key=lambda m: -st.mean(L[m][s]["ef"] for s in S)):
    ef = st.mean(L[m][s]["ef"] for s in S)
    hits = sum(L[m][s]["hits"] for s in S)
    rate = hits / (len(S) * BUDGET)
    if m in F:
        f2 = st.mean(F[m][s]["hits"] / ((BUDGET - F[m][s]["n_unknown"])
                                        * F[m][s]["total_hits"] / F[m][s]["n_genes"])
                     for s in S)
        f2h = sum(F[m][s]["hits"] for s in S)
        tail = f"{f2:7.3f} {f2h:8,d} {ef - f2:+7.3f}"
    else:
        tail = f"{'-':>7s} {'-':>8s} {'-':>7s}"
    print(f"{SHORT[m]:20s} {ef:8.3f} {hits:7,d} {rate:9.2%} | {tail}")

print(f"\n=== paired, vs {SHORT[a.base]} ===")
print(f"{'contrast':20s} {'dEF':>7s} {'publication 95% CI':>22s} {'d_hits':>7s} {'w/l':>8s}")
B = L[a.base]
for m in ARMS:
    if m == a.base:
        continue
    d = {s: L[m][s]["ef"] - B[s]["ef"] for s in S}
    dmean, lo, hi = boot(d)   # clustered mean, matching the interval
    dh = sum(L[m][s]["hits"] - B[s]["hits"] for s in S)
    star = "*" if lo > 0 or hi < 0 else " "
    print(f"{SHORT[m]:20s} {dmean:+7.3f} [{lo:+7.3f},{hi:+7.3f}]{star} "
          f"{dh:+7d} {sum(1 for v in d.values() if v > 0):3d}/{len(S):<4d}")
print("\n* = publication-level 95% interval excludes zero.")
print("Absolute EF is not comparable to the f2 columns -- different denominator (rule 4).")
