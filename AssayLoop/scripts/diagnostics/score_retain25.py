# The point estimate MUST be the clustered mean that `boot` returns, not a screen
# mean: the interval is built on publication means, so printing a screen mean beside
# it puts the centre and the interval on different units. They agree on balanced
# sets and diverge badly on unbalanced ones -- on the 146-screen test fold the dm.50
# contrast reads -0.659 by screen and -0.225 by publication, and only the second is
# consistent with its own [-0.769, +0.238]. `dscreen` is still printed, labelled.
"""Three LLM protocols x five arms on paper_test, in canonical adjusted EF.

The protocols differ only in what rounds 1-3 do with GLM-5's proposals:

  COLD     no LLM at all. All 10 rounds are the arm's own.
  WARM     GLM-5's recorded round, verbatim. The reply supplies ~97 of 100 slots
           and the shortfall is LOST, so the campaign spends 5,829 of 6,000 warm
           slots and ~944 of 1,000 overall. The LLM owns ~29% of the budget.
  RETAIN25 the teammate's protocol (`output/runs/external/teammate/README.md`,
           "Proposal handling"): keep the 25 best-scoring proposals of that round,
           the ARM fills the other 75. No slot is ever unfilled, the numerical
           method is live from round 1, and the LLM owns 7.5% of the budget.

Rounds 4-10 are identical in all three. So COLD -> RETAIN25 -> WARM is a dose
response on how much of the budget the LLM is allowed to spend: 0%, 7.5%, 29%.

Metric: canonical adjusted EF = hits / (eff * H/D), eff = budget - n2, with the
unfilled slots CHARGED (`score_glm5_handoff_hsa.py` has the derivation). This is
what makes the protocols comparable: WARM's lost slots are not forgiven.

Clustering is by PUBLICATION: paper_test's 20 screens span 15 publications.
"""
import argparse, statistics as st
import json
from collections import defaultdict

import numpy as np, yaml

ap = argparse.ArgumentParser()
ap.add_argument("--cold", nargs="+",
                default=["output/runs/external/paper_test_hsa.json",
                         "output/runs/external/paper_test_hsa_realmean.json"])
ap.add_argument("--warm", default="output/runs/external/glm5_handoff_hsa.json")
ap.add_argument("--retain", default="output/runs/external/glm5_retain25_hsa.json")
ap.add_argument("--manifest",
                default=".venv/lib/python3.11/site-packages/assaybench/data/"
                        "screen_sets/assayloop-test.yaml")
ap.add_argument("--base", default="AF + DepMap", help="reference arm for the paired block")
ap.add_argument("--n", type=int, default=20000)
a = ap.parse_args()

pub = {s["dataset_name"]: str(s["source_id"])
       for s in yaml.safe_load(open(a.manifest))["screens"]}
BUDGET = 1000
ARMS = ["AF alone", "AF + DepMap", "HSA only (ours)",
        "DepMap + HSA dm.25", "DepMap + HSA dm.50"]
SHORT = {"AF alone": "AF", "AF + DepMap": "AF + DepMap", "HSA only (ours)": "HSA",
         "DepMap + HSA dm.25": "DepMap + HSA .25",
         "DepMap + HSA dm.50": "DepMap + HSA .50"}


def adj(r):
    """eff = budget - n2. Unfilled slots are charged, so NOT n_measured."""
    assert r["n_measured"] + r["n_unknown"] == r["picks"], r["screen"]
    return r["hits"] / ((BUDGET - r["n_unknown"]) * r["total_hits"] / r["n_genes"])


def load(files):
    by = defaultdict(dict)
    for f in files:
        for r in json.load(open(f)):
            by[r["arm"]].setdefault(r["screen"], r)
    return by


P = {"COLD": load(a.cold), "WARM": load([a.warm]), "RETAIN25": load([a.retain])}
S = sorted(set.intersection(*(set(P[p][ARMS[0]]) for p in P)))
missing = [(p, m) for p in P for m in ARMS if m not in P[p]]
if missing:
    raise SystemExit(f"missing arm/protocol pairs: {missing}")

CL = sorted({pub[s] for s in S})


def boot(vals):
    """Publication-clustered bootstrap over a {screen: value} map."""
    cl = defaultdict(list)
    for s, v in vals.items():
        cl[pub[s]].append(v)
    k = list(cl)
    m = np.array([np.mean(cl[x]) for x in k])
    r = np.random.default_rng(0)
    draws = [m[r.integers(0, len(k), len(k))].mean() for _ in range(a.n)]
    return float(m.mean()), *np.percentile(draws, [2.5, 97.5])


print(f"{len(S)} screens over {len(CL)} publications, pool {P['COLD'][ARMS[0]][S[0]]['pool']:,}\n")
for p in ("WARM", "RETAIN25"):
    r0 = P[p][ARMS[0]]
    llm = sum(r0[s].get("llm_picks", r0[s].get("warm_picks", 0)) for s in S)
    fill = sum(r0[s].get("filled_picks", 0) for s in S)
    picks = sum(r0[s]["picks"] for s in S)
    print(f"{p:9s} LLM picks {llm:5,d}  numerically filled in rounds 1-3 {fill:5,d}  "
          f"slots used {picks:,d} of {len(S)*BUDGET:,d}  LLM share {llm/(len(S)*BUDGET):5.1%}")

print("\n=== canonical adjusted EF ===")
print(f"{'arm':20s} {'COLD':>8s} {'WARM':>8s} {'RETAIN25':>9s} "
      f"{'R25-COLD':>9s} {'R25-WARM':>9s}")
for m in ARMS:
    e = {p: st.mean(adj(P[p][m][s]) for s in S) for p in P}
    print(f"{SHORT[m]:20s} {e['COLD']:8.3f} {e['WARM']:8.3f} {e['RETAIN25']:9.3f} "
          f"{e['RETAIN25']-e['COLD']:+9.3f} {e['RETAIN25']-e['WARM']:+9.3f}")

print("\n=== total hits ===")
print(f"{'arm':20s} {'COLD':>8s} {'WARM':>8s} {'RETAIN25':>9s} "
      f"{'R25-COLD':>9s} {'R25-WARM':>9s}")
for m in ARMS:
    h = {p: sum(P[p][m][s]["hits"] for s in S) for p in P}
    print(f"{SHORT[m]:20s} {h['COLD']:8,d} {h['WARM']:8,d} {h['RETAIN25']:9,d} "
          f"{h['RETAIN25']-h['COLD']:+9,d} {h['RETAIN25']-h['WARM']:+9,d}")

print("\n=== (1) inside RETAIN25, vs " + SHORT[a.base] + " ===")
print(f"{'contrast':20s} {'dEF':>7s} {'publication 95% CI':>22s} {'d_hits':>7s} {'w/l':>8s}")
B = P["RETAIN25"][a.base]
for m in ARMS:
    if m == a.base:
        continue
    d = {s: adj(P["RETAIN25"][m][s]) - adj(B[s]) for s in S}
    dmean, lo, hi = boot(d)   # clustered mean, matching the interval
    dh = sum(P["RETAIN25"][m][s]["hits"] - B[s]["hits"] for s in S)
    star = "*" if lo > 0 or hi < 0 else " "
    print(f"{SHORT[m]:20s} {dmean:+7.3f} [{lo:+7.3f},{hi:+7.3f}]{star} "
          f"{dh:+7d} {sum(1 for v in d.values() if v > 0):3d}/{len(S):<4d}")

for label, p0 in (("(2) RETAIN25 minus COLD -- is a 7.5% LLM share worth having?", "COLD"),
                  ("(3) RETAIN25 minus WARM -- does capping the LLM help?", "WARM")):
    print(f"\n=== {label} ===")
    print(f"{'arm':20s} {'dEF':>7s} {'publication 95% CI':>22s} {'d_hits':>7s} {'w/l':>8s}")
    for m in ARMS:
        d = {s: adj(P["RETAIN25"][m][s]) - adj(P[p0][m][s]) for s in S}
        dmean, lo, hi = boot(d)   # clustered mean, matching the interval
        dh = sum(P["RETAIN25"][m][s]["hits"] - P[p0][m][s]["hits"] for s in S)
        star = "*" if lo > 0 or hi < 0 else " "
        print(f"{SHORT[m]:20s} {dmean:+7.3f} [{lo:+7.3f},{hi:+7.3f}]{star} "
              f"{dh:+7d} {sum(1 for v in d.values() if v > 0):3d}/{len(S):<4d}")

print("\n* = publication-level 95% interval excludes zero.")
print("Reference, same warm start, uncapped: AF 4.805/1,522 = sweep-glm5tu-pt-af; "
      "AF + DepMap 5.166/1,606 = sweep-dmcd-hyb-glm5tu-pt.")
