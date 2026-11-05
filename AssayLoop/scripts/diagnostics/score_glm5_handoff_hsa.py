# The point estimate MUST be the clustered mean that `boot` returns, not a screen
# mean: the interval is built on publication means, so printing a screen mean beside
# it puts the centre and the interval on different units. They agree on balanced
# sets and diverge badly on unbalanced ones -- on the 146-screen test fold the dm.50
# contrast reads -0.659 by screen and -0.225 by publication, and only the second is
# consistent with its own [-0.769, +0.238]. `dscreen` is still printed, labelled.
"""Score GLM-5 warm start -> HSA / DepMap + HSA against the cold arms and HANDOFF 4.61.

Two questions, and they need different baselines:
  (1) does HSA beat AF + DepMap *as the receiver of the same warm start*?
      -> paired contrasts inside this run.
  (2) is the warm start worth having at all for an HSA receiver?
      -> each arm here against the SAME arm run cold on paper_test.

Metric: canonical adjusted EF. The warm start under-supplies (291.4 of 300 picks), so
the shortfall is charged and `n_measured < 1000` -- the two EF conventions therefore do
NOT coincide here, unlike the external set. `eff = n_measured` only because every pick
comes from the f2 pool; asserted below.

Clustering is by PUBLICATION: paper_test's 20 screens span 15 publications.
"""
import argparse, json, statistics as st
from collections import defaultdict
import numpy as np, yaml

ap = argparse.ArgumentParser()
ap.add_argument("--warm", default="output/runs/external/glm5_handoff_hsa.json")
ap.add_argument("--cold", nargs="+",
                default=["output/runs/external/paper_test_hsa.json",
                         "output/runs/external/paper_test_hsa_realmean.json"])
ap.add_argument("--manifest",
                default=".venv/lib/python3.11/site-packages/assaybench/data/"
                        "screen_sets/assayloop-test.yaml")
ap.add_argument("--n", type=int, default=20000)
a = ap.parse_args()

pub = {s["dataset_name"]: str(s["source_id"])
       for s in yaml.safe_load(open(a.manifest))["screens"]}

BUDGET = 1000

def adj(r):
    """Canonical adjusted EF: hits / (eff * H/D), eff = budget - n2.

    n2 = picks off this screen's library but inside the pool -- FORGIVEN. Everything
    else is charged, INCLUDING THE UNFILLED SLOTS the warm start failed to supply
    (n3 = budget - n1 - n2), which is why `eff` is `budget - n_unknown` and NOT
    `n_measured`. Using n_measured here understates the budget and inflated every
    EF by about +0.07 on this run, because the warm start delivers only 944 of 1000
    picks on average. `metrics.hits_auc.adjusted_ef_value` is the reference.
    """
    assert r["n_measured"] + r["n_unknown"] == r["picks"], r["screen"]
    eff = BUDGET - r["n_unknown"]
    return r["hits"] / (eff * r["total_hits"] / r["n_genes"])

def load(files):
    by = defaultdict(dict)
    for f in files:
        for r in json.load(open(f)):
            by[r["arm"]].setdefault(r["screen"], r)
    return by

W = load([a.warm])
C = load(a.cold)
S = sorted(W[next(iter(W))])
print(f"{len(S)} screens over {len(set(pub[s] for s in S))} publications")
w0 = W[next(iter(W))]
print(f"warm start: rounds 1-{w0[S[0]]['n_warm']}, "
      f"{sum(w0[s]['warm_picks'] for s in S)} of {len(S)*300} picks supplied\n")

def boot(vals):
    cl = defaultdict(list)
    for s, v in vals.items(): cl[pub[s]].append(v)
    k = list(cl); m = np.array([np.mean(cl[x]) for x in k]); rng = np.random.default_rng(0)
    return float(m.mean()), *np.percentile(
        [m[rng.integers(0, len(k), len(k))].mean() for _ in range(a.n)], [2.5, 97.5])

print("=== WARM (GLM-5 rounds 1-3, then the arm) ===")
print(f"{'arm':24s} {'meanEF':>8s} {'hits':>7s} {'measured':>9s}")
for arm in sorted(W, key=lambda x: -st.mean(adj(W[x][s]) for s in S)):
    print(f"{arm:24s} {st.mean(adj(W[arm][s]) for s in S):8.3f} "
          f"{sum(W[arm][s]['hits'] for s in S):7,d} "
          f"{st.mean(W[arm][s]['n_measured'] for s in S):9.1f}")

B = "AF + DepMap"
print(f"\n=== (1) inside the warm run, vs {B} ===")
print(f"{'contrast':24s} {'dEF':>7s} {'publication 95% CI':>24s} {'d_hits':>7s} {'w/l':>7s}")
for arm in W:
    if arm == B: continue
    d = {s: adj(W[arm][s]) - adj(W[B][s]) for s in S}
    m, lo, hi = boot(d)
    dh = sum(W[arm][s]["hits"] - W[B][s]["hits"] for s in S)
    star = "*" if lo > 0 or hi < 0 else " "
    print(f"{arm:24s} {m:+7.3f} [{lo:+8.3f},{hi:+8.3f}]{star} "
          f"{dh:+7d} {sum(1 for v in d.values() if v>0):3d}/{len(S):<3d}")

print("\n=== (2) warm MINUS cold, same arm, same screens ===")
print(f"{'arm':24s} {'cold EF':>8s} {'warm EF':>8s} {'dEF':>7s} "
      f"{'publication 95% CI':>24s} {'d_hits':>7s}")
for arm in W:
    if arm not in C: 
        print(f"{arm:24s} {'(no cold twin)':>8s}"); continue
    d = {s: adj(W[arm][s]) - adj(C[arm][s]) for s in S}
    m, lo, hi = boot(d)
    dh = sum(W[arm][s]["hits"] - C[arm][s]["hits"] for s in S)
    star = "*" if lo > 0 or hi < 0 else " "
    print(f"{arm:24s} {st.mean(adj(C[arm][s]) for s in S):8.3f} "
          f"{st.mean(adj(W[arm][s]) for s in S):8.3f} {m:+7.3f} "
          f"[{lo:+8.3f},{hi:+8.3f}]{star} {dh:+7d}")
print("\n* = publication-level interval excludes zero.")
print("Controls, same warm start, canonical EF -- REPRODUCED EXACTLY:")
print("  AF alone    4.805 / 1,522 hits = sweep-glm5tu-pt-af")
print("  AF + DepMap 5.166 / 1,606 hits = sweep-dmcd-hyb-glm5tu-pt")
print("  (5.251 / 1,672 is sweep-dmne-hyb-glm5tu-pt, the TWO-TERM receiver: + screen-kNN)")
