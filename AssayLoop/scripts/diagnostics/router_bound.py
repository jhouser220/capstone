"""How much could routing each screen to AF + DepMap or HSA be worth? Read-only, no API calls.

218-screen validation fold, canonical adjusted EF. Routing signals tried:
  oracle        per-screen max (needs the answer key)
  total_hits    best threshold in-sample (needs the answer key: NOT deployable)
  n_genes       library size, known before the campaign (deployable)
  description   kNN on the 1536-D description embedding, leave-one-publication-out
"""
import json, sys
from collections import defaultdict
import numpy as np, yaml

R = "output/runs/external/validation_fold_hsa_REALCOMBINED.json"
MAN = "_reproduction/donor/validation_fold_218.yaml"
A, B = "AF + DepMap", sys.argv[1] if len(sys.argv) > 1 else "HSA only (ours)"

by = defaultdict(dict)
for r in json.load(open(R)):
    assert r["picks"] == 1000 and r["n_measured"] + r["n_unknown"] == 1000
    r["ef"] = r["hits"] / (r["n_measured"] * r["total_hits"] / r["n_genes"])
    by[r["arm"]][r["screen"]] = r
pub = {s["dataset_name"]: str(s["source_id"]) for s in yaml.safe_load(open(MAN))["screens"]}

# description embeddings from the cache only (EmbeddingCache overlays the user cache)
from assayloop.tasks.screen_sets import load_screens  # noqa: E402
from assayloop.amortized.text_embed import (EMBED_MODEL, EmbeddingCache, _text_key,
                                            screen_description_text)
cache = EmbeddingCache()
emb = {}
for s in load_screens(target_set="validation"):
    v = cache.get(_text_key(screen_description_text(s), model=EMBED_MODEL))
    if v is None:
        raise SystemExit(f"cache miss {s.dataset_name}; refusing to call the API")
    emb[s.dataset_name] = np.asarray(v, float) / np.linalg.norm(v)


def pubmean(S, f):
    g = defaultdict(list)
    for s in S:
        g[pub[s]].append(f(s))
    return np.mean([np.mean(v) for v in g.values()])


def report(S, label):
    ea = {s: by[A][s]["ef"] for s in S}
    eb = {s: by[B][s]["ef"] for s in S}
    d = {s: eb[s] - ea[s] for s in S}
    H = {s: by[A][s]["total_hits"] for s in S}
    N = {s: by[A][s]["n_genes"] for s in S}
    print(f"\n=== {label}: {len(S)} screens, {len(set(pub[s] for s in S))} pubs.  B = {B}")
    print(f"mean EF   A {np.mean(list(ea.values())):.3f}   B {np.mean(list(eb.values())):.3f}"
          f"   (per-pub A {pubmean(S, ea.get):.3f}  B {pubmean(S, eb.get):.3f})")
    print(f"B wins on {sum(v > 0 for v in d.values())}/{len(S)}")

    def gain(route):  # route[s] True -> use B
        e = {s: eb[s] if route[s] else ea[s] for s in S}
        return np.mean([e[s] - ea[s] for s in S]), pubmean(S, lambda s: e[s] - ea[s])

    g, gp = gain({s: d[s] > 0 for s in S})
    print(f"ORACLE router              +{g:.3f} /screen   +{gp:.3f} /pub")
    from scipy.stats import spearmanr
    for nm, x in (("total_hits", H), ("n_genes", N)):
        rho, p = spearmanr([x[s] for s in S], [d[s] for s in S])
        best = max(((gain({s: x[s] >= t for s in S}), t) for t in sorted(set(x.values()))),
                   key=lambda z: z[0][0])
        (g, gp), t = best
        print(f"{nm:<11} spearman(x, dEF) {rho:+.3f} p={p:.2g};  best in-sample rule "
              f"'B if {nm}>={t}': +{g:.3f} /screen +{gp:.3f} /pub")
    # n_genes threshold chosen leave-one-pub-out (honest)
    route = {}
    for s in S:
        tr = [u for u in S if pub[u] != pub[s]]
        def g_tr(t):
            return np.mean([(eb[u] if N[u] >= t else ea[u]) - ea[u] for u in tr])
        t = max(sorted(set(N[u] for u in tr)), key=g_tr)
        route[s] = N[s] >= t
    g, gp = gain(route)
    print(f"n_genes rule, threshold chosen leave-one-pub-out: {g:+.3f} /screen  {gp:+.3f} /pub")
    # description kNN, leave-one-pub-out
    for k in (5, 10, 25):
        route = {}
        for s in S:
            tr = [u for u in S if pub[u] != pub[s]]
            sims = np.array([emb[s] @ emb[u] for u in tr])
            nn = np.argsort(-sims)[:k]
            route[s] = np.mean([d[tr[i]] for i in nn]) > 0
        g, gp = gain(route)
        print(f"description kNN k={k:<3} leave-one-pub-out: {g:+.3f} /screen  {gp:+.3f} /pub"
              f"   (routes {sum(route.values())} to B)")
    # same, but WITHOUT holding out the publication -- shows the leakage channel
    route = {}
    for s in S:
        tr = [u for u in S if u != s]
        sims = np.array([emb[s] @ emb[u] for u in tr])
        nn = np.argsort(-sims)[:5]
        route[s] = np.mean([d[tr[i]] for i in nn]) > 0
    g, gp = gain(route)
    print(f"description kNN k=5, leave-one-SCREEN-out (LEAKY): {g:+.3f} /screen {gp:+.3f} /pub")


S_all = sorted(by[A])
report(S_all, "all 218")
report([s for s in S_all if by[A][s]["n_genes"] >= 1002], "genome-scale (n_genes>=1002)")
