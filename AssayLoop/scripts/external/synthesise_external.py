"""Pool every external directional screen and test what predicts the DepMap delta."""
import json, numpy as np
from scipy.stats import spearmanr, pearsonr

FILES = {
    "Schukken": ("output/runs/external/schukken_results.json", lambda s: True),
    "Abt":      ("output/runs/external/abt_results.json",      lambda s: True),
    "Wilson":   ("output/runs/external/wilson_results.json",   lambda s: not s.endswith("_merged")),
    "Nguyen":   ("output/runs/external/nguyen_results.json",
                 lambda s: not s.endswith("_merged")),
    "Shan":     ("output/runs/external/shan_results.json",     lambda s: True),
    "Organoid": ("output/runs/external/organoid_results.json", lambda s: True),
}
rec = []
for ds, (path, keep) in FILES.items():
    rows = json.loads(open(path).read())
    by = {(r["arm"], r["screen"]): r for r in rows}
    for s in sorted({r["screen"] for r in rows}):
        if not keep(s) or ("AF alone", s) not in by or ("AF + DepMap", s) not in by:
            continue
        b, t = by[("AF alone", s)], by[("AF + DepMap", s)]
        ceil = b["n_genes"] / max(b["total_hits"], 1000)
        rec.append(dict(ds=ds, screen=s, H=b["total_hits"], D=b["n_genes"],
                        frac=b["total_hits"]/b["n_genes"],
                        pct_ceil=100*b["ef"]/ceil, af_ef=b["ef"], af_rate=b["hit_rate"],
                        dEF=t["ef"]-b["ef"], dHits=t["hits"]-b["hits"]))
print(f"{len(rec)} external directional screens\n")
print(f"{'dataset':10s} {'screen':22s} {'H':>5} {'%ceil':>6} {'AF EF':>6} {'dEF':>7}")
for r in sorted(rec, key=lambda r: r["pct_ceil"]):
    print(f"{r['ds']:10s} {r['screen'][:22]:22s} {r['H']:5d} {r['pct_ceil']:5.0f}% {r['af_ef']:6.2f} {r['dEF']:+7.3f}")

y = np.array([r["dEF"] for r in rec])
print(f"\nmean dEF over all {len(rec)}: {y.mean():+.3f}   positive on {int((y>0).sum())}/{len(y)}")
print("\npredictor correlations with dEF:")
for nm, key in (("AF % of ceiling","pct_ceil"), ("hit count H","H"), ("hit fraction","frac"),
                ("AF-alone EF","af_ef"), ("AF-alone hit rate","af_rate")):
    x = np.array([r[key] for r in rec], dtype=float)
    pr, pp = pearsonr(x, y); sr, sp = spearmanr(x, y)
    print(f"  {nm:18s} pearson {pr:+.3f} (p={pp:.3f})   spearman {sr:+.3f} (p={sp:.3f})")

print("\nsplit by hit count:")
for lab, sel in (("sparse  H<500", [r for r in rec if r["H"]<500]),
                 ("dense   H>=500",[r for r in rec if r["H"]>=500])):
    v=np.array([r["dEF"] for r in sel])
    print(f"  {lab:15s} n={len(v):2d}  mean dEF {v.mean():+.3f}  positive {int((v>0).sum())}/{len(v)}")
print("\nsplit by % of ceiling (my claim last turn: helps below ~25-43%):")
for lab, sel in (("<35% of ceiling", [r for r in rec if r["pct_ceil"]<35]),
                 (">=35% of ceiling",[r for r in rec if r["pct_ceil"]>=35])):
    v=np.array([r["dEF"] for r in sel])
    print(f"  {lab:17s} n={len(v):2d}  mean dEF {v.mean():+.3f}  positive {int((v>0).sum())}/{len(v)}")

# ---- clustered bootstrap on the sparsity split -------------------------------
def clust(r):
    """Resampling unit: the experiment that screens share structure within."""
    if r["ds"] == "Schukken":
        return "SCH_" + r["screen"].replace("SCHUKKEN_", "").split("_")[0]  # cell line
    if r["ds"] == "Nguyen":
        return "NGU_" + r["screen"].split("_")[1]                            # drug
    if r["ds"] == "Shan":
        return "SHA_" + r["screen"].split("_")[1]                            # cell line
    if r["ds"] == "Organoid":
        # library suffix, not the sample id, is what must be stripped -- and
        # "v1_1" itself contains an underscore, so rsplit("_") mis-splits it.
        n = r["screen"].replace("ORGANOID_", "")
        for lib in ("_minLib", "_v1_1"):
            if n.endswith(lib):
                return "ORG_" + n[: -len(lib)]
        return "ORG_" + n
    return r["ds"]                                                           # one each

rng2 = np.random.default_rng(7)
def cb(sel, B=20000):
    from collections import defaultdict
    d = defaultdict(list)
    for r in sel:
        d[clust(r)].append(r["dEF"])
    ks = list(d)
    flat = [x for k in ks for x in d[k]]
    dr = np.array([np.mean([x for i in rng2.integers(0, len(ks), len(ks)) for x in d[ks[i]]])
                   for _ in range(B)])
    return np.mean(flat), np.percentile(dr, 2.5), np.percentile(dr, 97.5), len(ks)

print("\n=== clustered bootstrap (unit = cell line / drug / experiment) ===")
for lab, sel in (("ALL 23", rec),
                 ("SPARSE  H<500", [r for r in rec if r["H"] < 500]),
                 ("DENSE   H>=500", [r for r in rec if r["H"] >= 500])):
    o, lo, hi, nk = cb(sel)
    tag = "excludes 0" if (lo > 0 or hi < 0) else "crosses 0"
    print(f"  {lab:15s} n={len(sel):2d} screens, {nk:2d} clusters   dEF {o:+.3f}  [{lo:+.3f}, {hi:+.3f}]  {tag}")

print("\n  sparse screens, individually:")
for r in sorted([r for r in rec if r["H"] < 500], key=lambda r: r["H"]):
    print(f"    {r['ds']:9s} {r['screen'][:24]:24s} H={r['H']:4d}  dEF {r['dEF']:+.3f}")
