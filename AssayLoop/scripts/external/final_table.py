"""AF alone vs AF + DepMap on every screen run, all in ONE metric regime.

Every row here was produced by run_external.py with ``universe_genes`` set to
the screen's own library, so no pick is ever outside it, nothing is forgiven,
and the effective budget is exactly 1000 everywhere. The six external datasets
were always run this way; the AssayBench folds were re-run this way
(``benchlib_results.json``) so they sit on the same footing as the rest rather
than on their native fold-union pool, where 562 of every 1000 picks are
forgiven and the effective budget is ~437.

Resampling cluster per dataset: hit-set family for the benchmark (this
project's convention), cell line for Schukken and Shan, drug for Nguyen,
organoid for the biobank, and the whole dataset for the single-experiment sets
(Abt, Wilson, Tau).
``_merged`` arms are excluded: they are unions of their own directional
screens, not independent units.
"""
from __future__ import annotations
import json
from collections import defaultdict
import numpy as np

EXT = "output/runs/external"
rng = np.random.default_rng(0)


def boot(pairs, n=20000):
    g = defaultdict(list)
    for k, v in pairs:
        g[k].append(v)
    fm = np.array([np.mean(v) for v in g.values()])
    dr = fm[rng.integers(0, len(fm), size=(n, len(fm)))].mean(axis=1)
    return fm.mean(), np.percentile(dr, 2.5), np.percentile(dr, 97.5), len(fm)


def load(path, keep=lambda s: True, cluster=None):
    rows = json.loads(open(path).read())
    by = {(r["arm"], r["screen"]): r for r in rows}
    out = []
    for s in sorted({r["screen"] for r in rows}):
        if not keep(s) or ("AF alone", s) not in by:
            continue
        b, t = by[("AF alone", s)], by[("AF + DepMap", s)]
        out.append(dict(screen=s, fam=(cluster(s) if cluster else s),
                        H=b["total_hits"], af=b["ef"], dm=t["ef"],
                        d=t["ef"] - b["ef"], hb=b["hits"], ht=t["hits"]))
    return out


def org_cluster(s):
    n = s.replace("ORGANOID_", "")
    for lib in ("_minLib", "_v1_1"):
        if n.endswith(lib):
            return n[: -len(lib)]
    return n


nomerge = lambda s: not s.endswith("_merged")

# Benchmark families come from the stratification run, keyed by screen.
fam = {r["screen"]: f"{r['sset']}|{r['family']}"
       for r in json.loads(open(f"{EXT}/benchmark_stratified.json").read())}

DATASETS = [
    ("AssayBench folds (held out)", "benchlib_results.json",
     "curated, mixed phenotypes", lambda s: True, lambda s: fam.get(s, s)),
    ("Schukken 2026",  "schukken_results.json", "knockout fitness",
     nomerge, lambda s: s.replace("SCHUKKEN_", "").split("_")[0]),
    ("Abt 2026",       "abt_results.json",      "knockout fitness (SCPC)",
     nomerge, None),
    ("Wilson 2024",    "wilson_results.json",   "reporter, FACS-sorted",
     nomerge, lambda s: "wilson"),
    ("Nguyen 2025",    "nguyen_results.json",   "drug response (ADE)",
     nomerge, lambda s: s.split("_")[1]),
    ("Shan 2025",      "shan_results.json",     "drug / genetic interaction",
     nomerge, lambda s: s.split("_")[1]),
    ("Organoid 2026",  "organoid_results.json", "organoid dependency",
     nomerge, org_cluster),
    # The only CRISPRi dataset here, and the sparsest: one genome-wide screen
    # split into its two tails, so the cluster is the whole experiment.
    ("Tau CRISPRi 2026", "tau_results.json",    "reporter, FACS-sorted (CRISPRi)",
     nomerge, lambda s: "tau"),
]

print("AF alone  vs  AF + DepMap      (universe = library, b_eff = 1000, everywhere)\n")
hdr = (f"{'dataset':28s} {'n':>4} {'clus':>5} {'avail H':>8} | "
       f"{'AF/scr':>7} {'+DM/scr':>8} {'d/scr':>7} | "
       f"{'AF EF':>7} {'+DM EF':>7} {'dEF':>7} {'95% CI':>19} | {'up':>9}")
print(hdr); print("-" * len(hdr))

ext_pairs, ext_n = [], 0
ext_hb = ext_ht = ext_H = 0
for name, f, ph, keep, cl in DATASETS:
    rs = load(f"{EXT}/{f}", keep, cl)
    o, lo, hi, nf = boot([(r["fam"], r["d"]) for r in rs])
    hb, ht = sum(r["hb"] for r in rs), sum(r["ht"] for r in rs)
    up = sum(1 for r in rs if r["d"] > 0)
    if nf < 3:
        ci = f"{'(n too small)':>21s}"
    else:
        ci = f"[{lo:+6.3f},{hi:+6.3f}]" + ("*" if (lo > 0 or hi < 0) else " ")
    n = len(rs)
    print(f"{name:28s} {n:4d} {nf:5d} {np.mean([r['H'] for r in rs]):8.0f} | "
          f"{hb/n:7.1f} {ht/n:8.1f} {(ht-hb)/n:+7.1f} | "
          f"{np.mean([r['af'] for r in rs]):7.3f} {np.mean([r['dm'] for r in rs]):7.3f} "
          f"{o:+7.3f} {ci} | {up:4d}/{n:<4d}")
    if name.startswith("AssayBench"):
        continue
    ext_pairs += [(f"{name}|{r['fam']}", r["d"]) for r in rs]
    ext_n += len(rs); ext_hb += hb; ext_ht += ht
    ext_H += sum(r['H'] for r in rs)

print("-" * len(hdr))
o, lo, hi, nf = boot(ext_pairs)
up = sum(1 for _, v in ext_pairs if v > 0)
print(f"{'EXTERNAL POOLED (7 papers)':28s} {ext_n:4d} {nf:5d} {ext_H/ext_n:8.0f} | "
      f"{ext_hb/ext_n:7.1f} {ext_ht/ext_n:8.1f} {(ext_ht-ext_hb)/ext_n:+7.1f} | "
      f"{'':7s} {'':7s} {o:+7.3f} [{lo:+6.3f},{hi:+6.3f}]* | {up:4d}/{ext_n:<4d}")
print("\n* interval excludes zero")
print("avail H = mean hits available per screen; /scr columns are hits FOUND per screen in 1000 picks")
