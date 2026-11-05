"""Score the HSA arms on the teammate's 11 frozen external screens.

Pool is each screen's own library (intersected with the bank axis, uniformly across
every arm), so no pick can land outside it: n2 = n3 = 0, eff = 1000, and the project's
canonical adjusted EF is IDENTICALLY the runner's `n_hits_vs_random`. One convention
here, unlike every benchmark set (HANDOFF 24.29/24.33).

Clustering is by PAPER, not by screen: 12 rows come from 7 papers, Shan alone supplies 3
and the organoid study 3. `SCHUKKEN_PAIRED_WGL_{paper,rank194}` are two readings of ONE
screen whose hit list is not reproducible (SOURCES.md), so they are never both counted --
`--schukken` picks which, or `drop` excludes the screen.
"""
import argparse, json, statistics as st
from collections import defaultdict
import numpy as np

PAPER = {
    "ABT_PARCB_dependency": "abt",
    "NGUYEN_AraC_resistant": "nguyen",
    "ORGANOID_WTSI-COLIVM_031_minLib": "organoid",
    "ORGANOID_WTSI-OESO_005_minLib": "organoid",
    "ORGANOID_WTSI-OVAR_018_minLib": "organoid",
    "SCHUKKEN_PAIRED_WGL_paper": "schukken",
    "SCHUKKEN_PAIRED_WGL_rank194": "schukken",
    "SHAN_HeLa_LDLRko": "shan",
    "SHAN_HeLa_lovastatin": "shan",
    "SHAN_HepG2_lovastatin": "shan",
    "TAU_CRISPRi_dec": "tau",
    "WILSON_MERCS_inc": "wilson",
}

ap = argparse.ArgumentParser()
ap.add_argument("--results", default="output/runs/external/team11_hsa.json")
ap.add_argument("--base", nargs="+", default=["AF alone", "AF + DepMap"])
ap.add_argument("--schukken", choices=("paper", "rank194", "drop"), default="paper")
ap.add_argument("--n", type=int, default=20000)
a = ap.parse_args()

rows = json.load(open(a.results))
by = defaultdict(dict)
for r in rows:
    by[r["arm"]][r["screen"]] = r

drop = {"paper": {"SCHUKKEN_PAIRED_WGL_rank194"},
        "rank194": {"SCHUKKEN_PAIRED_WGL_paper"},
        "drop": {"SCHUKKEN_PAIRED_WGL_paper", "SCHUKKEN_PAIRED_WGL_rank194"}}[a.schukken]
S = [s for s in sorted(by[a.base[0]]) if s not in drop]
papers = sorted({PAPER[s] for s in S})
print(f"{len(S)} screens over {len(papers)} papers  (schukken reading: {a.schukken})")
print(f"papers: {', '.join(papers)}\n")

def ef(r):
    # identical to canonical adjusted EF here; asserted, not assumed
    assert r["picks"] == 1000 and r["n_measured"] == 1000, (r["arm"], r["screen"])
    return r["ef"]

def boot(vals):
    cl = defaultdict(list)
    for s, v in vals.items():
        cl[PAPER[s]].append(v)
    k = list(cl)
    m = np.array([np.mean(cl[x]) for x in k])
    rng = np.random.default_rng(0)
    d = [m[rng.integers(0, len(k), len(k))].mean() for _ in range(a.n)]
    return float(m.mean()), float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))

H = sum(by[a.base[0]][s]["total_hits"] for s in S)
print(f"{'arm':24s} {'meanEF':>8s} {'medEF':>7s} {'hits':>7s} {'% avail':>8s} {'per-paper EF':>13s}")
for arm in sorted(by, key=lambda x: -st.mean(ef(by[x][s]) for s in S)):
    e = [ef(by[arm][s]) for s in S]
    h = sum(by[arm][s]["hits"] for s in S)
    pp = st.mean(st.mean(ef(by[arm][s]) for s in S if PAPER[s] == p) for p in papers)
    print(f"{arm:24s} {st.mean(e):8.3f} {st.median(e):7.3f} {h:7,d} {h/H:7.1%} {pp:13.3f}")

for b in a.base:
    if b not in by:
        continue
    print(f"\n== vs {b} " + "=" * 62)
    print(f"{'contrast':24s} {'dEF':>7s} {'screen 95% CI':>22s} {'paper 95% CI':>22s} "
          f"{'d_hits':>7s} {'w/l':>7s}")
    for arm in by:
        if arm == b:
            continue
        d = {s: ef(by[arm][s]) - ef(by[b][s]) for s in S}
        v = list(d.values())
        rng = np.random.default_rng(0)
        bs = [np.mean(rng.choice(v, len(v))) for _ in range(a.n)]
        slo, shi = np.percentile(bs, [2.5, 97.5])
        m, plo, phi = boot(d)
        dh = sum(by[arm][s]["hits"] - by[b][s]["hits"] for s in S)
        w = sum(1 for x in v if x > 0)
        star = "*" if plo > 0 or phi < 0 else " "
        print(f"{arm:24s} {st.mean(v):+7.3f} [{slo:+7.3f},{shi:+7.3f}] "
              f"[{plo:+7.3f},{phi:+7.3f}]{star} {dh:+7d} {w:3d}/{len(S):<3d}")
print("\n* = PAPER-level interval excludes zero. Seven papers is a small cluster count;")
print("  read the paper-level interval, and treat either as indicative.")
