"""Paired bootstrap over screens for the paper_test HSA arms (HANDOFF 24.11)."""
import json, sys
import numpy as np

srcs = sys.argv[1:] or ["output/runs/external/paper_test_hsa.json"]
rows = [r for s in srcs for r in json.load(open(s))]
by = {}
for r in rows:
    by.setdefault(r["arm"], {})[r["screen"]] = r
S = sorted(by["AF alone"])
rng = np.random.default_rng(0)

def boot(d, n=20000):
    d = np.asarray(d, float)
    m = d[rng.integers(0, len(d), (n, len(d)))].mean(1)
    return d.mean(), np.percentile(m, 2.5), np.percentile(m, 97.5)

BASES = ["AF alone", "AF + DepMap", "HSA only (ours)"]
PAIRS = [(a, b) for b in BASES for a in sorted(by) if a != b and b in by
         and all(s in by[a] for s in S)]
last = None
print(f"{'contrast':<42}{'mean dEF':>10}{'95% CI':>22}{'d_hits':>9}{'wins/losses':>13}")
for a, b in PAIRS:
    if b != last:
        print(f"-- vs {b} " + "-"*(80-len(b)))
        last = b
    d = [by[a][s]["ef"] - by[b][s]["ef"] for s in S]
    dh = sum(by[a][s]["hits"] - by[b][s]["hits"] for s in S)
    m, lo, hi = boot(d)
    w = sum(1 for x in d if x > 1e-9); l = sum(1 for x in d if x < -1e-9)
    star = "*" if (lo > 0 or hi < 0) else " "
    print(f"{a+' - '+b:<42}{m:>+10.3f}  [{lo:>+7.3f},{hi:>+7.3f}]{star}{dh:>+9}{f'{w}/{l}':>13}")
print("\n* = 95% interval excludes zero.  Bootstrap is over the 20 screens, not hit-set")
print("families: paper_test is already a diversity-stratified subset.")
