"""Paired scoring of the paper_test HSA arms against AF alone."""
import json, statistics as st, sys
from pathlib import Path

srcs = sys.argv[1:] or ["output/runs/external/paper_test_hsa.json"]
rows = [r for s in srcs for r in json.load(open(s))]
by = {}
for r in rows:
    by.setdefault(r["arm"], {})[r["screen"]] = r
arms = list(by)
base = "AF alone"
screens = sorted(by[base])

print(f"{'screen':<16}" + "".join(f"{a[:14]:>15}" for a in arms))
print(f"{'':<16}" + "".join(f"{'EF':>15}" for _ in arms))
for s in screens:
    print(f"{s:<16}" + "".join(f"{by[a][s]['ef']:>15.3f}" if s in by[a] else f"{'-':>15}"
                               for a in arms))
print()
print(f"{'arm':<26}{'mean EF':>10}{'median':>10}{'hits':>8}{'dEF vs AF':>12}"
      f"{'median d':>11}{'wins':>7}{'losses':>8}")
for a in arms:
    efs = [by[a][s]["ef"] for s in screens if s in by[a]]
    hits = sum(by[a][s]["hits"] for s in screens if s in by[a])
    d = [by[a][s]["ef"] - by[base][s]["ef"] for s in screens if s in by[a]]
    w = sum(1 for x in d if x > 1e-9); l = sum(1 for x in d if x < -1e-9)
    print(f"{a:<26}{st.mean(efs):>10.3f}{st.median(efs):>10.3f}{hits:>8}"
          f"{st.mean(d):>+12.3f}{st.median(d):>+11.3f}{w:>7}{l:>8}")
print(f"\nn = {len(screens)} screens   pool = {rows[0]['pool']} genes (f2)")
print("EF here is the ADJUSTED EF: off-library picks are forgiven, so absolute levels")
print("are inflated relative to a library-pool run. Paired deltas are the quotable part.")
