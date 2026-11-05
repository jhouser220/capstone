"""mean=0 vs the reconstructed REAL mean, same arms, same screens (HANDOFF 24.21)."""
import argparse, json
from collections import defaultdict
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--zero", nargs="+", required=True, help="result file(s) built with mean=0")
ap.add_argument("--real", nargs="+", required=True, help="result file(s) built with the real mean")
ap.add_argument("--label", default="")
a = ap.parse_args()

def load(fs):
    by = defaultdict(dict)
    for f in fs:
        for r in json.load(open(f)):
            by[r["arm"]][r["screen"]] = r
    return by

Z, R = load(a.zero), load(a.real)
# the mean=0 paper_test sweep used different LABELS for identical weight settings
ALIAS = {"DepMap + HSA dm.50": "DepMap + HSA (no AF term)",   # dm=0.5, hsa=1, af=0
         "DepMap + HSA dm1.0": "HSA + DepMap eq"}             # dm=1.0, hsa=1, af=0
for new_name, old_name in ALIAS.items():
    if new_name in R and new_name not in Z and old_name in Z:
        Z[new_name] = Z[old_name]
        print(f"note: comparing real '{new_name}' against mean=0 '{old_name}' "
              f"(identical weights, different label)")
arms = [x for x in R if x in Z]
print(f"### {a.label}   arms re-run with the real mean: {len(arms)}")
print(f"{'arm':<24}{'screens':>8}{'hits(0)':>9}{'hits(real)':>11}{'dhits':>7}"
      f"{'EF(0)':>8}{'EF(real)':>9}{'dEF mean':>10}{'dEF med':>9}{'identical':>11}")
for arm in arms:
    S = sorted(set(Z[arm]) & set(R[arm]))
    hz = sum(Z[arm][s]["hits"] for s in S); hr = sum(R[arm][s]["hits"] for s in S)
    ez = np.array([Z[arm][s]["ef"] for s in S]); er = np.array([R[arm][s]["ef"] for s in S])
    d = er - ez
    same = sum(1 for s in S if Z[arm][s]["picks_sha"] == R[arm][s]["picks_sha"])
    print(f"{arm:<24}{len(S):>8}{hz:>9,}{hr:>11,}{hr-hz:>+7}"
          f"{ez.mean():>8.3f}{er.mean():>9.3f}{d.mean():>+10.4f}{np.median(d):>+9.4f}"
          f"{f'{same}/{len(S)}':>11}")
