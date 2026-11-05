"""The four own-library runs side by side, with AF + DepMap shown as a level, not a baseline.

`score_libpool.py` prints contrasts AGAINST `AF + DepMap`, so that arm appears only as the
zero line. This prints absolute EF and hits for every arm on every set, plus what each arm's
absolute EF did when the f2 refund was removed -- which is where AF + DepMap's own movement
lives (it is the arm the refund was helping most).

EF here is the one unambiguous own-library number: stored ef == canonical adjusted EF ==
hits/(1000*H/D), because no pick can land off-library (section 4.9). The f2 column is CANONICAL
adjusted EF on the same screens. Absolute values are not comparable ACROSS the two pools --
different denominator, section 7 rule 4 -- which is exactly why the `moved` column is the
quantity of interest rather than either level.
"""
import json
import statistics as st
from collections import defaultdict

ARMS = ["AF alone", "AF + DepMap", "HSA only (ours)",
        "DepMap + HSA dm.25", "DepMap + HSA dm.50"]
SHORT = {"AF alone": "AF", "AF + DepMap": "AF + DepMap", "HSA only (ours)": "HSA",
         "DepMap + HSA dm.25": "DepMap + HSA .25", "DepMap + HSA dm.50": "DepMap + HSA .50"}

SETS = {
    # REAL-MEAN FILE FIRST: `load` keeps the first file an arm appears in, and the
    # own-library runs use `banks_real_*`. Listing the mean=0 file first silently
    # compares a real-mean library run against a mean=0 f2 run for `HSA only` (and,
    # on paper_validation, for dm.25/dm.50 too, which the mean=0 file also carries).
    "paper_test": ("output/runs/external/paper_test_hsa_libpool.json",
                   ["output/runs/external/paper_test_hsa_realmean.json",
                    "output/runs/external/paper_test_hsa.json"]),
    "paper_validation": ("output/runs/external/paper_validation_hsa_libpool.json",
                         ["output/runs/external/paper_validation_hsa_realmean.json",
                          "output/runs/external/paper_validation_hsa.json"]),
    "val fold, genome-scale": ("output/runs/external/validation_fold_hsa_libpool.json",
                               ["output/runs/external/validation_fold_hsa_REALCOMBINED.json"]),
    # no f2 HSA run exists on the test fold, so its "f2" source is itself: the
    # refund-removal columns are zero by construction and carry no information
    "TEST fold, genome-scale": ("output/runs/external/test_fold_hsa_libpool.json",
                                ["output/runs/external/test_fold_hsa_libpool.json"]),
}


def load(fs):
    by = defaultdict(dict)
    for f in fs:
        for r in json.load(open(f)):
            by[r["arm"]].setdefault(r["screen"], r)
    return by


def canon(r):
    return r["hits"] / ((1000 - r["n_unknown"]) * r["total_hits"] / r["n_genes"])


DATA = {}
for name, (lib, f2) in SETS.items():
    L, F = load([lib]), load(f2)
    S = sorted(L["AF + DepMap"])
    DATA[name] = (L, F, S)

print("=== OWN-LIBRARY POOL: absolute EF / total hits (no refund, one metric) ===\n")
hdr = f"{'arm':20s}" + "".join(f"{n:>26s}" for n in SETS)
print(hdr)
print(f"{'':20s}" + "".join(f"{'EF':>9s}{'hits':>9s}{'rate':>8s}" for _ in SETS))
for m in ARMS:
    row = f"{SHORT[m]:20s}"
    for name in SETS:
        L, _, S = DATA[name]
        ef = st.mean(L[m][s]["ef"] for s in S)
        h = sum(L[m][s]["hits"] for s in S)
        row += f"{ef:9.3f}{h:9,d}{h / (len(S) * 1000):8.2%}"
    print(row)
for name in SETS:
    L, _, S = DATA[name]
    print(f"  {name}: {len(S)} screens, {sum(L['AF + DepMap'][s]['total_hits'] for s in S):,} "
          f"hits available")

print("\n=== WHAT EACH ARM'S ABSOLUTE EF DID WHEN THE REFUND WAS REMOVED ===")
print("(canonical f2 EF -> own-library EF, same screens; the f2 level is the inflated one)\n")
for name in SETS:
    L, F, S = DATA[name]
    print(f"{name}  ({len(S)} screens)")
    print(f"  {'arm':20s} {'f2 EF':>8s} {'lib EF':>8s} {'moved':>7s} "
          f"{'off-lib/1k':>11s} {'f2 hits':>9s} {'lib hits':>9s} {'d_hits':>7s}")
    for m in ARMS:
        f2 = st.mean(canon(F[m][s]) for s in S)
        lb = st.mean(L[m][s]["ef"] for s in S)
        off = st.mean(F[m][s]["n_unknown"] for s in S)
        fh = sum(F[m][s]["hits"] for s in S)
        lh = sum(L[m][s]["hits"] for s in S)
        print(f"  {SHORT[m]:20s} {f2:8.3f} {lb:8.3f} {lb - f2:+7.3f} {off:11.1f} "
              f"{fh:9,d} {lh:9,d} {lh - fh:+7d}")
    print()
print("The split is by FAMILY, not strictly by off-library count: on paper_test and the fold the")
print("two AF-based arms lose 0.19-0.28 EF while the three HSA arms lose <=0.11 or gain, matching")
print("their ~60-vs-~30 off-library picks. Within a family the order does not always follow the")
print("count (fold: AF 62.4 picks / -0.185 against AF + DepMap 58.3 / -0.256), and on")
print("paper_validation every arm moves -0.02 to -0.06 with no ordering at all -- there the f2")
print("pool is only 20,893 against a 19,283-gene mean library, so there is barely a refund to")
print("remove. Treat the family split as the finding and the per-arm ordering as noise.")
