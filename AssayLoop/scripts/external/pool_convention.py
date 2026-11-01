"""What the candidate-pool convention is worth, and the matched comparison.

Three conventions for the same 1000-pick budget:

  LIBRARY POOL (ours, default)   candidates = this screen's own library, so every
      pick has a measured outcome. 1000 picks, 1000 measured. No forgiveness
      terms, EF = hits / (1000 * H/D) with D = the library.
  FIXED POOL, CHARGED (theirs)   candidates = one cross-screen pool for every
      screen. A pick on a gene the screen never measured is revealed as a miss
      and still costs a budget slot, so 1000 picks buy only ~770-970 measured
      ones. `--universe` reproduces this exactly: GeneBatchTask already reveals
      off-library picks as ``hit=False, in_library=False``.
  FIXED POOL, FORGIVEN (AssayBench native)   same pool, but off-library picks
      advance neither axis. On the folds that forgives 562 of every 1000 picks
      (section 13.1), which is why the benchmark was re-run under the library
      pool before being compared with anything external.

**EF is not comparable across the first two.** ``hits_auc`` takes
``library_size = max(len(screen.genes), pool_size)`` (hits_auc.py:245), so
widening the pool widens the random baseline and EF can RISE while raw hits
FALL -- Abt goes 431 -> 349 hits but EF 4.458 -> 5.296. Compare **raw hits in
1000 picks**, or hits per 1000 measured picks. Not EF.

Their pool (20,893 genes) is not in their public package. The stand-in here is
``gene_universe(load_screens("validation"), min_screen_freq=2)`` = 21,213 genes,
which is the same kind of object, is independent of all seven external papers,
covers 98.7-100% of each library and leaves 98.3-100% of hits reachable.
"""
from __future__ import annotations
import json
from collections import defaultdict
import numpy as np

EXT = "output/runs/external"
TM = f"{EXT}/teammate"

OURS = {
    "E224_SCHUKKEN_PAIRED_WGL": ("schukken_paired_results.json", "SCHUKKEN_PAIRED_WGL_paper"),
    "E225_ABT_PARCB_DEPENDENCY": ("abt_results.json", "ABT_PARCB_dependency"),
    "E226_WILSON_MERCS_HIGH": ("wilson_results.json", "WILSON_MERCS_inc"),
    "E232_NGUYEN_K562_ARAC_KO_SENSITIZATION": ("nguyen_results.json", "NGUYEN_AraC_resistant"),
    "E236_SHAN_HEPG2_LOVASTATIN": ("shan_results.json", "SHAN_HepG2_lovastatin"),
    "E236_SHAN_HELA_LOVASTATIN": ("shan_results.json", "SHAN_HeLa_lovastatin"),
    "E236_SHAN_HELA_LDLR_KO": ("shan_results.json", "SHAN_HeLa_LDLRko"),
    "E239_ORGANOID_COLO_031": ("organoid_results.json", "ORGANOID_WTSI-COLIVM_031_minLib"),
    "E239_ORGANOID_OESO_005": ("organoid_results.json", "ORGANOID_WTSI-OESO_005_minLib"),
    "E239_ORGANOID_OVAR_018": ("organoid_results.json", "ORGANOID_WTSI-OVAR_018_minLib"),
    "E245_SAMELSON_T22_LOWERING": ("tau_results.json", "TAU_CRISPRi_dec"),
}

cache: dict[str, dict] = {}
def get(f, arm, screen):
    if f not in cache:
        cache[f] = {(x["arm"], x["screen"]): x for x in json.loads(open(f"{EXT}/{f}").read())}
    return cache[f][(arm, screen)]

camp = defaultdict(list)
for r in json.loads(open(f"{TM}/external_campaigns.json").read()):
    camp[(r["screen"], r["method"])].append(r)

FP = "team11_fixedpool_results.json"

hdr = (f"{'screen':38s} {'H':>5} | {'AF':>4} {'+DM':>4} {'d':>4} | "
       f"{'AF':>4} {'+DM':>4} {'d':>4} {'meas':>5} | "
       f"{'Nat':>6} {'HSA':>6} {'d':>6} {'meas':>5}")
print("                                             library pool  |  fixed pool (charged) |  theirs: fixed pool, LLM")
print(hdr); print("-" * len(hdr))

rows = []
for tid, (f, screen) in OURS.items():
    lb, lt = get(f, "AF alone", screen), get(f, "AF + DepMap", screen)
    pb, pt = get(FP, "AF alone", screen), get(FP, "AF + DepMap", screen)
    nat = camp[(tid, "Native")]; hsa = camp[(tid, "HSA")]
    nh = np.mean([x["hits"] for x in nat]); nm = np.mean([x["measured"] for x in nat])
    hh = np.mean([x["hits"] for x in hsa]); hm = np.mean([x["measured"] for x in hsa])
    print(f"{tid.replace('_', ' ')[:38]:38s} {lb['total_hits']:5d} | "
          f"{lb['hits']:4d} {lt['hits']:4d} {lt['hits'] - lb['hits']:+4d} | "
          f"{pb['hits']:4d} {pt['hits']:4d} {pt['hits'] - pb['hits']:+4d} "
          f"{(pb['n_measured'] + pt['n_measured']) // 2:5d} | "
          f"{nh:6.2f} {hh:6.2f} {hh - nh:+6.2f} {(nm + hm) / 2:5.0f}")
    rows.append([lb["hits"], lt["hits"], pb["hits"], pt["hits"],
                 pb["n_measured"], pt["n_measured"], nh, hh, nm, hm])

a = np.array(rows, dtype=float)
print("-" * len(hdr))
print(f"{'MEAN':38s} {'':5s} | "
      f"{a[:,0].mean():4.0f} {a[:,1].mean():4.0f} {(a[:,1]-a[:,0]).mean():+4.1f} | "
      f"{a[:,2].mean():4.0f} {a[:,3].mean():4.0f} {(a[:,3]-a[:,2]).mean():+4.1f} "
      f"{(a[:,4]+a[:,5]).mean()/2:5.0f} | "
      f"{a[:,6].mean():6.2f} {a[:,7].mean():6.2f} {(a[:,7]-a[:,6]).mean():+6.2f} "
      f"{(a[:,8]+a[:,9]).mean()/2:5.0f}")

print(f"""
what the convention alone costs, AF alone, same model and same seed:
  library pool  {a[:,0].mean():7.1f} hits per screen from 1000 measured picks
  fixed pool    {a[:,2].mean():7.1f} hits per screen from {(a[:,4]).mean():.0f} measured picks
                {100*(a[:,2].mean()/a[:,0].mean()-1):+.1f}% in raw hits, purely from pool accounting

sign of the delta, by convention:
  AF -> AF + DepMap, library pool   up {int((a[:,1]>a[:,0]).sum())}/11   mean {(a[:,1]-a[:,0]).mean():+.1f} hits
  AF -> AF + DepMap, fixed pool     up {int((a[:,3]>a[:,2]).sum())}/11   mean {(a[:,3]-a[:,2]).mean():+.1f} hits
  Native -> HSA,     their pool     up {int((a[:,7]>a[:,6]).sum())}/11   mean {(a[:,7]-a[:,6]).mean():+.1f} hits

matched-pool comparison (all three now on a fixed ~21k pool with unknowns charged):
  their Native  {a[:,6].mean():7.1f}   their HSA {a[:,7].mean():7.1f}
  our AF alone  {a[:,2].mean():7.1f}   our +DM   {a[:,3].mean():7.1f}
  AF alone beats their HSA on {int((a[:,2]>a[:,7]).sum())}/11 screens, their Native on {int((a[:,2]>a[:,6]).sum())}/11
""")
