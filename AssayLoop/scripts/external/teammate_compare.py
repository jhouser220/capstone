"""Our AF vs AF + DepMap next to the teammate's Native vs HSA, on the same 11 screens.

The two experiments are NOT the same comparison, and the numbers must not be read
as one leaderboard. Their records (vendored under
``output/runs/external/teammate/``, from github.com/jhouser220/capstone branch
``Odysseus``, ``experiments/odysseus/assayloop_comparison/recorded_results``):

  arms       ``Native`` = the original AssayLoop loop, an LLM proposes and
             pre-batch AssayFormer fills short batches.
             ``HSA`` = Historical Subspace Adaptation, their method: 32 historical
             response directions appended to the source's 10 features, task
             coefficients fitted on measured labels, ranked by availability times
             integrated hit probability.
  actors     google.gemma-4-31b and qwen.qwen3-32b, 2 draws each -> 4 campaigns
             per screen per arm. Both arms are LLM-in-the-loop.
  pool       a fixed 20 893-gene candidate pool, NOT the screen's library, and
             picks with no measured outcome are charged to the budget. So of
             1000 accepted picks only ~770-950 land on a measured gene.

Ours: no LLM anywhere, ``universe_genes`` = the screen's own library, so all 1000
picks are measured. ``hits/1k meas`` rescales both sides to hits per 1000
MEASURED picks, which is the closest thing to a common footing; it still does not
equalise the proposal mechanism, the budget semantics or the baseline.

The honest reading: their intervention is tested against a weaker baseline
(LLM-driven AssayLoop) than ours (AssayFormer's own ranking), and on several of
these screens our AF-alone control already finds more hits than their HSA arm.
Their Native-to-HSA delta is positive on all 7 publications; our AF-to-AF+DepMap
delta on the same screens is ~0. Both can be true at once.
"""
from __future__ import annotations
import json
from collections import defaultdict
import numpy as np

EXT = "output/runs/external"
TM = f"{EXT}/teammate"

# team screen id -> (our results file, our screen name)
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

camp = defaultdict(list)
for r in json.loads(open(f"{TM}/external_campaigns.json").read()):
    camp[(r["screen"], r["method"])].append(r)

cache: dict[str, dict] = {}
def ours(f, screen, arm):
    if f not in cache:
        cache[f] = {(x["arm"], x["screen"]): x for x in json.loads(open(f"{EXT}/{f}").read())}
    return cache[f][(arm, screen)]

hdr = (f"{'screen':40s} {'H':>5} | {'Native':>7} {'HSA':>7} {'d':>6} {'meas':>6} | "
       f"{'AF':>5} {'+DM':>5} {'d':>5} | {'Nat/1k':>7} {'HSA/1k':>7} {'AF/1k':>6}")
print(__doc__.strip().splitlines()[0])
print()
print(hdr); print("-" * len(hdr))

rows = []
for tid, (f, screen) in OURS.items():
    nat = camp[(tid, "Native")]; hsa = camp[(tid, "HSA")]
    nh = np.mean([x["hits"] for x in nat]); nm = np.mean([x["measured"] for x in nat])
    hh = np.mean([x["hits"] for x in hsa]); hm = np.mean([x["measured"] for x in hsa])
    b = ours(f, screen, "AF alone"); t = ours(f, screen, "AF + DepMap")
    print(f"{tid.replace('_', ' ')[:40]:40s} {b['total_hits']:5d} | "
          f"{nh:7.2f} {hh:7.2f} {hh - nh:+6.2f} {(nm + hm) / 2:6.0f} | "
          f"{b['hits']:5d} {t['hits']:5d} {t['hits'] - b['hits']:+5d} | "
          f"{1000 * nh / nm:7.1f} {1000 * hh / hm:7.1f} {b['hits']:6d}")
    rows.append((nh, hh, b["hits"], t["hits"], nm, hm))

print("-" * len(hdr))
a = np.array(rows)
print(f"{'MEAN over the 11 screens':40s} {'':5s} | "
      f"{a[:, 0].mean():7.2f} {a[:, 1].mean():7.2f} {(a[:, 1] - a[:, 0]).mean():+6.2f} "
      f"{(a[:, 4] + a[:, 5]).mean() / 2:6.0f} | "
      f"{a[:, 2].mean():5.1f} {a[:, 3].mean():5.1f} {(a[:, 3] - a[:, 2]).mean():+5.1f} | "
      f"{1000 * a[:, 0].sum() / a[:, 4].sum():7.1f} {1000 * a[:, 1].sum() / a[:, 5].sum():7.1f} "
      f"{a[:, 2].mean():6.1f}")
print(f"\narms up: HSA over Native {int((a[:, 1] > a[:, 0]).sum())}/11   "
      f"AF+DepMap over AF alone {int((a[:, 3] > a[:, 2]).sum())}/11")
print(f"AF alone beats HSA on {int((a[:, 2] > a[:, 1]).sum())}/11 screens "
      f"(raw hits, unequal budgets -- see module docstring)")
print("\nE224 is shown with the paper's own hit rule (188 hits). Their 194-hit "
      "definition is not in the public package.")
