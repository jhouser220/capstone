"""AF alone vs AF + DepMap on exactly the 11 screens in ``external_sources.json``.

The team's file names screens by internal experiment id (``E224``...``E245``) and
carries no gene list, only ``source_genes``, ``source_hits`` and a PMID/DOI. Each
entry was matched to a screen already built here by that pair of counts, which is
unambiguous: 10 of the 11 agree on hits exactly and on gene count to within 7.

``source_genes`` differs by a handful on three of them (Wilson 18458 vs 18464,
Nguyen 19111 vs 19114, Samelson 18818 vs 18825) because their duplicate-symbol
collapse differs from ours by 3-7 genes out of ~19 000. That is a 0.04% library
difference and cannot move EF meaningfully; it is reported, not patched.

E224 is the one entry that had no counterpart and had to be built:
``build_schukken_paired.py``, the paper's PAIRED aneuploid-vs-euploid contrast on
the quantile-normalised ``WGL`` sheet. **Its hit rule is not recoverable** -- see
that module -- so it is carried in both readings, and they disagree in sign.
"""
from __future__ import annotations
import json
from collections import defaultdict
import numpy as np

EXT = "output/runs/external"

# team id -> (results file, our screen name, paper cluster, note)
MAP = [
    ("E224_SCHUKKEN_PAIRED_WGL", "schukken_paired_results.json",
     "SCHUKKEN_PAIRED_WGL_paper", "schukken", "rule = paper's p<0.05 & beta<-0.2"),
    ("E224_SCHUKKEN_PAIRED_WGL", "schukken_paired_results.json",
     "SCHUKKEN_PAIRED_WGL_rank194", "schukken", "rule = top 194 by AneuploidRank"),
    ("E225_ABT_PARCB_DEPENDENCY", "abt_results.json",
     "ABT_PARCB_dependency", "abt", ""),
    ("E226_WILSON_MERCS_HIGH", "wilson_results.json",
     "WILSON_MERCS_inc", "wilson", ""),
    ("E232_NGUYEN_K562_ARAC_KO_SENSITIZATION", "nguyen_results.json",
     "NGUYEN_AraC_resistant", "nguyen", "their label says sensitization; 240 is the resistant arm"),
    ("E236_SHAN_HEPG2_LOVASTATIN", "shan_results.json",
     "SHAN_HepG2_lovastatin", "shan", ""),
    ("E236_SHAN_HELA_LOVASTATIN", "shan_results.json",
     "SHAN_HeLa_lovastatin", "shan", ""),
    ("E236_SHAN_HELA_LDLR_KO", "shan_results.json",
     "SHAN_HeLa_LDLRko", "shan", ""),
    ("E239_ORGANOID_COLO_031", "organoid_results.json",
     "ORGANOID_WTSI-COLIVM_031_minLib", "organoid", ""),
    ("E239_ORGANOID_OESO_005", "organoid_results.json",
     "ORGANOID_WTSI-OESO_005_minLib", "organoid", ""),
    ("E239_ORGANOID_OVAR_018", "organoid_results.json",
     "ORGANOID_WTSI-OVAR_018_minLib", "organoid", ""),
    ("E245_SAMELSON_T22_LOWERING", "tau_results.json",
     "TAU_CRISPRi_dec", "samelson", ""),
]

team = json.loads(open("external_sources.json").read())
cache: dict[str, dict] = {}


def row(f, screen):
    if f not in cache:
        cache[f] = {(r["arm"], r["screen"]): r for r in json.loads(open(f"{EXT}/{f}").read())}
    by = cache[f]
    return by[("AF alone", screen)], by[("AF + DepMap", screen)]


hdr = (f"{'team id':42s} {'our screen':32s} | {'genes':>6} {'H':>5} | "
       f"{'AF':>5} {'+DM':>5} {'d':>5} | {'AF EF':>7} {'+DM EF':>7} {'dEF':>7}")
print("AF alone  vs  AF + DepMap   on the 11 screens of external_sources.json")
print("(universe = library, b_eff = 1000, 10 rounds x 100 picks)\n")
print(hdr)
print("-" * len(hdr))

rows, notes = [], []
for tid, f, screen, clus, note in MAP:
    b, t = row(f, screen)
    want = team[tid]
    flags = []
    if want["source_hits"] != b["total_hits"]:
        flags.append(f"hits {want['source_hits']} vs {b['total_hits']}")
    if want["source_genes"] != b["n_genes"]:
        flags.append(f"genes {want['source_genes']} vs {b['n_genes']}")
    if flags or note:
        notes.append((tid, screen, "; ".join(flags + ([note] if note else []))))
    print(f"{tid:42s} {screen:32s} | {b['n_genes']:6d} {b['total_hits']:5d} | "
          f"{b['hits']:5d} {t['hits']:5d} {t['hits'] - b['hits']:+5d} | "
          f"{b['ef']:7.3f} {t['ef']:7.3f} {t['ef'] - b['ef']:+7.3f}"
          + ("   <-- see note" if flags or note else ""))
    rows.append(dict(tid=tid, clus=clus, d=t["ef"] - b["ef"],
                     hb=b["hits"], ht=t["hits"], H=b["total_hits"],
                     alt=screen.endswith("_rank194")))

print("-" * len(hdr))

# The two E224 readings are the same screen under two hit rules, so only one can
# enter a pooled number. The paper's own rule is the primary.
for label, keep in (("primary (E224 = paper's rule) ", lambda r: not r["alt"]),
                    ("alt     (E224 = rank194)      ", lambda r: r["tid"] != "E224_SCHUKKEN_PAIRED_WGL" or r["alt"])):
    sub = [r for r in rows if keep(r)]
    d = np.array([r["d"] for r in sub])
    g = defaultdict(list)
    for r in sub:
        g[r["clus"]].append(r["d"])
    fm = np.array([np.mean(v) for v in g.values()])
    rng = np.random.default_rng(0)
    dr = fm[rng.integers(0, len(fm), size=(20000, len(fm)))].mean(axis=1)
    lo, hi = np.percentile(dr, [2.5, 97.5])
    print(f"{label} n={len(sub):2d}  papers={len(fm)}  "
          f"hits/scr {np.mean([r['hb'] for r in sub]):6.1f} -> "
          f"{np.mean([r['ht'] for r in sub]):6.1f} "
          f"({np.mean([r['ht'] - r['hb'] for r in sub]):+5.1f})  "
          f"mean dEF {d.mean():+.3f}  paper-clustered 95% CI [{lo:+.3f}, {hi:+.3f}]"
          f"{'*' if (lo > 0 or hi < 0) else ''}  up {int((d > 0).sum())}/{len(d)}")

if notes:
    print("\nnotes")
    for tid, screen, n in notes:
        print(f"  {tid} / {screen}: {n}")
