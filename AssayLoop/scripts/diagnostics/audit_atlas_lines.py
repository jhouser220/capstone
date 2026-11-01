"""Does the ablated atlas really contain no evaluation screen's own cell line?

`build_noeval_atlas.py` drops every Achilles line that `screen_depmap.resolve` matches to
an evaluation screen's `cell_line`. That is a NAME match, so the claim it supports -- "no
screen's own cell line is in the atlas by construction" -- is only as good as the
resolver. This audits it from the other direction: take every screen whose cell line does
NOT resolve, and ask whether DepMap nonetheless holds that line under a name the resolver
cannot reach.

Three passes, cheapest first:
  1. resolved / resolved-but-not-in-Achilles / unresolved counts
  2. strip engineering decorations (ACE2, SunTag, Cas9, hTERT prefix, trailing clone
     digits) from each unresolved name and retry against lines that SURVIVE the ablation
  3. prefix search: any surviving Achilles line whose name extends the screen's, which is
     how clonal subclones are named (hTERT-RPE1 -> RPE1-ss48)

Pass 3 is the one that found the real gap; passes 1-2 came back clean.

    python scripts/diagnostics/audit_atlas_lines.py --depmap-dir ~/depmap
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "src")
sys.path.insert(0, "scripts/diagnostics")
from assayloop.tasks import load_screens                       # noqa: E402
from screen_depmap import line_map, norm, resolve              # noqa: E402

SETS = ("validation", "test")
# decorations that make an engineered derivative's name miss its parent
SUFFIXES = ["ACE2", "TMPRSS2", "SUNTAG", "CAS9", "GFP110", "GFP", "DELTASTAT1",
            "IR", "CELL", "CELLS", "KO", "SCLYDEFICIENT", "RELAMNEONGREEN"]


def variants(k: str) -> set[str]:
    out = {k}
    if k.startswith("HTERT"):
        out.add(k[5:])
    for s in SUFFIXES:
        for base in list(out):
            if base.endswith(s) and len(base) > len(s) + 2:
                out.add(base[: -len(s)])
    for base in list(out):          # HUH751 -> HUH75 -> HUH7
        m = base
        while re.search(r"\d$", m) and len(m) > 3:
            m = m[:-1]
            out.add(m)
    return {x for x in out if len(x) >= 3}


def main() -> None:
    import logging
    logging.disable(logging.CRITICAL)
    ap = argparse.ArgumentParser()
    ap.add_argument("--depmap-dir", default="/home/jeff/depmap")
    ap.add_argument("--atlas", default=None, help="default: <depmap-dir>/achilles_noeval.npz")
    a = ap.parse_args()
    D = Path(a.depmap_dir)
    atlas = Path(a.atlas) if a.atlas else D / "achilles_noeval.npz"

    by = line_map(D / "Model.csv")
    full = set(np.load(D / "achilles.npz", allow_pickle=True)["rows"])
    kept = set(np.load(atlas, allow_pickle=True)["rows"])
    mod = pd.read_csv(D / "Model.csv", low_memory=False)
    name_of = {r["ModelID"]: r["StrippedCellLineName"] for _, r in mod.iterrows()}
    # normalised name -> ModelID, for lines that SURVIVE the ablation
    surviving: dict[str, str] = {}
    for _, r in mod.iterrows():
        if r["ModelID"] in kept:
            for c in ("StrippedCellLineName", "CellLineName"):
                k = norm(r[c])
                if k:
                    surviving.setdefault(k, r["ModelID"])

    screens, lines = {}, Counter()
    for sset in SETS:
        for s in load_screens(target_set=sset):
            screens.setdefault(s.dataset_name, (sset, s.cell_line))
    for sset, cl in screens.values():
        lines[cl] += 1
    print(f"atlas {atlas.name}: {len(kept)} of {len(full)} Achilles lines kept")
    print(f"{len(screens)} unique screens over {SETS}, {len(lines)} distinct cell_line strings\n")

    unresolved = Counter()
    n_in = n_out = 0
    for cl, n in lines.items():
        mid = resolve(cl, by)
        if mid is None:
            unresolved[cl] = n
        elif mid in full:
            n_in += n
        else:
            n_out += n
    print(f"  pass 1  resolve -> in Achilles      {n_in:4d} screens (their lines were dropped)")
    print(f"          resolve -> not in Achilles  {n_out:4d} screens (nothing to drop)")
    print(f"          DO NOT resolve              {sum(unresolved.values()):4d} screens  <- audited below\n")

    bad: dict[str, list] = defaultdict(list)
    for cl, n in unresolved.items():
        for cand in sorted(variants(norm(cl)), key=len, reverse=True):
            if cand in surviving:
                bad["pass 2: stripped form matches a surviving line"].append(
                    (cl, n, cand, surviving[cand]))
                break
    for cl, n in unresolved.items():
        k = norm(cl)
        base = k[5:] if k.startswith("HTERT") else k
        if len(base) < 4:
            continue
        hits = sorted({m for kk, m in surviving.items() if kk.startswith(base) and kk != base})
        if hits:
            bad["pass 3: a surviving line's name EXTENDS this one (clonal subclone)"].append(
                (cl, n, base, hits))

    if not bad:
        print("CLEAN: no unresolved evaluation cell line has a counterpart surviving in the atlas.")
        return
    total = 0
    for k, rows in bad.items():
        print(f"!! {k}")
        for cl, n, form, mids in rows:
            mids = mids if isinstance(mids, list) else [mids]
            total += n
            print(f"     {str(cl)[:28]:28s} {n:3d} screens  via {form:14s} -> "
                  + ", ".join(f"{m} {name_of.get(m,'?')}" for m in mids))
    print(f"\n{total} screens are NOT protected by construction. "
          f"Add their ModelIDs to the builder's drop set and rebuild.")
    raise SystemExit(1)


if __name__ == "__main__":
    main()
