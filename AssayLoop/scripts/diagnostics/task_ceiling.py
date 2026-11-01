"""How much room is there in the TASK, as opposed to inside AssayFormer's û?

Two different ceilings get confused easily:

  oracle-û ceiling   the best ranking expressible as a(g) = û·v_g + b_g for ANY
                     û. Bounds every intervention that works by improving the
                     screen latent. Measured in oracle_u.py: +0.41 / +0.84.

  perfect-play ceiling   rank every true hit first. Bounds nothing about the
                     architecture; it is what the benchmark itself allows.
                     With a budget of 1000 picks out of a D-gene library holding
                     H hits, random expects 1000H/D hits, so
                         EF_max = min(1000, H) / (1000H/D)
                     which is D/1000 (about 19) whenever H <= 1000.

The gap between the two is the part of the task that is unreachable *in this
representation* rather than unreachable in principle.
"""
from __future__ import annotations

import argparse, json, sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "src")
from assayloop.scripts.full_genome_table import _adj_ef_from_run   # noqa: E402
from assayloop.tasks import gene_universe, load_screens            # noqa: E402
from run_paths import run_dir  # noqa: E402

BUDGET = 1000


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--screen-set", default="paper_validation")
    ap.add_argument("--sweep", default="sweep-fd2017be")
    ap.add_argument("--dir", default="output/runs/baseline/gemini")
    a = ap.parse_args()

    screens = load_screens(target_set=a.screen_set)
    uni = set(gene_universe(screens, min_screen_freq=2))
    rows = []
    for i, s in enumerate(screens):
        p = run_dir(a.dir, a.sweep, i, s.dataset_name) / "result.json"
        if not p.is_file():
            continue
        lib, H, D = set(s.genes), s.total_hits, len(s.genes)
        ef = _adj_ef_from_run(p, lib, uni, H, D, BUDGET)
        found = sum(h for st in json.loads(p.read_text())["steps"] for h in (st["hits"] or []))
        rand = BUDGET * H / D
        ef_max = min(BUDGET, H) / rand
        rows.append(dict(name=s.dataset_name, H=H, D=D, found=found, ef=ef,
                         ef_max=ef_max, frac_ceiling=ef / ef_max,
                         frac_hits=found / H, reachable=min(BUDGET, H)))

    rows.sort(key=lambda r: r["ef"])
    hdr = (f"{'screen':16s} {'hits':>6s} {'found':>6s} {'% of hits':>9s} "
           f"{'AF EF':>7s} {'perfect EF':>10s} {'% of ceiling':>12s}")
    print(f"{a.screen_set}: budget {BUDGET}\n")
    print(hdr); print("-" * len(hdr))
    for r in rows:
        print(f"{r['name']:16s} {r['H']:6d} {r['found']:6d} {100*r['frac_hits']:8.1f}% "
              f"{r['ef']:7.2f} {r['ef_max']:10.2f} {100*r['frac_ceiling']:11.1f}%")
    A = lambda k: np.array([r[k] for r in rows])
    print(f"\nmean   AF EF {A('ef').mean():.2f}   perfect-play EF {A('ef_max').mean():.2f}   "
          f"AF reaches {100*A('frac_ceiling').mean():.1f}% of the task ceiling")
    hard = [r for r in rows if r["ef"] < 5]
    easy = [r for r in rows if r["ef"] >= 5]
    for lab, grp in (("EF < 5 (hard)", hard), ("EF >= 5 (easy)", easy)):
        if not grp:
            continue
        g = lambda k: np.array([r[k] for r in grp])
        print(f"  {lab:16s} n={len(grp):2d}  AF {g('ef').mean():5.2f}  "
              f"perfect {g('ef_max').mean():5.2f}  "
              f"-> {100*g('frac_ceiling').mean():4.1f}% of ceiling, "
              f"{100*g('frac_hits').mean():4.1f}% of hits found")


if __name__ == "__main__":
    main()
