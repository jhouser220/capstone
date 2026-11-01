"""Score `min_seeds=1` against the deployed arm, split by whether the feature could
fire at all under the default.

§4.23 measured that the deployed configuration never fires on 42 of 198 fold screens
because they never reach 2 confirmed hits (mean total hits 11, AssayFormer early hit
rate 0.8%). Lowering `min_seeds` to 1 is the only change that enlarges the method's
REACH rather than adjusting its dose, so the question is not just "is the overall
number better" but "did the gain come from the screens that previously got nothing".

Those two can come apart in both directions, which is why the split matters:
  * a positive overall number driven by the ALREADY-FIRING screens would mean the
    one-seed direction helps where a two-seed one already worked -- interesting, but
    not the hypothesis, and it would mean `min_seeds` is not what changed anything;
  * a positive number on the never-fired group with a negative overall would argue
    for firing at 1 seed ONLY when 2 are unreachable, a different rule.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "src")
sys.path.insert(0, "scripts/diagnostics")
from assayloop.scripts.full_genome_table import _adj_ef_from_run   # noqa: E402
from assayloop.tasks import gene_universe, load_screens            # noqa: E402
from score_ctx import boot_ci, families                            # noqa: E402
from run_paths import run_dir  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--screen-set", default="validation")
    ap.add_argument("--exclude-set", default="paper_validation")
    ap.add_argument("--dir", default="output/runs")
    ap.add_argument("--baseline", default="sweep-dm-validation-w0")
    ap.add_argument("--deployed", default="sweep-dm-validation-dm05knn")
    ap.add_argument("--arm", default="sweep-ms1-validation")
    a = ap.parse_args()

    screens = load_screens(target_set=a.screen_set)
    uni = set(gene_universe(screens, min_screen_freq=2))
    drop = {x.dataset_name for x in load_screens(target_set=a.exclude_set)} \
        if a.exclude_set else set()

    rows: dict[str, dict] = {}
    for i, s in enumerate(screens):
        if s.dataset_name in drop:
            continue
        p = {k: run_dir(a.dir, k, i, s.dataset_name) / "result.json"
             for k in (a.baseline, a.deployed, a.arm)}
        if not all(v.is_file() for v in p.values()):
            continue
        lib = set(s.genes)
        ef = {k: _adj_ef_from_run(v, lib, uni, s.total_hits, len(lib), 1000)
              for k, v in p.items()}
        # did the DEPLOYED arm (min_seeds=2) ever fire on this screen?
        tr = [st.get("acquisition_trace") or {}
              for st in json.loads(p[a.deployed].read_text())["steps"]]
        fired2 = sum(1 for t in tr if t.get("codep_active") or t.get("knn_active"))
        tr1 = [st.get("acquisition_trace") or {}
               for st in json.loads(p[a.arm].read_text())["steps"]]
        fired1 = sum(1 for t in tr1 if t.get("codep_active") or t.get("knn_active"))
        rows[s.dataset_name] = dict(
            vs_base=ef[a.arm] - ef[a.baseline],
            dep_vs_base=ef[a.deployed] - ef[a.baseline],
            vs_dep=ef[a.arm] - ef[a.deployed],
            fired2=fired2, fired1=fired1, total=s.total_hits)

    kept = [s for s in screens if s.dataset_name in rows]
    fam = families(kept)
    print(f"{len(rows)} fresh screens, {len(set(fam.values()))} families\n")

    f1 = np.array([v["fired1"] for v in rows.values()])
    f2 = np.array([v["fired2"] for v in rows.values()])
    print(f"rounds the feature fired, mean of 10:  min_seeds=2 {f2.mean():.2f}  ->  "
          f"min_seeds=1 {f1.mean():.2f}")
    print(f"screens that never fired:              min_seeds=2 {(f2 == 0).sum()}  ->  "
          f"min_seeds=1 {(f1 == 0).sum()}\n")

    def show(label, sel):
        if len(sel) < 3:
            print(f"{label:34s} n={len(sel):3d}  (too few)")
            return
        f = {k: fam[k] for k in sel}
        m1, (l1, h1), _ = boot_ci({k: rows[k]["vs_base"] for k in sel}, f)
        m0, _, _ = boot_ci({k: rows[k]["dep_vs_base"] for k in sel}, f)
        md, (ld, hd), _ = boot_ci({k: rows[k]["vs_dep"] for k in sel}, f)
        up = sum(1 for k in sel if rows[k]["vs_dep"] > 1e-9)
        print(f"{label:34s} n={len(sel):3d}  deployed {m0:+.3f}  ms1 {m1:+.3f}  "
              f"diff {md:+.3f} [{ld:+.3f}, {hd:+.3f}]  better on {up}/{len(sel)}")

    print("per-FAMILY means. 'diff' is min_seeds=1 minus the deployed arm:")
    show("ALL screens", list(rows))
    show("never fired at min_seeds=2", [k for k, v in rows.items() if v["fired2"] == 0])
    show("did fire at min_seeds=2", [k for k, v in rows.items() if v["fired2"] > 0])
    show("  of those: <= 50 total hits",
         [k for k, v in rows.items() if v["fired2"] > 0 and v["total"] <= 50])
    show("  of those: > 50 total hits",
         [k for k, v in rows.items() if v["fired2"] > 0 and v["total"] > 50])

    newly = [k for k, v in rows.items() if v["fired2"] == 0 and v["fired1"] > 0]
    print(f"\nscreens NEWLY activated by min_seeds=1: {len(newly)}")
    if newly:
        d = [rows[k]["vs_dep"] for k in newly]
        print(f"  their mean dEF change {np.mean(d):+.3f}, "
              f"{sum(1 for x in d if x > 0)} up / {sum(1 for x in d if x < 0)} down, "
              f"mean total hits {np.mean([rows[k]['total'] for k in newly]):.0f}")


if __name__ == "__main__":
    main()
