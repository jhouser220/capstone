"""What would a different min_seeds / hit-rate gate have suppressed?

min_seeds cannot be scored exactly from finished runs: the campaign is sequential,
so round 5 of the feature arm depends on rounds 1-4 of the feature arm, and you
cannot splice rounds from two different campaigns. But the acquisition records
`n_seeds` per round, so you CAN read off exactly which rounds a higher threshold
would have suppressed, and check whether those are the rounds on the screens that
lost. That predicts the rerun's outcome for free.

Both candidate gates use ONLY revealed labels, so neither leaks:
  n_seeds        how many genes have come back as hits so far
  observed rate  hits / genes tested so far
"""
from __future__ import annotations

import argparse, json, sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "src")
from assayloop.scripts.full_genome_table import _adj_ef_from_run   # noqa: E402
from assayloop.tasks import gene_universe, load_screens            # noqa: E402
from run_paths import run_dir  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--screen-set", default="validation")
    ap.add_argument("--off", default="sweep-cf-validation-w0")
    ap.add_argument("--on", default="sweep-cf-validation-ess-w05")
    ap.add_argument("--exclude-set", default="paper_validation")
    ap.add_argument("--dir", default="output/runs")
    a = ap.parse_args()

    screens = load_screens(target_set=a.screen_set)
    uni = set(gene_universe(screens, min_screen_freq=2))
    drop = {x.dataset_name for x in load_screens(target_set=a.exclude_set)} if a.exclude_set else set()

    rows = []
    for i, s in enumerate(screens):
        if s.dataset_name in drop:
            continue
        po = run_dir(a.dir, a.off, i, s.dataset_name) / "result.json"
        pn = run_dir(a.dir, a.on, i, s.dataset_name) / "result.json"
        if not (po.is_file() and pn.is_file()):
            continue
        lib = set(s.genes)
        d = (_adj_ef_from_run(pn, lib, uni, s.total_hits, len(lib), 1000)
             - _adj_ef_from_run(po, lib, uni, s.total_hits, len(lib), 1000))
        B = json.loads(pn.read_text())
        seeds, rate, active = [], [], 0
        cum_h = cum_n = 0
        for st in B["steps"]:
            tr = st.get("acquisition_trace") or {}
            if tr.get("codep_active"):
                active += 1
            seeds.append(tr.get("n_seeds", 0) or 0)
            rate.append(cum_h / cum_n if cum_n else 0.0)
            cum_h += sum(st.get("hits") or []); cum_n += len(st.get("acquired_batch") or [])
        rows.append(dict(name=s.dataset_name, d=d, total_hits=s.total_hits,
                         active=active, seeds=np.array(seeds), rate=np.array(rate),
                         final_seeds=seeds[-1] if seeds else 0))

    print(f"{len(rows)} screens with both arms\n")
    # how many rounds survive each candidate threshold, and what the delta is on those screens
    print(f"{'rule':34s} {'screens it fires on':>20s} {'mean dEF on those':>18s} {'mean dEF overall':>17s}")
    print("-" * 94)
    base = np.mean([r["d"] for r in rows])
    for label, pred in [
        ("current: min_seeds 2", lambda r: (r["seeds"] >= 2).any()),
        ("min_seeds 5", lambda r: (r["seeds"] >= 5).any()),
        ("min_seeds 10", lambda r: (r["seeds"] >= 10).any()),
        ("min_seeds 20", lambda r: (r["seeds"] >= 20).any()),
        ("min_seeds 10 + rate <= 0.40", lambda r: ((r["seeds"] >= 10) & (r["rate"] <= 0.40)).any()),
        ("min_seeds 10 + rate <= 0.25", lambda r: ((r["seeds"] >= 10) & (r["rate"] <= 0.25)).any()),
        ("min_seeds 20 + rate <= 0.25", lambda r: ((r["seeds"] >= 20) & (r["rate"] <= 0.25)).any()),
    ]:
        fires = [r for r in rows if pred(r)]
        on = np.mean([r["d"] for r in fires]) if fires else float("nan")
        # screens the rule silences revert to the baseline, i.e. delta 0
        overall = np.mean([r["d"] if pred(r) else 0.0 for r in rows])
        print(f"{label:34s} {len(fires):11d} / {len(rows):<6d} {on:+18.3f} {overall:+17.3f}")

    print(f"\n(overall assumes a silenced screen scores exactly its weight-0 run, which is true for a\n"
          f" per-screen gate but only APPROXIMATE for min_seeds: raising it changes which genes are\n"
          f" picked in the rounds that still fire, so the rerun will not match this exactly.)")

    low = [r for r in rows if r["final_seeds"] < 10]
    print(f"\nscreens that never reached 10 confirmed hits: {len(low)}")
    if low:
        print(f"   their mean dEF {np.mean([r['d'] for r in low]):+.3f}  "
              f"({sum(r['d'] > 0 for r in low)} up / {sum(r['d'] < 0 for r in low)} down)")
    hi = [r for r in rows if r["total_hits"] > 1000]
    if hi:
        print(f"budget-saturated screens: {len(hi)}, mean dEF {np.mean([r['d'] for r in hi]):+.3f}, "
              f"mean observed rate by round 3 {100*np.mean([r['rate'][2] for r in hi]):.1f}%")


if __name__ == "__main__":
    main()
