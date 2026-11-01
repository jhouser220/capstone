"""codep_seed with a DEPLOYABLE gate: switch it off on saturated screens.

The campaigns show codep_seed helps where the 1000-pick budget leaves headroom
and hurts where it does not: on validation the 8 budget-saturated screens
(total hits > 1000, the near-duplicate fitness family) lost 158 hits and improved
on 0/8, while the 12 screens with headroom gained 65 and improved on 8/12.

Total hits is not observable during a campaign, but the early hit RATE is, and it
separates the two groups cleanly -- AssayFormer hits at 71-95% in rounds 1-3 on
the saturated screens and 2-23% on the rest. So "disable the feature if the first
three rounds hit above `--thresh`" is a rule a real campaign could apply. This
scores that policy by splicing each screen's run from whichever arm the gate picks,
using only round 1-3 information to decide.
"""
from __future__ import annotations

import argparse, json, sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "src")
sys.path.insert(0, "scripts/diagnostics")
from assayloop.scripts.full_genome_table import _adj_ef_from_run   # noqa: E402
from assayloop.tasks import gene_universe, load_screens            # noqa: E402
from score_ctx import boot_ci, families                            # noqa: E402
from run_paths import run_dir  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--screen-set", default="paper_validation")
    ap.add_argument("--off", default="sweep-cdv-w0", help="the arm the gate falls back to")
    ap.add_argument("--on", default="sweep-cdv-w05", help="the arm the gate selects")
    ap.add_argument("--dir", default="output/runs")
    ap.add_argument("--thresh", type=float, default=0.40,
                    help="rounds 1-3 hit rate above which the feature is disabled")
    ap.add_argument("--exclude-set", default=None,
                    help="drop every screen belonging to this screen set, e.g. paper_test. Use it "
                         "when confirming on a whole fold: the curated 20 screens inside that fold "
                         "are the ones the weight was chosen on, so they are not held out.")
    ap.add_argument("--quiet", action="store_true", help="summary only, no per-screen table")
    a = ap.parse_args()

    screens = load_screens(target_set=a.screen_set)
    uni = set(gene_universe(screens, min_screen_freq=2))   # the fold's own universe
    drop: set[str] = set()
    if a.exclude_set:
        drop = {x.dataset_name for x in load_screens(target_set=a.exclude_set)}
        print(f"excluding {len(drop)} screens that belong to {a.exclude_set}")
    keep = [x for x in screens if x.dataset_name not in drop]
    fam = families(keep)
    rows, gated = {}, {}
    for i, s in enumerate(screens):
        # run dirs are indexed by position in the FULL fold, so enumerate that and
        # filter afterwards rather than re-indexing
        if s.dataset_name in drop:
            continue
        po = run_dir(a.dir, a.off, i, s.dataset_name) / "result.json"
        pn = run_dir(a.dir, a.on, i, s.dataset_name) / "result.json"
        if not (po.is_file() and pn.is_file()):
            continue
        lib = set(s.genes)
        ef_o = _adj_ef_from_run(po, lib, uni, s.total_hits, len(lib), 1000)
        ef_n = _adj_ef_from_run(pn, lib, uni, s.total_hits, len(lib), 1000)
        # the gate's only input: AssayFormer's own hit rate in rounds 1-3
        steps = json.loads(po.read_text())["steps"][:3]
        r13 = np.mean([h for st in steps for h in (st["hits"] or [])])
        use_on = r13 <= a.thresh
        rows[s.dataset_name] = (ef_o, ef_n, r13, use_on)
        gated[s.dataset_name] = (ef_n if use_on else ef_o) - ef_o

    ungated = {k: v[1] - v[0] for k, v in rows.items()}
    n_on = sum(1 for v in rows.values() if v[3])
    gm, (glo, ghi), nf = boot_ci(gated, fam)
    um, (ulo, uhi), _ = boot_ci(ungated, fam)
    print(f"{a.screen_set}: {len(rows)} screens, {nf} families, "
          f"gate ON for {n_on} of {len(rows)} (rounds 1-3 hit rate <= {a.thresh:.0%})\n")
    print(f"  ungated  dEF {um:+.3f} [{ulo:+.3f}, {uhi:+.3f}]  "
          f"better on {sum(v > 0 for v in ungated.values())}/{len(ungated)}")
    print(f"  GATED    dEF {gm:+.3f} [{glo:+.3f}, {ghi:+.3f}]  "
          f"better on {sum(v > 0 for v in gated.values())}/{len(gated)}")
    if a.quiet:
        return
    print(f"\n{'screen':16s} {'r1-3 hits':>10s} {'gate':>6s} {'off EF':>8s} {'on EF':>8s} {'used':>8s}")
    for n, (eo, en, r13, on) in sorted(rows.items(), key=lambda kv: kv[1][2]):
        print(f"{n:16s} {100*r13:9.1f}% {'ON' if on else 'off':>6s} "
              f"{eo:8.2f} {en:8.2f} {gated[n]:+8.3f}")


if __name__ == "__main__":
    main()
