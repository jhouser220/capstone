"""Re-check every discounted candidate on the HARD screens, using campaigns already run.

Why. §4.28 found that BioPlex PPI's apparent value came entirely from screens
AssayFormer already handles (+0.024 on the 28 hard ones), and that DepMap itself falls
from +0.618 to +0.166 there. Benefit is therefore NOT distributed evenly across
difficulty, and several candidates were discounted on a pooled average that could have
hidden a hard-screen effect. This re-checks them, with DEPLOYED campaign data rather
than offline estimates -- no new runs needed.

Two filters applied throughout:
  * the 20 curated screens are dropped (they chose the settings), and
  * metric-degenerate screens are dropped: 98 of 198 fold screens have a library
    SMALLER than the 1000-pick budget, which makes perfect-play EF < 1 (§4.12
    correction). Averages over them are not interpretable.

Difficulty is bucketed by AssayFormer's own EF on the shared weight-0 baseline, so
every arm is bucketed identically regardless of which family it belongs to.

Each arm is paired with ITS OWN family's weight-0 baseline, because the centroid,
two-feature and diffmeans families each ran their own nested baseline.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path


sys.path.insert(0, "src")
sys.path.insert(0, "scripts/diagnostics")
from assayloop.scripts.full_genome_table import _adj_ef_from_run   # noqa: E402
from assayloop.tasks import gene_universe, load_screens            # noqa: E402
from score_ctx import boot_ci, families                           # noqa: E402
from run_paths import run_dir  # noqa: E402

# label = arm sweep : its own nested baseline sweep
ARMS = [
    "deployed diffmeans+kNN=sweep-dm-validation-dm05knn:sweep-dm-validation-w0",
    "codep alone (diffmeans)=sweep-dm-validation-dm05:sweep-dm-validation-w0",
    "codep weight 1.0=sweep-dm-validation-dm10:sweep-dm-validation-w0",
    "centroid alone=sweep-cf-validation-ess-w05:sweep-cf-validation-w0",
    "EXPRESSION kernel=sweep-cf-validation-expr-w05:sweep-cf-validation-w0",
    "screen-kNN alone=sweep-tf-validation-knn05:sweep-tf-validation-w0",
    "centroid+kNN=sweep-tf-validation-both05:sweep-tf-validation-w0",
    "co-hit added (0.25)=sweep-ch-validation-h025:sweep-dm-validation-w0",
    "co-hit added (0.5)=sweep-ch-validation-h05:sweep-dm-validation-w0",
    "all three at 1/3=sweep-ch-validation-eq:sweep-dm-validation-w0",
    "min_seeds 1 (not 2)=sweep-ms1-validation:sweep-dm-validation-w0",
    "knn weight 0.25=sweep-wt-validation-c05k025:sweep-dm-validation-w0",
    "knn weight 0.75=sweep-wt-validation-c05k075:sweep-dm-validation-w0",
    "knn weight 1.0=sweep-wt-validation-c05k10:sweep-dm-validation-w0",
    "codep weight 0.25=sweep-wt-validation-c025k05:sweep-dm-validation-w0",
    "codep weight 0.75=sweep-wt-validation-c075k05:sweep-dm-validation-w0",
]
BUCKETS = ((0.0, 2.0, "hard"), (2.0, 5.0, "mid"), (5.0, 1e9, "easy"))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--screen-set", default="validation")
    ap.add_argument("--exclude-set", default="paper_validation")
    ap.add_argument("--dir", default="output/runs")
    ap.add_argument("--ref-baseline", default="sweep-dm-validation-w0",
                    help="the arm whose EF defines the difficulty buckets for everyone")
    ap.add_argument("--min-library", type=int, default=1001)
    ap.add_argument("--arm", action="append", default=None,
                    help="LABEL=arm_sweep:baseline_sweep, repeatable")
    a = ap.parse_args()

    screens = load_screens(target_set=a.screen_set)
    uni = set(gene_universe(screens, min_screen_freq=2))
    drop = {x.dataset_name for x in load_screens(target_set=a.exclude_set)} \
        if a.exclude_set else set()
    keepable = [(i, s) for i, s in enumerate(screens)
                if s.dataset_name not in drop and len(s.genes) > a.min_library]
    print(f"{len(screens)} screens -> {len(keepable)} after dropping {len(drop)} curated "
          f"and every library <= {a.min_library} genes\n")

    def ef_of(sweep: str) -> dict[str, float]:
        out = {}
        for i, s in keepable:
            p = run_dir(a.dir, sweep, i, s.dataset_name) / "result.json"
            if p.is_file():
                lib = set(s.genes)
                out[s.dataset_name] = _adj_ef_from_run(
                    p, lib, uni, s.total_hits, len(lib), 1000)
        return out

    cache: dict[str, dict[str, float]] = {}

    def ef_cached(sweep: str) -> dict[str, float]:
        if sweep not in cache:
            cache[sweep] = ef_of(sweep)
        return cache[sweep]

    ref = ef_cached(a.ref_baseline)
    if not ref:
        raise SystemExit(f"no runs for the reference baseline {a.ref_baseline}")
    fam = families([s for _i, s in keepable if s.dataset_name in ref])
    bucket = {n: lab for n, v in ref.items()
              for lo, hi, lab in BUCKETS if lo <= v < hi}
    counts = {lab: sum(1 for b in bucket.values() if b == lab) for _l, _h, lab in BUCKETS}
    print("difficulty buckets by AssayFormer's own EF: "
          + ", ".join(f"{lab} {counts[lab]}" for _l, _h, lab in BUCKETS))
    print("(per-FAMILY mean dEF vs each arm's own nested baseline; "
          "CI shown for the hard bucket only)\n")

    hdr = (f"{'arm':26s} {'n':>4s} {'ALL':>8s} {'hard':>8s} {'95% CI (hard)':>18s} "
           f"{'mid':>8s} {'easy':>8s}")
    print(hdr)
    print("-" * len(hdr))
    for spec in (a.arm or ARMS):
        label, _, rest = spec.rpartition("=")[0], "=", spec.rpartition("=")[2]
        arm, _, base = rest.partition(":")
        ea, eb = ef_cached(arm), ef_cached(base)
        common = sorted(set(ea) & set(eb) & set(bucket))
        if len(common) < 10:
            print(f"{label:26s} {len(common):4d}  (runs missing -- skipped)")
            continue
        d = {n: ea[n] - eb[n] for n in common}
        cells = {}
        for _lo, _hi, lab in BUCKETS:
            sel = [n for n in common if bucket[n] == lab]
            if len(sel) < 3:
                cells[lab] = (float("nan"), (float("nan"),) * 2)
                continue
            mm, ci, _ = boot_ci({n: d[n] for n in sel}, {n: fam[n] for n in sel})
            cells[lab] = (mm, ci)
        allm, _, _ = boot_ci(d, {n: fam[n] for n in common})
        hm, (hlo, hhi) = cells["hard"]
        print(f"{label:26s} {len(common):4d} {allm:+8.3f} {hm:+8.3f} "
              f"[{hlo:+6.3f}, {hhi:+6.3f}] {cells['mid'][0]:+8.3f} {cells['easy'][0]:+8.3f}")

    # PAIRED contrasts on the hard bucket only. The table above compares arms through
    # their own baselines, so its intervals are wide and small orderings in it are not
    # evidence. Differencing two arms per screen cancels the screen variance.
    dep = ef_cached("sweep-dm-validation-dm05knn")
    hard = [n for n in bucket if bucket[n] == "hard" and n in dep]
    print(f"\nPAIRED against the deployed arm, HARD screens only (n={len(hard)})")
    print(f"  {'arm':26s} {'delta':>9s} {'95% CI':>18s} {'better on':>11s}")
    print("  " + "-" * 66)
    for spec in (a.arm or ARMS):
        label, rest = spec.rpartition("=")[0], spec.rpartition("=")[2]
        arm = rest.partition(":")[0]
        if arm == "sweep-dm-validation-dm05knn":
            continue
        ea = ef_cached(arm)
        sel = [n for n in hard if n in ea]
        if len(sel) < 5:
            continue
        d = {n: ea[n] - dep[n] for n in sel}
        mm, (lo, hi), _ = boot_ci(d, {n: fam[n] for n in sel})
        print(f"  {label:26s} {mm:+9.3f} [{lo:+6.3f}, {hi:+6.3f}] "
              f"{sum(1 for v in d.values() if v > 1e-9):5d}/{len(sel):<5d}")

    print("\nRead the `hard` column against `all`. A candidate discounted on its pooled\n"
          "average is worth revisiting ONLY if its hard-bucket figure is both positive and\n"
          "larger than the deployed arm's. §4.28 is the precedent: PPI looked worth +0.190\n"
          "overall and was +0.024 where the headroom is.")


if __name__ == "__main__":
    main()
