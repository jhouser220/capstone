"""Score a joint weight sweep on a held-out fold, with PAIRED contrasts.

Two things this does that `codep_gated.py` does not.

1. **Paired contrasts against a reference arm, not just against the baseline.**
   The question is whether (0.5, 0.75) beats (0.5, 0.5), and both arms are measured
   on the same screens. Differencing them per screen cancels the screen-to-screen
   variance that dominates an absolute delta, which is where the power comes from:
   the absolute deltas on this fold have CIs ~0.22 wide, far too coarse to resolve
   the +0.02-0.05 differences at stake.

2. **It excludes the curated screens by name**, so the reported numbers are on the
   fresh screens of the fold only (198 of 218 for `validation`).

Both per-screen and per-family means are reported. Family weighting is the correct
default -- this benchmark contains near-duplicate screens at hit-set Jaccard up to
1.00 (§6) -- but it roughly doubles the per-screen figure here because the three
largest families are negative, so quote both.

    uv run python scripts/diagnostics/score_weight_sweep.py \
        --baseline sweep-dm-validation-w0 \
        --reference sweep-dm-validation-dm05knn \
        --arm 'codep 0.5 / knn 0.25=sweep-wt-validation-c05k025' \
        --arm 'codep 0.5 / knn 0.75=sweep-wt-validation-c05k075'
"""
from __future__ import annotations

import argparse
import sys

import numpy as np

sys.path.insert(0, "src")
sys.path.insert(0, "scripts/diagnostics")
from assayloop.tasks import gene_universe, load_screens   # noqa: E402
from score_ctx import boot_ci, ef_by_screen, families     # noqa: E402


def boot_screens(deltas: dict[str, float], n=20000, seed=0):
    """Bootstrap the PER-SCREEN mean by resampling screens.

    `score_ctx.boot_ci` resamples FAMILIES and returns the family-weighted mean, so
    its interval belongs beside the per-family figure. Printing it beside the
    per-screen mean produces an interval that does not contain the number it sits
    next to -- which is exactly what the first version of this script did.
    """
    v = np.array(list(deltas.values()))
    rng = np.random.default_rng(seed)
    draws = v[rng.integers(0, len(v), size=(n, len(v)))].mean(axis=1)
    return float(v.mean()), tuple(np.percentile(draws, [2.5, 97.5]))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--screen-set", default="validation")
    ap.add_argument("--exclude-set", default="paper_validation",
                    help="screens to drop because they chose a setting earlier")
    ap.add_argument("--runs-dir", default="output/runs")
    ap.add_argument("--baseline", required=True, help="the weight-0 nested baseline")
    ap.add_argument("--reference", default=None,
                    help="arm to form paired contrasts against (the deployed config)")
    ap.add_argument("--ref-label", default="DEPLOYED 0.5/0.5")
    ap.add_argument("--arm", action="append", default=[], help="LABEL=sweep-id, repeatable")
    a = ap.parse_args()

    screens = load_screens(target_set=a.screen_set)
    uni = set(gene_universe(screens, min_screen_freq=2))
    drop = {x.dataset_name for x in load_screens(target_set=a.exclude_set)} \
        if a.exclude_set else set()

    def load(sweep):
        ef = ef_by_screen(sweep, screens, uni, root=a.runs_dir)
        return {k: v for k, v in ef.items() if k not in drop}

    base = load(a.baseline)
    if not base:
        raise SystemExit(f"no runs found for baseline {a.baseline} in {a.runs_dir}")
    arms: list[tuple[str, dict[str, float]]] = []
    if a.reference:
        arms.append((a.ref_label, load(a.reference)))
    for spec in a.arm:
        label, _, sid = spec.partition("=")
        if not sid:
            raise SystemExit(f"--arm wants LABEL=sweep-id, got {spec!r}")
        got = load(sid)
        if not got:
            print(f"!! no runs for {sid} -- skipped")
            continue
        arms.append((label, got))

    kept = [s for s in screens if s.dataset_name in base]
    fam = families(kept)
    print(f"{a.screen_set}: {len(base)} fresh screens "
          f"({len(drop)} excluded as {a.exclude_set}), "
          f"{len(set(fam.values()))} families")
    print(f"baseline {a.baseline}: mean adjusted EF {np.mean(list(base.values())):.3f}\n")

    hdr = (f"{'arm':24s} {'n':>4s} {'dEF/screen':>11s} {'95% CI (screens)':>18s} "
           f"{'dEF/family':>11s} {'95% CI (fams)':>18s} {'up/dn/same':>12s}")
    print(hdr)
    print("-" * len(hdr))
    for label, ef in arms:
        common = sorted(set(ef) & set(base))
        d = {k: ef[k] - base[k] for k in common}
        f2 = {k: fam[k] for k in common}
        # two means, each with ITS OWN interval. The family mean weights each family
        # once, which is correct given near-duplicates, and runs roughly 2x the
        # per-screen number here because the largest families are negative.
        sm, (slo, shi) = boot_screens(d)
        fm, (flo, fhi), _ = boot_ci(d, f2)
        up = sum(1 for v in d.values() if v > 1e-9)
        dn = sum(1 for v in d.values() if v < -1e-9)
        print(f"{label:24s} {len(d):4d} {sm:+11.3f} [{slo:+6.3f}, {shi:+6.3f}] "
              f"{fm:+11.3f} [{flo:+6.3f}, {fhi:+6.3f}] "
              f"{f'{up}/{dn}/{len(d) - up - dn}':>12s}")

    if not a.reference:
        return
    ref = dict(arms[0][1])
    print(f"\nPAIRED contrasts against {a.ref_label} (same screens, variance cancels)")
    h2 = (f"{'arm vs reference':24s} {'n':>4s} {'per screen':>11s} {'95% CI':>18s} "
          f"{'per family':>11s} {'95% CI':>18s} {'wins':>9s}")
    print(h2)
    print("-" * len(h2))
    for label, ef in arms[1:]:
        common = sorted(set(ef) & set(ref))
        d = {k: ef[k] - ref[k] for k in common}
        f2 = {k: fam[k] for k in common}
        sm, (slo, shi) = boot_screens(d)
        fm, (flo, fhi), _ = boot_ci(d, f2)
        print(f"{label:24s} {len(d):4d} {sm:+11.3f} [{slo:+6.3f}, {shi:+6.3f}] "
              f"{fm:+11.3f} [{flo:+6.3f}, {fhi:+6.3f}] "
              f"{sum(1 for v in d.values() if v > 1e-9):4d}/{len(d):<4d}")
    print("\nA contrast whose CI includes zero means this fold cannot separate that arm\n"
          "from the deployed configuration -- which for differences of ~0.02-0.05 is the\n"
          "expected outcome even if the effect is real. Do not promote on a point estimate.")


if __name__ == "__main__":
    main()
