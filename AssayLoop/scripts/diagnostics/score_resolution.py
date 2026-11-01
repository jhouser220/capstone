"""Measure the SAME campaigns with a better instrument, to buy resolution for free.

Why. Three candidates have now been closed on unmeasurability rather than on a null
(§4.18 Perturb-seq, §4.21 co-hit, §4.28 BioPlex), and the best remaining candidate is
expected at ~+0.06 deployed against a fold that resolves to about +-0.1. When three of
your own verdicts are "below resolution", resolution is the binding constraint, not the
candidates. Nothing here changes a campaign; it re-reads runs already on disk.

Two independent gains, both standard and both free:

  1. nAUC over the budget curve instead of EF at one budget. Adjusted EF at budget 1000
     is a single endpoint -- the state after all ten rounds. The per-round curve has ten
     readings of the same campaign, and averaging them cuts endpoint noise without
     changing what is being estimated. NOTE `_adj_ef_from_run`'s `budget` argument does
     NOT truncate the pick list; it only rescales the denominator, so the curve has to
     be built from pick PREFIXES, which is what `ef_curve` does here. EF at the last
     round boundary reproduces `_adj_ef_from_run` exactly, and that is asserted.

  2. ANCOVA on the paired delta. The paired design already removes each screen's main
     effect, but the deltas still vary with things known BEFORE the intervention --
     library size, hit count, selection direction. Regressing the delta on centred
     covariates leaves the point estimate algebraically unchanged (OLS intercept on
     centred covariates IS the mean) and shrinks its interval by whatever those
     covariates explain. The bootstrap still resamples FAMILIES, with the fit redone
     inside each replicate, so the clustering §6 requires is preserved.

Two covariate sets, because one of them is not strictly pre-treatment:

  strict  log library size, log total hits, selection-direction dummies. Nothing from
          either run -- these are properties of the screen as posed.
  full    strict + the BASELINE arm's EF and its seed count at round 3. Legitimate for
          a paired delta (a baseline covariate), and the intercept is still the mean,
          but it is derived from the control arm rather than from the screen, so it is
          reported separately and never alone.

    uv run python scripts/diagnostics/score_resolution.py \
        --baseline sweep-dm-validation-w0 \
        --arm "DEPLOYED 0.5/0.5=sweep-dm-validation-dm05knn"
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "src")
sys.path.insert(0, "scripts/diagnostics")
from assayloop.metrics.hits_auc import adjusted_ef_value           # noqa: E402
from assayloop.screen_direction import classify_selection          # noqa: E402
from assayloop.scripts.full_genome_table import _adj_ef_from_run   # noqa: E402
from assayloop.tasks import gene_universe, load_screens            # noqa: E402
from score_ctx import families                                     # noqa: E402
from run_paths import run_dir  # noqa: E402


def ef_curve(path: Path, lib: set[str], uni: set[str], total_hits: int,
             domain: int) -> list[float]:
    """Adjusted EF after each round, from pick PREFIXES.

    Round k is scored at budget 100*k against the first 100*k picks, which is the
    quantity "EF at budget b" is supposed to mean. `_adj_ef_from_run` counts every pick
    whatever budget it is passed, so it cannot produce this curve.
    """
    rd = json.loads(path.read_text())
    out, n1, n2, hits = [], 0, 0, 0
    steps = rd.get("steps", [])
    for k, step in enumerate(steps, start=1):
        for g, h in zip(step.get("acquired_batch", []) or [],
                        step.get("hits", []) or []):
            if g in lib:
                n1 += 1
                hits += bool(h)
            elif g in uni:
                n2 += 1
        # the INTENDED budget at this round, not the picks delivered. Identical for the
        # codep arms (every slot filled with a valid in-universe gene, so n1+n2 == 100k,
        # which the endpoint assertion confirms) but NOT for LLM arms, which leave slots
        # unfilled and emit invalid symbols -- a shortfall the benchmark charges and
        # `n1 + n2` would forgive. Kept correct here so this scorer is safe on any arm.
        budget_k = int(round(1000 * k / max(len(steps), 1)))
        out.append(adjusted_ef_value(hits, n1, n2, domain, total_hits, budget_k))
    return out


def reg_adjusted(y: np.ndarray, Z: np.ndarray | None,
                 centre: np.ndarray | None) -> float:
    """The regression (ANCOVA) estimator of the mean delta: the replicate's own mean,
    corrected for how far the replicate's covariates drifted from the FULL sample's.

        est = y_rep.mean() + beta' (x_full - x_rep.mean())

    The fixed `centre` is what makes this work. Re-centring Z inside the replicate makes
    the intercept identically y_rep.mean(), so the adjustment cancels and the interval
    does not move -- that was the first version of this function, and it silently
    reported "0% width change" for every covariate set.
    """
    if Z is None or Z.shape[1] == 0 or centre is None:
        return float(y.mean())
    X = np.column_stack([np.ones(len(y)), Z])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return float(y.mean() + beta[1:] @ (centre - Z.mean(axis=0)))


def boot(y: np.ndarray, Z: np.ndarray | None, n: int, seed: int):
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(y), size=(n, len(y)))
    centre = None if Z is None else Z.mean(axis=0)
    v = np.array([reg_adjusted(y[i], None if Z is None else Z[i], centre)
                  for i in idx])
    return float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--screen-set", default="validation")
    ap.add_argument("--exclude-set", default="paper_validation")
    ap.add_argument("--runs-dir", default="output/runs")
    ap.add_argument("--baseline", required=True)
    ap.add_argument("--arm", action="append", default=[], help="LABEL=sweep-id, repeatable")
    ap.add_argument("--min-library", type=int, default=0,
                    help="drop screens with library <= this (1001 excludes the "
                         "metric-degenerate ones, §4.29)")
    ap.add_argument("--boot", type=int, default=4000)
    a = ap.parse_args()

    screens = load_screens(target_set=a.screen_set)
    uni = set(gene_universe(screens, min_screen_freq=2))
    drop = ({x.dataset_name for x in load_screens(target_set=a.exclude_set)}
            if a.exclude_set else set())

    def curves(sweep: str) -> dict[str, list[float]]:
        out = {}
        for i, s in enumerate(screens):
            if s.dataset_name in drop or len(s.genes) <= a.min_library:
                continue
            p = run_dir(a.runs_dir, sweep, i, s.dataset_name) / "result.json"
            if not p.is_file():
                continue
            lib = set(s.genes)
            c = ef_curve(p, lib, uni, s.total_hits, len(lib))
            # the last point must agree with the project's own EF function
            ref = _adj_ef_from_run(p, lib, uni, s.total_hits, len(lib), 1000)
            if abs(c[-1] - ref) > 1e-6:
                raise SystemExit(f"curve endpoint {c[-1]:.6f} != _adj_ef_from_run "
                                 f"{ref:.6f} for {s.dataset_name}")
            out[s.dataset_name] = c
        return out

    base = curves(a.baseline)
    if not base:
        raise SystemExit(f"no runs for baseline {a.baseline}")
    by_name = {s.dataset_name: s for s in screens}
    kept = [by_name[k] for k in base]
    fam = families(kept)
    print(f"{a.screen_set}: {len(base)} screens, {len(set(fam.values()))} families, "
          f"min-library {a.min_library}, baseline {a.baseline}")
    print(f"endpoint check passed: per-round curve agrees with _adj_ef_from_run "
          f"on all {len(base)} screens\n")

    # covariates, per screen, in the order `names` gives
    def covars(name: str) -> tuple[list[float], list[float]]:
        s = by_name[name]
        d = classify_selection(getattr(s, "phenotype", "") or "")[0]
        strict = [np.log(len(s.genes)), np.log(max(s.total_hits, 1)),
                  float(d == "dropout"), float(d == "enrichment")]
        c = base[name]
        seeds3 = c[2] if len(c) > 2 else 0.0   # baseline progress by round 3
        return strict, strict + [c[-1], seeds3]

    for spec in a.arm:
        label, _, sid = spec.partition("=")
        arm = curves(sid)
        common = [k for k in base if k in arm]
        if not common:
            print(f"!! no runs for {sid} -- skipped")
            continue
        print(f"=== {label}  ({sid}), {len(common)} paired screens")
        rounds = min(min(len(base[k]) for k in common), min(len(arm[k]) for k in common))
        d_end = {k: arm[k][rounds - 1] - base[k][rounds - 1] for k in common}
        d_auc = {k: float(np.mean(arm[k][:rounds]) - np.mean(base[k][:rounds]))
                 for k in common}

        # per-round deltas, so the curve itself is visible
        print("    dEF by round boundary: " + "  ".join(
            f"r{r + 1} {np.mean([arm[k][r] - base[k][r] for k in common]):+.3f}"
            for r in range(rounds)))

        for mlabel, dd in (("dEF at budget 1000", d_end),
                           (f"dnAUC over {rounds} rounds", d_auc)):
            print(f"    {mlabel}")
            for unit in ("screens", "families"):
                if unit == "screens":
                    keys = common
                    y = np.array([dd[k] for k in keys])
                    C = [covars(k) for k in keys]
                else:
                    groups: dict[int, list[str]] = {}
                    for k in common:
                        groups.setdefault(fam[k], []).append(k)
                    keys = list(groups)
                    y = np.array([np.mean([dd[k] for k in groups[g]]) for g in keys])
                    C = [([np.mean([covars(k)[0][j] for k in groups[g]])
                           for j in range(4)],
                          [np.mean([covars(k)[1][j] for k in groups[g]])
                           for j in range(6)]) for g in keys]
                Zs = np.array([c[0] for c in C], dtype=float)
                Zf = np.array([c[1] for c in C], dtype=float)
                base_lo, base_hi = boot(y, None, a.boot, 0)
                w0 = base_hi - base_lo
                print(f"      {unit:8s} n={len(y):3d}  mean {y.mean():+.3f}")
                for aname, Z in (("unadjusted", None), ("strict", Zs), ("full", Zf)):
                    lo, hi = boot(y, Z, a.boot, 0)
                    r2 = ""
                    if Z is not None:
                        Zc = Z - Z.mean(axis=0, keepdims=True)
                        X = np.column_stack([np.ones(len(y)), Zc])
                        bh = np.linalg.lstsq(X, y, rcond=None)[0]
                        ss = 1 - ((y - X @ bh) ** 2).sum() / ((y - y.mean()) ** 2).sum()
                        r2 = f"  R2 {ss:.3f}"
                    tag = "" if aname == "unadjusted" else \
                        f"  width {(hi - lo) / w0 - 1:+.0%}"
                    snr = y.mean() / ((hi - lo) / 3.92) if hi > lo else float("nan")
                    print(f"        {aname:10s} [{lo:+.3f}, {hi:+.3f}]"
                          f"  width {hi - lo:.3f}{tag}  SNR {snr:4.2f}{r2}")
        print()


if __name__ == "__main__":
    main()
