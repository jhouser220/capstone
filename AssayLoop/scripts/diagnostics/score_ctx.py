"""Score the test-time context-intervention sweeps against the nested baseline.

Reports, for each sweep: mean adjusted EF, the paired delta vs the baseline with
a 95% bootstrap CI over screen FAMILIES (the validation set contains a
near-duplicate fitness family; bootstrapping over screens would badly overstate
n), the per-screen deltas, and whether a deployable gate -- selection direction,
read from the phenotype text -- picks out the screens where the variant wins.
"""
from __future__ import annotations

import argparse, sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "src")
from assayloop.screen_direction import classify_selection                  # noqa: E402
from assayloop.scripts.full_genome_table import _adj_ef_from_run           # noqa: E402
from assayloop.tasks import gene_universe, load_screens                    # noqa: E402
from run_paths import run_dir  # noqa: E402

JACCARD = 0.3


def families(screens) -> dict[str, int]:
    """Single-linkage clusters of screens whose hit sets overlap at Jaccard >= 0.3."""
    names = [s.dataset_name for s in screens]
    sets = [{g for g, h in zip(s.genes, s.hits) if h} for s in screens]
    parent = list(range(len(names)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]; i = parent[i]
        return i

    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            u = len(sets[i] | sets[j])
            if u and len(sets[i] & sets[j]) / u >= JACCARD:
                parent[find(i)] = find(j)
    roots = {}
    return {n: roots.setdefault(find(i), len(roots)) for i, n in enumerate(names)}


def ef_by_screen(sweep: str, screens, universe, root="output/runs") -> dict[str, float]:
    out = {}
    for i, s in enumerate(screens):
        for d in (run_dir(root, sweep, i, s.dataset_name),):
            if (d / "result.json").is_file():
                lib = set(s.genes)
                out[s.dataset_name] = _adj_ef_from_run(
                    d / "result.json", lib, universe, s.total_hits, len(lib), 1000)
    return out


def boot_ci(deltas: dict[str, float], fam: dict[str, int], n=20000, seed=0):
    """Paired bootstrap of the mean delta, resampling families with replacement.

    A family contributes the mean of its screens, so an 8-screen family counts
    once -- the same weighting the point estimate uses.
    """
    groups: dict[int, list[float]] = {}
    for k, v in deltas.items():
        groups.setdefault(fam[k], []).append(v)
    fam_means = np.array([np.mean(v) for v in groups.values()])
    rng = np.random.default_rng(seed)
    draws = fam_means[rng.integers(0, len(fam_means), size=(n, len(fam_means)))].mean(axis=1)
    return float(fam_means.mean()), tuple(np.percentile(draws, [2.5, 97.5])), len(fam_means)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--poor-below", type=float, default=5.0,
                    help="baseline adjusted EF under which a screen counts as poorly served")
    ap.add_argument("--screen-set", default="paper_validation")
    ap.add_argument("--base", default="sweep-ctx-base", help="the nested-baseline sweep id")
    ap.add_argument("--base-dir", default="output/runs",
                    help="where the baseline sweep lives (runs move to output/runs/baseline/<topic>/)")
    ap.add_argument("--runs-dir", default="output/runs", help="where the variant sweeps live")
    ap.add_argument("--glob", default="sweep-ctx-*", help="which variant sweeps to score")
    a = ap.parse_args()

    screens = load_screens(target_set=a.screen_set)
    universe = set(gene_universe(screens, min_screen_freq=2))
    fam = families(screens)
    direction = {s.dataset_name: classify_selection(s.phenotype)[0] for s in screens}

    sweeps = sorted({p.name.rsplit("-", 2)[0]
                     for p in Path(a.runs_dir).glob(f"*/{a.glob}-00-*")}
                    | {p.name.rsplit("-", 2)[0]
                       for p in Path(a.runs_dir).glob(f"{a.glob}-00-*")})
    base = ef_by_screen(a.base, screens, universe, a.base_dir)
    if not base:
        sys.exit(f"baseline {a.base} not found under {a.base_dir}")
    ef = {sw: ef_by_screen(sw, screens, universe, a.runs_dir) for sw in sweeps}
    ef[a.base] = base
    sweeps = [sw for sw in sweeps if sw != a.base]
    print(f"{len(set(fam.values()))} families over {len(screens)} screens; "
          f"baseline mean adjusted EF {np.mean(list(base.values())):.3f}\n")

    hdr = f"{'sweep':26s} {'EF':>6s} {'dEF':>7s} {'95% CI':>18s} {'wins':>6s} {'poor-screen dEF':>16s}"
    print(hdr); print("-" * len(hdr))
    poor = [n for n, v in base.items() if v < a.poor_below]
    for sw in sweeps:
        d = {n: ef[sw][n] - base[n] for n in base if n in ef[sw]}
        if len(d) < len(base):
            print(f"{sw:26s} incomplete ({len(d)}/{len(base)} screens)"); continue
        m, (lo, hi), nf = boot_ci(d, fam)
        pm = np.mean([d[n] for n in poor]) if poor else float("nan")
        print(f"{sw:26s} {np.mean(list(ef[sw].values())):6.3f} {m:+7.3f} "
              f"[{lo:+6.3f}, {hi:+6.3f}] {sum(v > 0 for v in d.values()):3d}/{len(d):<2d} {pm:+16.3f}")

    print(f"\n(poor-screen column = mean delta over the {len(poor)} screens with "
          f"baseline EF < {a.poor_below}: {', '.join(sorted(poor))})")

    # -- per-screen detail + the deployable-gate question ---------------------
    other = list(sweeps)
    print("\nPer-screen adjusted EF (delta vs baseline):\n")
    print(f"{'screen':16s} {'fam':>3s} {'dir':>11s} {'base':>6s} " +
          " ".join(f"{sw.replace('sweep-ctx-',''):>12s}" for sw in other))
    for n in sorted(base, key=lambda n: base[n]):
        row = " ".join(f"{ef[sw].get(n, float('nan')) - base[n]:+12.3f}" for sw in other)
        print(f"{n:16s} {fam[n]:3d} {direction[n]:>11s} {base[n]:6.2f} {row}")

    # Direction is read from the phenotype sentence, so gating on it needs no
    # labels and is deployable. Report the gated subset with its own CI.
    print("\nGated on direction == enrichment (deployable: read from the phenotype text)\n")
    enr = [n for n in base if direction[n] == "enrichment"]
    print(f"{'sweep':26s} {'dEF (enrichment only)':>22s} {'95% CI':>18s} {'wins':>6s}")
    for sw in other:
        d = {n: ef[sw][n] - base[n] for n in enr if n in ef[sw]}
        m, (lo, hi), nf = boot_ci(d, fam)
        print(f"{sw:26s} {m:+22.3f} [{lo:+6.3f}, {hi:+6.3f}] {sum(v > 0 for v in d.values()):3d}/{len(d):<2d}")
    print(f"  ({len(enr)} enrichment screens in {len({fam[n] for n in enr})} families)")

    print("\nIs selection direction a usable gate? (mean delta by direction)\n")
    print(f"{'sweep':26s} " + " ".join(f"{d:>12s}" for d in ("enrichment", "dropout", "both", "unknown")))
    for sw in other:
        cells = []
        for dd in ("enrichment", "dropout", "both", "unknown"):
            vals = [ef[sw][n] - base[n] for n in base if direction[n] == dd and n in ef[sw]]
            cells.append(f"{np.mean(vals):+12.3f}" if vals else f"{'-':>12s}")
        print(f"{sw:26s} " + " ".join(cells))


if __name__ == "__main__":
    main()
