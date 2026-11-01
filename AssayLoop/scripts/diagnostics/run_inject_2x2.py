"""IN-CAMPAIGN INJECTION, as a 2x2: would a perfect NOMINATOR help `codep_seed`? (§4.50)

The better-posed version of §4.49. That test handed the scorers seeds from outside the
campaign, which forced those genes OUT of the candidate pool so they could not be banked
as free hits -- and that in turn denied the arm any credit for finding them itself, a
handicap growing with the dose. Its 12-seed arm was clean enough to conclude from, but
the dose-response was confounded.

THIS DESIGN HAS NO SUCH PROBLEM. The injected genes enter the campaign as REAL PICKS in
the early rounds, taking the lowest-ranked slots of each batch:

  * the budget stays exactly 10 x 100 = 1000, so EF's denominator is untouched,
  * nothing is confiscated, and
  * it is what a nominator actually does -- occupy a few early slots with genes the
    model would never choose.

The injected hits DO credit the score, so the experiment is only interpretable as a 2x2:

  kernel benefit, ordinary seeds  =  deployed        - baseline          (on disk, +0.175)
  kernel benefit, injected seeds  =  inject+kernel   - inject+no-kernel
  SYNERGY                         =  the second minus the first

The free-hit credit is identical in both injected arms and cancels exactly in that
difference of differences. No accounting correction is needed, which is the whole point.

Injected genes are true hits that NEITHER reference arm ever acquired -- i.e. genes
AssayFormer demonstrably cannot reach -- and the SAME genes go to both arms of a pair, or
the credit would not cancel.

DOSES: 4 per round x 3 rounds = 12 per screen matches the LLM's own measured new-hit
rate (~12 per screen, §2e). 12 per round x 3 = 36 is the generous arm.

    uv run python scripts/diagnostics/run_inject_2x2.py --dose 4 --dose 12
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "src")
sys.path.insert(0, "scripts/diagnostics")
from assayloop.acquisitions.codep_seed_acq import CodepSeedAcquisition  # noqa: E402
from assayloop.experiment.runner import RunConfig, run_sweep           # noqa: E402
from assayloop.tasks import gene_universe, load_screens                # noqa: E402
from run_paths import run_dir, sweep_root  # noqa: E402

CKPT = ("/home/jeff/.cache/huggingface/hub/models--Genentech--assayformer/"
        "snapshots/4927703e5903984ae0254a9fe8e93ec0b29c27e4")


def acquired(run_dir: Path) -> set[str]:
    p = run_dir / "result.json"
    if not p.is_file():
        return set()
    rd = json.loads(p.read_text())
    return {str(g) for st in rd.get("steps", []) for g in (st.get("acquired_batch") or [])}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--screen-set", default="validation")
    ap.add_argument("--dose", action="append", type=int, default=None,
                    help="injected genes PER ROUND (rounds 1-3); repeatable")
    ap.add_argument("--rounds", type=int, default=3)
    ap.add_argument("--runs-dir", default="output/runs")
    ap.add_argument("--base", default="sweep-dm-validation-w0")
    ap.add_argument("--ref", default="sweep-dm-validation-dm05knn")
    ap.add_argument("--depmap", default="/home/jeff/depmap/achilles_noeval.npz",
                    help="DepMap gene-effect matrix. Default since 2026-10-28 is the "
                         "ABLATED atlas: 768 of 808 Achilles lines, every line used by "
                         "any evaluation set dropped. Pass achilles.npz for all 808.")
    ap.add_argument("--train-hits", default="output/analysis/train_hits")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--only", default=None, help="comma-separated names, smoke test")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    doses = a.dose or [4, 12]

    screens = load_screens(target_set=a.screen_set)
    universe = gene_universe(load_screens(target_set=a.screen_set), min_screen_freq=2)
    rng = np.random.default_rng(a.seed)

    # per screen: true hits neither reference arm ever acquired, shuffled ONCE so both
    # arms of a pair inject the identical genes in the identical order
    inject: dict[str, list[str]] = {}
    for i, s in enumerate(screens):
        seen = set()
        for sweep in (a.base, a.ref):
            seen |= acquired(run_dir(a.runs_dir, sweep, i, s.dataset_name))
        hit_of = dict(zip(s.genes, s.hits))
        pool = sorted(g for g in s.genes if hit_of.get(g) and g not in seen)
        rng.shuffle(pool)
        inject[s.dataset_name] = pool
    sizes = np.array([len(v) for v in inject.values()])
    print(f"{a.screen_set}: {len(screens)} screens, universe {len(universe)} genes")
    print(f"unreachable true hits per screen: median {np.median(sizes):.0f}, "
          f"min {sizes.min()}, max {sizes.max()}, {(sizes == 0).sum()} with none")
    for d in doses:
        need = d * a.rounds
        print(f"  dose {d}/round x {a.rounds} = {need}: "
              f"{(sizes >= need).sum()}/{len(sizes)} screens have enough")
    if a.dry_run:
        return

    for d in doses:
        for kernel_on in (False, True):
            tag = "k1" if kernel_on else "k0"
            sweep_id = f"sweep-inj-{a.screen_set}-d{d}{tag}"
            if list(sweep_root(a.runs_dir, sweep_id).glob(f"{sweep_id}-[0-9]*")):
                print(f"SKIP {sweep_id}")
                continue

            def factory(screen, d=d, kernel_on=kernel_on):
                return CodepSeedAcquisition(
                    profile_npz=a.depmap, train_hits=a.train_hits,
                    weight=0.5 if kernel_on else 0.0,
                    knn_weight=0.5 if kernel_on else 0.0,
                    query="diffmeans", fuse="rrf",
                    inject_genes=inject[screen.dataset_name],
                    inject_per_round=d, inject_rounds=a.rounds)

            cfg = RunConfig(
                model="amortized_ranker", acq="codep_seed", screen_set=a.screen_set,
                n_steps=10, batch_size=100,
                model_overrides={"checkpoint": CKPT, "ckpt_file": "model_last.pt",
                                 "device": "cpu"},
                acq_overrides={}, seed=a.seed, tag=sweep_id,
                universe_genes=universe,
                metrics=["hits_auc", "andcg", "batch_hits", "batch_diversity"],
                dataset_names=([x.strip() for x in a.only.split(",")] if a.only else None),
            )
            print(f"\n=== {sweep_id}: {d}/round x {a.rounds} rounds, "
                  f"kernel {'ON' if kernel_on else 'OFF'}")
            run_sweep(cfg, sweep_id=sweep_id, verbose=False, acq_factory=factory)
            n = len(list(sweep_root(a.runs_dir, sweep_id).glob(f"{sweep_id}-[0-9]*")))
            print(f"  {n} run dirs")


if __name__ == "__main__":
    main()
