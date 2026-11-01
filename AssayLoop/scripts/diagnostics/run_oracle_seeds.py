"""ORACLE-SEED BOUND: would a better NOMINATOR help `codep_seed`? (§4.48 -> §4.49)

The question. §4 killed the LLM as a picker (thirteen variants) and §4.11 killed it as
an additive feature (zero incremental over span(V)). Neither covers the LLM as a
NOMINATOR feeding the kernel: it changes which genes get tried, the kernel changes how
candidates are ranked, so the pair is not redundant by construction. The synergy story
is that the LLM surfaces ~12 hits per screen AssayFormer never reaches (§2e), those
become seeds, and the co-essentiality direction is pointed somewhere new.

Rather than pay for LLM calls, bound the whole family with an oracle -- the pattern that
settled §4.7, §4.9 and §4.11. Hand the scorers genes the ANSWER KEY says are hits and
measure the deployed gain. If even a perfect nominator is flat, no real one can help and
the family closes for free.

TWO PROPERTIES MAKE THIS A BOUND AND NOT A CHEAT, and both are enforced:

  1. Oracle seeds are REMOVED FROM THE CANDIDATE POOL (`CodepSeedAcquisition` does this),
     so the arm cannot bank them as hits. Otherwise the result is just "we gave it N
     free hits", which is trivially positive and means nothing.
  2. The seeds are drawn only from hits that NEITHER comparison arm ever acquired --
     neither the weight-0 baseline nor the deployed configuration. Those arms score zero
     from those genes, so denying them to the oracle arm costs it nothing relative to
     them, and the contrast stays fair in both directions.

ARMS, pre-specified, a dose-response plus one control:
  orc12     12 unreachable hits  -- the LLM's own measured new-hit rate (§2e)
  orc50     50
  orc200   200                   -- far beyond any real nominator; the upper bound
  rnd50     50 random NON-hit genes -- the control. If a nonsense direction scores like
            a true one, the instrument is measuring perturbation, not information.

    uv run python scripts/diagnostics/run_oracle_seeds.py --arm orc12 --arm rnd50
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
from assayloop.tasks import load_screens                               # noqa: E402
from run_paths import run_dir, sweep_root  # noqa: E402

CKPT = ("/home/jeff/.cache/huggingface/hub/models--Genentech--assayformer/"
        "snapshots/4927703e5903984ae0254a9fe8e93ec0b29c27e4")
ARMS = {"orc12": (12, False), "orc50": (50, False),
        "orc200": (200, False), "rnd50": (50, True)}


def acquired(run_dir: Path) -> set[str]:
    p = run_dir / "result.json"
    if not p.is_file():
        return set()
    rd = json.loads(p.read_text())
    return {str(g) for st in rd.get("steps", []) for g in (st.get("acquired_batch") or [])}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--screen-set", default="validation")
    ap.add_argument("--arm", action="append", required=True, choices=sorted(ARMS))
    ap.add_argument("--runs-dir", default="output/runs")
    ap.add_argument("--base", default="sweep-dm-validation-w0")
    ap.add_argument("--ref", default="sweep-dm-validation-dm05knn")
    ap.add_argument("--depmap", default="/home/jeff/depmap/achilles_noeval.npz",
                    help="DepMap gene-effect matrix. Default since 2026-10-28 is the "
                         "ABLATED atlas: 768 of 808 Achilles lines, every line used by "
                         "any evaluation set dropped. Pass achilles.npz for all 808.")
    ap.add_argument("--train-hits", default="output/analysis/train_hits")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--only", default=None,
                    help="comma-separated dataset_names, for a smoke test")
    a = ap.parse_args()

    screens = load_screens(target_set=a.screen_set)
    rng = np.random.default_rng(a.seed)
    from assayloop.tasks import gene_universe
    # cli.py computes this over the SCREEN SET, not over all screens. Using the wrong
    # one silently changes the candidate pool and voids the contrast against the arms
    # already on disk, which were run at 21213 genes for `validation`.
    universe = gene_universe(load_screens(target_set=a.screen_set), min_screen_freq=2)
    print(f"candidate universe: {len(universe)} genes (f2 over {a.screen_set}); "
          f"the arms on disk used 21213")

    # per-screen pool of TRUE hits that neither comparison arm ever touched
    pools: dict[str, list[str]] = {}
    nonhits: dict[str, list[str]] = {}
    for i, s in enumerate(screens):
        seen = set()
        for sweep in (a.base, a.ref):
            seen |= acquired(run_dir(a.runs_dir, sweep, i, s.dataset_name))
        hit_of = dict(zip(s.genes, s.hits))
        pools[s.dataset_name] = sorted(g for g in s.genes
                                       if hit_of.get(g) and g not in seen)
        nonhits[s.dataset_name] = sorted(g for g in s.genes if not hit_of.get(g))
    sizes = np.array([len(v) for v in pools.values()])
    print(f"{a.screen_set}: {len(screens)} screens. Unreachable true hits per screen "
          f"(never acquired by {a.base} or {a.ref}):")
    print(f"  median {np.median(sizes):.0f}, quartiles {np.percentile(sizes,25):.0f}/"
          f"{np.percentile(sizes,75):.0f}, min {sizes.min()}, max {sizes.max()}; "
          f"{(sizes == 0).sum()} screens have none")
    if a.dry_run:
        return

    for arm in a.arm:
        n_seed, use_random = ARMS[arm]
        sweep_id = f"sweep-orc-{a.screen_set}-{arm}"
        if list(sweep_root(a.runs_dir, sweep_id).glob(f"{sweep_id}-[0-9]*")):
            print(f"SKIP {sweep_id} (run dirs already present)")
            continue

        def factory(screen, n_seed=n_seed, use_random=use_random):
            pool = nonhits[screen.dataset_name] if use_random \
                else pools[screen.dataset_name]
            k = min(n_seed, len(pool))
            picked = [] if k == 0 else list(
                rng.choice(np.array(pool, dtype=object), size=k, replace=False))
            return CodepSeedAcquisition(
                profile_npz=a.depmap, train_hits=a.train_hits,
                weight=0.5, knn_weight=0.5, query="diffmeans", fuse="rrf",
                oracle_seeds=[str(g) for g in picked])

        # `--full-genome --min-screen-freq 2` is a CLI concept: it resolves to
        # `universe_genes`, the paper's f2 universe, computed over ALL screens (not just
        # this set) exactly as cli.py does it, so the candidate pool matches every other
        # arm on disk. Getting this wrong would change the pool and void the contrast.
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
        print(f"\n=== {sweep_id}: {n_seed} "
              f"{'RANDOM NON-HIT' if use_random else 'unreachable true hit'} seeds")
        res = run_sweep(cfg, sweep_id=sweep_id, verbose=False, acq_factory=factory)
        print(f"  done: {getattr(res, 'n_ok', '?')} ok / "
              f"{getattr(res, 'n_failed', '?')} failed")


if __name__ == "__main__":
    main()
