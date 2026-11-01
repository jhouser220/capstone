"""LLM warm start -> AssayFormer handoff, evaluated the way the paper's table script does it.

`assayloop eval-ranker-handoff` (the documented CLI) differs from scripts/full_genome_table.py in two ways:
  1. it replays each warm-start round from the *library-filtered* acquired batch, whereas the table
     script re-parses the raw LLM answer against the shared f2 universe (21,147 genes for the
     curated 20-screen sets);
  2. it lets AssayFormer pick from each screen's own library, whereas the table script runs the
     ranker over the same shared universe the LLM used.

This script does (1) and (2) for any screen set and any warm-start sweep, by adding `universe_genes`
to the same RunConfig the CLI builds. It also refuses to run if any screen has no warm-start trace,
because such a screen would otherwise silently become a cold-start AssayFormer run.

    uv run python handoff_fullgenome.py --checkpoint "$CKPT" \\
        --screen-set paper_validation --warm-dir output/runs --warm-prefix sweep-gem35val- --n 3

Warm-start run folders must be in the standard layout: <warm-dir>/<prefix>*/result.json, with
llm_calls.jsonl beside it (used for the raw-answer replay).
"""
from __future__ import annotations

import argparse

from assayloop.acquisitions.glm_handoff_acq import GlmHandoffAcquisition
from assayloop.acquisitions.codep_seed_acq import CodepSeedAcquisition
from assayloop.acquisitions.greedy_from_model import GreedyFromModel
from assayloop.amortized.warmstart import WarmStart
from assayloop.experiment.runner import RunConfig, run_sweep
from assayloop.screen_direction import DROPOUT, ENRICHMENT, classify_selection, should_use_llm
from assayloop.tasks import gene_universe, load_screens


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--checkpoint", required=True, help="AssayFormer checkpoint directory")
    ap.add_argument("--ckpt-file", default="model_last.pt", help="weights file inside --checkpoint (paper: model_last.pt)")
    ap.add_argument("--screen-set", default="paper_validation")
    ap.add_argument("--screen", default="", help="comma-separated dataset names (subset of --screen-set)")
    ap.add_argument("--selection-mode", choices=["all", ENRICHMENT, DROPOUT], default="all",
                    help="run only screens of this selection direction (see assayloop/screen_direction.py)")
    ap.add_argument("--include-ambiguous", action="store_true",
                    help="with --selection-mode, also run screens whose direction is 'both' or 'unknown'")
    ap.add_argument("--warm-dir", default="output/runs")
    ap.add_argument("--warm-prefix", required=True, help="e.g. sweep-gem35val-")
    ap.add_argument("--n", type=int, default=3, help="warm-start rounds taken from the LLM (paper: 3 for Gemini/GPT)")
    ap.add_argument("--n-steps", type=int, default=10)
    ap.add_argument("--batch-size", type=int, default=100, help="must match the warm-start runs")
    ap.add_argument("--min-screen-freq", type=int, default=2, help="2 = the paper's f2 universe")
    # What runs AFTER the warm-start rounds. "greedy" is the paper's handoff and the
    # default, so this flag changes nothing unless it is passed. "codep_seed" hands the
    # LLM's confirmed hits to the fused ranker instead, which is the combined method:
    # the LLM seeds rounds 1-n, our kernel drives rounds n+1..10 off those seeds.
    ap.add_argument("--inner-acq", choices=["greedy", "codep_seed"], default="greedy")
    ap.add_argument("--profile-npz", default="/home/jeff/depmap/achilles_noeval.npz",
                    help="dependency atlas, for --inner-acq codep_seed")
    ap.add_argument("--train-hits", default="output/analysis/train_hits",
                    help="training hit matrix, for --inner-acq codep_seed")
    ap.add_argument("--weight", type=float, default=0.5)
    ap.add_argument("--knn-weight", type=float, default=0.5)
    ap.add_argument("--tag", default="")
    ap.add_argument("--device", default="auto", help="auto|cpu|cuda")
    ap.add_argument("--parallel", type=int, default=1)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--sweep-id", default=None, help="fixed sweep id (needed for --resume)")
    ap.add_argument("--resume", action="store_true", help="skip screens whose result.json already exists")
    args = ap.parse_args(argv)

    all_screens = load_screens(target_set=args.screen_set)
    # The universe always comes from the whole screen set, so a subset run stays comparable to a full one.
    universe = gene_universe(all_screens, min_screen_freq=args.min_screen_freq)
    names = [x for x in args.screen.split(",") if x]
    screens = [s for s in all_screens if not names or s.dataset_name in names]
    if args.selection_mode != "all":
        keep = [s for s in screens
                if should_use_llm(classify_selection(s.phenotype)[0], args.selection_mode, args.include_ambiguous)]
        for s in screens:
            d = classify_selection(s.phenotype)[0]
            print(f"  {s.dataset_name:15s} {d:10s} {'run' if s in keep else 'skipped'}")
        screens = keep
    print(f"{len(all_screens)} screens in {args.screen_set!r}; running {len(screens)}; "
          f"shared universe = {len(universe)} genes")
    if not screens:
        raise SystemExit("no screens selected")

    ws = WarmStart.from_run_dirs(
        args.warm_dir, args.warm_prefix, universe_genes=universe, batch_size=args.batch_size
    )
    missing = [s.dataset_name for s in screens if not ws.traces.get(s.dataset_name)]
    if missing:
        raise SystemExit(
            f"No warm-start trace for {len(missing)}/{len(screens)} screen(s): {missing}\n"
            f"Looked for {args.warm_dir}/{args.warm_prefix}*/result.json. Refusing to continue: those "
            "screens would silently become cold-start AssayFormer runs under a handoff label."
        )

    device = args.device
    if device == "auto":
        import torch
        device = "cuda" if torch.cuda.is_available() else "cpu"

    def make_inner():
        if args.inner_acq == "greedy":
            return GreedyFromModel(seed=args.seed)
        # One instance per screen: CodepSeedAcquisition keeps per-campaign state, so a
        # shared one would carry a previous screen's injection/trace into the next.
        return CodepSeedAcquisition(
            profile_npz=args.profile_npz, train_hits=args.train_hits,
            weight=args.weight, knn_weight=args.knn_weight,
            query="diffmeans", fuse="rrf", seed=args.seed)

    def acq_factory(screen):
        rounds = ws.traces[screen.dataset_name][0]          # first trace for this screen
        return GlmHandoffAcquisition(rounds, args.n, make_inner())

    cfg = RunConfig(
        screen_set=args.screen_set,
        dataset_names=[s.dataset_name for s in screens] if len(screens) != len(all_screens) else None,
        model="amortized_ranker",
        acq="greedy",
        batch_size=args.batch_size,
        n_steps=args.n_steps,
        seed=args.seed,
        model_overrides={"checkpoint": args.checkpoint, "device": device, "ckpt_file": args.ckpt_file},
        metrics=["hits_auc", "andcg", "batch_hits"],
        persist=True,
        parallel=args.parallel,
        tag=args.tag or f"handoff-fullgenome-n{args.n}-{args.inner_acq}",
        universe_genes=universe,
    )
    sweep = run_sweep(cfg, sweep_id=args.sweep_id, verbose=True, acq_factory=acq_factory, resume=args.resume)
    agg = sweep.aggregate or {}
    print(f"\nHandoff sweep {sweep.sweep_id}: mean_n_hits_vs_random={agg.get('mean_n_hits_vs_random')}, "
          f"mean_hits_auc={agg.get('mean_hits_auc')}")
    return sweep


if __name__ == "__main__":
    main()
