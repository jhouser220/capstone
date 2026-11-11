from __future__ import annotations

import argparse, json, hashlib
from pathlib import Path
import numpy as np
import pandas as pd

from assaybench.data.screens import gene_universe
from assayloop.tasks.screen_sets import load_screens
from assayloop.experiment.runner import RunConfig, make_model, make_acquisition, run_one_screen
from assayloop.acquisitions.greedy_from_model import GreedyFromModel
from assayloop.acquisitions.glm_handoff_acq import GlmHandoffAcquisition
from assayloop.amortized.text_embed import EMBED_MODEL, EmbeddingCache, _text_key, screen_description_text

from experiments.odysseus.assayloop_comparison.adapter import Recipient, SourceState
from experiments.odysseus.assayloop_comparison.assets import HistoricalBank

# Compatibility patch:
# Jeff's newer AssayLoop runner calls model.coverage() before starting a run.
# Odysseus's HSA Recipient predates that interface addition, so delegate
# coverage reporting to the underlying AssayFormer source model.
if not hasattr(Recipient, "coverage"):
    def _recipient_coverage(self, candidates):
        return self.source.coverage(candidates)
    Recipient.coverage = _recipient_coverage


def cached_description(screen):
    cache = EmbeddingCache()
    key = _text_key(screen_description_text(screen), model=EMBED_MODEL)
    v = cache.get(key)
    if v is None:
        raise RuntimeError(
            f"{screen.dataset_name}: description embedding missing from shipped cache. "
            "Refusing to call an embedding API."
        )
    return np.asarray(v, dtype=np.float32)


def source_model(checkpoint: Path, sid: str, desc, device: str):
    model = make_model(
        "amortized_ranker",
        overrides={
            "checkpoint": str(checkpoint),
            "ckpt_file": "model_last.pt",
            "device": device,
            "desc_emb_by_name": {sid: desc},
        },
        seed=0,
    )
    for p in model.net.parameters():
        p.requires_grad_(False)
    return model


def run_arm(rec, universe, rounds, checkpoint, bank_dir, arm, device, atlas=None):
    sid = rec.dataset_name
    desc = cached_description(rec)
    src = source_model(checkpoint, sid, desc, device)

    if arm == "AssayFormer":
        model_obj = src
        base = GreedyFromModel(seed=0, epsilon=0)

    elif arm == "HSA":
        bank = HistoricalBank.load(bank_dir / f"{sid}.npz")
        model_obj = Recipient(
            "historical_only",
            src,
            SourceState(src, desc),
            bank,
            sid,
            draw=0,
        )
        base = GreedyFromModel(seed=0, epsilon=0)

    elif arm == "DepMap":
        if atlas is None:
            raise ValueError("DepMap arm requested without an atlas.")
        model_obj = src
        base = make_acquisition(
            "codep_seed",
            overrides={
                "profile_npz": str(atlas),
                "weight": 0.5,
                "knn_weight": 0.0,
                "query": "diffmeans",
                "fuse": "rrf",
            },
            seed=0,
        )

    else:
        raise ValueError(arm)

    acq = GlmHandoffAcquisition(
        rounds=rounds,
        n_warm=3,
        base=base,
        retain=None,
        rank="model",
    )

    cfg = RunConfig(
        model="amortized_ranker",
        acq="greedy",
        batch_size=100,
        n_steps=10,
        warm_start_size=0,
        seed=0,
        model_overrides={},
        acq_overrides={},
        metrics=["hits_auc", "batch_hits"],
        universe_genes=list(universe),
        max_shortfall_frac=None,
        persist=False,
    )

    res = run_one_screen(
        rec,
        cfg,
        run_id=f"chris-gpt56-{arm}-{sid}",
        sweep_id="chris_gpt56_handoff",
        verbose=False,
        model_obj=model_obj,
        acq_obj=acq,
    )

    picks = [g for st in res.history for g in (st.acquired_batch or [])]
    assert len(picks) == len(set(picks)), f"{sid}/{arm}: duplicate picks"

    lib = set(rec.genes)
    measured = sum(g in lib for g in picks)
    unknown = len(picks) - measured
    m = res.history[-1].metrics
    hits = int(m["n_hits"])

    round_hits = []
    for st in res.history:
        sm = st.metrics or {}
        round_hits.append(int(sm.get("n_hits", 0)))

    return {
        "screen": sid,
        "arm": arm,
        "n_genes": int(rec.num_genes),
        "total_hits": int(rec.total_hits),
        "hit_density": float(rec.total_hits / rec.num_genes),
        "picks": len(picks),
        "n_measured": int(measured),
        "n_unknown": int(unknown),
        "hits": hits,
        "recall": float(hits / rec.total_hits) if rec.total_hits else np.nan,
        "ef": float(m.get("n_hits_vs_random", np.nan)),
        "round_cumulative_hits": json.dumps(round_hits),
        "picks_sha": hashlib.sha256("\n".join(picks).encode()).hexdigest()[:16],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rounds", required=True, type=Path)
    ap.add_argument("--banks", required=True, type=Path)
    ap.add_argument("--checkpoint", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--atlas", default="", type=str)
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--screen", default="")
    args = ap.parse_args()

    warm = json.loads(args.rounds.read_text())
    screens = load_screens(target_set="paper_test")
    universe = sorted(gene_universe(screens, min_screen_freq=2))
    screen_by = {s.dataset_name: s for s in screens}

    usable = [sid for sid in warm if sid in screen_by]
    if args.screen:
        usable = [sid for sid in usable if sid == args.screen]

    arms = ["AssayFormer", "HSA"]
    atlas = Path(args.atlas) if args.atlas else None
    if atlas is not None:
        arms.append("DepMap")

    rows = []
    args.out.parent.mkdir(parents=True, exist_ok=True)

    print(f"Usable screens: {len(usable)} | pool={len(universe)} | arms={arms}")
    for i, sid in enumerate(sorted(usable), 1):
        print(f"\n[{i}/{len(usable)}] {sid}")
        rec = screen_by[sid]
        for arm in arms:
            row = run_arm(
                rec=rec,
                universe=universe,
                rounds=warm[sid],
                checkpoint=args.checkpoint,
                bank_dir=args.banks,
                arm=arm,
                device=args.device,
                atlas=atlas,
            )
            rows.append(row)
            print(
                f"  {arm:12s} hits={row['hits']:4d} "
                f"recall={row['recall']:.4f} ef={row['ef']:.3f}"
            )

        pd.DataFrame(rows).to_csv(args.out, index=False)

    df = pd.DataFrame(rows)
    df.to_csv(args.out, index=False)
    print(f"\nSaved {len(df)} rows to {args.out}")


if __name__ == "__main__":
    main()
