"""Run AF alone and AF + DepMap over externally supplied ScreenRecords.

``run_sweep`` resolves its screens through ``load_screens``, and that loader
only ever selects rows of the published AssayBench dataset -- a YAML manifest
names the public dataset in its ``source`` block, so there is no route through
it for a screen that is not in the benchmark. ``run_one_screen`` takes a
``ScreenRecord`` directly, so that is the injection point this driver uses, and
the rest is the ordinary per-screen loop.

Both arms are ``--acq codep_seed`` and differ ONLY in ``weight``: 0.0 disables
the DepMap fusion and leaves AssayFormer's own ranking (bit-identical to
``--acq greedy``, which ``tests/test_codep_seed.py`` pins), 0.5 is the standard
method. Running the control through the same acquisition rather than through
``greedy`` keeps every other code path identical, so the difference between the
two columns is the DepMap term and nothing else.

``universe_genes`` is set to the library. All 14 screens share one Brunello
library, so the candidate pool IS each screen's own library: no pick can land
outside it, which makes the adjusted-EF forgiveness terms (``N_U`` for genes
outside this library, ``N_3`` for invalid symbols and unfilled slots) all zero.
Adjusted EF therefore collapses to ``hits / (budget * H/D)``, which is exactly
the ``n_hits_vs_random`` the runner already records -- so EF here needs no
reconstruction and cannot disagree with the notebook's convention.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import pickle
from pathlib import Path

from assayloop.experiment.runner import RunConfig, run_one_screen

log = logging.getLogger("run_external")

CKPT = ("/home/jeff/.cache/huggingface/hub/models--Genentech--assayformer/"
        "snapshots/4927703e5903984ae0254a9fe8e93ec0b29c27e4")

# Exactly the acq_overrides of sweep-dmne-cd-* (the standard AF + DepMap arm),
# minus train_hits, which codep_seed loads only when knn_weight != 0.
# "AF alone (greedy)" is the independent check that weight=0.0 really is
# AssayFormer untouched: it runs the plain greedy acquisition, which never loads
# the atlas at all. If the two AF-alone columns agree, the weight=0 control is
# sound and the DepMap delta is attributable.
ARMS = {
    "AF alone":          dict(acq="codep_seed", weight=0.0, knn_weight=0.0,
                              fuse="rrf", query="diffmeans"),
    "AF + DepMap":       dict(acq="codep_seed", weight=0.5, knn_weight=0.0,
                              fuse="rrf", query="diffmeans"),
    "AF alone (greedy)": dict(acq="greedy"),
}


def run_arm(screens, arm: str, atlas: str, device: str, tag: str,
            batch_size: int, n_steps: int, seed: int,
            universe: list[str] | None = None) -> list[dict]:
    rows = []
    overrides = dict(ARMS[arm])
    acq = overrides.pop("acq")
    if acq == "codep_seed":
        # Required by the constructor whatever the weight is; at weight=0.0 the
        # profiles are loaded and then contribute nothing, which is how
        # sweep-dm-validation-w0 ran its own AF-alone control.
        overrides["profile_npz"] = atlas
    for i, s in enumerate(screens):
        cfg = RunConfig(
            model="amortized_ranker",
            acq=acq,
            batch_size=batch_size,
            n_steps=n_steps,
            warm_start_size=0,
            seed=seed,
            model_overrides={"checkpoint": CKPT, "ckpt_file": "model_last.pt",
                             "device": device},
            acq_overrides=overrides,
            metrics=["hits_auc", "batch_hits"],
            # Default: the candidate pool is this library and nothing else --
            # see module docstring on why that makes adjusted EF unambiguous.
            # With --universe it is a FIXED cross-screen pool instead, so picks
            # can land on genes this screen never measured. GeneBatchTask reveals
            # those as ``hit=False, in_library=False`` and they still consume a
            # budget slot: the "unknown picks are charged" convention.
            universe_genes=list(universe) if universe else list(s.genes),
            max_shortfall_frac=None,
            persist=False,
        )
        run_id = f"{tag}-{arm.replace(' ', '').replace('+', '')}-{i:02d}-{s.dataset_name}"
        res = run_one_screen(s, cfg, run_id=run_id, sweep_id=tag, verbose=False)
        per_round = []
        for st in res.history:
            m = st.metrics
            per_round.append({
                "round": int(m["step"]),
                "picks": int(m["n_acquired"]),
                "cum_hits": int(m["n_hits"]),
                "round_hits": int(m.get("batch_n_hits", 0)),
                "ef": float(m.get("n_hits_vs_random", float("nan"))),
            })
        final = per_round[-1]
        lib = set(s.genes)
        n_unknown = sum(1 for st in res.history
                        for g in (st.acquired_batch or []) if g not in lib)
        rows.append({
            "arm": arm,
            "screen": s.dataset_name,
            "n_genes": s.num_genes,
            "n_unknown": n_unknown,
            "n_measured": final["picks"] - n_unknown,
            "total_hits": s.total_hits,
            "hits": final["cum_hits"],
            "picks": final["picks"],
            "hit_rate": final["cum_hits"] / max(final["picks"], 1),
            "recall": final["cum_hits"] / max(s.total_hits, 1),
            "ef": final["ef"],
            "per_round": per_round,
        })
        log.info("%-13s %-26s hits=%4d  rate=%.3f  EF=%.3f  measured=%d",
                 arm, s.dataset_name, final["cum_hits"],
                 rows[-1]["hit_rate"], final["ef"], rows[-1]["n_measured"])
    return rows


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    ap = argparse.ArgumentParser()
    ap.add_argument("--screens", required=True, type=Path, help="pickle of ScreenRecords")
    ap.add_argument("--atlas", default="/home/jeff/depmap/achilles_noeval.npz")
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--tag", default="schukken")
    ap.add_argument("--batch-size", type=int, default=100)
    ap.add_argument("--n-steps", type=int, default=10)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--arms", default="AF alone,AF + DepMap")
    ap.add_argument("--limit", type=int, default=0, help="first N screens only (smoke test)")
    ap.add_argument("--universe", type=Path, default=None,
                    help="JSON list of genes used as a FIXED candidate pool for every "
                         "screen, instead of each screen's own library. Picks landing "
                         "outside a screen's library are revealed as misses and still "
                         "cost a budget slot.")
    ap.add_argument("--screen-names", default="",
                    help="comma-separated dataset_names to keep")
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()

    screens = pickle.load(open(a.screens, "rb"))
    if a.screen_names:
        keep = {x.strip() for x in a.screen_names.split(",") if x.strip()}
        screens = [s for s in screens if s.dataset_name in keep]
        missing = keep - {s.dataset_name for s in screens}
        if missing:
            raise SystemExit(f"not in {a.screens}: {sorted(missing)}")
    if a.limit:
        screens = screens[: a.limit]
    universe = json.loads(a.universe.read_text()) if a.universe else None
    if universe:
        log.info("fixed candidate pool: %d genes", len(universe))
    log.info("%d screens, %d genes, atlas=%s", len(screens),
             screens[0].num_genes, Path(a.atlas).name)

    rows = []
    for arm in [x.strip() for x in a.arms.split(",") if x.strip()]:
        rows += run_arm(screens, arm, a.atlas, a.device, a.tag,
                        a.batch_size, a.n_steps, a.seed, universe)

    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(rows, indent=2))
    with open(a.out.with_suffix(".csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=[k for k in rows[0] if k != "per_round"])
        w.writeheader()
        for r in rows:
            w.writerow({k: v for k, v in r.items() if k != "per_round"})
    print(f"\nwrote {a.out} and {a.out.with_suffix('.csv')}")


if __name__ == "__main__":
    main()
