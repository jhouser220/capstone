"""Measure what an atlas revision changes, on exactly the screens it can touch.

§21 found four hTERT-RPE1 subclones surviving in `achilles_noeval.npz` that 61 evaluation
screens are grown from, and `achilles_noeval2.npz` drops them. Only those 61 screens can
move: for every other screen the two atlases are identical, so re-running them would be
wasted compute and a wider panel would dilute the effect with guaranteed zeros.

Re-runs ONE arm (AF + DepMap, `weight=0.5`) against the new atlas with the stored run's
exact settings -- same screen set, same f2 pool, same seed, `device=cpu`, `train_hits`,
`query=diffmeans`, `fuse=rrf` -- and pairs each result against the recorded v1 run. The
AF-alone baseline needs no re-run: it carries `weight=0.0`, so the atlas is loaded and
contributes nothing, and the two atlases give bit-identical AssayFormer rankings.

EF is recomputed for BOTH arms here, by the notebook's prefix accounting, so the two
columns cannot disagree through a metric difference.

    python scripts/diagnostics/score_atlas_revision.py \
        --atlas ~/depmap/achilles_noeval2.npz --out output/runs/external/atlas_rev.json
"""
from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "src")
sys.path.insert(0, "scripts/diagnostics")
from assayloop.experiment.runner import RunConfig, run_one_screen      # noqa: E402
from assayloop.metrics.hits_auc import adjusted_ef_value               # noqa: E402
from assayloop.tasks import gene_universe, load_screens                # noqa: E402
from run_paths import run_dir                                          # noqa: E402
from screen_depmap import line_map, norm, resolve                      # noqa: E402

log = logging.getLogger("score_atlas_revision")
CKPT = ("/home/jeff/.cache/huggingface/hub/models--Genentech--assayformer/"
        "snapshots/4927703e5903984ae0254a9fe8e93ec0b29c27e4")
# the stored arms this reproduces, per fold
ARM = {"validation": "sweep-dmne-cd-validation", "test": "sweep-dmne-cd-test"}
BASE = {"validation": "sweep-dm-validation-w0", "test": "sweep-tf-test-w0"}
BUDGET = 1000


def ef_from_run(path: Path, lib: set, uni: set, total_hits: int, domain: int):
    rd = json.loads(path.read_text())
    n1 = n2 = hits = 0
    for st in rd.get("steps", []):
        for g, h in zip(st.get("acquired_batch") or [], st.get("hits") or []):
            if g in lib:
                n1 += 1
                hits += bool(h)
            elif g in uni:
                n2 += 1
    return hits, n1, n2, adjusted_ef_value(hits, n1, n2, domain, total_hits, BUDGET)


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    logging.getLogger("httpx").setLevel(logging.WARNING)
    ap = argparse.ArgumentParser()
    ap.add_argument("--atlas", required=True)
    ap.add_argument("--depmap-dir", default="/home/jeff/depmap")
    ap.add_argument("--train-hits", default="output/analysis/train_hits")
    ap.add_argument("--runs", default="output/runs")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()

    by = line_map(Path(a.depmap_dir) / "Model.csv")
    rows = []
    for fold in ("validation", "test"):
        screens = load_screens(target_set=fold)
        uni = set(gene_universe(screens, min_screen_freq=2))
        # the screens the revision can touch: unresolved cell line mentioning RPE1
        todo = [(i, s) for i, s in enumerate(screens)
                if resolve(s.cell_line, by) is None and "RPE1" in norm(s.cell_line)]
        log.info("%s: %d affected screens of %d", fold, len(todo), len(screens))
        if a.limit:
            todo = todo[: a.limit]
        for i, s in todo:
            lib = set(s.genes)
            old = run_dir(a.runs, ARM[fold], i, s.dataset_name) / "result.json"
            if not old.is_file():
                log.warning("no v1 run for %s/%s, skipping", fold, s.dataset_name)
                continue
            h0, n10, n20, ef0 = ef_from_run(old, lib, uni, s.total_hits, len(lib))
            cfg = RunConfig(
                model="amortized_ranker", acq="codep_seed",
                batch_size=100, n_steps=10, warm_start_size=0, seed=0,
                model_overrides={"checkpoint": CKPT, "ckpt_file": "model_last.pt",
                                 "device": "cpu"},
                acq_overrides={"profile_npz": a.atlas, "fuse": "rrf",
                               "train_hits": a.train_hits, "query": "diffmeans",
                               "weight": 0.5, "knn_weight": 0.0},
                metrics=["hits_auc", "batch_hits"],
                universe_genes=sorted(uni), max_shortfall_frac=None, persist=False,
            )
            res = run_one_screen(s, cfg, run_id=f"atlasrev-{fold}-{i:02d}-{s.dataset_name}",
                                 sweep_id="atlas-revision", verbose=False)
            n1 = n2 = hits = 0
            for st in res.history:
                # StepRecord carries Observations, not a parallel hits list
                for obs in st.new_observations or []:
                    g = obs.candidate
                    if g in lib:
                        n1 += 1
                        lbl = obs.label
                        hits += bool(lbl.get("hit") if isinstance(lbl, dict) else lbl)
                    elif g in uni:
                        n2 += 1
            ef1 = adjusted_ef_value(hits, n1, n2, len(lib), s.total_hits, BUDGET)
            rows.append(dict(fold=fold, screen=s.dataset_name, cell_line=s.cell_line,
                             total_hits=s.total_hits, n_genes=len(lib),
                             hits_v1=h0, hits_v2=hits, ef_v1=ef0, ef_v2=ef1,
                             d_hits=hits - h0, d_ef=ef1 - ef0))
            log.info("%-14s %-16s hits %4d -> %4d (%+d)   EF %7.3f -> %7.3f (%+.3f)",
                     fold, s.dataset_name, h0, hits, hits - h0, ef0, ef1, ef1 - ef0)

    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(rows, indent=1))
    d_ef = np.array([r["d_ef"] for r in rows])
    d_h = np.array([r["d_hits"] for r in rows])
    rng = np.random.default_rng(0)
    bs = d_ef[rng.integers(0, len(d_ef), size=(20000, len(d_ef)))].mean(axis=1)
    print(f"\n{len(rows)} affected screens re-run against {Path(a.atlas).name}")
    print(f"  dEF   mean {d_ef.mean():+.4f}  [{np.percentile(bs,2.5):+.4f}, "
          f"{np.percentile(bs,97.5):+.4f}]   changed on {int((d_ef != 0).sum())}/{len(rows)}")
    print(f"  dhits mean {d_h.mean():+.3f}  min {d_h.min():+d}  max {d_h.max():+d}")
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
