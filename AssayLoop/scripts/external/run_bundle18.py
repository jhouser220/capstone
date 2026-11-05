"""AF alone and AF + DepMap on the 18 screens of the HSA reproduction bundle.

Puts our standard method on EXACTLY the footing the HSA comparison used, so the
three columns differ only in the acquisition:

  * same screens        -- `bundle/screens/<id>.json`, their frozen ScreenRecords
  * same candidate pool -- `bundle/UNIVERSE_{train,validation,test}.json`, picked per
    screen by its manifest `fold`. These are gene-for-gene identical to
    `gene_universe(paper_validation|paper_test|public_train, min_screen_freq=2)`
    (HANDOFF section 15.1), but note they are NOT the pools our own
    `sweep-dmne-cd-{validation,test}` sweeps used -- those ran on the 218- and
    334-screen FOLD universes (21,213 / 22,656). Hence this re-run.
  * same description vector -- `bundle/screens/<id>.npy`, 1,536-D L2-normalised,
    injected through `desc_emb_by_name` so no text embedding is recomputed and no
    API call is made.
  * same checkpoint     -- `bundle/checkpoint`, byte-identical to our cached
    `Genentech/assayformer` revision 4927703e.
  * same budget         -- 10 rounds x 100, seed 0.
  * same metric         -- unique MEASURED hits in 1000 selections. Off-library picks
    are revealed `hit=False, in_library=False` and still consume a slot, which is the
    "unknowns are charged" convention the HSA numbers use.

`AF alone` here is `codep_seed` at weight 0, which is the nested control our sweeps
use. `AF alone (greedy)` is the independent check that weight 0 really is AssayFormer
untouched. Both should reproduce the bundle's own `source` arm exactly; if they do,
the bridge is sound and the `AF + DepMap` column is attributable to the DepMap term.

DepMap atlas: `achilles_noeval.npz` (reduced), per this project's standing rule.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
from pathlib import Path

import numpy as np

from assaybench.data.screens import ScreenRecord
from assayloop.experiment.runner import RunConfig, run_one_screen

log = logging.getLogger("run_bundle18")

# Every codep_seed arm shares one recipe and differs only in the three ensemble
# weights, so the columns are attributable term by term:
#   af_weight   AssayFormer's own rank          (1.0 = historical behaviour)
#   weight      DepMap co-essentiality kernel
#   hsa_weight  the 42-D historical-subspace scorer, supplied as a callable
# HSA already contains AssayFormer internally (its prior mean is the screen latent
# `z` and its first ten design columns are `W[ids]`), so arms that keep af_weight=1
# alongside hsa_weight>0 double-count it; `(no AF term)` is the variant that does not.
_CD = dict(acq="codep_seed", fuse="rrf", query="diffmeans", knn_weight=0.0)
ARMS = {
    "AF alone":                  dict(_CD, weight=0.0),
    "AF + DepMap":               dict(_CD, weight=0.5),
    "AF + HSA":                  dict(_CD, weight=0.0, hsa_weight=0.5),
    "AF + DepMap + HSA":         dict(_CD, weight=0.5, hsa_weight=0.5),
    "AF + DepMap + HSA .25":     dict(_CD, weight=0.5, hsa_weight=0.25),
    "HSA only (ours)":           dict(_CD, weight=0.0, hsa_weight=1.0, af_weight=0.0),
    "DepMap + HSA (no AF term)": dict(_CD, weight=0.5, hsa_weight=1.0, af_weight=0.0),
    # Equal-weight variants: at dm=0.5 / hsa=1.0 the HSA term outweighs DepMap 2:1 in
    # the RRF sum, so DepMap can only reorder the tail. These give it parity and more.
    "HSA + DepMap eq":           dict(_CD, weight=1.0, hsa_weight=1.0, af_weight=0.0),
    "HSA + DepMap dm2":          dict(_CD, weight=2.0, hsa_weight=1.0, af_weight=0.0),
    "AF alone (greedy)":         dict(acq="greedy"),
    # DepMap weight sweep against hsa=1.0 with NO duplicate AF rank term. HSA already
    # contains AssayFormer (prior centre `z`, design columns 0-9 = W[ids]), so af_weight=0
    # is the non-double-counting form. The argument is structural; the measured value is
    # small -- +0.033 EF on paper_test in the project's canonical adjusted EF. (HANDOFF 24.12
    # quoted +0.305, which is `n_hits_vs_random` units and overstates it 9x; see 24.29.)
    "DepMap + HSA dm.25":        dict(_CD, weight=0.25, hsa_weight=1.0, af_weight=0.0),
    "DepMap + HSA dm.50":        dict(_CD, weight=0.50, hsa_weight=1.0, af_weight=0.0),
    "DepMap + HSA dm.75":        dict(_CD, weight=0.75, hsa_weight=1.0, af_weight=0.0),
    "DepMap + HSA dm1.0":        dict(_CD, weight=1.00, hsa_weight=1.0, af_weight=0.0),
}
NEEDS_BANK = {k for k, v in ARMS.items() if v.get("hsa_weight")}


def make_hsa_scorer(bundle: Path, sid: str, model, desc, draw: int = 0,
                    zero_mean: bool = False):
    """Wrap the reference `historical_only` recipient as a per-round score callable.

    Imported lazily and from the caller's own path: their package is not a dependency
    of ours, and only these arms need it. `Recipient.predict` counts its own rounds,
    so a FRESH one is built per (screen, arm).
    """
    from experiments.odysseus.assayloop_comparison.adapter import Recipient, SourceState
    from experiments.odysseus.assayloop_comparison.assets import HistoricalBank
    bank = HistoricalBank.load(bundle / "banks" / f"{sid}.npz")
    if zero_mean:
        # Ablate the only per-target array. `mean` is a ~0.01-logit offset inside a
        # prior 47 units wide (HANDOFF section 24.10), so this asks whether the
        # frozen `W_text32` is needed at all to run HSA on a new screen.
        bank.mean = np.zeros_like(bank.mean)
    rec = Recipient("historical_only", model, SourceState(model, desc), bank, sid, draw)

    def score(observations, candidates, task_context):
        return rec.predict(observations, candidates, task_context).scores

    return score


def load_bundle(bundle: Path):
    """[(screen_id, fold, cohort, ScreenRecord, desc)], one per screen."""
    folds = {}
    for e in json.loads((bundle / "manifest.json").read_text()):
        folds[e["cell"]["screen_id"]] = (e["fold"], e["cohort"])
    out = []
    for p in sorted((bundle / "screens").glob("*.json")):
        sid = p.stem
        rec = ScreenRecord(**json.loads(p.read_text()))
        desc = np.load(bundle / "screens" / f"{sid}.npy", allow_pickle=False)
        fold, cohort = folds[sid]
        out.append((sid, fold, cohort, rec, desc))
    return out


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    ap = argparse.ArgumentParser()
    ap.add_argument("--bundle", required=True, type=Path)
    ap.add_argument("--atlas", default="/home/jeff/depmap/achilles_noeval.npz")
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--arms", default="AF alone,AF + DepMap,AF alone (greedy)")
    ap.add_argument("--draw", type=int, default=0, choices=(0, 1),
                    help="repetition index handed to the HSA scorer (it pins 0/1)")
    ap.add_argument("--screen", default="", help="one screen id (smoke test)")
    ap.add_argument("--zero-mean", action="store_true",
                    help="zero bank.mean in every HSA arm (ablates the text prior)")
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()

    screens = load_bundle(a.bundle)
    if a.screen:
        screens = [s for s in screens if s[0] == a.screen]
    universes = {f: json.loads((a.bundle / f"UNIVERSE_{f}.json").read_text())
                 for f in sorted({s[1] for s in screens})}
    for f, u in universes.items():
        log.info("universe %s: %d genes", f, len(u))

    from assayloop.experiment.runner import make_acquisition, make_model

    arms = [x.strip() for x in a.arms.split(",") if x.strip()]
    rows = []
    for i, (sid, fold, cohort, rec, desc) in enumerate(screens):
        # one model per screen, reused by every arm: it is the expensive object and
        # it depends only on the screen's description vector
        model = make_model("amortized_ranker", overrides={
            "checkpoint": str(a.bundle / "checkpoint"),
            "ckpt_file": "model_last.pt", "device": a.device,
            "desc_emb_by_name": {sid: desc}}, seed=0)
        for p in model.net.parameters():
            p.requires_grad_(False)
        for arm in arms:
            over = dict(ARMS[arm])
            acq_name = over.pop("acq")
            if acq_name == "codep_seed":
                over["profile_npz"] = a.atlas
            if arm in NEEDS_BANK:
                over["hsa_scorer"] = make_hsa_scorer(a.bundle, sid, model, desc, a.draw,
                                                     zero_mean=a.zero_mean)
            acq = make_acquisition(acq_name, overrides=over, seed=0)
            cfg = RunConfig(
                model="amortized_ranker",
                acq=acq_name,
                batch_size=100,
                n_steps=10,
                warm_start_size=0,
                seed=0,
                model_overrides={},
                acq_overrides={k: v for k, v in over.items() if k != "hsa_scorer"},
                metrics=["hits_auc", "batch_hits"],
                universe_genes=list(universes[fold]),
                max_shortfall_frac=None,
                persist=False,
            )
            res = run_one_screen(rec, cfg, run_id=f"b18-{i:02d}-{sid}",
                                 sweep_id="bundle18", verbose=False,
                                 model_obj=model, acq_obj=acq)
            lib = set(rec.genes)
            picks = [g for st in res.history for g in (st.acquired_batch or [])]
            assert len(set(picks)) == len(picks), f"{sid}: duplicate picks"
            unknown = sum(1 for g in picks if g not in lib)
            m = res.history[-1].metrics
            rows.append(dict(arm=arm, screen=sid, fold=fold, cohort=cohort,
                             pool=len(universes[fold]), n_genes=rec.num_genes,
                             total_hits=rec.total_hits,
                             picks=len(picks), n_unknown=unknown,
                             n_measured=len(picks) - unknown,
                             hits=int(m["n_hits"]),
                             ef=float(m.get("n_hits_vs_random", float("nan"))),
                             zero_mean=bool(a.zero_mean),
                             picks_sha=hashlib.sha256(
                                 "\n".join(picks).encode()).hexdigest()[:16],
                             pick_list=picks))
            log.info("%-18s %-16s hits=%4d measured=%4d ef=%.3f",
                     arm, sid, rows[-1]["hits"], rows[-1]["n_measured"], rows[-1]["ef"])

    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(rows, indent=1))
    with open(a.out.with_suffix(".csv"), "w", newline="") as fh:
        cols = [c for c in rows[0] if c != "pick_list"]
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {a.out} and {a.out.with_suffix('.csv')}  ({len(rows)} rows)")


if __name__ == "__main__":
    main()
