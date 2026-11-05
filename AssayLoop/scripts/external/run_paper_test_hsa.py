"""HSA arms on the paper's full 20-screen test set, using locally built banks.

This is the run that `W_text32` was thought to block. It does not: the donor side is
shared and frozen, and `mean=zeros(32)` is numerically immaterial (HANDOFF 24.10).
Banks come from `scripts/external/build_paper_test_banks.py`.

Settings match the project's standard for a paper_test comparison:
  * pool      f2 universe, `gene_universe(paper_test, min_screen_freq=2)` = 21,147 genes
  * budget    10 rounds x 100, seed 0
  * atlas     achilles_noeval.npz (reduced), per the standing rule
  * descriptions served from the shipped embedding cache -- a MISS ABORTS rather than
    calling the API, so this run costs nothing.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
from pathlib import Path

import numpy as np

from assaybench.data.screens import gene_universe
from assayloop.experiment.runner import RunConfig, run_one_screen
from assayloop.tasks.screen_sets import load_screens

from run_bundle18 import ARMS, NEEDS_BANK   # noqa: E402  (same directory)

log = logging.getLogger("run_paper_test_hsa")


def cached_descriptions(screens) -> dict[str, np.ndarray]:
    """{dataset_name: 1536-D vector} from the shipped cache. Never calls the API."""
    from assayloop.amortized.text_embed import (
        EMBED_MODEL, EmbeddingCache, _text_key, screen_description_text)
    cache = EmbeddingCache()
    out, missing = {}, []
    for s in screens:
        v = cache.get(_text_key(screen_description_text(s), model=EMBED_MODEL))
        if v is None:
            missing.append(s.dataset_name)
        else:
            out[s.dataset_name] = np.asarray(v, dtype=np.float32)
    if missing:
        raise SystemExit(
            f"{len(missing)} screen description(s) not in the shipped cache: {missing}. "
            "Refusing to call the embedding API. Embed them deliberately first.")
    dims = {v.shape for v in out.values()}
    if dims != {(1536,)}:
        raise SystemExit(f"unexpected embedding shapes: {dims}")
    return out


def make_hsa_scorer(bank_dir: Path, sid: str, model, desc, draw: int = 0):
    from experiments.odysseus.assayloop_comparison.adapter import Recipient, SourceState
    from experiments.odysseus.assayloop_comparison.assets import HistoricalBank
    bank = HistoricalBank.load(bank_dir / f"{sid}.npz")
    rec = Recipient("historical_only", model, SourceState(model, desc), bank, sid, draw)

    def score(observations, candidates, task_context):
        return rec.predict(observations, candidates, task_context).scores

    return score


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    ap = argparse.ArgumentParser()
    ap.add_argument("--banks", required=True, type=Path)
    ap.add_argument("--checkpoint", required=True, type=Path)
    ap.add_argument("--atlas", default="/home/jeff/depmap/achilles_noeval.npz")
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--arms", default="AF alone,AF + DepMap,HSA only (ours),"
                                     "AF + HSA,AF + DepMap + HSA")
    ap.add_argument("--screen-set", default="paper_test")
    ap.add_argument("--screens", default="", type=str,
                    help="pickle of ScreenRecords instead of a benchmark screen set. "
                         "`load_screens` can only return rows of the published AssayBench "
                         "dataset (see run_external.py), so externally built screens have to "
                         "come in this way. The pool then becomes EACH SCREEN'S OWN LIBRARY, "
                         "which is the external convention: no pick can land outside it, so "
                         "the forgiveness terms are zero and adjusted EF is identically "
                         "`n_hits_vs_random` -- one metric, no re-scoring (HANDOFF 24.33).")
    ap.add_argument("--per-screen-pool", action="store_true",
                    help="keep the benchmark screen set but make the pool EACH SCREEN'S OWN "
                         "LIBRARY instead of the f2 universe. The f2 pool exists only to put "
                         "model arms on the LLM's footing (`full_genome_table.py` docstring: "
                         "'Re-evaluate baselines with full-genome candidate pool for fair LLM "
                         "comparison'); with no LLM in the table it buys nothing and costs the "
                         "forgiveness asymmetry, which is what reorders AF + DepMap and dm.50 "
                         "(SUMMARY section 7 rule 1). With this flag no pick can land off-library, "
                         "so n2 = n3 = 0 and the two EF conventions are identically equal. "
                         "Absolute EF is NOT comparable to an f2 run -- different denominator, "
                         "section 7 rule 4. Combine with --restrict-to-bank, as the team11 run did.")
    ap.add_argument("--draw", type=int, default=0, choices=(0, 1))
    ap.add_argument("--screen", default="")
    ap.add_argument("--min-genes", type=int, default=0,
                    help="skip screens whose library is this small or smaller. REQUIRED with "
                         "--per-screen-pool on a whole fold: the 218-screen validation fold has "
                         "median library 2,318 genes and 98 screens at or under 1,001, where a "
                         "1,000-pick budget exhausts the library and EF is meaningless "
                         "(SUMMARY section 7 rule 6). 18000 keeps the genome-scale screens.")
    ap.add_argument("--restrict-to-bank", action="store_true",
                    help="intersect the pool with the bank gene axis, UNIFORMLY across every "
                         "arm. Needed for the 218-screen validation fold, whose f2 universe "
                         "holds 36 pseudogenes/lncRNAs absent from the bank (and from "
                         "AssayFormer's vocabulary). Dropping them from only the HSA arms "
                         "would invalidate the comparison, so it applies to all.")
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()

    axis = None
    if a.restrict_to_bank:
        import numpy as _np
        probe = sorted(a.banks.glob("*.npz"))[0]
        with _np.load(probe, allow_pickle=False) as _b:
            axis = set(_b["genes"].astype(str).tolist())

    per_screen_pool = bool(a.screens) or a.per_screen_pool
    universe = None              # per-screen pools resolve it inside the loop
    if a.screens:
        import pickle
        with open(a.screens, "rb") as fh:
            screens = pickle.load(fh)
        log.info("%s: %d externally built screens; pool = each screen's own library",
                 a.screens, len(screens))
    else:
        screens = load_screens(target_set=a.screen_set)
        if per_screen_pool:
            log.info("%s: %d screens; pool = EACH SCREEN'S OWN LIBRARY, no f2 universe",
                     a.screen_set, len(screens))
        else:
            universe = sorted(gene_universe(screens, min_screen_freq=2))
            log.info("%s: %d screens, f2 universe %d genes",
                     a.screen_set, len(screens), len(universe))
        if universe is not None and axis is not None:
            dropped = [g for g in universe if g not in axis]
            universe = [g for g in universe if g in axis]
            log.info("restrict-to-bank: dropped %d gene(s) -> pool %d (uniform across all arms)",
                     len(dropped), len(universe))
            log.info("  dropped: %s",
                     ", ".join(dropped[:40]) + (" ..." if len(dropped) > 40 else ""))
    if a.min_genes:
        keep = [s for s in screens if s.num_genes > a.min_genes]
        log.info("--min-genes %d: %d of %d screens kept (dropped %d with libraries too small "
                 "for a 1,000-pick budget)",
                 a.min_genes, len(keep), len(screens), len(screens) - len(keep))
        screens = keep
    desc_by = cached_descriptions(screens)
    log.info("all %d descriptions served from cache (0 API calls)", len(desc_by))
    if a.screen:
        screens = [s for s in screens if s.dataset_name == a.screen]

    from assayloop.experiment.runner import make_acquisition, make_model

    arms = [x.strip() for x in a.arms.split(",") if x.strip()]
    rows = []
    for i, rec in enumerate(screens):
        sid = rec.dataset_name
        desc = desc_by[sid]
        if per_screen_pool:
            # the screen's own library, intersected with the bank axis when asked --
            # applied to EVERY arm, never to the HSA arms alone
            lib_pool = set(rec.genes)
            if axis is not None:
                lib_pool &= axis
            universe = sorted(lib_pool)
            log.info("%-16s pool %d of %d library genes", sid, len(universe), rec.num_genes)
        model = make_model("amortized_ranker", overrides={
            "checkpoint": str(a.checkpoint), "ckpt_file": "model_last.pt",
            "device": a.device, "desc_emb_by_name": {sid: desc}}, seed=0)
        for p in model.net.parameters():
            p.requires_grad_(False)
        for arm in arms:
            over = dict(ARMS[arm])
            acq_name = over.pop("acq")
            if acq_name == "codep_seed":
                over["profile_npz"] = a.atlas
            if arm in NEEDS_BANK:
                over["hsa_scorer"] = make_hsa_scorer(a.banks, sid, model, desc, a.draw)
            acq = make_acquisition(acq_name, overrides=over, seed=0)
            cfg = RunConfig(
                model="amortized_ranker", acq=acq_name, batch_size=100, n_steps=10,
                warm_start_size=0, seed=0, model_overrides={},
                acq_overrides={k: v for k, v in over.items() if k != "hsa_scorer"},
                metrics=["hits_auc", "batch_hits"], universe_genes=list(universe),
                max_shortfall_frac=None, persist=False)
            res = run_one_screen(rec, cfg, run_id=f"pt-{i:02d}-{sid}",
                                 sweep_id="paper_test_hsa", verbose=False,
                                 model_obj=model, acq_obj=acq)
            lib = set(rec.genes)
            picks = [g for st in res.history for g in (st.acquired_batch or [])]
            assert len(set(picks)) == len(picks), f"{sid}: duplicate picks"
            unknown = sum(1 for g in picks if g not in lib)
            if per_screen_pool:
                # pool subset of library, so nothing can land off-library: both EF
                # conventions must coincide, and a non-zero count means they do not
                assert unknown == 0, f"{sid}/{arm}: {unknown} off-library picks from a library pool"
            m = res.history[-1].metrics
            rows.append(dict(arm=arm, screen=sid, pool=len(universe),
                             pool_mode="library" if per_screen_pool else "f2",
                             restricted=bool(a.restrict_to_bank),
                             n_genes=rec.num_genes, total_hits=rec.total_hits,
                             picks=len(picks), n_unknown=unknown,
                             n_measured=len(picks) - unknown, hits=int(m["n_hits"]),
                             ef=float(m.get("n_hits_vs_random", float("nan"))),
                             picks_sha=hashlib.sha256("\n".join(picks).encode()).hexdigest()[:16]))
            log.info("%-20s %-16s hits=%4d measured=%4d ef=%.3f",
                     arm, sid, rows[-1]["hits"], rows[-1]["n_measured"], rows[-1]["ef"])
        # flush after every screen: this run is long enough that losing it all to a
        # killed process would be the expensive failure
        a.out.parent.mkdir(parents=True, exist_ok=True)
        a.out.write_text(json.dumps(rows, indent=1))

    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps(rows, indent=1))
    with open(a.out.with_suffix(".csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {a.out} and {a.out.with_suffix('.csv')}  ({len(rows)} rows)")


if __name__ == "__main__":
    main()
