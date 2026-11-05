"""GLM-5 warm start for rounds 1-3, then HSA / DepMap + HSA for the final 7.

Extends HANDOFF 4.61, which handed the same warm start to `AF alone` and to
`AF + DepMap`. Those two arms are re-run here as controls: if they reproduce 4.805
and 5.251, the replay is sound and the HSA columns are attributable.

THE WARM START IS REPLAYED FROM DISK -- rounds 1-3 are read verbatim out of
`sweep-glm5tu-pt-af-*/result.json`, so NO LLM IS CALLED and the run costs nothing.
That is asserted at the end by checking no `llm_calls.jsonl` appeared.

No driver survived for the original `sweep-dmne-hyb-glm5tu-pt` sweeps (their recorded
`config.json` says `acq: greedy` with no overrides, because the acquisition was built
in the driver and passed as `acq_obj`), so this rebuilds the mechanism on
`GlmHandoffAcquisition`, which is what it was written for and had never been run
(HANDOFF 24.31).

Pool, budget and atlas match the hybrid sweeps exactly: f2 universe of paper_test
(21,147 genes), 10 rounds x 100, seed 0, reduced atlas.

`--retain K` switches to the teammate's proposal-capped protocol: a warm round keeps
only its K best-scoring recorded proposals and the arm itself fills the other 100-K
slots, so the numerical method is live from round 1 and no slot goes unfilled. Their
package uses K=25. Still a replay -- the capped subset is a subset of proposals
already on disk, so nothing is billed either way.
"""
from __future__ import annotations

import argparse
import csv
import glob
import hashlib
import json
import logging
from pathlib import Path

from assaybench.data.screens import gene_universe
from assayloop.acquisitions.glm_handoff_acq import GlmHandoffAcquisition
from assayloop.experiment.runner import RunConfig, make_acquisition, make_model, run_one_screen
from assayloop.tasks.screen_sets import load_screens

from run_bundle18 import ARMS, NEEDS_BANK            # noqa: E402
from run_paper_test_hsa import cached_descriptions, make_hsa_scorer   # noqa: E402

log = logging.getLogger("glm5_handoff_hsa")


def warm_rounds(sweep: str, runs_dir: Path, n_warm: int) -> dict[str, list[list[str]]]:
    """{screen: [round1, round2, ...]} straight out of the recorded handoff arm.

    Rounds 1..n_warm of that arm ARE the LLM's own picks: its handoff is at k=3, so
    round 4 onwards is already AssayFormer's and must not be taken.
    """
    out = {}
    for d in sorted(glob.glob(str(runs_dir / "*" / f"{sweep}-*-*"))):
        d = Path(d)
        name = json.loads((d / "config.json").read_text())["screen_name"]
        steps = json.loads((d / "result.json").read_text())["steps"]
        out[name] = [list(s.get("acquired_batch") or []) for s in steps[:n_warm]]
    if not out:
        raise SystemExit(f"no runs matched {sweep} under {runs_dir}")
    supplied = sum(len(r) for rs in out.values() for r in rs)
    log.info("%s: %d screens, %d of %d warm-start picks supplied (%.1f per screen)",
             sweep, len(out), supplied, len(out) * n_warm * 100, supplied / len(out))
    return out


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    ap = argparse.ArgumentParser()
    ap.add_argument("--banks", required=True, type=Path)
    ap.add_argument("--checkpoint", required=True, type=Path)
    ap.add_argument("--warm-sweep", default="sweep-glm5tu-pt-af",
                    help="the recorded handoff arm whose rounds 1..n_warm are the LLM's")
    ap.add_argument("--runs-dir", default=Path("output/runs"), type=Path)
    ap.add_argument("--screen-set", default="paper_test")
    ap.add_argument("--n-warm", type=int, default=3)
    ap.add_argument("--retain", type=int, default=None,
                    help="keep at most this many proposals per warm round and let the "
                         "arm fill the rest (teammate's protocol uses 25); default: "
                         "keep them all and leave the round short")
    ap.add_argument("--rank", default="model", choices=["model", "llm"],
                    help="order the proposal pool by the model's score (teammate's "
                         "source/MAP ranking) or by the recorded reply order")
    ap.add_argument("--atlas", default="/home/jeff/depmap/achilles_noeval.npz")
    ap.add_argument("--device", default="cpu")
    ap.add_argument("--arms", default="AF alone,AF + DepMap,HSA only (ours),"
                                     "DepMap + HSA dm.25,DepMap + HSA dm.50")
    ap.add_argument("--screen", default="")
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()

    screens = load_screens(target_set=a.screen_set)
    universe = sorted(gene_universe(screens, min_screen_freq=2))
    log.info("%s: %d screens, pool %d genes", a.screen_set, len(screens), len(universe))
    warm = warm_rounds(a.warm_sweep, a.runs_dir, a.n_warm)
    desc_by = cached_descriptions(screens)
    log.info("all %d descriptions from cache (0 API calls)", len(desc_by))

    screens = [s for s in screens if s.dataset_name in warm]
    if a.screen:
        screens = [s for s in screens if s.dataset_name == a.screen]
    before = set(glob.glob("**/llm_calls.jsonl", recursive=True))

    arms = [x.strip() for x in a.arms.split(",") if x.strip()]
    rows = []
    for i, rec in enumerate(screens):
        sid = rec.dataset_name
        desc = desc_by[sid]
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
                over["hsa_scorer"] = make_hsa_scorer(a.banks, sid, model, desc, 0)
            base = make_acquisition(acq_name, overrides=over, seed=0)
            acq = GlmHandoffAcquisition(rounds=warm[sid], n_warm=a.n_warm, base=base,
                                        retain=a.retain, rank=a.rank)
            cfg = RunConfig(
                model="amortized_ranker", acq=acq_name, batch_size=100, n_steps=10,
                warm_start_size=0, seed=0, model_overrides={},
                acq_overrides={k: v for k, v in over.items() if k != "hsa_scorer"},
                metrics=["hits_auc", "batch_hits"], universe_genes=list(universe),
                max_shortfall_frac=None, persist=False)
            res = run_one_screen(rec, cfg, run_id=f"glmhsa-{i:02d}-{sid}",
                                 sweep_id="glm5_handoff_hsa", verbose=False,
                                 model_obj=model, acq_obj=acq)
            lib = set(rec.genes)
            picks = [g for st in res.history for g in (st.acquired_batch or [])]
            assert len(set(picks)) == len(picks), f"{sid}: duplicate picks"
            n_warm_picks = sum(len(r) for r in warm[sid])
            # what the warm rounds actually did, read back off the traces
            tr = [st.acquisition_trace or {} for st in res.history[:a.n_warm]]
            n_retained = sum(int(t.get("n_retained", t.get("n_picked", 0))) for t in tr)
            n_filled = sum(int(t.get("n_filled", 0)) for t in tr)
            m = res.history[-1].metrics
            rows.append(dict(arm=arm, screen=sid, pool=len(universe),
                             n_genes=rec.num_genes, total_hits=rec.total_hits,
                             n_warm=a.n_warm, warm_picks=n_warm_picks,
                             retain=a.retain if a.retain is not None else -1,
                             rank=a.rank if a.retain is not None else "",
                             llm_picks=n_retained, filled_picks=n_filled,
                             picks=len(picks),
                             n_unknown=sum(1 for g in picks if g not in lib),
                             n_measured=sum(1 for g in picks if g in lib),
                             hits=int(m["n_hits"]),
                             ef=float(m.get("n_hits_vs_random", float("nan"))),
                             picks_sha=hashlib.sha256(
                                 "\n".join(picks).encode()).hexdigest()[:16]))
            log.info("%-26s %-16s llm=%3d fill=%3d picks=%4d hits=%4d ef=%.3f",
                     arm, sid, n_retained, n_filled, len(picks),
                     rows[-1]["hits"], rows[-1]["ef"])
        a.out.parent.mkdir(parents=True, exist_ok=True)
        a.out.write_text(json.dumps(rows, indent=1))

    after = set(glob.glob("**/llm_calls.jsonl", recursive=True))
    assert after == before, f"LLM CALLS WERE MADE: {sorted(after - before)}"
    log.info("no llm_calls.jsonl written: the warm start was replayed, nothing was billed")

    a.out.write_text(json.dumps(rows, indent=1))
    with open(a.out.with_suffix(".csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {a.out} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
