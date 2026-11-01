"""OFFLINE SCREEN: are the DepMap kernel and HSA's 32 historical directions redundant?

The gate before building a combined arm (HANDOFF section 17 scope). On a few bundle
screens it replays the RECORDED HSA trajectory so both scorers see identical state, and
at every round computes three candidate score vectors on the same candidate set:

    af   AssayFormer's own score            (model.predict)
    cd   our DepMap co-essentiality kernel  (codep_seed's own _P / _row / _seed_query)
    hsa  HSA's posterior score q*prob       (THEIR Recipient, arm='historical_only')

`hsa` comes from their own code, not a reimplementation, so the comparison cannot be
wrong about what HSA scores. `cd` is read out of a live CodepSeedAcquisition's loaded
state, so it is the same kernel the deployed arm uses.

What decides the build: if `fused(af, cd)` already selects the genes `hsa` selects, the
two signals are near-substitutes and a combined arm cannot add much.
"""
import gzip, json, sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.stats import spearmanr

import argparse

_ap = argparse.ArgumentParser(description=__doc__)
_ap.add_argument("bundle", help="path to _reproduction/bundle")
_ap.add_argument("atlas", help="DepMap co-essentiality .npz")
_ap.add_argument("screens", help="comma-separated screen ids")
_ap.add_argument("--out", type=Path, default=None,
                 help="write the per-round rows as JSON (the printed table is a summary)")
_a = _ap.parse_args()
B = Path(_a.bundle)
ATLAS = _a.atlas
SCREENS = _a.screens.split(",")
RRF_K = 60.0

from assaybench.data.screens import ScreenRecord
from assayloop.models.amortized_ranker import AmortizedRankerModel
from assayloop.acquisitions.codep_seed_acq import CodepSeedAcquisition
from assayloop.tasks.gene_batch import AssayBenchGeneBatchTask
from experiments.odysseus.assayloop_comparison.adapter import SourceState, Recipient
from experiments.odysseus.assayloop_comparison.assets import HistoricalBank

folds = {e["cell"]["screen_id"]: e["fold"]
         for e in json.loads((B / "manifest.json").read_text())}


def ranks(order, n):
    r = np.full(n, np.inf)
    for pos, i in enumerate(order):
        r[i] = pos
    return r


def rrf_term(v, n):
    ok = np.isfinite(v)
    if not ok.any():
        return np.zeros(n)
    order = [int(i) for i in np.argsort(-np.where(ok, v, -np.inf)) if ok[i]]
    r = ranks(order, n)
    t = np.where(np.isfinite(r), 1.0 / (RRF_K + r), np.nan)
    t[~np.isfinite(t)] = np.nanmean(t)
    return t


def top(v, k=100):
    ok = np.isfinite(v)
    idx = np.argsort(-np.where(ok, v, -np.inf))
    return set(int(i) for i in idx[:k])


records = []
print(f"{'screen':14s} {'rd':>3s} {'cand':>6s} {'cov':>6s} | "
      f"{'rho(cd,hsa)':>12s} {'rho(af,hsa)':>12s} {'rho(af,cd)':>11s} | "
      f"{'top100 cd^hsa':>13s} {'fused^hsa':>10s} {'af^hsa':>7s}")
print("-" * 118)
summary = {}
for sid in SCREENS:
    rec = ScreenRecord(**json.loads((B / "screens" / f"{sid}.json").read_text()))
    desc = np.load(B / "screens" / f"{sid}.npy", allow_pickle=False)
    uni = json.loads((B / f"UNIVERSE_{folds[sid]}.json").read_text())
    bank = HistoricalBank.load(B / "banks" / f"{sid}.npz")
    traj = json.loads(gzip.decompress((B / "recorded" / f"{sid}__scalar__none__0.json.gz").read_bytes()))

    model = AmortizedRankerModel(checkpoint=str(B / "checkpoint"), ckpt_file="model_last.pt",
                                 device="cpu", desc_emb_by_name={sid: desc})
    for p in model.net.parameters():
        p.requires_grad_(False)
    acq = CodepSeedAcquisition(profile_npz=ATLAS, weight=0.5, knn_weight=0.0,
                               fuse="rrf", query="diffmeans", seed=0)
    task = AssayBenchGeneBatchTask(rec, universe_genes=list(uni), seed=0)
    ctx = task.context()
    recip = Recipient("historical_only", model, SourceState(model, desc), bank, sid, 0)

    obs_all, acquired, hist, rows = [], set(), [], []
    for rd in traj["rounds"]:
        cands = [g for g in uni if g not in acquired]
        af_d = model.predict(obs_all, cands, ctx).scores
        af = np.array([af_d.get(c, -np.inf) for c in cands], float)
        q, nseed = acq._seed_query(hist)
        cd = np.full(len(cands), np.nan)
        if q is not None:
            have = np.array([acq._row.get(str(c), -1) for c in cands])
            cov = have >= 0
            cd[cov] = acq._P[have[cov]] @ q
        hsa_d = recip.predict(obs_all, cands, ctx).scores
        hsa = np.array([hsa_d[c] for c in cands], float)

        fused = 1.0 / (RRF_K + ranks(sorted(range(len(cands)),
                       key=lambda i: (-af[i], str(cands[i]))), len(cands)))
        if q is not None:
            fused = fused + 0.5 * rrf_term(cd, len(cands))

        ok = np.isfinite(cd)
        r_ch = spearmanr(cd[ok], hsa[ok])[0] if ok.sum() > 2 else float("nan")
        r_ah = spearmanr(af, hsa)[0]
        r_ac = spearmanr(af[ok], cd[ok])[0] if ok.sum() > 2 else float("nan")
        t_ch = len(top(cd) & top(hsa)) if ok.sum() > 2 else -1
        t_fh = len(top(fused) & top(hsa))
        t_ah = len(top(af) & top(hsa))
        rows.append((r_ch, r_ah, r_ac, t_ch, t_fh, t_ah))
        records.append(dict(screen=sid, fold=folds[sid], round=rd["round"],
                            n_candidates=len(cands), n_covered=int(ok.sum()),
                            rho_cd_hsa=None if np.isnan(r_ch) else float(r_ch),
                            rho_af_hsa=float(r_ah),
                            rho_af_cd=None if np.isnan(r_ac) else float(r_ac),
                            top100_cd_hsa=int(t_ch), top100_fused_hsa=int(t_fh),
                            top100_af_hsa=int(t_ah)))
        print(f"{sid:14s} {rd['round']:3d} {len(cands):6d} {int(ok.sum()):6d} | "
              f"{r_ch:12.3f} {r_ah:12.3f} {r_ac:11.3f} | {t_ch:13d} {t_fh:10d} {t_ah:7d}")

        for o in rd["observations"]:
            obs_all.append(SimpleNamespace(candidate=o["candidate"], label=o["label"],
                                           metadata=o["metadata"]))
        hist.append(SimpleNamespace(new_observations=obs_all[-len(rd["observations"]):]))
        acquired.update(rd["accepted"])
    a = np.array([[x if x is not None else np.nan for x in r] for r in rows], float)
    summary[sid] = np.nanmean(a, axis=0)
    print("-" * 118)

print(f"\n{'screen':14s} {'rho(cd,hsa)':>12s} {'rho(af,hsa)':>12s} {'rho(af,cd)':>11s} | "
      f"{'cd^hsa':>7s} {'fused^hsa':>10s} {'af^hsa':>7s}   (means over 10 rounds)")
for sid, m in summary.items():
    print(f"{sid:14s} {m[0]:12.3f} {m[1]:12.3f} {m[2]:11.3f} | "
          f"{m[3]:7.1f} {m[4]:10.1f} {m[5]:7.1f}")

if _a.out:
    _a.out.parent.mkdir(parents=True, exist_ok=True)
    _a.out.write_text(json.dumps(records, indent=1))
    print(f"\nwrote {_a.out}  ({len(records)} screen-rounds)")
