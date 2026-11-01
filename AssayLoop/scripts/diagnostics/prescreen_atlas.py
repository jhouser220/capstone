"""Generic three-gate pre-screen for any candidate atlas in the shared .npz layout.

Generalises `prescreen_ppi.py` so a new atlas costs one command. Runs, in increasing
cost, the gates §4.27 established:

  0. POSITIVE CONTROL on known protein complexes. Do co-complex genes cluster in this
     space more than random pairs? This exists because the Perturb-seq atlas once
     scored AUROC 0.49 -- a clean-looking null -- which turned out to be 6832
     infinities and `np.nan_to_num`'s default mapping them to 3.4e38 (§4.17). A null
     is only interpretable once the space is known to carry *any* functional signal.
  1. FAME BIAS. corr(feature, cross-screen hit frequency) must be near zero. "This
     gene is often a hit" is what `v_g` already encodes, so a correlated feature is
     probably inside span(V) -- the §4.11 LLM-score failure, stated mechanically.
     Reference values: PRISM ~0.00, BioPlex +0.203, condition-matching +0.595.
  2. COVERAGE of the hits that matter. §4.25: the hard screens' hits appear in a
     median of 9 of 1349 training screens. A feature defined only for well-studied
     genes cannot reach them whatever its margin. Reference: DepMap 83%, BioPlex 61%,
     Perturb-seq 46%, LINCS OE 13.5% (killed on this).
  3. INCREMENTAL VALUE at the oracle, split by difficulty, plus the decisive paired
     contrast against DepMap alone. An upper bound and a filter only: six offline
     margins have now overstated deployed results, one with a sign flip (§4.24).

BOTH CENTRING MODES are always reported. `gene` subtracts each gene's mean across
dimensions, right for a profile whose dimensions are comparable measurements of one
quantity (knockout effect per cell line). `dim` standardises each dimension across
genes, right for a learned or heterogeneous feature space (morphology features,
embedding coordinates). Getting it backwards produces a null for the wrong reason.

    uv run python scripts/diagnostics/prescreen_atlas.py \
        --atlas output/analysis/jump_orf.npz --label "JUMP ORF"
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import torch
from scipy.sparse import load_npz

sys.path.insert(0, "src")
sys.path.insert(0, "scripts/diagnostics")
from assayloop.metrics.hits_auc import adjusted_ef_value            # noqa: E402
from assayloop.models.amortized_ranker import AmortizedRankerModel  # noqa: E402
from assayloop.scripts.full_genome_table import _adj_ef_from_run    # noqa: E402
from assayloop.tasks import gene_universe, load_screens             # noqa: E402
from augment_u import fit_logistic                                  # noqa: E402
from score_ctx import boot_ci, families                             # noqa: E402
from screen_kernels import auroc                                    # noqa: E402
from screen_seed_aggs import CKPT, unit_profiles                    # noqa: E402
from run_paths import run_dir  # noqa: E402

BUCKETS = ((0.0, 2.0, "hard"), (2.0, 5.0, "mid"), (5.0, 1e9, "easy"))
COMPLEXES = {
    "ribosome": ["RPL3", "RPL4", "RPL5", "RPL6", "RPL7", "RPL8"],
    "RNA pol II": ["POLR2A", "POLR2B", "POLR2C", "POLR2D", "POLR2E"],
    "eIF3": ["EIF3A", "EIF3B", "EIF3C", "EIF3D", "EIF3E"],
    "proteasome": ["PSMA1", "PSMA2", "PSMA3", "PSMB1", "PSMB2"],
    "SF3B": ["SF3B1", "SF3B2", "SF3B3", "SF3B5"],
    "mediator": ["MED11", "MED17", "MED18", "MED20", "MED27"],
}


def load(path: Path, center: str, uni: list[str]):
    z = np.load(path, allow_pickle=True)
    vals = np.asarray(z["values"], dtype=np.float32)
    cols = {str(c): j for j, c in enumerate(z["cols"])}
    idx = np.array([cols.get(g, -1) for g in uni])
    ok = idx >= 0
    X = np.zeros((len(uni), vals.shape[0]), dtype=np.float64)
    # posinf/neginf -> 0, never numpy's 3.4e38 default (§4.17)
    X[ok] = np.nan_to_num(vals[:, idx[ok]].T, posinf=0.0, neginf=0.0)
    if center == "gene":
        X -= X.mean(axis=1, keepdims=True)
    else:
        mu = X[ok].mean(axis=0, keepdims=True)
        sd = X[ok].std(axis=0, keepdims=True)
        X[ok] = (X[ok] - mu) / np.where(sd > 0, sd, 1.0)
    n = np.linalg.norm(X, axis=1, keepdims=True)
    return X / np.where(n > 0, n, 1.0), ok


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--atlas", required=True)
    ap.add_argument("--label", default="candidate")
    ap.add_argument("--screen-set", default="validation")
    ap.add_argument("--exclude-set", default="paper_validation")
    ap.add_argument("--base", default="sweep-dm-validation-w0")
    ap.add_argument("--dir", default="output/runs")
    ap.add_argument("--train-hits", default="output/analysis/train_hits")
    ap.add_argument("--profile-npz", default="/home/jeff/depmap/achilles.npz")
    ap.add_argument("--min-library", type=int, default=1001)
    ap.add_argument("--l2", type=float, default=1.0)
    ap.add_argument("--noise-draws", type=int, default=5)
    a = ap.parse_args()

    m = AmortizedRankerModel(checkpoint=CKPT, ckpt_file="model_last.pt", device="cpu")
    screens = load_screens(target_set=a.screen_set)
    uni = sorted(gene_universe(screens, min_screen_freq=2))
    ids = torch.tensor([m.vocab.to_idx(g) for g in uni])
    V = m.net.gene_emb.weight[ids].detach().numpy().astype(np.float64)
    b = m.net.gene_bias[ids].detach().numpy().astype(np.float64)
    pos = {g: i for i, g in enumerate(uni)}
    drop = {x.dataset_name for x in load_screens(target_set=a.exclude_set)}

    P = {}
    for c in ("gene", "dim"):
        P[f"{a.label} [{c}-centred]"], okm = load(Path(a.atlas), c, uni)
    Pd, _ = unit_profiles(Path(a.profile_npz), uni)
    Pd = Pd.astype(np.float64)
    print(f"{a.label}: {int(okm.sum())} of {len(uni)} universe genes covered "
          f"({100 * okm.mean():.1f}%), {P[list(P)[0]].shape[1]} dims\n")

    print("GATE 0  positive control: do co-complex genes cluster?")
    rng = np.random.default_rng(0)
    for lab, X in P.items():
        ridx = rng.choice(np.flatnonzero(okm), 4000)
        rnd = (X[ridx[:2000]] * X[ridx[2000:]]).sum(axis=1)
        zs = []
        for nm, gs in COMPLEXES.items():
            gi = [pos[g] for g in gs if g in pos and okm[pos[g]]]
            if len(gi) < 2:
                continue
            C = X[gi] @ X[gi].T
            off = C[~np.eye(len(gi), dtype=bool)]
            zs.append((nm, (off.mean() - rnd.mean()) / (rnd.std() or 1)))
        print(f"  {lab:28s} random pair cosine {rnd.mean():+.3f} (sd {rnd.std():.3f}) | "
              + "  ".join(f"{n} {z:+.1f}sd" for n, z in zs))

    K = load_npz(Path(a.train_hits).with_suffix(".mat.npz")).tocsc()
    kmeta = np.load(Path(a.train_hits).with_suffix(".npz"), allow_pickle=True)
    kby = {str(g): j for j, g in enumerate(kmeta["gene_names"])}
    freq = np.asarray(K.sum(axis=0)).ravel()
    fame = np.array([freq[kby[g]] if g in kby else 0.0 for g in uni])

    acc, af, cov_hits = {}, {}, []
    for i, s in enumerate(screens):
        if s.dataset_name in drop or len(s.genes) <= a.min_library:
            continue
        run = run_dir(a.dir, a.base, i, s.dataset_name) / "result.json"
        if not run.is_file():
            continue
        hit_of = dict(zip(s.genes, s.hits))
        lib = set(s.genes)
        in_lib = np.array([g in lib for g in uni])
        y = np.array([1.0 if hit_of.get(g) else 0.0 for g in uni])
        steps = json.loads(run.read_text())["steps"][:3]
        seeds = [pos[g] for st in steps
                 for g, h in zip(st["acquired_batch"] or [], st["hits"] or []) if h and g in pos]
        miss = [pos[g] for st in steps
                for g, h in zip(st["acquired_batch"] or [], st["hits"] or []) if not h and g in pos]
        if len(seeds) < 2 or len(miss) < 10:
            continue

        def ef(sc):
            o = np.argsort(-sc)[:1000]
            pk = [uni[j] for j in o]
            n1 = sum(1 for g in pk if g in lib)
            hh = sum(1 for g in pk if hit_of.get(g))
            return adjusted_ef_value(hh, n1, len(pk) - n1, len(lib), s.total_hits, 1000)

        def fit(c):
            c = np.asarray(c, dtype=np.float64).reshape(len(uni), -1)
            sd = c.std(axis=0, keepdims=True)
            X = np.column_stack([V, c / np.where(sd > 0, sd, 1.0)])
            w, _ = fit_logistic(X[in_lib], y[in_lib], a.l2)
            return ef(X @ w + b)

        def direction(M):
            w = M[seeds].mean(axis=0) - M[miss].mean(axis=0)
            n = np.linalg.norm(w)
            return M @ (w / n) if n > 0 else np.zeros(len(uni))

        r = np.random.default_rng(abs(hash(s.dataset_name)) % 2**32)
        n1 = float(np.mean([fit(r.standard_normal(len(uni))) for _ in range(a.noise_draws)]))
        n2 = float(np.mean([fit(r.standard_normal((len(uni), 2))) for _ in range(a.noise_draws)]))
        dcol = direction(Pd)
        acc.setdefault("DepMap", {})[s.dataset_name] = fit(dcol) - n1
        for lab, X in P.items():
            col = direction(X)
            acc.setdefault(lab, {})[s.dataset_name] = fit(col) - n1
            acc.setdefault(f"DepMap + {lab}", {})[s.dataset_name] = \
                fit(np.column_stack([dcol, col])) - n2
        af[s.dataset_name] = _adj_ef_from_run(run, lib, set(uni), s.total_hits, len(lib), 1000)
        hm = in_lib & (y > 0)
        cov_hits.append((s.dataset_name, af[s.dataset_name], 100.0 * okm[hm].mean(),
                         auroc(P[list(P)[0]][in_lib] @ np.zeros(P[list(P)[0]].shape[1]) * 0
                               + direction(P[list(P)[0]])[in_lib], y[in_lib])))

    fam = families([s for s in screens if s.dataset_name in af])
    bucket = {n: lb for n, v in af.items() for lo, hi, lb in BUCKETS if lo <= v < hi}
    cnt = {lb: sum(1 for x in bucket.values() if x == lb) for _l, _h, lb in BUCKETS}
    hard = [n for n in af if bucket[n] == "hard"]
    print(f"\nGATE 2  coverage of HITS: median "
          f"{np.median([c[2] for c in cov_hits]):.1f}%  (hard screens median "
          f"{np.median([c[2] for c in cov_hits if c[1] < 2]):.1f}%)")
    print(f"\n{len(af)} non-degenerate held-out screens: "
          + ", ".join(f"{lb} {cnt[lb]}" for _l, _h, lb in BUCKETS))

    print("\nGATE 1  fame bias")
    for lab, X in P.items():
        col = X @ X[[pos[g] for g in uni[:1]]].mean(axis=0)   # arbitrary probe direction
        print(f"  {lab:28s} corr(|profile norm|, hit frequency) = "
              f"{np.corrcoef(np.linalg.norm(X[okm], axis=1), fame[okm])[0, 1]:+.3f}")

    print("\nGATE 3  oracle margin over width-matched noise (per-family)")
    hdr = f"{'kernel':30s} {'ALL':>8s} {'HARD':>8s} {'95% CI (hard)':>18s} {'mid':>8s} {'easy':>8s}"
    print(hdr)
    print("-" * len(hdr))
    for lab in ["DepMap", *P]:
        d = acc[lab]
        cells = {}
        for _lo, _hi, bl in BUCKETS:
            sl = [n for n in d if bucket.get(n) == bl]
            cells[bl] = boot_ci({n: d[n] for n in sl}, {n: fam[n] for n in sl})[:2] \
                if len(sl) >= 3 else (float("nan"), (float("nan"),) * 2)
        allm, _, _ = boot_ci(d, {n: fam[n] for n in d})
        hm, (hlo, hhi) = cells["hard"]
        print(f"{lab:30s} {allm:+8.3f} {hm:+8.3f} [{hlo:+6.3f}, {hhi:+6.3f}] "
              f"{cells['mid'][0]:+8.3f} {cells['easy'][0]:+8.3f}")

    print("\nDECISIVE: adding each to DepMap, paired, HARD screens only")
    for lab in P:
        pl = f"DepMap + {lab}"
        sl = [n for n in hard if n in acc[pl]]
        dd = {n: acc[pl][n] - acc["DepMap"][n] for n in sl}
        mm, (lo, hi), _ = boot_ci(dd, {n: fam[n] for n in sl})
        print(f"  {lab:30s} {mm:+.3f} [{lo:+.3f}, {hi:+.3f}]  "
              f"up {sum(1 for v in dd.values() if v > 0)}/{len(sl)}")


if __name__ == "__main__":
    main()
