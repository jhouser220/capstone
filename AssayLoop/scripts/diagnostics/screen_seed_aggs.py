"""Compare seed-aggregation rules and screen-kNN, all as single added features.

Two ideas, screened before any campaign is built:

A. The deployed `codep_seed` scores a gene by its correlation with the MEAN of the
   confirmed hits' DepMap profiles. That centroid is a poor query when a screen's
   hits span several modules: the average of two unrelated profiles resembles
   neither. `max` and `top3` keep the multi-module structure.

C. Screen-kNN. AssayFormer's v_g is a 10-number compression of 1349 training
   screens; §4.9 found only 3.5% of a screen's hit structure survives it. This
   scores a gene by how often it is a hit in the training screens whose hit sets
   most resemble what this screen has revealed -- the raw matrix, uncompressed.
   `models/screen_knn.py` does this as a standalone POLICY, where it is a weak
   baseline; nothing has tested it as an ADDITIVE feature, and this project has
   twice seen a feature that loses as a replacement still add as a term.

Every candidate is one column, so capacity is not a confound, but each is still
reported against its own width-matched Gaussian control.
"""
from __future__ import annotations

import argparse, json, sys
from pathlib import Path

import numpy as np
import torch
from scipy.sparse import load_npz

sys.path.insert(0, "src")
sys.path.insert(0, "scripts/diagnostics")
from assayloop.metrics.hits_auc import adjusted_ef_value            # noqa: E402
from assayloop.models.amortized_ranker import AmortizedRankerModel  # noqa: E402
from assayloop.tasks import gene_universe, load_screens             # noqa: E402
from augment_u import fit_logistic                                  # noqa: E402
from score_ctx import boot_ci, families                             # noqa: E402
from run_paths import run_dir  # noqa: E402

CKPT = ("/home/jeff/.cache/huggingface/hub/models--Genentech--assayformer/"
        "snapshots/4927703e5903984ae0254a9fe8e93ec0b29c27e4")
MAX_SEEDS = 150      # cap for the pairwise rules; a full matmul over 900 seeds is too slow


def unit_profiles(npz: Path, uni):
    z = np.load(npz, allow_pickle=True)
    vals, cols = np.asarray(z["values"], dtype=np.float32), [str(c) for c in z["cols"]]
    by = {g: j for j, g in enumerate(cols)}
    idx = np.array([by.get(g, -1) for g in uni])
    ok = idx >= 0
    P = np.zeros((len(uni), vals.shape[0]), dtype=np.float32)
    # posinf/neginf -> 0, not numpy's default 3.4e38, which would make the norm of
    # any row holding one infinity astronomical and its correlation with every query
    # ~0. DepMap has none; the Replogle Perturb-seq atlas has ~6.8k.
    P[ok] = np.nan_to_num(vals[:, idx[ok]].T, posinf=0.0, neginf=0.0)
    P -= P.mean(axis=1, keepdims=True)
    n = np.linalg.norm(P, axis=1, keepdims=True)
    return P / np.where(n > 0, n, 1.0), ok


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile-npz", required=True)
    ap.add_argument("--train-hits", default="output/analysis/train_hits",
                    help="prefix written by build_train_hits.py (expects <p>.npz and <p>.mat.npz)")
    ap.add_argument("--screen-set", default="paper_validation")
    ap.add_argument("--base-sweep", default="sweep-fd2017be")
    ap.add_argument("--base-dir", default="output/runs/baseline/gemini")
    ap.add_argument("--seed-rounds", type=int, default=3)
    ap.add_argument("--knn-k", type=int, default=50)
    ap.add_argument("--l2", type=float, default=1.0)
    ap.add_argument("--noise-draws", type=int, default=5,
                    help="Gaussian control blocks to average per feature. ONE draw is not enough: "
                         "a single control is itself a noisy estimate of the capacity baseline, and "
                         "changing only the RNG consumption order moved these margins by ~0.1 EF.")
    a = ap.parse_args()

    m = AmortizedRankerModel(checkpoint=CKPT, ckpt_file="model_last.pt", device="cpu")
    screens = load_screens(target_set=a.screen_set)
    uni = sorted(gene_universe(screens, min_screen_freq=2))
    ids = torch.tensor([m.vocab.to_idx(g) for g in uni])
    V = m.net.gene_emb.weight[ids].detach().numpy().astype(np.float64)
    b = m.net.gene_bias[ids].detach().numpy().astype(np.float64)
    pos = {g: i for i, g in enumerate(uni)}
    P, _covered = unit_profiles(Path(a.profile_npz), uni)

    # screen-kNN inputs: the raw training hit matrix, aligned to this universe
    K = Kcols = None
    tp = Path(a.train_hits)
    if tp.with_suffix(".mat.npz").is_file():
        K = load_npz(tp.with_suffix(".mat.npz")).tocsc()
        meta = np.load(tp.with_suffix(".npz"), allow_pickle=True)
        kby = {str(g): j for j, g in enumerate(meta["gene_names"])}
        Kcols = np.array([kby.get(g, -1) for g in uni])
        print(f"screen-kNN: {K.shape[0]} training screens x {K.shape[1]} genes; "
              f"{int((Kcols >= 0).sum())} of {len(uni)} universe genes present")
    else:
        print(f"screen-kNN skipped: run build_train_hits.py first ({tp}.mat.npz missing)")

    acc: dict[str, dict[str, float]] = {}
    for i, s in enumerate(screens):
        hit_of = dict(zip(s.genes, s.hits))
        lib = set(s.genes)
        in_lib = np.array([g in lib for g in uni])
        y = np.array([1.0 if hit_of.get(g) else 0.0 for g in uni])
        run = run_dir(a.base_dir, a.base_sweep, i, s.dataset_name) / "result.json"
        if not run.is_file():
            continue
        steps = json.loads(run.read_text())["steps"][:a.seed_rounds]
        seeds = [pos[g] for st in steps
                 for g, h in zip(st["acquired_batch"] or [], st["hits"] or [])
                 if h and g in pos]
        nonseeds = [pos[g] for st in steps
                    for g, h in zip(st["acquired_batch"] or [], st["hits"] or [])
                    if not h and g in pos]
        if len(seeds) < 2:
            continue
        use = seeds[:MAX_SEEDS]

        feats: dict[str, np.ndarray] = {}
        q = P[use].mean(axis=0)
        nq = np.linalg.norm(q)
        feats["A centroid (deployed)"] = P @ (q / nq) if nq else np.zeros(len(uni))
        S = P @ P[use].T                                 # (genes, seeds) similarity
        feats["A max-over-seeds"] = S.max(axis=1)
        k3 = min(3, S.shape[1])
        feats["A mean-top3-seeds"] = np.sort(S, axis=1)[:, -k3:].mean(axis=1)

        # --- idea 1: difference of means. The centroid above is a nearest-centroid
        # --- classifier using only the positive class; this is the standard two-class
        # --- centroid difference, the Fisher direction under isotropic covariance.
        # --- It uses ~300 labels instead of ~20 and still needs no fitting.
        def _proj(w):
            n = np.linalg.norm(w)
            return P @ (w / n) if n > 0 else np.zeros(len(uni))

        if len(nonseeds) >= 10:
            m1 = P[use].mean(axis=0)
            m0 = P[nonseeds[:MAX_SEEDS * 4]].mean(axis=0)
            feats["1 diff-of-means"] = _proj(m1 - m0)

            # --- idea 3: scale each cell-line dimension by its variance among the
            # --- labelled genes (nearest shrunken centroid / diagonal LDA), which
            # --- down-weights directions that vary a lot regardless of hit status.
            lab = np.concatenate([P[use], P[nonseeds[:MAX_SEEDS * 4]]])
            var = lab.var(axis=0)
            feats["3 diag-shrunk diff"] = _proj((m1 - m0) / (var + np.median(var)))

        # --- idea 5: the seeds' own dominant mode of VARIATION, after centring them.
        # --- If a screen's hits span two modules the centroid points between them
        # --- while this direction separates them -- the multi-module hypothesis that
        # --- max-over-seeds failed to capture, tested properly this time.
        if len(use) >= 3:
            Q = P[use] - P[use].mean(axis=0)
            try:
                _u, _sv, vt = np.linalg.svd(Q, full_matrices=False)
                feats["5 seed-PC1"] = _proj(vt[0])
            except np.linalg.LinAlgError:
                pass

        # --- idea 5 combined: several fixed directions, no per-screen fitting.
        # --- seed-PC1 measured NEGATIVE on both sets, so the 3-way combo has to be
        # --- compared against the 2-way ones to see whether it wins because of PC1
        # --- or in spite of it.
        def _combo(name, keys):
            cols = [feats[k] for k in keys if k in feats]
            if len(cols) == len(keys) >= 2:
                feats[name] = np.column_stack(cols)

        _combo("5 combo: cen+diff", ("A centroid (deployed)", "1 diff-of-means"))
        _combo("5 combo: cen+shrunk", ("A centroid (deployed)", "3 diag-shrunk diff"))
        _combo("5 combo: diff+shrunk", ("1 diff-of-means", "3 diag-shrunk diff"))
        _combo("5 combo: 3 directions",
               ("A centroid (deployed)", "1 diff-of-means", "5 seed-PC1"))
        _combo("5 combo: cen+diff+shrunk",
               ("A centroid (deployed)", "1 diff-of-means", "3 diag-shrunk diff"))

        if K is not None:
            seed_genes = [uni[j] for j in seeds]
            kidx = np.array([kby[g] for g in seed_genes if g in kby])
            if len(kidx) >= 2:
                overlap = np.asarray(K[:, kidx].sum(axis=1)).ravel()
                sz = np.asarray(K.sum(axis=1)).ravel()
                sim = overlap / np.sqrt(np.maximum(sz, 1))      # overlap, size-corrected
                top = np.argsort(-sim)[:a.knn_k]
                w = sim[top]
                if w.sum() > 0:
                    prof = np.asarray(K[top].T.dot(w)).ravel() / w.sum()
                    v = np.zeros(len(uni))
                    ok = Kcols >= 0
                    v[ok] = prof[Kcols[ok]]
                    feats["C screen-kNN"] = v

        # the decisive question for C: does screen-kNN ADD to the centroid, or just
        # restate it? A 2-column block beside a 2-column noise control answers it.
        if "C screen-kNN" in feats:
            feats["A+C centroid & kNN"] = np.column_stack(
                [feats["A centroid (deployed)"], feats["C screen-kNN"]])

        def ef(score):
            order = np.argsort(-score)[:1000]
            picks = [uni[j] for j in order]
            n1 = sum(1 for g in picks if g in lib)
            hits = sum(1 for g in picks if hit_of.get(g))
            return adjusted_ef_value(hits, n1, len(picks) - n1, len(lib), s.total_hits, 1000)

        w0, _ = fit_logistic(V[in_lib], y[in_lib], a.l2)
        acc.setdefault("plain", {})[s.dataset_name] = ef(V @ w0 + b)
        for name, phi in feats.items():
            phi = np.asarray(phi, dtype=np.float64)
            phi = phi[:, None] if phi.ndim == 1 else phi

            def _fit(col):
                c = np.asarray(col, dtype=np.float64)
                sd = c.std(axis=0, keepdims=True)
                X = np.column_stack([V, c / np.where(sd > 0, sd, 1.0)])
                w, _ = fit_logistic(X[in_lib], y[in_lib], a.l2)
                return ef(X @ w + b)

            acc.setdefault(name, {})[s.dataset_name] = _fit(phi)
            # a fresh generator per (screen, feature) so the draws do not depend on
            # how many features ran before, which previously moved results by ~0.1
            nrng = np.random.default_rng(abs(hash((s.dataset_name, name, phi.shape[1]))) % 2**32)
            acc.setdefault(f"{name} [noise]", {})[s.dataset_name] = float(np.mean(
                [_fit(nrng.standard_normal(phi.shape)) for _ in range(a.noise_draws)]))

    kept = [s for s in screens if s.dataset_name in acc.get("plain", {})]
    fam = families(kept)
    plain = acc["plain"]
    print(f"\n{a.screen_set}: {len(plain)} screens, {len(set(fam.values()))} families, "
          f"plain oracle {np.mean(list(plain.values())):.3f}\n")
    hdr = f"{'added feature':24s} {'NET of noise':>13s} {'95% CI':>18s} {'wins':>7s}"
    print(hdr); print("-" * len(hdr))
    for name in [k for k in acc if k not in ("plain",) and not k.endswith("[noise]")]:
        net = {n: acc[name][n] - acc[f"{name} [noise]"][n] for n in acc[name]}
        mm, (lo, hi), _ = boot_ci(net, fam)
        print(f"{name:24s} {mm:+13.3f} [{lo:+6.3f}, {hi:+6.3f}] "
              f"{sum(v > 0 for v in net.values()):3d}/{len(net):<3d}")


if __name__ == "__main__":
    main()
