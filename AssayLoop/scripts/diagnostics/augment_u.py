"""Can capacity be ADDED to AssayFormer's score without displacing what it knows?

The released score is  a(g) = u.v_g + b_g,  u in R^10. Over a universe of N genes
the reachable score vectors form a 10-dimensional affine subspace
    S = { V u + b : u in R^10 },   V = gene_emb (N x 10).
The oracle ceiling (oracle_u.py) is the best ranking in S, so it bounds anything
that works by improving u. Appending features is a different move:
    a(g) = u.v_g + b_g + sum_k w_k phi_k(g)
which lives in a (10+k)-dimensional subspace STRICTLY CONTAINING S, so it is NOT
bounded by that ceiling -- w = 0 recovers the original exactly.

Whether it pays depends only on whether phi has a component orthogonal to
span(V, b, 1) that correlates with the labels. This script measures that two ways:
  * R^2 of the label vector regressed on span(V, b, 1)  -- how much of "which
    genes are hits" is even expressible in the current basis;
  * the AUGMENTED oracle ceiling with phi appended, against the plain one.
"""
from __future__ import annotations

import argparse, json, sys
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, "src")
from assayloop.metrics.hits_auc import adjusted_ef_value            # noqa: E402
from assayloop.models.amortized_ranker import AmortizedRankerModel  # noqa: E402
from assayloop.tasks import gene_universe, load_screens             # noqa: E402
from run_paths import run_dir  # noqa: E402

CKPT = ("/home/jeff/.cache/huggingface/hub/models--Genentech--assayformer/"
        "snapshots/4927703e5903984ae0254a9fe8e93ec0b29c27e4")
GMT = Path("src/assayloop/data/gene_sets/ReactomePathways.gmt")


def modules(lo=2, hi=150):
    out = []
    for line in GMT.read_text().splitlines():
        p = line.rstrip("\n").split("\t")
        gs = {g for g in p[2:] if g}
        if lo <= len(gs) <= hi:
            out.append(gs)
    return out


def pathway_feature(mods, seeds, weights, uni):
    """score(g) = sum over modules with >=2 seeds of (seed evidence)/|module|.

    ``weights`` maps a seed gene to its revealed relevance; pass all-ones for the
    unweighted version. This is the strength-aware variant of the propagation
    scorer, now entering as an ADDITIVE term rather than a competing ranking.
    """
    sset = set(seeds)
    sc = {}
    for gs in mods:
        hit = sset & gs
        if len(hit) < 2:
            continue
        w = sum(weights.get(g, 1.0) for g in hit) / len(gs)
        for g in gs:
            sc[g] = sc.get(g, 0.0) + w
    v = np.array([sc.get(g, 0.0) for g in uni])
    return v / (v.std() or 1.0)


def fit_logistic(X, y, l2=1.0, iters=400):
    X = torch.as_tensor(X, dtype=torch.float32)
    y = torch.as_tensor(y, dtype=torch.float32)
    w = torch.zeros(X.shape[1], requires_grad=True)
    c = torch.zeros(1, requires_grad=True)
    opt = torch.optim.LBFGS([w, c], max_iter=iters, line_search_fn="strong_wolfe")
    lf = torch.nn.BCEWithLogitsLoss()

    def closure():
        opt.zero_grad()
        loss = lf(X @ w + c, y) + l2 * (w @ w) / len(y)
        loss.backward()
        return loss
    opt.step(closure)
    return w.detach().numpy(), float(c.detach())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--screen-set", default="paper_validation")
    ap.add_argument("--base-sweep", default="sweep-ctx-base")
    ap.add_argument("--base-dir", default="output/runs")
    ap.add_argument("--seed-rounds", type=int, default=3)
    ap.add_argument("--l2", type=float, default=1.0)
    a = ap.parse_args()

    m = AmortizedRankerModel(checkpoint=CKPT, ckpt_file="model_last.pt", device="cpu")
    screens = load_screens(target_set=a.screen_set)
    uni = sorted(gene_universe(screens, min_screen_freq=2))
    ids = torch.tensor([m.vocab.to_idx(g) for g in uni])
    V = m.net.gene_emb.weight[ids].detach().numpy().astype(np.float64)
    b = m.net.gene_bias[ids].detach().numpy().astype(np.float64)
    mods = modules()

    # span(V, b, 1): the subspace every AssayFormer score lives in
    Base = np.column_stack([V, b, np.ones(len(uni))])
    Q, _ = np.linalg.qr(Base)

    names = ["oracle relevance", "pathway", "pathway x strength", "gene_bias", "random"]
    rows = {k: [] for k in ["plain", *names]}
    r2s = []
    rng = np.random.default_rng(0)

    for i, s in enumerate(screens):
        hit_of = dict(zip(s.genes, s.hits))
        rel_of = dict(zip(s.genes, s.relevance_scores))
        lib = set(s.genes)
        in_lib = np.array([g in lib for g in uni])
        y = np.array([1.0 if hit_of.get(g) else 0.0 for g in uni])

        # how much of the label vector is expressible in the current basis?
        yl = y[in_lib] - y[in_lib].mean()
        proj = Q[in_lib] @ (Q[in_lib].T @ yl)
        r2s.append(1 - ((yl - proj) ** 2).sum() / (yl ** 2).sum())

        # features
        unmask = (np.array(s.unmasked_relevance_scores) if len(s.unmasked_relevance_scores)
                  else np.array(s.relevance_scores))
        umap = dict(zip(s.genes, unmask))
        feats = {"oracle relevance": np.array([umap.get(g, 0.0) for g in uni])}
        run = run_dir(a.base_dir, a.base_sweep, i, s.dataset_name) / "result.json"
        seeds, wts = [], {}
        if run.is_file():
            steps = json.loads(run.read_text())["steps"][:a.seed_rounds]
            for st in steps:
                for g, h in zip(st["acquired_batch"] or [], st["hits"] or []):
                    if h:
                        seeds.append(g); wts[g] = float(rel_of.get(g, 1.0))
        feats["pathway"] = pathway_feature(mods, seeds, {}, uni) if len(seeds) > 1 else np.zeros(len(uni))
        feats["pathway x strength"] = (pathway_feature(mods, seeds, wts, uni)
                                       if len(seeds) > 1 else np.zeros(len(uni)))
        feats["gene_bias"] = b.copy()
        feats["random"] = rng.standard_normal(len(uni))

        def ef(score):
            order = np.argsort(-score)[:1000]
            picks = [uni[j] for j in order]
            n1 = sum(1 for g in picks if g in lib)
            hits = sum(1 for g in picks if hit_of.get(g))
            return adjusted_ef_value(hits, n1, len(picks) - n1, len(lib), s.total_hits, 1000)

        w, c = fit_logistic(V[in_lib], y[in_lib], a.l2)
        rows["plain"].append(ef(V @ w + b))
        for k, phi in feats.items():
            sd = phi.std() or 1.0
            X = np.column_stack([V, phi / sd])
            w2, _ = fit_logistic(X[in_lib], y[in_lib], a.l2)
            rows[k].append(ef(np.column_stack([V, phi / sd]) @ w2 + b))

    print(f"{a.screen_set}: {len(uni)} genes, u in R^{V.shape[1]}\n")
    print(f"mean R^2 of the hit vector on span(V, b, 1): {np.mean(r2s):.4f} "
          f"(median {np.median(r2s):.4f}, max {np.max(r2s):.4f})")
    print(f"  -> {100*np.mean(r2s):.1f}% of 'which genes are hits' is expressible in the basis\n")
    plain = np.mean(rows["plain"])
    print(f"{'added feature':22s} {'oracle EF':>10s} {'vs plain oracle':>16s}")
    print(f"{'(none)':22s} {plain:10.3f} {0.0:+16.3f}")
    for k in names:
        v = np.mean(rows[k])
        print(f"{k:22s} {v:10.3f} {v - plain:+16.3f}")


if __name__ == "__main__":
    main()
