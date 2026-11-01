"""How good can ANY score of the form a(g) = u.v_g + b_g be?

AssayFormer's whole screen-specific channel is a 10-dimensional vector u, which
its transformer predicts from the screen description plus the revealed labels.
Fitting u to a screen's FULL label vector by maximum likelihood gives the best
member of that family -- an oracle, not a policy -- and therefore the ceiling on
every intervention that works by giving AssayFormer a better u. If the ceiling is
close to the released model, the family is exhausted and no such work can pay.

Also reports a deployable mid-point: u fitted to the labels revealed by round k
only, ranking the rest. That one IS a policy.
"""
from __future__ import annotations

import argparse, sys
from pathlib import Path

import numpy as np
import torch

sys.path.insert(0, "src")
from assayloop.metrics.hits_auc import adjusted_ef_value          # noqa: E402
from assayloop.models.amortized_ranker import AmortizedRankerModel  # noqa: E402
from assayloop.tasks import gene_universe, load_screens           # noqa: E402
from run_paths import run_dir  # noqa: E402

CKPT = ("/home/jeff/.cache/huggingface/hub/models--Genentech--assayformer/"
        "snapshots/4927703e5903984ae0254a9fe8e93ec0b29c27e4")


def fit_u_ridge(V, y, l2=1.0):
    """Least-squares u for a CONTINUOUS target (strength-aware oracles)."""
    V = np.asarray(V, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    A = V.T @ V + l2 * np.eye(V.shape[1])
    return np.linalg.solve(A, V.T @ y), 0.0


def build_target(kind, hit_mask, rel):
    """The oracle's regression target.

    binary    1 / 0 -- the label the model actually consumes
    revealed  the relevance_score the task reveals (hits only; non-hits are
              masked to 0), i.e. strictly MORE information than binary
    rank      hits ranked by relevance and mapped to (0, 1], non-hits 0 --
              scale-free, so the compression of relevance near 1.0 cannot hide
              a real ordering
    """
    y = np.zeros(len(hit_mask))
    if kind == "binary":
        y[hit_mask] = 1.0
        return y
    if kind == "revealed":
        y[hit_mask] = rel[hit_mask]
        return y
    if kind == "rank":
        r = rel[hit_mask]
        order = np.argsort(np.argsort(r))
        y[hit_mask] = (order + 1) / len(r)
        return y
    raise ValueError(kind)


def fit_u(V, y, b, l2=1.0, iters=400):
    """Logistic MLE for u in sigmoid(u.v_g + b_g + c), b_g as a fixed offset."""
    V = torch.as_tensor(V, dtype=torch.float32)
    y = torch.as_tensor(y, dtype=torch.float32)
    b = torch.as_tensor(b, dtype=torch.float32)
    u = torch.zeros(V.shape[1], requires_grad=True)
    c = torch.zeros(1, requires_grad=True)
    opt = torch.optim.LBFGS([u, c], max_iter=iters, line_search_fn="strong_wolfe")
    lossfn = torch.nn.BCEWithLogitsLoss()

    def closure():
        opt.zero_grad()
        loss = lossfn(V @ u + b + c, y) + l2 * (u @ u) / len(y)
        loss.backward()
        return loss
    opt.step(closure)
    return u.detach().numpy(), float(c.detach())


def ef_of_ranking(order, lib_set, uni_list, hit_of, total_hits, domain, budget=1000):
    picks = [uni_list[i] for i in order[:budget]]
    n1 = sum(1 for g in picks if g in lib_set)
    n2 = len(picks) - n1
    hits = sum(1 for g in picks if hit_of.get(g))
    return adjusted_ef_value(hits, n1, n2, domain, total_hits, budget)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--screen-set", default="paper_validation")
    ap.add_argument("--reveal-rounds", type=int, default=3,
                    help="for the deployable variant: fit u on the first k*100 picks of a run")
    ap.add_argument("--base-sweep", default="sweep-ctx-base")
    ap.add_argument("--base-dir", default="output/runs")
    ap.add_argument("--l2", type=float, default=1.0)
    ap.add_argument("--target", default="binary", choices=["binary", "revealed", "rank"],
                    help="what the oracle fits u to; see build_target()")
    ap.add_argument("--estimator", default="auto", choices=["auto", "logistic", "ridge"],
                    help="auto = logistic for binary, ridge otherwise. Pass it explicitly to "
                         "compare targets with the estimator held fixed -- otherwise a target "
                         "change and a loss change are confounded.")
    a = ap.parse_args()

    m = AmortizedRankerModel(checkpoint=CKPT, ckpt_file="model_last.pt", device="cpu")
    screens = load_screens(target_set=a.screen_set)
    uni = sorted(gene_universe(screens, min_screen_freq=2))
    ids = torch.tensor([m.vocab.to_idx(g) for g in uni])
    V = m.net.gene_emb.weight[ids].detach().numpy()
    b = m.net.gene_bias[ids].detach().numpy()
    pos = {g: i for i, g in enumerate(uni)}

    print(f"{a.screen_set}: {len(uni)} genes, u in R^{V.shape[1]}\n")
    print(f"target = {a.target}, estimator = {a.estimator}\n")
    print(f"{'screen':16s} {'AF':>7s} {'oracle u':>9s} {'fit on r1-' + str(a.reveal_rounds):>12s}")
    rows = []
    for i, s in enumerate(screens):
        hit_of = dict(zip(s.genes, s.hits))
        rel_of = dict(zip(s.genes, s.relevance_scores))
        lib = set(s.genes)
        y_all = np.array([1.0 if hit_of.get(g) else 0.0 for g in uni])
        in_lib = np.array([g in lib for g in uni])

        # oracle: fit on every labelled gene in the library
        rel = np.array([rel_of.get(g, 0.0) for g in uni])
        tgt = build_target(a.target, in_lib & (y_all > 0), rel)
        est = a.estimator if a.estimator != "auto" else (
            "logistic" if a.target == "binary" else "ridge")
        if est == "logistic":
            u, c = fit_u(V[in_lib], tgt[in_lib], b[in_lib], a.l2)
        else:
            u, c = fit_u_ridge(V[in_lib], tgt[in_lib], a.l2)
        sc = V @ u + b
        oracle = ef_of_ranking(np.argsort(-sc), lib, uni, hit_of, s.total_hits, len(lib))

        # deployable: fit on what a real run had revealed after k rounds
        run = run_dir(a.base_dir, a.base_sweep, i, s.dataset_name) / "result.json"
        policy = float("nan"); af = float("nan")
        if run.is_file():
            import json
            steps = json.loads(run.read_text())["steps"]
            seen = [g for st in steps[:a.reveal_rounds] for g in (st["acquired_batch"] or [])]
            seen = [g for g in seen if g in pos and g in lib]
            if len(seen) > 20:
                idx = np.array([pos[g] for g in seen])
                u2, _ = fit_u(V[idx], y_all[idx], b[idx], a.l2)
                sc2 = V @ u2 + b
                sc2[idx] = -np.inf                       # already tested
                order = list(idx) + list(np.argsort(-sc2))
                policy = ef_of_ranking(order, lib, uni, hit_of, s.total_hits, len(lib))
            af = ef_of_ranking(
                np.array([pos[g] for st in steps for g in (st["acquired_batch"] or [])
                          if g in pos]), lib, uni, hit_of, s.total_hits, len(lib))
        rows.append((af, oracle, policy))
        print(f"{s.dataset_name:16s} {af:7.2f} {oracle:9.2f} {policy:12.2f}")

    A = np.array(rows, dtype=float)
    print(f"\nmean   AssayFormer {np.nanmean(A[:,0]):.3f}   "
          f"oracle u {np.nanmean(A[:,1]):.3f}   fit-on-r1-{a.reveal_rounds} {np.nanmean(A[:,2]):.3f}")
    print(f"headroom inside the 10-d family: {np.nanmean(A[:,1]) - np.nanmean(A[:,0]):+.3f} EF")


if __name__ == "__main__":
    main()
