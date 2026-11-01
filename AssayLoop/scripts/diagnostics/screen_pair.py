"""Screen COMBINATIONS of profile atlases and screen-kNN as added feature blocks.

The open question from §4.17. DepMap co-essentiality (fitness across 808 cell lines)
is worth +0.175 EF deployed. Perturb-seq (transcriptional response in one cell line)
screens real but 2-4x weaker and does not clear the bar alone. They measure different
axes, so the pair may ADD -- which is exactly what happened with DepMap co-essentiality
and screen-kNN (+0.106 + 0.033 offline, +0.145 measured deployed, §4.14), the one time
in this project an offline additivity prediction survived deployment.

Every block is read against a WIDTH-MATCHED Gaussian control averaged over
`--noise-draws` draws. This is not optional here: the whole question is whether a
second column adds signal, and a second column adds fitting capacity whatever is in
it. §4.17's own pre-screen fell into exactly this trap from the other direction, and
§4.15 found that multi-direction blocks buy capacity rather than signal once the
control is honest. The `NET of noise` column is therefore the only interpretable one,
and the paired-contrast table below it is the decisive output: it asks whether adding
an atlas beats adding noise of the same width, per screen, bootstrapped over families.

All directions are difference-of-means (§4.15), the deployed choice.

    uv run python scripts/diagnostics/screen_pair.py \
        --atlas depmap=/home/jeff/depmap/achilles.npz \
        --atlas pseq=/home/jeff/perturbseq/k562_gwps.npz \
        --screen-set paper_validation
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
from assayloop.tasks import gene_universe, load_screens             # noqa: E402
from augment_u import fit_logistic                                  # noqa: E402
from screen_seed_aggs import CKPT, unit_profiles                    # noqa: E402
from score_ctx import boot_ci, families                            # noqa: E402
from run_paths import run_dir  # noqa: E402

MAX_LABELS = 600          # matches the acquisition's own cap


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--atlas", action="append", default=[],
                    help="LABEL=PATH, repeatable")
    ap.add_argument("--train-hits", default="output/analysis/train_hits")
    ap.add_argument("--screen-set", default="paper_validation")
    ap.add_argument("--base-sweep", default="sweep-fd2017be")
    ap.add_argument("--base-dir", default="output/runs/baseline/gemini")
    ap.add_argument("--seed-rounds", type=int, default=3)
    ap.add_argument("--knn-k", type=int, default=50)
    ap.add_argument("--l2", type=float, default=1.0)
    ap.add_argument("--noise-draws", type=int, default=5)
    a = ap.parse_args()

    atlases = []
    for spec in a.atlas:
        label, _, path = spec.partition("=")
        if not path:
            raise SystemExit(f"--atlas wants LABEL=PATH, got {spec!r}")
        atlases.append((label, Path(path)))

    m = AmortizedRankerModel(checkpoint=CKPT, ckpt_file="model_last.pt", device="cpu")
    screens = load_screens(target_set=a.screen_set)
    uni = sorted(gene_universe(screens, min_screen_freq=2))
    ids = torch.tensor([m.vocab.to_idx(g) for g in uni])
    V = m.net.gene_emb.weight[ids].detach().numpy().astype(np.float64)
    b = m.net.gene_bias[ids].detach().numpy().astype(np.float64)
    pos = {g: i for i, g in enumerate(uni)}

    P: dict[str, np.ndarray] = {}
    for label, path in atlases:
        X, ok = unit_profiles(path, uni)
        P[label] = X.astype(np.float64)
        print(f"{label:10s} {X.shape[1]:6d} dims, {int(ok.sum()):6d} of {len(uni)} "
              f"universe genes covered")

    K = load_npz(Path(a.train_hits).with_suffix(".mat.npz")).tocsc()
    meta = np.load(Path(a.train_hits).with_suffix(".npz"), allow_pickle=True)
    kby = {str(g): j for j, g in enumerate(meta["gene_names"])}
    Kcols = np.array([kby.get(g, -1) for g in uni])
    k_size = np.sqrt(np.maximum(np.asarray(K.sum(axis=1)).ravel(), 1.0))
    print(f"screen-kNN {K.shape[0]} training screens, "
          f"{int((Kcols >= 0).sum())} of {len(uni)} universe genes present")

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
        if len(seeds) < 2 or len(nonseeds) < 10:
            continue

        cols: dict[str, np.ndarray] = {}
        for label, X in P.items():
            w = (X[seeds[:MAX_LABELS]].mean(axis=0)
                 - X[nonseeds[:MAX_LABELS * 4]].mean(axis=0))
            n = np.linalg.norm(w)
            cols[label] = X @ (w / n) if n > 0 else np.zeros(len(uni))

        kidx = np.array([kby[uni[j]] for j in seeds if uni[j] in kby])
        if len(kidx) >= 2:
            sim = np.asarray(K[:, kidx].sum(axis=1)).ravel() / k_size
            top = np.argsort(-sim)[:a.knn_k]
            w = sim[top]
            if w.sum() > 0:
                prof = np.asarray(K[top].T.dot(w)).ravel() / w.sum()
                v = np.zeros(len(uni))
                okc = Kcols >= 0
                v[okc] = prof[Kcols[okc]]
                cols["knn"] = v

        names = [lbl for lbl, _ in atlases if lbl in cols]
        blocks: dict[str, list[str]] = {k: [k] for k in cols}
        if len(names) >= 2:
            blocks["+".join(names)] = list(names)
        if "knn" in cols:
            for lbl in names:
                blocks[f"{lbl}+knn"] = [lbl, "knn"]
            if len(names) >= 2:
                blocks["+".join(names) + "+knn"] = [*names, "knn"]

        def ef(score):
            order = np.argsort(-score)[:1000]
            picks = [uni[j] for j in order]
            n1 = sum(1 for g in picks if g in lib)
            hits = sum(1 for g in picks if hit_of.get(g))
            return adjusted_ef_value(hits, n1, len(picks) - n1, len(lib), s.total_hits, 1000)

        def _fit(col):
            c = np.asarray(col, dtype=np.float64).reshape(len(uni), -1)
            sd = c.std(axis=0, keepdims=True)
            X = np.column_stack([V, c / np.where(sd > 0, sd, 1.0)])
            w, _ = fit_logistic(X[in_lib], y[in_lib], a.l2)
            return ef(X @ w + b)

        # Mechanism check: two components that correlate are substitutes, and no
        # amount of fusion weight makes a substitute additive. Measured on library
        # genes only, which is where the budget is actually spent.
        for x, z in [(p, q) for p in cols for q in cols if p < q]:
            u, v2 = cols[x][in_lib], cols[z][in_lib]
            if u.std() > 0 and v2.std() > 0:
                acc.setdefault(f"~corr {x}|{z}", {})[s.dataset_name] = float(
                    np.corrcoef(u, v2)[0, 1])

        w0, _ = fit_logistic(V[in_lib], y[in_lib], a.l2)
        acc.setdefault("plain", {})[s.dataset_name] = ef(V @ w0 + b)
        for name, keys in blocks.items():
            phi = np.column_stack([cols[k] for k in keys])
            acc.setdefault(name, {})[s.dataset_name] = _fit(phi)
            nrng = np.random.default_rng(
                abs(hash((s.dataset_name, name, phi.shape[1]))) % 2**32)
            acc.setdefault(f"{name} [noise]", {})[s.dataset_name] = float(np.mean(
                [_fit(nrng.standard_normal(phi.shape)) for _ in range(a.noise_draws)]))

    kept = [s for s in screens if s.dataset_name in acc.get("plain", {})]
    fam = families(kept)
    plain = acc["plain"]
    print(f"\n{a.screen_set}: {len(plain)} screens, {len(set(fam.values()))} families, "
          f"plain oracle {np.mean(list(plain.values())):.3f}\n")

    for name in sorted(k for k in acc if k.startswith("~corr ")):
        v = list(acc[name].values())
        print(f"{name:26s} mean r {np.mean(v):+.3f}  (|r| > ~0.3 means substitutes, "
              f"not complements)" if v else "")
    print()

    nets: dict[str, dict[str, float]] = {}
    hdr = f"{'feature block':26s} {'cols':>4s} {'NET of noise':>13s} {'95% CI':>18s} {'wins':>8s}"
    print(hdr)
    print("-" * len(hdr))
    for name in [k for k in acc if k != "plain" and not k.endswith("[noise]")
                 and not k.startswith("~corr ")]:
        net = {n: acc[name][n] - acc[f"{name} [noise]"][n] for n in acc[name]}
        nets[name] = net
        mm, (lo, hi), _ = boot_ci(net, fam)
        print(f"{name:26s} {name.count('+') + 1:4d} {mm:+13.3f} [{lo:+6.3f}, {hi:+6.3f}] "
              f"{sum(v > 0 for v in net.values()):3d}/{len(net):<3d}")

    # The decisive table. Each row asks: does adding this component beat adding
    # nothing, with both sides already net of their own width-matched noise? A
    # positive block mean means little if the single-atlas block is just as good.
    print("\npaired contrasts (both sides net of their own width-matched noise)")
    pairs = [(b, s_) for b in nets for s_ in nets
             if b != s_ and set(s_.split("+")) < set(b.split("+"))]
    print(f"{'contrast':38s} {'delta':>9s} {'95% CI':>18s} {'wins':>8s}")
    print("-" * 76)
    for big, small in sorted(pairs, key=lambda t: (len(t[0]), t[0])):
        common = set(nets[big]) & set(nets[small])
        d = {n: nets[big][n] - nets[small][n] for n in common}
        if not d:
            continue
        mm, (lo, hi), _ = boot_ci(d, fam)
        print(f"{big + ' vs ' + small:38s} {mm:+9.3f} [{lo:+6.3f}, {hi:+6.3f}] "
              f"{sum(v > 0 for v in d.values()):3d}/{len(d):<3d}")


if __name__ == "__main__":
    main()
