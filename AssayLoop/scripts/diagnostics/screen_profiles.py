"""Screen PROFILE features -- blocks of k columns, not single numbers.

Why this exists. The earlier feature screens (screen_depmap.py, augment_u.py) all
collapsed a rich object into one scalar: mean essentiality over 808 cell lines,
expression in one line, one pathway-propagation score. Each of those is a
hand-chosen projection that answers roughly "how likely is this gene to be a
hit", which is exactly what û·v_g already computes -- so finding them redundant
says little. This appends the PROFILE instead and lets the fitted weights pick
which directions matter:

    a(g) = û·v_g + b_g + sum_j w_j * z_j(g)        z = k PCA components

CRITICAL: a k-column block raises any oracle fit through capacity alone, so every
block is reported beside a WIDTH-MATCHED Gaussian noise block. The honest read is
block minus its own noise control, never block minus the plain baseline.

Blocks:
  ess-profile    each gene's knockout effect across 808 DepMap cell lines, PCA'd.
                 Genes that co-vary across lines are functionally coupled; this is
                 how complexes are normally recovered from DepMap.
  expr-profile   each gene's expression across 1673 lines, PCA'd (co-expression).
  pathway-basis  the gene x pathway membership matrix from the local GMTs, reduced
                 by truncated SVD. Unlike pathway propagation this fixes no
                 aggregation: the fit chooses which pathway directions matter.
  codep-seed     SCREEN-CONDITIONAL, 2 columns: correlation of each gene's
                 essentiality / expression profile with the mean profile of the
                 hits this screen revealed in rounds 1-3. Pathway propagation with
                 a data-derived graph instead of a curated one.
"""
from __future__ import annotations

import argparse, json, sys
from pathlib import Path

import numpy as np
import torch

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
GMTS = ["ReactomePathways.gmt", "c5.go.bp.v2023.2.Hs.symbols.gmt",
        "c2.cp.v2023.2.Hs.symbols.gmt", "h.all.v2023.2.Hs.symbols.gmt"]


def load_npz(path: Path):
    z = np.load(path, allow_pickle=True)
    return z["values"], list(z["rows"]), list(z["cols"])


def profile_block(path: Path, uni, k: int):
    """Genes x cell-lines profile -> (k PCA components, n_present, raw matrix, present mask).

    Genes absent from the atlas are mean-imputed so they sit at the centroid
    rather than at zero, which for a signed quantity like gene effect would be a
    meaningful value. That imputation creates a confound: "is this gene in DepMap
    at all" tracks how well studied it is, which tracks being a hit. The
    `depmap-presence` control block exists to measure that separately.
    """
    from sklearn.decomposition import PCA
    vals, _rows, cols = load_npz(path)
    byname = {g: j for j, g in enumerate(cols)}
    idx = np.array([byname.get(g, -1) for g in uni])
    M = np.zeros((len(uni), vals.shape[0]), dtype=np.float32)
    ok = idx >= 0
    M[ok] = vals[:, idx[ok]].T
    M = np.nan_to_num(M)
    # mean-impute missing genes so they sit at the origin rather than at zero,
    # which for a signed quantity like gene effect is a meaningful value
    M[~ok] = M[ok].mean(axis=0)
    Z = PCA(n_components=k, random_state=0).fit_transform(M)
    return Z, int(ok.sum()), M, ok


def pathway_block(gene_sets_dir: Path, uni, k: int):
    from scipy.sparse import csr_matrix
    from sklearn.decomposition import TruncatedSVD
    pos = {g: i for i, g in enumerate(uni)}
    rows, cols, n_sets = [], [], 0
    for fn in GMTS:
        p = gene_sets_dir / fn
        if not p.is_file():
            continue
        for line in p.read_text().splitlines():
            parts = line.rstrip("\n").split("\t")
            members = [pos[g] for g in parts[2:] if g in pos]
            if 2 <= len(members) <= 500:
                rows.extend(members)
                cols.extend([n_sets] * len(members))
                n_sets += 1
    M = csr_matrix((np.ones(len(rows), dtype=np.float32), (rows, cols)),
                   shape=(len(uni), n_sets))
    Z = TruncatedSVD(n_components=k, random_state=0).fit_transform(M)
    return Z, n_sets


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--depmap-dir", required=True)
    ap.add_argument("--screen-set", default="paper_validation")
    ap.add_argument("--base-sweep", default="sweep-fd2017be")
    ap.add_argument("--base-dir", default="output/runs/baseline/gemini")
    ap.add_argument("--k", type=int, default=50, help="components per block")
    ap.add_argument("--seed-rounds", type=int, default=3)
    ap.add_argument("--l2", type=float, default=1.0)
    ap.add_argument("--per-screen", action="store_true",
                    help="print the net margin per screen (block minus its noise control)")
    a = ap.parse_args()
    D = Path(a.depmap_dir)

    m = AmortizedRankerModel(checkpoint=CKPT, ckpt_file="model_last.pt", device="cpu")
    screens = load_screens(target_set=a.screen_set)
    uni = sorted(gene_universe(screens, min_screen_freq=2))
    ids = torch.tensor([m.vocab.to_idx(g) for g in uni])
    V = m.net.gene_emb.weight[ids].detach().numpy().astype(np.float64)
    b = m.net.gene_bias[ids].detach().numpy().astype(np.float64)
    pos = {g: i for i, g in enumerate(uni)}

    blocks: dict[str, np.ndarray] = {}
    ess_M = expr_M = None
    if (D / "achilles.npz").is_file():
        Z, n, ess_M, ess_ok = profile_block(D / "achilles.npz", uni, a.k)
        blocks["ess-profile"] = Z
        # CONFOUND CONTROL: a single indicator for "this gene is in the atlas".
        # If this alone captures the ess-profile margin, the profile is a proxy
        # for gene coverage rather than for biology.
        blocks["depmap-presence"] = ess_ok.astype(float)[:, None]
        print(f"ess-profile    {Z.shape}  ({n} of {len(uni)} genes present)")
    if (D / "Expression.npz").is_file():
        Z, n, expr_M, _ = profile_block(D / "Expression.npz", uni, a.k)
        blocks["expr-profile"] = Z
        print(f"expr-profile   {Z.shape}  ({n} of {len(uni)} genes present)")
    gs = Path("src/assayloop/data/gene_sets")
    Z, n_sets = pathway_block(gs, uni, a.k)
    blocks["pathway-basis"] = Z
    print(f"pathway-basis  {Z.shape}  ({n_sets} gene sets)")

    def zs(X):
        X = np.asarray(X, dtype=np.float64)
        sd = X.std(axis=0, keepdims=True)
        return (X - X.mean(axis=0, keepdims=True)) / np.where(sd > 0, sd, 1.0)

    blocks = {k: zs(v) for k, v in blocks.items()}
    rng = np.random.default_rng(0)
    acc: dict[str, list[float]] = {}

    for i, s in enumerate(screens):
        hit_of = dict(zip(s.genes, s.hits))
        lib = set(s.genes)
        in_lib = np.array([g in lib for g in uni])
        y = np.array([1.0 if hit_of.get(g) else 0.0 for g in uni])

        def ef(score):
            order = np.argsort(-score)[:1000]
            picks = [uni[j] for j in order]
            n1 = sum(1 for g in picks if g in lib)
            hits = sum(1 for g in picks if hit_of.get(g))
            return adjusted_ef_value(hits, n1, len(picks) - n1, len(lib), s.total_hits, 1000)

        w0, _ = fit_logistic(V[in_lib], y[in_lib], a.l2)
        acc.setdefault("plain", []).append(ef(V @ w0 + b))

        here = dict(blocks)
        # screen-conditional: correlation with the profile of THIS screen's seeds
        run = run_dir(a.base_dir, a.base_sweep, i, s.dataset_name) / "result.json"
        seeds: list[int] = []
        if run.is_file():
            steps = json.loads(run.read_text())["steps"][:a.seed_rounds]
            seeds = [pos[g] for st in steps
                     for g, h in zip(st["acquired_batch"] or [], st["hits"] or [])
                     if h and g in pos]
        cols = []
        for M in (ess_M, expr_M):
            if M is None or len(seeds) < 2:
                continue
            P = M - M.mean(axis=1, keepdims=True)
            nrm = np.linalg.norm(P, axis=1, keepdims=True)
            P = P / np.where(nrm > 0, nrm, 1.0)
            q = P[seeds].mean(axis=0)
            q = q / (np.linalg.norm(q) or 1.0)
            cols.append(P @ q)
        if cols:
            here["codep-seed"] = zs(np.column_stack(cols))

        def deployable(X):
            """Fit the weights on what rounds 1..seed_rounds revealed, then rank
            the rest. This is the only number a real campaign could achieve."""
            if len(seeds) < 20:
                return float("nan")
            idx = np.array(seeds + [pos[g] for st in
                                    json.loads(run.read_text())["steps"][:a.seed_rounds]
                                    for g, h in zip(st["acquired_batch"] or [],
                                                    st["hits"] or []) if not h and g in pos])
            w, _ = fit_logistic(X[idx], y[idx], a.l2)
            sc = X @ w + b
            sc[idx] = -np.inf
            order = list(idx) + list(np.argsort(-sc))
            picks = [uni[j] for j in order[:1000]]
            n1 = sum(1 for g in picks if g in lib)
            hits = sum(1 for g in picks if hit_of.get(g))
            return adjusted_ef_value(hits, n1, len(picks) - n1, len(lib),
                                     s.total_hits, 1000)

        for name, Zb in here.items():
            for tag, Xb in ((name, Zb), (f"{name} [noise]",
                                         rng.standard_normal(Zb.shape))):
                X = np.column_stack([V, Xb])
                w, _ = fit_logistic(X[in_lib], y[in_lib], a.l2)
                acc.setdefault(tag, []).append(ef(X @ w + b))
                acc.setdefault(f"{tag} [deploy]", []).append(deployable(X))

    plain = float(np.mean(acc["plain"]))
    print(f"\n{a.screen_set}: plain oracle {plain:.3f} over {len(acc['plain'])} screens\n")
    fam = families(screens)
    sname = [x.dataset_name for x in screens]
    hdr = (f"{'block':22s} {'cols':>5s} {'NET of noise':>13s} {'95% CI':>18s} "
           f"{'wins':>6s} {'deploy NET':>11s} {'95% CI':>18s}")
    print(hdr); print("-" * len(hdr))
    names = [k for k in acc if k not in ("plain",) and not k.endswith("[noise]")
             and not k.endswith("[deploy]")]
    for name in names:
        cols = (blocks.get(name) if name in blocks else np.zeros((1, 2))).shape[1]
        net = {sname[j]: acc[name][j] - acc[f"{name} [noise]"][j]
               for j in range(len(acc[name]))}
        m, (lo, hi), _ = boot_ci(net, fam)
        dnet = {sname[j]: acc[f"{name} [deploy]"][j] - acc[f"{name} [noise] [deploy]"][j]
                for j in range(len(acc[f"{name} [deploy]"]))
                if acc[f"{name} [deploy]"][j] == acc[f"{name} [deploy]"][j]}
        if dnet:
            dm, (dlo, dhi), _ = boot_ci(dnet, fam)
        else:
            dm = dlo = dhi = float("nan")
        print(f"{name:22s} {cols:5d} {m:+13.3f} [{lo:+6.3f}, {hi:+6.3f}] "
              f"{sum(v > 0 for v in net.values()):3d}/{len(net):<2d} "
              f"{dm:+11.3f} [{dlo:+6.3f}, {dhi:+6.3f}]")
    if a.per_screen:
        print("\nNet margin per screen (block minus its own width-matched noise control):")
        print(f"{'screen':16s} " + " ".join(f"{n[:13]:>14s}" for n in names))
        for j, sc in enumerate(screens):
            cells = []
            for n in names:
                if j < len(acc[n]) and j < len(acc[f"{n} [noise]"]):
                    cells.append(f"{acc[n][j] - acc[f'{n} [noise]'][j]:+14.3f}")
                else:
                    cells.append(f"{'-':>14s}")
            print(f"{sc.dataset_name:16s} " + " ".join(cells))
    print(f"\nplain oracle {plain:.3f}. Both columns are the margin over a WIDTH-MATCHED noise\n"
          "block, bootstrapped over screen families. A k-column block lifts any fit through\n"
          "capacity alone, so the raw 'vs plain' number is not evidence of anything.")


if __name__ == "__main__":
    main()
