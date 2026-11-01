"""Screen DepMap features for AssayFormer's additive slot.

Asks the only question that matters before building anything: appended to
a(g) = u.v_g + b_g as an extra term, how much could this feature be worth if its
weight were fitted to the full answer key? Read every number against the two
reference points the script prints -- Gaussian noise (floor) and the screen's own
relevance scores (oracle ceiling).

DepMap inputs (public, free; download with the figshare file-download API):
  achilles.csv   gene effect, rows = DepMap_ID, cols = "SYMBOL (entrez)"
  Model.csv      ModelID -> cell line name / lineage
  Expression.csv log2(TPM+1), rows = DepMap_ID, cols = "SYMBOL (entrez)"   [optional]

Essentiality is signed so that MORE NEGATIVE = more essential; the fitted weight
is free to choose the sign, so no orientation decision is baked in here.
"""
from __future__ import annotations

import argparse, re, sys
from pathlib import Path

import numpy as np
import pandas as pd
import torch

sys.path.insert(0, "src")
sys.path.insert(0, "scripts/diagnostics")
from assayloop.metrics.hits_auc import adjusted_ef_value            # noqa: E402
from assayloop.models.amortized_ranker import AmortizedRankerModel  # noqa: E402
from assayloop.tasks import gene_universe, load_screens             # noqa: E402
from augment_u import fit_logistic                                  # noqa: E402

CKPT = ("/home/jeff/.cache/huggingface/hub/models--Genentech--assayformer/"
        "snapshots/4927703e5903984ae0254a9fe8e93ec0b29c27e4")

# Engineered derivatives and subclones that DepMap lists under the parent line.
ALIAS = {
    "HUH751": "HUH7",            # Huh-7.5.1 is a Huh-7 subclone
    "HELARELAMNEONGREEN": "HELA",
    "CALU1ACE2": "CALU1",
    "HEK293TACE2": "HEK293T",
}


def norm(x):
    return re.sub(r"[^A-Z0-9]", "", str(x).upper())


def load_matrix(path: Path):
    """DepMap wide matrix -> DataFrame indexed by DepMap_ID, columns = symbol.

    Prefers a sibling .npz built by depmap_to_npz.py: parsing an 18k-column CSV
    takes minutes and dominates the runtime, the fits take seconds.
    """
    npz = path.with_suffix(".npz")
    if npz.is_file():
        z = np.load(npz, allow_pickle=True)
        return pd.DataFrame(z["values"], index=list(z["rows"]), columns=list(z["cols"]))
    df = pd.read_csv(path, index_col=0, low_memory=False)
    df.columns = [c.split(" (")[0] for c in df.columns]
    df = df.loc[:, ~df.columns.duplicated()]
    return df


def line_map(model_csv: Path):
    mod = pd.read_csv(model_csv, low_memory=False)
    by = {}
    for _, r in mod.iterrows():
        for col in ("StrippedCellLineName", "CellLineName"):
            k = norm(r[col])
            if k and k not in by:
                by[k] = r["ModelID"]
    return by


def resolve(cell_line, by):
    for cand in (cell_line, str(cell_line).split("(")[0], str(cell_line).split("/")[0]):
        k = norm(cand)
        k = ALIAS.get(k, k)
        if k in by:
            return by[k]
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--depmap-dir", required=True)
    ap.add_argument("--screen-set", default="paper_validation")
    ap.add_argument("--l2", type=float, default=1.0)
    ap.add_argument("--no-expression", action="store_true")
    a = ap.parse_args()
    D = Path(a.depmap_dir)

    m = AmortizedRankerModel(checkpoint=CKPT, ckpt_file="model_last.pt", device="cpu")
    screens = load_screens(target_set=a.screen_set)
    uni = sorted(gene_universe(screens, min_screen_freq=2))
    ids = torch.tensor([m.vocab.to_idx(g) for g in uni])
    V = m.net.gene_emb.weight[ids].detach().numpy().astype(np.float64)
    b = m.net.gene_bias[ids].detach().numpy().astype(np.float64)

    by = line_map(D / "Model.csv")
    ess = load_matrix(D / "achilles.csv")
    expr = None
    if not a.no_expression and (D / "Expression.csv").is_file():
        expr = load_matrix(D / "Expression.csv")
    print(f"essentiality {ess.shape[0]} lines x {ess.shape[1]} genes"
          + (f";  expression {expr.shape[0]} x {expr.shape[1]}" if expr is not None else ""))

    # gene-only summaries (no cell-line matching needed)
    ess_mean = ess.mean(axis=0)
    ess_frac = (ess < -0.5).mean(axis=0)
    col = lambda sr: np.array([sr.get(g, np.nan) for g in uni], dtype=float)

    def fill(v):
        v = v.copy()
        v[~np.isfinite(v)] = np.nanmean(v) if np.isfinite(v).any() else 0.0
        return v

    F_mean, F_frac = fill(col(ess_mean)), fill(col(ess_frac))
    rng = np.random.default_rng(0)
    names = ["ess: pan mean", "ess: frac essential", "ess: THIS cell line",
             "expr: THIS cell line", "expr: not expressed flag",
             "oracle relevance", "random noise"]
    # keyed by screen, so a feature that covers only some screens is compared
    # against the baseline on exactly those screens and not on all 20
    acc = {k: {} for k in ["plain", *names]}
    matched = 0

    for s in screens:
        hit_of = dict(zip(s.genes, s.hits))
        lib = set(s.genes)
        in_lib = np.array([g in lib for g in uni])
        y = np.array([1.0 if hit_of.get(g) else 0.0 for g in uni])
        mid = resolve(s.cell_line, by)
        has_ess = mid in ess.index
        matched += bool(has_ess)

        feats = {
            "ess: pan mean": F_mean,
            "ess: frac essential": F_frac,
            "ess: THIS cell line": fill(col(ess.loc[mid])) if has_ess else F_mean,
            "oracle relevance": np.array([
                dict(zip(s.genes, s.unmasked_relevance_scores or s.relevance_scores)).get(g, 0.0)
                for g in uni]),
            "random noise": rng.standard_normal(len(uni)),
        }
        if expr is not None and mid in expr.index:
            e = fill(col(expr.loc[mid]))
            feats["expr: THIS cell line"] = e
            feats["expr: not expressed flag"] = (e < 1.0).astype(float)

        def ef(score):
            order = np.argsort(-score)[:1000]
            picks = [uni[j] for j in order]
            n1 = sum(1 for g in picks if g in lib)
            hits = sum(1 for g in picks if hit_of.get(g))
            return adjusted_ef_value(hits, n1, len(picks) - n1, len(lib), s.total_hits, 1000)

        w, _ = fit_logistic(V[in_lib], y[in_lib], a.l2)
        acc["plain"][s.dataset_name] = ef(V @ w + b)
        for k, phi in feats.items():
            sd = phi.std() or 1.0
            X = np.column_stack([V, phi / sd])
            w2, _ = fit_logistic(X[in_lib], y[in_lib], a.l2)
            acc[k][s.dataset_name] = ef(X @ w2 + b)

    print(f"cell line matched in DepMap for {matched}/{len(screens)} screens\n")
    print(f"{'added feature':26s} {'oracle EF':>10s} {'vs plain':>10s} {'screens':>8s}")
    print(f"{'(none)':26s} {np.mean(list(acc['plain'].values())):10.3f} "
          f"{0.0:+10.3f} {len(acc['plain']):8d}")
    for k in names:
        if not acc[k]:
            print(f"{k:26s} {'-- no data --':>21s}")
            continue
        shared = sorted(acc[k])                      # the screens this feature covers
        v = np.mean([acc[k][n] for n in shared])
        base = np.mean([acc["plain"][n] for n in shared])
        print(f"{k:26s} {v:10.3f} {v - base:+10.3f} {len(shared):8d}"
              + ("" if len(shared) == len(acc["plain"]) else "   (vs baseline on the same screens)"))


if __name__ == "__main__":
    main()
