"""Kernel pre-screen: how well does a candidate gene space propagate confirmed hits?

*** KNOWN LIMITATION -- READ BEFORE TRUSTING A NUMBER FROM THIS SCRIPT. ***
This measures STANDALONE ranking ability: the AUROC of similarity-to-confirmed-hits
inside a candidate space. It does NOT measure incremental value over `span(V)`, and
this project has already proved those two can diverge completely -- the LLM gene
score scored AUROC 0.553-0.584 standalone and exactly 0.00 incremental (§4.11).
Including `v_g` as a reference row does not repair this: it compares standalone
ability against standalone ability. Worse, the DepMap row falsifies the obvious
decision rule, because DepMap scores BELOW `v_g` here (0.552 vs 0.634 on validation)
and is nonetheless worth +0.175 EF deployed -- a kernel does not need to beat `v_g`,
it needs to be complementary to it, which is a different quantity that this script
never computes.

Measured against the calibrated instrument the two disagree: Perturb-seq ranks ABOVE
DepMap here on both coverage conventions, and BELOW it by 2-4x on the EF oracle
(§4.17). **The EF oracle is the one to believe** -- `augment_u.py` and
`screen_seed_aggs.py` fit the added column jointly with û and so measure the
increment, and their scale is calibrated to a deployed result.

What this script is still good for: a cheap sanity check that a new space carries
gene-level information at all, read against its own `[shuf]` control. That is a
genuine use -- it is what caught 6832 infinities silently zeroing the Perturb-seq
kernel -- but it is a smoke test, not a screen.

Original intent, for the record:

Why this exists. Four additive features have been screened and all landed at ~0
(pathway, DepMap essentiality, DepMap expression, LLM score -- §4.10, §4.11), and
the one thing that worked did so by a different route: not a column, but a KERNEL.
`codep_seed` scores an untested gene by its similarity to the hits this screen has
already confirmed, measured in DepMap profile space, and that is worth +0.175 EF
held out. So the question for any new gene embedding -- ESM-2, Perturb-seq,
Geneformer, co-expression -- is not "is it predictive" but "is its geometry better
than `v_g`'s for propagating hits".

That question is answerable in seconds per space, with no campaign and no fitting,
which is the point: the four previous nulls each cost a build.

THE FLOOR IS `v_g`, NOT RANDOM. AssayFormer's own gene embedding defines a
similarity too, and the model already has it. A candidate space that beats random
but not `v_g` tells us nothing we can use -- that was precisely the LLM-score
trap, where AUROC 0.553-0.584 against random looked like signal and the
incremental value was 0.00. DepMap is included as the known-good reference: it is
the only space here with a deployed result, so it calibrates what a USEFUL number
on this scale looks like.

Three measurements per space, because "predictive" is not one thing:

* `AUROC all`    -- over every untested library gene. Comparable across spaces.
* `AUROC tail`   -- restricted to genes AssayFormer ranks OUTSIDE its top 1000,
                    i.e. the hits it actually misses. §4.12 located essentially
                    all the task's headroom here, so this is the decision-relevant
                    column and it is where a merely-redundant space collapses.
* `tail top1k`   -- hit rate among the kernel's own top 1000 of that tail, as a
                    practical read in the units the campaign cares about.

Two controls:

* `[shuf]` permutes which gene owns which profile row, keeping each space's
  marginal distribution and destroying only the gene correspondence. Any space
  must beat its own shuffle or its geometry is not carrying gene information.
* every space is scored on the SAME gene set (the intersection of all spaces'
  coverage) unless `--no-intersect`. Coverage differs a lot -- the DepMap atlas
  holds 18.1k genes, a Perturb-seq library ~9.9k -- and a space evaluated on a
  smaller, better-measured subset would otherwise look better for free.

    uv run python scripts/diagnostics/screen_kernels.py \
        --space depmap:/home/jeff/depmap/achilles.npz \
        --space perturbseq:/home/jeff/perturbseq/k562_gwps.npz:dim \
        --screen-set paper_validation
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from run_paths import run_dir  # noqa: E402
import torch

sys.path.insert(0, "src")
sys.path.insert(0, "scripts/diagnostics")
from assayloop.models.amortized_ranker import AmortizedRankerModel  # noqa: E402
from assayloop.tasks import gene_universe, load_screens             # noqa: E402
from screen_seed_aggs import CKPT                                   # noqa: E402
from score_ctx import boot_ci, families                             # noqa: E402


def auroc(score: np.ndarray, y: np.ndarray) -> float:
    """Mann-Whitney AUROC, ties averaged. nan if either class is empty."""
    pos, neg = y > 0, y <= 0
    if not pos.any() or not neg.any():
        return float("nan")
    r = np.empty(len(score), dtype=np.float64)
    order = np.argsort(score)
    s = score[order]
    i = 0
    while i < len(s):
        j = i
        while j + 1 < len(s) and s[j + 1] == s[i]:
            j += 1
        r[order[i:j + 1]] = 0.5 * (i + j) + 1.0
        i = j + 1
    n1, n0 = int(pos.sum()), int(neg.sum())
    return float((r[pos].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


def load_space(path: Path, center: str) -> tuple[np.ndarray, dict[str, int]]:
    """A gene space as unit-normalised rows, so a dot product is a cosine.

    `center="gene"` subtracts each GENE's own mean across dimensions. That is right
    for a profile whose dimensions are comparable measurements of one quantity (a
    knockout effect per cell line), and it is what `codep_seed_acq.py` ships.

    `center="dim"` standardises each DIMENSION across genes instead. That is right
    for a learned embedding, where dimension 7 has no units in common with
    dimension 400 and a per-gene mean is meaningless. Getting this backwards
    produces a null for the wrong reason, so it is explicit per space rather than
    global.
    """
    z = np.load(path, allow_pickle=True)
    vals = np.asarray(z["values"], dtype=np.float32)        # (dims, genes)
    cols = [str(c) for c in z["cols"]]
    # float32 throughout: a Perturb-seq atlas is 8248 dims x ~10k genes, and the
    # aligned copies would be 1.4 GB each in float64 on a 14 GB machine. Only the
    # RANKING matters here, which float32 cosines preserve.
    # posinf/neginf must be 0, not the default 3.4e38: one infinity in a row makes
    # that row's norm astronomical and its cosine with everything ~0, which reads as
    # a clean null. The Replogle atlas has ~6.8k of them.
    X = np.nan_to_num(vals.T, posinf=0.0, neginf=0.0)       # (genes, dims)
    if center == "gene":
        X -= X.mean(axis=1, keepdims=True)
    elif center == "dim":
        X -= X.mean(axis=0, keepdims=True)
        sd = X.std(axis=0, keepdims=True)
        X /= np.where(sd > 0, sd, 1.0)
    else:
        raise ValueError(f"center must be 'gene' or 'dim', got {center!r}")
    n = np.linalg.norm(X, axis=1, keepdims=True)
    return X / np.where(n > 0, n, 1.0), {g: i for i, g in enumerate(cols)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--space", action="append", default=[],
                    help="LABEL:PATH[:gene|dim] -- repeatable. Default centring is 'gene'.")
    ap.add_argument("--screen-set", default="paper_validation")
    ap.add_argument("--base-sweep", default="sweep-fd2017be")
    ap.add_argument("--base-dir", default="output/runs/baseline/gemini")
    ap.add_argument("--seed-rounds", type=int, default=3)
    ap.add_argument("--budget", type=int, default=1000)
    ap.add_argument("--no-intersect", action="store_true",
                    help="score each space on its own coverage instead of the shared set. "
                         "Reports are then NOT comparable between spaces; use for coverage "
                         "diagnosis only.")
    a = ap.parse_args()

    specs = []
    for sp in a.space:
        parts = sp.split(":")
        if len(parts) == 2:
            parts.append("gene")
        if len(parts) != 3:
            raise SystemExit(f"--space wants LABEL:PATH[:gene|dim], got {sp!r}")
        specs.append((parts[0], Path(parts[1]), parts[2]))

    m = AmortizedRankerModel(checkpoint=CKPT, ckpt_file="model_last.pt", device="cpu")
    screens = load_screens(target_set=a.screen_set)
    uni = sorted(gene_universe(screens, min_screen_freq=2))
    pos = {g: i for i, g in enumerate(uni)}

    # the FLOOR: AssayFormer's own gene geometry, unit-normalised the same way
    ids = torch.tensor([m.vocab.to_idx(g) for g in uni])
    V = m.net.gene_emb.weight[ids].detach().numpy().astype(np.float64)
    Vn = V - V.mean(axis=0, keepdims=True)
    Vn /= np.where(np.linalg.norm(Vn, axis=1, keepdims=True) > 0,
                   np.linalg.norm(Vn, axis=1, keepdims=True), 1.0)

    # label -> (aligned rows, covered mask, row order). The order is how the
    # shuffled control is built: permuting an index array rather than copying a
    # 1.4 GB matrix.
    spaces: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    ident = np.arange(len(uni))
    spaces["v_g (FLOOR)"] = (Vn.astype(np.float32), np.ones(len(uni), dtype=bool), ident)
    for label, path, center in specs:
        X, by = load_space(path, center)
        idx = np.array([by.get(g, -1) for g in uni])
        ok = idx >= 0
        A = np.zeros((len(uni), X.shape[1]), dtype=np.float32)
        A[ok] = X[idx[ok]]
        del X
        spaces[label] = (A, ok, ident)
        print(f"{label:16s} {A.shape[1]:6d} dims, centring={center}, "
              f"{int(ok.sum()):6d} of {len(uni)} universe genes covered")

    shared = np.ones(len(uni), dtype=bool)
    if not a.no_intersect:
        for _X, ok, _o in spaces.values():
            shared &= ok
        print(f"\nshared gene set: {int(shared.sum())} genes "
              f"(every space scored on these, so the columns are comparable)")

    # the shuffled control: same rows, permuted gene ownership
    rng = np.random.default_rng(0)
    for label in list(spaces):
        if label.startswith("v_g"):
            continue
        A, ok, _o = spaces[label]
        spaces[f"{label} [shuf]"] = (A, ok, rng.permutation(len(uni)))

    res: dict[str, dict[str, dict[str, float]]] = {}
    for i, s in enumerate(screens):
        hit_of = dict(zip(s.genes, s.hits))
        run = run_dir(a.base_dir, a.base_sweep, i, s.dataset_name) / "result.json"
        if not run.is_file():
            continue
        steps = json.loads(run.read_text())["steps"]
        tested = {g for st in steps[:a.seed_rounds] for g in (st["acquired_batch"] or [])}
        seeds = [pos[g] for st in steps[:a.seed_rounds]
                 for g, h in zip(st["acquired_batch"] or [], st["hits"] or [])
                 if h and g in pos]
        if len(seeds) < 2:
            continue

        # AssayFormer's own ranking of the untested pool, read off the SAME run, so
        # the tail column means "genes this policy really did leave outside its top
        # 1000" rather than a re-scored approximation
        af_rank = {g: r for r, g in enumerate(
            [g for st in steps for g in (st["acquired_batch"] or [])])}

        lib = set(s.genes)
        cand = np.array([g in lib and g not in tested and shared[j]
                         for j, g in enumerate(uni)])
        if cand.sum() < 100:
            continue
        y = np.array([1.0 if hit_of.get(g) else 0.0 for g in uni])[cand]
        if y.sum() < 5:
            continue
        cg = [uni[j] for j in np.flatnonzero(cand)]
        # the tail: outside AssayFormer's own budget, where §4.12 put the headroom
        tail = np.array([af_rank.get(g, 10**9) >= a.budget for g in cg])

        ci = np.flatnonzero(cand)
        for label, (X, _ok, order) in spaces.items():
            sim = (X[order[ci]] @ X[order[seeds]].mean(axis=0)).astype(np.float64)
            d = res.setdefault(label, {})
            d.setdefault("all", {})[s.dataset_name] = auroc(sim, y)
            if tail.sum() > 50 and y[tail].sum() >= 3:
                d.setdefault("tail", {})[s.dataset_name] = auroc(sim[tail], y[tail])
                k = min(a.budget, int(tail.sum()))
                top = np.argsort(-sim[tail])[:k]
                d.setdefault("top1k", {})[s.dataset_name] = float(
                    100.0 * y[tail][top].mean())
                d.setdefault("base", {})[s.dataset_name] = float(100.0 * y[tail].mean())

    kept = [s for s in screens if s.dataset_name in res["v_g (FLOOR)"]["all"]]
    fam = families(kept)
    print(f"\n{a.screen_set}: {len(kept)} screens, {len(set(fam.values()))} families, "
          f"seeds from rounds 1-{a.seed_rounds}")
    base = res["v_g (FLOOR)"].get("base", {})
    if base:
        print(f"tail base hit rate (what random picking scores there): "
              f"{np.mean(list(base.values())):.2f}%\n")

    hdr = (f"{'gene space':24s} {'AUROC all':>19s} {'AUROC tail':>19s} "
           f"{'tail top1k %':>13s}")
    print(hdr)
    print("-" * len(hdr))
    for label, d in res.items():
        cells = []
        for key in ("all", "tail"):
            v = d.get(key, {})
            if not v:
                cells.append(f"{'-':>19s}")
                continue
            mm, (lo, hi), _ = boot_ci(v, fam)
            cells.append(f"{mm:7.3f} [{lo:.3f},{hi:.3f}]")
        t = d.get("top1k", {})
        tt = f"{np.mean(list(t.values())):13.2f}" if t else f"{'-':>13s}"
        print(f"{label:24s} {cells[0]:>19s} {cells[1]:>19s} {tt}")

    print("\nRead the tail column, against v_g (the floor) and depmap (the only space\n"
          "here with a deployed result, +0.175 EF). A space that beats random but not\n"
          "v_g is the LLM-score trap repeating: AUROC 0.58 standalone, 0.00 incremental.")


if __name__ == "__main__":
    main()
