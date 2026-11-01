"""GATE 0, re-run with a positive set appropriate to a GAIN-OF-FUNCTION atlas.

§4.39 failed JUMP ORF because co-complex genes cluster at ~0 sd, then explained the
failure by stoichiometry -- overexpressing one subunit of a complex creates an
IMBALANCE, so subunits should diverge. That mechanism PREDICTS the failure, so the
control cannot distinguish "this atlas has no functional signal" from "co-complex
membership is the wrong ground truth for overexpression data". Everything §4.39
concluded about gain-of-function rests on that one control.

This scores three ground truths against every atlas at once, so the comparison is a
DOUBLE DISSOCIATION rather than a single null (see `build_gate0_sets.py` for how the
groups are built and why these two predictions are opposite):

                      knockout atlas (DepMap)      overexpression atlas (JUMP ORF)
    complexes         converge (shared machine)    diverge (stoichiometric imbalance)
    paralogues        weak (mutual buffering)      converge (same substrates, no imbalance)

If JUMP shows paralogues > complexes while DepMap shows complexes > paralogues, the
overexpression atlas carries real functional signal of the kind it SHOULD carry, and
§4.39's conclusion ("gain-of-function may not support co-functional clustering at all")
does not hold -- which reopens a 63-66%-coverage atlas aimed at the CRISPRa screens
that hold 72-80% of the hard screens' hits (§4.26b, §4.30). If JUMP is flat on ALL
THREE, the null is structural and becomes a finding rather than a dead end.

Statistics. The headline is AUROC of within-set pair cosines against random covered
pairs: scale-free, so it compares across atlases whose cosine distributions differ in
width, which a raw sd gap does not. The sd gap §4.39 quoted is reported alongside for
continuity. Random pairs are drawn only from genes the atlas covers, so coverage cannot
drive it. Both centring modes are always shown, because getting centring backwards
produces a null for the wrong reason -- the §4.17 lesson.

    uv run python scripts/diagnostics/gate0_gof.py \
        --sets output/analysis/gate0_sets.json \
        --atlas "JUMP ORF=output/analysis/jump_orf.npz" \
        --atlas "DepMap=/home/jeff/depmap/achilles.npz"
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "src")
sys.path.insert(0, "scripts/diagnostics")
from assayloop.tasks import gene_universe, load_screens   # noqa: E402
from prescreen_atlas import load                          # noqa: E402  (one centring impl)

N_RANDOM = 200_000
CHUNK = 20_000   # X[i]*X[j] over 287k pathway pairs x 790 dims is 1.8 GB per temporary


def auroc_vs(w: np.ndarray, rnd_sorted: np.ndarray) -> float:
    """P(within-set pair cosine > random pair cosine), ties at half credit. The same
    quantity as `screen_kernels.auroc`, but O(n log m) against a pre-sorted null rather
    than an argsort of the pooled array, which is what makes 400 bootstrap replicates
    affordable on the 151k-pair pathway group."""
    lo = np.searchsorted(rnd_sorted, w, side="left")
    hi = np.searchsorted(rnd_sorted, w, side="right")
    return float(((lo + hi) / 2).mean() / len(rnd_sorted))


def boot_sets(w, fam, rnd_sorted, n, rng):
    """95% CI for the AUROC, resampling FAMILIES rather than pairs -- pairs inside one
    family share genes and are not independent, the same reason §6 bootstraps screen
    families rather than screens."""
    nf = int(fam.max()) + 1
    by = [np.flatnonzero(fam == k) for k in range(nf)]
    v = [auroc_vs(w[np.concatenate([by[k] for k in rng.integers(0, nf, nf)])], rnd_sorted)
         for _ in range(n)]
    return float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))


def pair_dots(X: np.ndarray, pr: np.ndarray) -> np.ndarray:
    """Cosines for index pairs, chunked -- the unchunked form OOMs on the pathway group."""
    out = np.empty(len(pr), dtype=np.float64)
    for s in range(0, len(pr), CHUNK):
        e = min(s + CHUNK, len(pr))
        out[s:e] = np.einsum("ij,ij->i", X[pr[s:e, 0]], X[pr[s:e, 1]])
    return out


def pairs_of(sets: dict[str, list[str]], pos: dict[str, int], ok: np.ndarray):
    """Within-set index pairs, keeping only genes the atlas covers."""
    out, fam, kept = [], [], 0
    for genes in sets.values():
        gi = [pos[g] for g in genes if g in pos and ok[pos[g]]]
        if len(gi) < 2:
            continue
        for x in range(len(gi)):
            for y in range(x + 1, len(gi)):
                out.append((gi[x], gi[y]))
                fam.append(kept)
        kept += 1
    return (np.array(out, dtype=np.int64).reshape(-1, 2),
            np.array(fam, dtype=np.int64), kept)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sets", default="output/analysis/gate0_sets.json")
    ap.add_argument("--atlas", action="append", required=True,
                    help="LABEL=PATH, repeatable")
    ap.add_argument("--screen-set", default="validation")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--boot", type=int, default=400)
    a = ap.parse_args()

    spec = json.loads(Path(a.sets).read_text())
    groups, spotlight = spec["groups"], spec["spotlight"]
    uni = sorted(gene_universe(load_screens(target_set=a.screen_set), min_screen_freq=2))
    pos = {g: i for i, g in enumerate(uni)}
    rng = np.random.default_rng(a.seed)

    rows = []
    for item in a.atlas:
        label, _, path = item.partition("=")
        for center in ("gene", "dim"):
            X, ok = load(Path(path), center, uni)
            cov = np.flatnonzero(ok)
            r = rng.choice(cov, 2 * N_RANDOM)
            rnd = pair_dots(X, np.column_stack([r[:N_RANDOM], r[N_RANDOM:]]))
            print(f"\n=== {label} [{center}-centred] === {len(cov)} of {len(uni)} "
                  f"universe genes covered ({100 * ok.mean():.1f}%), {X.shape[1]} dims")
            print(f"    random covered pair cosine {rnd.mean():+.4f} (sd {rnd.std():.4f})")
            print(f"    {'ground truth':10s} {'sets':>5s} {'pairs':>7s} {'mean cos':>9s} "
                  f"{'sd gap':>8s} {'AUROC':>7s}")
            rs, held = np.sort(rnd), {}
            for gname, sets in groups.items():
                pr, fam, kept = pairs_of(sets, pos, ok)
                if len(pr) < 10:
                    continue
                w = pair_dots(X, pr)
                z = (w.mean() - rnd.mean()) / (rnd.std() or 1)
                au = auroc_vs(w, rs)
                lo, hi = boot_sets(w, fam, rs, a.boot, rng)
                held[gname] = (w, fam)
                rows.append((label, center, gname, kept, len(pr), w.mean(), z, au, lo, hi))
                print(f"    {gname:10s} {kept:5d} {len(pr):7d} {w.mean():+9.4f} "
                      f"{z:+8.2f} {au:7.3f}  [{lo:.3f}, {hi:.3f}]")
            if "paralog" in held and "complex" in held:
                # The decisive statistic: is the atlas BETTER on the ground truth that
                # suits its perturbation direction than on the one that does not?
                wp, fp = held["paralog"]
                wc, fc = held["complex"]
                npf, ncf = int(fp.max()) + 1, int(fc.max()) + 1
                byp = [np.flatnonzero(fp == k) for k in range(npf)]
                byc = [np.flatnonzero(fc == k) for k in range(ncf)]
                d = np.array([
                    auroc_vs(wp[np.concatenate([byp[k] for k in rng.integers(0, npf, npf)])], rs)
                    - auroc_vs(wc[np.concatenate([byc[k] for k in rng.integers(0, ncf, ncf)])], rs)
                    for _ in range(a.boot)])
                print(f"    DISSOCIATION  paralog - complex = "
                      f"{auroc_vs(wp, rs) - auroc_vs(wc, rs):+.3f} AUROC "
                      f"[{np.percentile(d, 2.5):+.3f}, {np.percentile(d, 97.5):+.3f}]  "
                      f"P(>0) = {(d > 0).mean():.2f}")

            spot = {k: groups["paralog"][k] for k in spotlight if k in groups["paralog"]}
            if spot and center == "gene":
                print("    spotlight paralogue families (convergent when overexpressed):")
                for nm, gs in sorted(spot.items()):
                    pr, _, _ = pairs_of({nm: gs}, pos, ok)
                    if len(pr) < 1:
                        continue
                    w = pair_dots(X, pr)
                    print(f"      {nm[:48]:48s} {len(pr):3d} pairs "
                          f"{w.mean():+.4f}  {(w.mean() - rnd.mean()) / (rnd.std() or 1):+6.2f} sd")

    print("\n" + "=" * 78)
    print("DOUBLE DISSOCIATION  (AUROC of within-set vs random pairs; 0.500 = no signal)")
    labs = sorted({r[0] for r in rows}, key=lambda L: [r[0] for r in rows].index(L))
    for center in ("gene", "dim"):
        print(f"\n  [{center}-centred]")
        print(f"    {'ground truth':10s} " + "  ".join(f"{L:>22s}" for L in labs))
        for gname in groups:
            cells = []
            for L in labs:
                m = [r for r in rows if r[0] == L and r[1] == center and r[2] == gname]
                cells.append(f"{m[0][7]:.3f} [{m[0][8]:.3f},{m[0][9]:.3f}]" if m else "-")
            print(f"    {gname:10s} " + "  ".join(f"{c:>22s}" for c in cells))
    print("\n  READ: paralog > complex in the overexpression atlas AND complex > paralog in\n"
          "  the knockout atlas is the dissociation that reopens the atlas. Flat on all three\n"
          "  rows means the null is structural (see this file's docstring).")


if __name__ == "__main__":
    main()
