"""Does a second atlas recover gene-gene relationships the first one MISSES?

§4.18 established that low inter-column correlation (r ~ 0.11) is necessary but not
sufficient for additivity, and left the question there. §4.41 supplies something
stronger than a correlation: a MECHANISTIC reason for orthogonality. DepMap ranks
co-complex pairs far above paralogue pairs (AUROC 0.750 vs 0.590) and JUMP ORF ranks
them the other way (0.541 vs 0.555), because knockout converges on shared machines
while overexpression converges on shared substrates. Two atlases that disagree about
WHICH relationships they can see may be complementary even if the weaker one is weak.

This measures that directly, in similarity space, before anything is built. For pairs
covered by both atlases it asks: restricted to pairs the FIRST atlas cannot see -- its
cosine no higher than a random pair's -- does the second atlas still separate genuinely related
pairs from random ones? The random baseline is conditioned identically, so the
comparison is not confounded by the conditioning itself.

A conditional AUROC near 0.5 means the second atlas adds nothing where the first is
blind, which is the strongest cheap evidence against building the pair. Clearly above
0.5 means it sees something the first does not, which is a reason to spend a campaign.

    uv run python scripts/diagnostics/atlas_complement.py \
        --sets output/analysis/gate0_sets.json \
        --atlas "DepMap=/home/jeff/depmap/achilles.npz" \
        --atlas "JUMP ORF=output/analysis/jump_orf.npz"
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
from gate0_gof import N_RANDOM, auroc_vs, pair_dots, pairs_of   # noqa: E402
from prescreen_atlas import load                          # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sets", default="output/analysis/gate0_sets.json")
    ap.add_argument("--atlas", action="append", required=True, help="LABEL=PATH, exactly 2")
    ap.add_argument("--screen-set", default="validation")
    ap.add_argument("--center", default="gene")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    if len(a.atlas) != 2:
        raise SystemExit("need exactly two --atlas")

    groups = json.loads(Path(a.sets).read_text())["groups"]
    uni = sorted(gene_universe(load_screens(target_set=a.screen_set), min_screen_freq=2))
    pos = {g: i for i, g in enumerate(uni)}
    rng = np.random.default_rng(a.seed)

    labs, X, OK = [], {}, {}
    for item in a.atlas:
        L, _, path = item.partition("=")
        labs.append(L)
        X[L], OK[L] = load(Path(path), a.center, uni)
    both = OK[labs[0]] & OK[labs[1]]
    print(f"{labs[0]} and {labs[1]}, {a.center}-centred: "
          f"{int(both.sum())} of {len(uni)} genes covered by BOTH "
          f"({100 * both.mean():.1f}%)")

    cov = np.flatnonzero(both)
    r = rng.choice(cov, 2 * N_RANDOM)
    rpr = np.column_stack([r[:N_RANDOM], r[N_RANDOM:]])
    rnd = {L: pair_dots(X[L], rpr) for L in labs}

    print("\nCORRELATION between the two atlases' cosines on the same pairs")
    for gname, sets in groups.items():
        pr, _, _ = pairs_of(sets, pos, both)
        if len(pr) < 10:
            continue
        v = {L: pair_dots(X[L], pr) for L in labs}
        print(f"  {gname:8s} {len(pr):7d} pairs   pearson "
              f"{np.corrcoef(v[labs[0]], v[labs[1]])[0, 1]:+.3f}")
    print(f"  {'random':8s} {N_RANDOM:7d} pairs   pearson "
          f"{np.corrcoef(rnd[labs[0]], rnd[labs[1]])[0, 1]:+.3f}")

    print("\nCONDITIONAL AUROC: can atlas B see pairs atlas A is BLIND to?")
    print("  'blind' = A's cosine for that pair is at or below the median random pair,")
    print("  so A ranks it no higher than two unrelated genes. The random comparison set")
    print("  is filtered the same way, so the conditioning cannot create the effect.")
    for A, B in ((labs[0], labs[1]), (labs[1], labs[0])):
        thr = float(np.median(rnd[A]))
        rmask = rnd[A] <= thr
        rsB = np.sort(rnd[B][rmask])
        print(f"\n  A = {A} (blind threshold cos <= {thr:+.4f}), B = {B}"
              f"   | {rmask.sum()} random pairs qualify")
        print(f"    {'ground truth':10s} {'pairs':>7s} {'blind to A':>11s} "
              f"{'B AUROC all':>12s} {'B AUROC | A blind':>18s}")
        for gname, sets in groups.items():
            pr, _, _ = pairs_of(sets, pos, both)
            if len(pr) < 10:
                continue
            vA, vB = pair_dots(X[A], pr), pair_dots(X[B], pr)
            m = vA <= thr
            if m.sum() < 20:
                continue
            print(f"    {gname:10s} {len(pr):7d} {100 * m.mean():10.1f}% "
                  f"{auroc_vs(vB, np.sort(rnd[B])):12.3f} "
                  f"{auroc_vs(vB[m], rsB):18.3f}")
    print("\n  READ: the last column is the one that matters. Near 0.500 means B adds\n"
          "  nothing where A is blind. Above it, B sees relationships A cannot.")


if __name__ == "__main__":
    main()
