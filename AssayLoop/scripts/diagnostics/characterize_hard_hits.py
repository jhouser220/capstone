"""What KIND of gene is a hit on the screens nothing can touch?

The question this answers. §4.23 found 21 of 198 held-out screens where AssayFormer
returns ZERO hits in all 1000 picks, so no seed-based method can ever act on them,
and that bounds the whole family of methods in this project. But it does not say WHY
those screens are unreachable, and the two candidate explanations imply completely
different programmes of work:

  DATA COVERAGE -- their hits are genes that barely appear as hits anywhere in the
      1349 training screens, so the training signal for them does not exist. The fix
      is more screens (BioGRID ORCS), not a better model.

  REPRESENTATION -- their hits ARE well represented in training and still missed. The
      fix is the gene representation itself: the rank-10 BPMF embedding, i.e.
      retraining, i.e. a different project.

Method. Group the fold's screens by how AssayFormer does on them, pool each group's
hit genes, and compare those genes on two axes measured OUTSIDE the screen:

  * training-label coverage -- in how many of the 1349 training screens is this gene
    a hit? Zero means the training matrix says nothing about it.
  * DepMap profile amplitude -- the spread of its knockout effect across 808 cell
    lines. A near-flat profile carries no co-dependency signal even when present, so
    this is the ceiling on what the working method could ever extract for that gene.

The comparison that matters is BETWEEN GROUPS OF SCREENS, not against random genes:
any gene that is a hit somewhere is already unusual, so a random-gene baseline would
confound "is a hit" with "is a hit on a hard screen". A random-gene row is printed
anyway, as a scale reference.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.sparse import load_npz

sys.path.insert(0, "src")
sys.path.insert(0, "scripts/diagnostics")
from assayloop.scripts.full_genome_table import _adj_ef_from_run   # noqa: E402
from assayloop.tasks import gene_universe, load_screens            # noqa: E402
from run_paths import run_dir  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--screen-set", default="validation")
    ap.add_argument("--exclude-set", default="paper_validation")
    ap.add_argument("--dir", default="output/runs")
    ap.add_argument("--base", default="sweep-dm-validation-w0",
                    help="AssayFormer-alone arm; defines which screens yield no hits")
    ap.add_argument("--train-hits", default="output/analysis/train_hits")
    ap.add_argument("--profile-npz", default="/home/jeff/depmap/achilles.npz")
    a = ap.parse_args()

    screens = load_screens(target_set=a.screen_set)
    uni = set(gene_universe(screens, min_screen_freq=2))
    drop = {x.dataset_name for x in load_screens(target_set=a.exclude_set)} \
        if a.exclude_set else set()

    # training-label coverage: how many of the 1349 training screens call each gene a hit
    K = load_npz(Path(a.train_hits).with_suffix(".mat.npz")).tocsc()
    meta = np.load(Path(a.train_hits).with_suffix(".npz"), allow_pickle=True)
    kby = {str(g): j for j, g in enumerate(meta["gene_names"])}
    col_freq = np.asarray(K.sum(axis=0)).ravel()      # per gene, over 1349 screens

    # DepMap amplitude: spread of the knockout effect across cell lines
    z = np.load(Path(a.profile_npz), allow_pickle=True)
    vals = np.nan_to_num(np.asarray(z["values"], dtype=np.float32),
                         posinf=0.0, neginf=0.0)       # (lines, genes)
    dby = {str(c): i for i, c in enumerate(z["cols"])}
    amp = vals.std(axis=0)                             # per gene

    def stats(genes: list[str]) -> dict:
        if not genes:
            return {}
        freq = np.array([col_freq[kby[g]] if g in kby else 0.0 for g in genes])
        inatlas = np.array([g in dby for g in genes])
        amps = np.array([amp[dby[g]] for g in genes if g in dby])
        return dict(n=len(genes),
                    never=100.0 * float((freq == 0).mean()),
                    med_freq=float(np.median(freq)),
                    mean_freq=float(freq.mean()),
                    atlas=100.0 * float(inatlas.mean()),
                    med_amp=float(np.median(amps)) if len(amps) else float("nan"))

    groups: dict[str, list[str]] = {}
    per_screen: list[tuple] = []
    for i, s in enumerate(screens):
        if s.dataset_name in drop:
            continue
        p = run_dir(a.dir, a.base, i, s.dataset_name) / "result.json"
        if not p.is_file():
            continue
        r = json.loads(p.read_text())
        found = sum(int(h) for st in r["steps"] for h in (st["hits"] or []))
        lib = set(s.genes)
        ef = _adj_ef_from_run(p, lib, uni, s.total_hits, len(lib), 1000)
        hits = [g for g, h in zip(s.genes, s.hits) if h]
        if found == 0:
            key = "A. AssayFormer finds ZERO hits"
        elif ef < 2.0:
            key = "B. very hard (EF < 2, some hits)"
        elif ef < 5.0:
            key = "C. middling (2 <= EF < 5)"
        else:
            key = "D. AssayFormer does well (EF >= 5)"
        groups.setdefault(key, []).extend(hits)
        per_screen.append((key, s.dataset_name, found, ef, s.total_hits, len(lib)))

    print(f"{len(per_screen)} fresh screens from {a.screen_set}\n")
    hdr = (f"{'group':34s} {'screens':>7s} {'hits':>7s} {'% hits NEVER a hit in':>22s} "
           f"{'median train':>13s} {'% in':>6s} {'median DepMap':>14s}")
    print(hdr)
    print(f"{'':34s} {'':>7s} {'':>7s} {'any training screen':>22s} "
          f"{'screens hit':>13s} {'atlas':>6s} {'amplitude':>14s}")
    print("-" * len(hdr))
    for key in sorted(groups):
        st = stats(groups[key])
        ns = sum(1 for x in per_screen if x[0] == key)
        print(f"{key:34s} {ns:7d} {st['n']:7d} {st['never']:21.1f}% "
              f"{st['med_freq']:13.0f} {st['atlas']:5.1f}% {st['med_amp']:14.3f}")

    rng = np.random.default_rng(0)
    allg = sorted(uni)
    ref = stats([allg[i] for i in rng.choice(len(allg), 4000, replace=False)])
    print(f"{'(reference) 4000 random genes':34s} {'-':>7s} {ref['n']:7d} "
          f"{ref['never']:21.1f}% {ref['med_freq']:13.0f} {ref['atlas']:5.1f}% "
          f"{ref['med_amp']:14.3f}")

    print("\nthe zero-hit screens themselves:")
    print(f"{'screen':24s} {'total hits':>10s} {'library':>8s} {'hit rate':>9s}")
    for key, nm, found, ef, tot, nlib in sorted(per_screen, key=lambda t: t[4])[:40]:
        if key.startswith("A."):
            print(f"{nm:24s} {tot:10d} {nlib:8d} {100.0 * tot / nlib:8.3f}%")


if __name__ == "__main__":
    main()
