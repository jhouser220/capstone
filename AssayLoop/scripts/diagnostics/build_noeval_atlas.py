"""Build the ABLATED DepMap atlas -- the project's default since 2026-10-28.

achilles_noeval.npz is achilles.npz with every Achilles cell line used by ANY of the four
evaluation screen sets removed. 40 of 808 lines go (5%), touching 154 screens. The point
is not the score -- the paired effect is a few hundredths, see HANDOFF.md section 4.54 --
it is the CLAIM. With these lines gone, "no screen's own cell line is in the atlas" stops
being something we audited and becomes true by construction, so the leakage question
cannot be reopened against any arm that uses this file.

This file had no builder in the repo while it was already the default that every campaign
reads, which is a reproducibility hole: the default atlas could not be rebuilt or checked.
Running with --verify reproduces the shipped file and asserts it matches byte for byte.

    python scripts/diagnostics/build_noeval_atlas.py --depmap-dir ~/depmap --verify

The resolver is the one section 4.35 used, imported from screen_depmap.py rather than
restated, so the two cannot drift: Model.csv's StrippedCellLineName and CellLineName are
normalised to alphanumerics-uppercase, a small ALIAS table covers engineered subclones,
and a screen's cell_line is tried whole, then cut at "(", then cut at "/".
"""
from __future__ import annotations

import argparse, sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "src")
sys.path.insert(0, "scripts/diagnostics")
from assayloop.tasks import load_screens                  # noqa: E402
from screen_depmap import line_map, resolve               # noqa: E402

SETS = ("paper_validation", "paper_test", "validation", "test")

# Lines that `resolve()` CANNOT reach from any evaluation screen's `cell_line`, found by
# scripts/diagnostics/audit_atlas_lines.py and recorded in HANDOFF.md section 21.
#
# DepMap carries four single-cell subclones of hTERT RPE-1 -- RPE1-ss6/ss48/ss51/ss77 --
# and 61 evaluation screens (6 validation, 55 test) use hTERT-RPE1. norm("hTERT-RPE1") is
# "HTERTRPE1" and norm("RPE1-ss48") is "RPE1SS48": no shared form, and the mapping is
# one-to-many, which resolve() cannot express since it returns a single ModelID. Without
# these the "no screen's own cell line is in the atlas" claim is false for 11% of the
# folds. Add to this set, never to ALIAS, when the audit finds a one-to-many case.
EXTRA_DROP = {
    "ACH-002462",   # RPE1-ss48
    "ACH-002463",   # RPE1-ss77
    "ACH-002464",   # RPE1-ss6
    "ACH-002467",   # RPE1-ss51
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--depmap-dir", default="/home/jeff/depmap")
    ap.add_argument("--out", default=None,
                    help="default: <depmap-dir>/achilles_noeval2.npz, or "
                         "achilles_noeval.npz with --legacy")
    ap.add_argument("--legacy", action="store_true",
                    help="skip EXTRA_DROP and reproduce the original (superseded) "
                         "achilles_noeval.npz, for --verify of the historical file")
    ap.add_argument("--verify", action="store_true",
                    help="compare against the file already on disk and fail on any difference")
    a = ap.parse_args()
    D = Path(a.depmap_dir)
    out = Path(a.out) if a.out else D / (
        "achilles_noeval.npz" if a.legacy else "achilles_noeval2.npz")

    z = np.load(D / "achilles.npz", allow_pickle=True)
    values, rows, cols = z["values"], list(z["rows"]), list(z["cols"])
    print(f"full atlas {values.shape[0]} lines x {values.shape[1]} genes")

    by = line_map(D / "Model.csv")
    drop, hit_screens = set(), 0
    present = set(rows)
    for sset in SETS:
        n_set = 0
        for s in load_screens(target_set=sset):
            mid = resolve(s.cell_line, by)
            if mid in present:
                drop.add(mid); n_set += 1
        hit_screens += n_set
        print(f"  {sset:17s} {n_set:4d} screens resolve to a line in the atlas")
    print(f"name-resolved: {len(drop)} lines, touching {hit_screens} screen slots")
    if not a.legacy:
        extra = {m for m in EXTRA_DROP if m in present} - drop
        drop |= extra
        print(f"EXTRA_DROP (section 21, unreachable by resolve): +{len(extra)} lines {sorted(extra)}")
    print(f"dropping {len(drop)} of {len(rows)} lines in total")

    keep = [i for i, r in enumerate(rows) if r not in drop]
    V = np.ascontiguousarray(values[keep])
    R = np.array([rows[i] for i in keep], dtype=object)
    C = np.array(cols, dtype=object)

    # A gene whose profile is now constant has no direction left to contribute. None are
    # expected -- 18119 genes survived the original ablation -- so this is a tripwire.
    # nanstd, not std: achilles.npz carries 28922 NaNs (0.198%, at most 30 of 808 lines
    # for any one gene) and the acquisition nan_to_num's them at load. With plain std
    # every gene holding a single NaN reads as constant and the tripwire fires on 1936.
    live = int((np.nanstd(V, axis=0) > 0).sum())
    print(f"ablated atlas {V.shape[0]} lines x {V.shape[1]} genes, "
          f"{live} with a non-degenerate profile")
    if live != V.shape[1]:
        print(f"!! {V.shape[1] - live} genes went constant")

    if a.verify:
        if not out.is_file():
            sys.exit(f"--verify: {out} does not exist")
        old = np.load(out, allow_pickle=True)
        # equal_nan: the matrix carries NaNs, and NaN != NaN makes a plain array_equal
        # report DIFFERS on a file that is in fact identical.
        ok = (list(old["rows"]) == list(R) and list(old["cols"]) == list(C)
              and np.array_equal(old["values"], V, equal_nan=True))
        print(f"verify against {out}: {'IDENTICAL' if ok else 'DIFFERS'}")
        if not ok:
            print(f"  rows {len(old['rows'])} vs {len(R)}, "
                  f"shape {old['values'].shape} vs {V.shape}")
            sys.exit(1)
        return

    np.savez(out, values=V, rows=R, cols=C)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
