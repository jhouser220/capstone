"""Build HistoricalBank NPZs for the paper's 20-screen test set from assets in hand.

Rationale (HANDOFF section 24.8 / 24.10). A bank's donor side is frozen and shared:
the 15 bundle banks whose target is not a donor carry byte-identical `residual`,
`availability` and `scalar_c`. The only per-target array is `mean`, which needs the
unavailable `W_text32` -- and zeroing it is numerically immaterial (17 of 18 screens
at dEF exactly 0.0000). So every field is obtainable:

  genes, residual, availability, scalar_c, donor_publications  <- a 134-donor template
  query_publication                                            <- manifest `source_id`
  mean                                                         <- zeros(32)

The donor-exclusion invariant is enforced, not assumed: the build aborts if any target
publication appears among the template's donors.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import yaml

TEMPLATE_DEFAULT = "1566"   # any of the 15 with 134 donors; all byte-identical


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--bundle", required=True, type=Path)
    ap.add_argument("--manifest", required=True, type=Path,
                    help="assaybench assayloop-test.yaml (or assayloop-validation.yaml)")
    ap.add_argument("--template", default=TEMPLATE_DEFAULT,
                    help="bundle bank id to take the donor side from")
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()

    tpl = a.bundle / "banks" / f"{a.template}.npz"
    with np.load(tpl, allow_pickle=False) as t:
        genes = t["genes"].astype(str)
        residual = t["residual"]
        availability = t["availability"]
        scalar_c = float(t["scalar_c"])
        donors = t["donor_publications"].astype(str)
    if len(donors) != 134:
        raise SystemExit(f"{a.template}: expected a 134-donor template, got {len(donors)}")
    print(f"template {a.template}: {len(genes)} genes, {len(donors)} donors, "
          f"scalar_c={scalar_c:.5f}")

    # every 134-donor bank must agree with the template, or "shared donor side" is wrong
    import hashlib
    h = lambda x: hashlib.sha256(np.ascontiguousarray(x).tobytes()).hexdigest()
    ref = (h(residual), h(availability))
    n_checked = 0
    for p in sorted((a.bundle / "banks").glob("*.npz")):
        with np.load(p, allow_pickle=False) as b:
            if len(b["donor_publications"]) != 134:
                continue
            if (h(b["residual"]), h(b["availability"])) != ref:
                raise SystemExit(f"{p.name}: donor side differs from template")
            if float(b["scalar_c"]) != scalar_c:
                raise SystemExit(f"{p.name}: scalar_c differs from template")
            n_checked += 1
    print(f"donor side verified identical across {n_checked} bundle banks")

    screens = yaml.safe_load(a.manifest.read_text())["screens"]
    donor_set = set(donors.tolist())
    leaks = [(s["dataset_name"], str(s["source_id"])) for s in screens
             if str(s["source_id"]) in donor_set]
    if leaks:
        raise SystemExit("LEAKAGE: target publication(s) among the donors: " + repr(leaks))
    print(f"{len(screens)} targets, 0 publications among the 134 donors")

    a.out.mkdir(parents=True, exist_ok=True)
    mean = np.zeros(32, dtype=np.float64)
    written = []
    for s in screens:
        sid, pub = s["dataset_name"], str(s["source_id"])
        np.savez_compressed(
            a.out / f"{sid}.npz",
            genes=genes, residual=residual, mean=mean, availability=availability,
            scalar_c=np.float64(scalar_c), query_publication=np.array(pub),
            donor_publications=donors,
        )
        written.append(sid)
    print(f"wrote {len(written)} banks to {a.out}")

    # each must round-trip through their own loader and its validation
    import sys
    sys.path.insert(0, "_reproduction")
    from experiments.odysseus.assayloop_comparison.assets import HistoricalBank
    for sid in written:
        bk = HistoricalBank.load(a.out / f"{sid}.npz")
        assert bk.residual.shape == (len(genes), 32), sid
        assert bk.mean.shape == (32,) and not bk.mean.any(), sid
        assert bk.query_publication not in bk.donor_publications, sid
    print(f"all {len(written)} validated by HistoricalBank.load")
    (a.out / "PROVENANCE.json").write_text(json.dumps({
        "template_bank": a.template, "manifest": str(a.manifest),
        "mean": "zeros(32) -- see HANDOFF section 24.10",
        "donor_publications": len(donors), "scalar_c": scalar_c,
        "genes": len(genes), "screens": written,
    }, indent=1))


if __name__ == "__main__":
    main()
