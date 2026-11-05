"""Rebuild REAL HistoricalBanks via their own `build_bank`, from reconstructed donor arrays.

Everything is an actual value, not an estimate:
  M, Y, pubs, genes   reconstructed from Genentech/assaybench (HANDOFF 24.14); the resulting
                      `availability` matches a shipped bank to 3.3e-16, which validates them
  W (source_weights)  AssayFormer gene_emb.weight, rows aligned to the axis
  E (donor desc)      shipped text-embedding cache, 1,349 train screens, 1536-D
  query_description   shipped cache (bitwise identical to the bundle's own, HANDOFF 24.8)

`--verify <id>` rebuilds a screen that ships a frozen bank and compares, reporting the
SIGN-INVARIANT quantity `residual @ mean` (the round-1 prior contribution) alongside
`scalar_c` and `availability`, because SVD component sign/order is arbitrary.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np


def descriptions(names):
    from assayloop.amortized.text_embed import (
        EMBED_MODEL, EmbeddingCache, _text_key, screen_description_text)
    from assayloop.tasks.screen_sets import load_screens
    cache = EmbeddingCache()
    recs = {r.dataset_name: r for r in load_screens(target_set="train")}
    out, missing = [], []
    for n in names:
        r = recs.get(n)
        v = None if r is None else cache.get(
            _text_key(screen_description_text(r), model=EMBED_MODEL))
        if v is None:
            missing.append(n)
            out.append(np.zeros(1536, dtype=np.float32))
        else:
            out.append(np.asarray(v, dtype=np.float32))
    if missing:
        raise SystemExit(f"{len(missing)} donor description(s) not cached: {missing[:10]} "
                         "-- refusing to call the embedding API.")
    return np.stack(out)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--donor", required=True, type=Path)
    ap.add_argument("--manifest", help="target manifest yaml (dataset_name + source_id)")
    ap.add_argument("--target-set", help="assayloop screen-set name, for the query descriptions")
    ap.add_argument("--verify", help="bundle screen id to rebuild and compare")
    ap.add_argument("--bundle", type=Path, default=Path("_reproduction/bundle"))
    ap.add_argument("--screens", type=Path,
                    help="pickle of externally built ScreenRecords instead of --manifest/"
                         "--target-set. No external paper is one of the 134 donor "
                         "publications, so the excluded-publication set is EMPTY for every "
                         "target: the donor side is shared and only `mean` varies per screen "
                         "(HANDOFF 24.33). query_publication is written as EXTERNAL.")
    ap.add_argument("--skip-existing", action="store_true",
                    help="leave banks already in --out alone. build_bank holds several "
                         "1349x23742 float64 copies, about 1.5-2 GB peak, and a 12-target "
                         "run was OOM-killed at the ninth on a 14 GB box -- so this exists "
                         "to finish a partial directory a few targets at a time.")
    ap.add_argument("--limit", type=int, default=0,
                    help="build at most N targets this invocation (with --skip-existing, "
                         "the next N still missing). Keeps peak memory to one process.")
    ap.add_argument("--out", type=Path)
    a = ap.parse_args()

    sys.path.insert(0, "_reproduction")
    from experiments.odysseus.assayloop_comparison.assets import build_bank

    z = np.load(a.donor, allow_pickle=False)
    M, Y = z["M"].astype(np.float64), z["Y"].astype(np.float64)
    pubs, genes, W = z["pubs"].astype(str), z["genes"].astype(str), z["W"].astype(np.float64)
    names = z["names"].astype(str)
    print(f"donors {M.shape}, {len(set(pubs.tolist()))} publications, axis {len(genes)}")
    E = descriptions(names.tolist())
    print(f"donor descriptions {E.shape} (all from cache, 0 API calls)")

    from assayloop.amortized.text_embed import (
        EMBED_MODEL, EmbeddingCache, _text_key, screen_description_text)
    cache = EmbeddingCache()

    def build(qdesc, qpub):
        return build_bank(genes=list(genes), source_weights=W, donor_mask=M, donor_hits=Y,
                          donor_publications=pubs, donor_descriptions=E,
                          query_description=qdesc, query_publication=str(qpub))

    if a.verify:
        with np.load(a.bundle / "banks" / f"{a.verify}.npz", allow_pickle=False) as b:
            ref = {k: b[k] for k in b.files}
        qd = np.load(a.bundle / "screens" / f"{a.verify}.npy", allow_pickle=False).astype(float)
        qpub = str(ref["query_publication"].item())
        print(f"\nrebuilding {a.verify} (query publication {qpub}) ...")
        bk = build(qd, qpub)
        rr, rm = ref["residual"].astype(float), ref["mean"].astype(float)
        print(f"  donors            ours {len(bk.donor_publications)}  ref {len(ref['donor_publications'])}")
        print(f"  scalar_c          ours {bk.scalar_c:.5f}  ref {float(ref['scalar_c']):.5f}  "
              f"ratio {bk.scalar_c/float(ref['scalar_c']):.4f}")
        print(f"  availability      max abs diff {np.abs(bk.availability-ref['availability'].astype(float)).max():.3e}")
        # sign-invariant: the round-1 prior contribution per gene
        ours, theirs = bk.residual @ bk.mean, rr @ rm
        c = np.corrcoef(ours, theirs)[0, 1]
        print(f"  residual @ mean   pearson r = {c:+.4f}   spearman-ish rank r = "
              f"{np.corrcoef(np.argsort(np.argsort(ours)), np.argsort(np.argsort(theirs)))[0,1]:+.4f}")
        print(f"                    ||ours|| {np.linalg.norm(ours):.5f}  ||theirs|| {np.linalg.norm(theirs):.5f}")
        print(f"  ||mean||          ours {np.linalg.norm(bk.mean):.5f}  ref {np.linalg.norm(rm):.5f}")
        # subspace agreement, order/sign free
        qa, _ = np.linalg.qr(bk.residual); qb, _ = np.linalg.qr(rr)
        sv = np.linalg.svd(qa.T @ qb, compute_uv=False)
        print(f"  32-D subspace     principal cosines: mean {sv.mean():.4f}, "
              f"min {sv.min():.4f}, max {sv.max():.4f}")
        return

    if a.screens:
        import pickle
        with open(a.screens, "rb") as fh:
            ext = pickle.load(fh)
        # "EXTERNAL" matches no donor publication, so build_bank excludes nothing --
        # which is correct here and is asserted below.
        scr = [dict(dataset_name=r.dataset_name, source_id="EXTERNAL") for r in ext]
        recs = {r.dataset_name: r for r in ext}
        print(f"{len(scr)} external targets from {a.screens}; "
              "no donor publication is excluded for any of them")
    else:
        import yaml
        scr = yaml.safe_load(Path(a.manifest).read_text())["screens"]
        from assayloop.tasks.screen_sets import load_screens
        recs = {r.dataset_name: r for r in load_screens(target_set=a.target_set)}
    n_donor_pubs = len(set(pubs.tolist()))
    a.out.mkdir(parents=True, exist_ok=True)
    if a.skip_existing:
        n0 = len(scr)
        scr = [s for s in scr if not (a.out / f"{s['dataset_name']}.npz").exists()]
        print(f"skip-existing: {n0 - len(scr)} already built, {len(scr)} to go")
    if a.limit:
        scr = scr[:a.limit]
        print(f"limit: building {len(scr)} this invocation")
    if not scr:
        print("nothing to build")
        return
    meta = {}
    for i, s in enumerate(scr):
        sid, qpub = s["dataset_name"], str(s["source_id"])
        r = recs[sid]
        qd = np.asarray(cache.get(_text_key(screen_description_text(r), model=EMBED_MODEL)),
                        dtype=np.float64)
        bk = build(qd, qpub)
        if a.screens:
            # every donor retained: the whole point of an out-of-benchmark target
            assert len(bk.donor_publications) == n_donor_pubs, (
                f"{sid}: {len(bk.donor_publications)} donors, expected all {n_donor_pubs} "
                "-- an external target must not collide with a donor publication")
        np.savez_compressed(a.out / f"{sid}.npz", genes=np.array(bk.genes),
                            residual=bk.residual, mean=bk.mean, availability=bk.availability,
                            scalar_c=np.float64(bk.scalar_c),
                            query_publication=np.array(qpub),
                            donor_publications=np.array(bk.donor_publications))
        meta[sid] = dict(scalar_c=float(bk.scalar_c), mean_norm=float(np.linalg.norm(bk.mean)),
                         donors=len(bk.donor_publications))
        print(f"  [{i+1}/{len(scr)}] {sid:<16} donors={len(bk.donor_publications)} "
              f"scalar_c={bk.scalar_c:.4f} ||mean||={np.linalg.norm(bk.mean):.5f}")
    prov_path = a.out / "PROVENANCE.json"
    if prov_path.exists():
        try:
            meta = {**json.loads(prov_path.read_text()).get("banks", {}), **meta}
        except Exception:
            pass
    prov_path.write_text(json.dumps(
        {"source": "build_bank on reconstructed donor arrays (HANDOFF 24.14/24.20)",
         "donor_file": str(a.donor),
         "target_set": a.target_set or f"external:{a.screens}", "banks": meta}, indent=1))
    print(f"\nwrote {len(scr)} REAL banks to {a.out}")


if __name__ == "__main__":
    main()
