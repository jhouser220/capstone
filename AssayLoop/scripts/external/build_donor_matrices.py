"""Extract the HSA donor-side matrices (M, Y, pubs, genes) from the public HF dataset.

HANDOFF 24.14: the HSA donor set IS the AssayBench train fold -- 1,349 screens over 134
publications, set-identical to any bank's `donor_publications`. So `train_1349_archive.npz`
is reconstructible from `Genentech/assaybench`, config `biogrid`, `yearfold0=='train'`:

  M[i,j] = 1  if genes[j] was measured in donor i   (`relevance_genes`)
  Y[i,j] = 1  if it was called a hit                 (`hit`, parallel to relevance_genes)

The gene axis is intersected with AssayFormer's vocabulary, because `build_bank` requires
`source_weights` aligned to `genes` and the adapter rejects any gene absent from the source
vocab. Axis ORDER is taken from a shipped bank when one is given, so residual/availability
compare element-wise.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--checkpoint", required=True, type=Path)
    ap.add_argument("--reference-bank", type=Path,
                    help="take the gene axis ORDER from this bank, for element-wise comparison")
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()

    import datasets
    from assayloop.experiment.runner import make_model

    d = datasets.load_dataset("Genentech/assaybench", "biogrid", split="train")
    keep = [i for i, f in enumerate(d["yearfold0"]) if f == "train"]
    d = d.select(keep)
    print(f"train fold: {len(d)} screens")

    names = d["dataset_name"]
    pubs = np.array([str(x) for x in d["source_id"]])
    print(f"distinct publications: {len(set(pubs.tolist()))}")

    model = make_model("amortized_ranker", overrides={
        "checkpoint": str(a.checkpoint), "ckpt_file": "model_last.pt",
        "device": "cpu", "desc_emb_by_name": {}}, seed=0)
    vocab = model.vocab
    W_all = model.net.gene_emb.weight.detach().cpu().numpy().copy()
    print(f"AssayFormer vocab {len(vocab)}, gene_emb {W_all.shape}")

    rg = d["relevance_genes"]
    union = set()
    for g in rg:
        union.update(g)
    print(f"union of measured genes across the fold: {len(union)}")
    in_vocab = {g for g in union if vocab.to_idx(str(g)) != 0}
    print(f"of those, present in AssayFormer vocab: {len(in_vocab)}")

    if a.reference_bank:
        with np.load(a.reference_bank, allow_pickle=False) as b:
            ref = b["genes"].astype(str).tolist()
        print(f"reference bank axis: {len(ref)} genes")
        print(f"  reference == in-vocab union as SETS: {set(ref) == in_vocab}")
        print(f"  in reference not in union: {len(set(ref) - in_vocab)}")
        print(f"  in union not in reference: {len(in_vocab - set(ref))}")
        genes = ref
    else:
        genes = sorted(in_vocab)

    idx = {g: j for j, g in enumerate(genes)}
    n, p = len(d), len(genes)
    M = np.zeros((n, p), dtype=np.uint8)
    Y = np.zeros((n, p), dtype=np.uint8)
    hits = d["hit"]
    for i in range(n):
        gl, hl = rg[i], hits[i]
        cols = [idx[g] for g in gl if g in idx]
        sel = [k for k, g in enumerate(gl) if g in idx]
        M[i, cols] = 1
        hv = np.fromiter((bool(hl[k]) for k in sel), dtype=bool, count=len(sel))
        if hv.any():
            Y[i, np.asarray(cols)[hv]] = 1
        if (i + 1) % 250 == 0:
            print(f"  {i+1}/{n} rows")
    print(f"M {M.shape} density {M.mean():.4f}   Y hits total {int(Y.sum()):,}")

    vi = np.array([vocab.to_idx(str(g)) for g in genes], dtype=int)
    assert (vi != 0).all(), "gene absent from source vocabulary"
    W = W_all[vi]
    print(f"source_weights aligned to axis: {W.shape}")

    # AVAILABILITY CHECK -- depends only on M and pubs, so a match proves both are right.
    pm = np.asarray([M[pubs == q].mean(0) for q in sorted(set(pubs.tolist()))])
    avail = (1 + pm.sum(0)) / (2 + len(pm))
    print(f"\navailability: {len(pm)} publication masks, min {avail.min():.10f} "
          f"(= 1/(2+{len(pm)}) = {1/(2+len(pm)):.10f})")
    if a.reference_bank:
        with np.load(a.reference_bank, allow_pickle=False) as b:
            ra = b["availability"].astype(float)
        print(f"  vs reference bank: max abs diff {np.abs(avail - ra).max():.3e}")
        print(f"  BITWISE IDENTICAL: {np.array_equal(avail, ra)}")
        print(f"  allclose(1e-12): {np.allclose(avail, ra, atol=1e-12)}")

    a.out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(a.out, M=M, Y=Y, pubs=pubs, genes=np.array(genes),
                        W=W.astype(np.float32), names=np.array(names))
    print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
