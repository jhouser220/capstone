"""Rebuild the unlabeled Achilles features used by the numerical pilots."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def build(source, output):
    with source.open(newline="") as f:
        columns = next(csv.reader(f))[1:]
    genes = [c.split(" (")[0].upper() for c in columns]
    if len(genes) != len(set(genes)):
        raise ValueError("Ambiguous gene symbols in Achilles source")
    frame = pd.read_csv(source, usecols=columns).dropna(axis=0)
    values = frame.to_numpy(dtype=np.float32).T
    norms = np.linalg.norm(values, axis=1, keepdims=True)
    if np.any(norms == 0) or not np.isfinite(values).all():
        raise ValueError("Invalid feature values")
    values /= norms
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, genes=np.asarray(genes), values=values)
    record = dict(source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                  feature_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),
                  shape=list(values.shape), dtype=str(values.dtype),
                  source_url="https://figshare.com/ndownloader/files/49843176",
                  transformation="Original column order; uppercase symbols; drop incomplete cell-line rows; float32; normalize each gene to unit L2 norm")
    output.with_suffix('.manifest.json').write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--csv', type=Path, required=True)
    parser.add_argument('--output', type=Path,
                        default=Path(__file__).resolve().parent/'data/achilles_normalized.npz')
    args = parser.parse_args()
    build(args.csv, args.output)
