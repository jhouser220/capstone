"""One-off: cache the two DepMap matrices as .npz so repeated screens reload in
seconds instead of minutes. Wide CSVs (~18k columns) are the bottleneck, not the fits.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

D = Path(sys.argv[1])
for name in ("achilles", "Expression"):
    src, dst = D / f"{name}.csv", D / f"{name}.npz"
    if dst.is_file() or not src.is_file():
        print(f"skip {name}")
        continue
    df = pd.read_csv(src, index_col=0, low_memory=False)
    df.columns = [c.split(" (")[0] for c in df.columns]
    df = df.loc[:, ~df.columns.duplicated()]
    np.savez_compressed(dst, values=df.to_numpy(dtype=np.float32),
                        rows=np.array(df.index, dtype=object),
                        cols=np.array(df.columns, dtype=object))
    print(f"{name}: {df.shape} -> {dst.name}")
