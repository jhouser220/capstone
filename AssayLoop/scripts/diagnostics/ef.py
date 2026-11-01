import json, pathlib, sys
import numpy as np
sys.path.insert(0, "src")
from assayloop.scripts.full_genome_table import _adj_ef_from_run
from assayloop.tasks import gene_universe, load_screens

screens = load_screens(target_set="paper_validation")
universe = set(gene_universe(screens, min_screen_freq=2))
by = {s.dataset_name: s for s in screens}

efs = []
for d in sorted(pathlib.Path("output/runs").glob("*/sweep-afalone-val-*")):
    name = json.loads((d / "config.json").read_text())["screen_name"]
    s = by[name]
    lib = set(s.genes)
    efs.append(_adj_ef_from_run(d / "result.json", lib, universe,
                                s.total_hits, len(lib), 1000))
print(f"n={len(efs)}  mean adjusted EF = {np.mean(efs):.3f}   (HANDOFF records 6.10 for sweep-fd2017be)")
