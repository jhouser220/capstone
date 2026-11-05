"""Build `results/main_results.csv` — the one table behind every DepMap/HSA comparison.

Reads the per-run JSONs under `output/runs/external/` (gitignored, machine-local) and emits
one tidy, committable CSV: one row per (set, pool, protocol, arm).

THREE AXES, and a row is only comparable to another row that shares the first three:

  set       which screens. paper_test / paper_validation are the curated 20-screen panels
            (`num_hits_min: 50`, so no sparse screens); `validation fold` and `TEST fold`
            are whole AssayBench folds, genome-scale subset only; `team11 external` is the
            teammate's 11 out-of-benchmark screens.
  pool      what the policy could pick from. `f2` = the shared 21,147-gene universe the
            paper uses "for fair LLM comparison"; `library` = each screen's own genes.
            Under `f2` some picks land outside the screen's library and are REFUNDED from
            the EF denominator, and the refund is per-arm, so it reorders arms. Under
            `library` nothing can land outside, so all EF conventions agree.
  protocol  what rounds 1-3 did. `cold` = no LLM. `retain25` = keep GLM-5's 25 best-scoring
            proposals and let the arm fill the other 75. `warm` = the whole recorded round,
            shortfall lost.

TWO EF COLUMNS, because the project has two conventions and they rank arms differently:

  ef_screen / ef_pub   CANONICAL adjusted EF = hits / (eff * H/D), D = the screen's own
                       library, eff = 1000 - n2 (off-library-but-in-pool picks refunded;
                       unfilled slots stay charged). This is the paper's metric -- our
                       AssayFormer reads 4.828 against their published 4.83. `_screen` is a
                       plain screen mean, `_pub` averages publications. They can differ by
                       more than 1.0 where one publication dominates a set.
  ef_nvr               `n_hits_vs_random` = hits / (1000 * H/pool). What every result.json
                       stores and what `compare_hsa_depmap.ipynb` cells 0-37 plot. NOT the
                       paper's metric; do not put it beside 4.83 / 5.27 / 5.67.

Intervals are publication-clustered (screens from one paper are near-duplicates).
"""
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

BUDGET = 1000
ARMS = ["AF alone", "AF + DepMap", "HSA only (ours)",
        "DepMap + HSA dm.25", "DepMap + HSA dm.50"]
SHORT = {"AF alone": "AF", "AF + DepMap": "AF + DepMap", "HSA only (ours)": "HSA",
         "DepMap + HSA dm.25": "DepMap + HSA .25", "DepMap + HSA dm.50": "DepMap + HSA .50"}
BASE = "AF + DepMap"

# (set, pool, protocol, llm budget share, [files])
# FILE ORDER MATTERS: the first file an arm appears in wins, and the REAL-MEAN file must
# come first -- every own-library run uses `banks_real_*`, so putting the mean=0 file first
# silently compares different bank variants.
SOURCES = [
    ("paper_test",       "f2",      "cold",     0.000,
     ["paper_test_hsa_realmean.json", "paper_test_hsa.json"]),
    ("paper_test",       "library", "cold",     0.000, ["paper_test_hsa_libpool.json"]),
    ("paper_test",       "f2",      "retain25", 0.075, ["glm5_retain25_hsa.json"]),
    ("paper_test",       "f2",      "warm",     0.291, ["glm5_handoff_hsa.json"]),
    ("paper_validation", "f2",      "cold",     0.000,
     ["paper_validation_hsa_realmean.json", "paper_validation_hsa.json"]),
    ("paper_validation", "library", "cold",     0.000, ["paper_validation_hsa_libpool.json"]),
    ("validation fold",  "f2",      "cold",     0.000, ["validation_fold_hsa_REALCOMBINED.json"]),
    ("validation fold",  "library", "cold",     0.000, ["validation_fold_hsa_libpool.json"]),
    ("TEST fold",        "library", "cold",     0.000, ["test_fold_hsa_libpool.json"]),
    ("team11 external",  "library", "cold",     0.000, ["team11_hsa.json"]),
]

TEAM11_PAPER = {
    "ABT_PARCB_dependency": "abt", "NGUYEN_AraC_resistant": "nguyen",
    "ORGANOID_WTSI-COLIVM_031_minLib": "organoid", "ORGANOID_WTSI-OESO_005_minLib": "organoid",
    "ORGANOID_WTSI-OVAR_018_minLib": "organoid", "SCHUKKEN_PAIRED_WGL_paper": "schukken",
    "SCHUKKEN_PAIRED_WGL_rank194": "schukken", "SHAN_HeLa_LDLRko": "shan",
    "SHAN_HeLa_lovastatin": "shan", "SHAN_HepG2_lovastatin": "shan",
    "TAU_CRISPRi_dec": "tau", "WILSON_MERCS_inc": "wilson",
}


def manifests(root: Path, venv_site: Path) -> dict[str, Path]:
    return {
        "paper_test":       venv_site / "assaybench/data/screen_sets/assayloop-test.yaml",
        "paper_validation": venv_site / "assaybench/data/screen_sets/assayloop-validation.yaml",
        "validation fold":  root / "_reproduction/donor/validation_fold_218.yaml",
        "TEST fold":        root / "_reproduction/donor/test_fold_334.yaml",
    }


def pubmap(sset: str, mf: dict[str, Path]) -> dict[str, str]:
    if sset == "team11 external":
        return dict(TEAM11_PAPER)
    p = mf.get(sset)
    if p is None or not Path(p).exists():
        return {}
    return {s["dataset_name"]: str(s["source_id"])
            for s in yaml.safe_load(Path(p).read_text())["screens"]}


def boot(by_pub: dict[str, list[float]], n: int, seed: int = 0):
    k = list(by_pub)
    if not k:
        return np.nan, np.nan, np.nan
    m = np.array([np.mean(by_pub[x]) for x in k])
    if len(k) < 2:
        return float(m.mean()), np.nan, np.nan
    rng = np.random.default_rng(seed)
    bs = m[rng.integers(0, len(k), (n, len(k)))].mean(1)
    return float(m.mean()), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", type=Path, default=Path("output/runs/external"))
    ap.add_argument("--root", type=Path, default=Path("."))
    ap.add_argument("--venv-site", type=Path,
                    default=Path(".venv/lib/python3.11/site-packages"))
    ap.add_argument("--out", type=Path, default=Path("results/main_results.csv"))
    ap.add_argument("--n", type=int, default=20000)
    a = ap.parse_args()

    mf = manifests(a.root, a.venv_site)
    rows, missing = [], []
    for sset, pool, proto, share, files in SOURCES:
        pm = pubmap(sset, mf)
        seen, recs, src = set(), [], {}
        for fn in files:
            p = a.runs / fn
            if not p.exists():
                continue
            for r in json.loads(p.read_text()):
                if r["arm"] not in ARMS or (r["arm"], r["screen"]) in seen:
                    continue
                seen.add((r["arm"], r["screen"]))
                src.setdefault(r["arm"], fn)
                picks = int(r.get("picks", BUDGET))
                n_meas = int(r.get("n_measured", picks))
                eff = BUDGET - (picks - n_meas)        # unfilled slots stay charged
                prev = r["total_hits"] / r["n_genes"]
                recs.append(dict(arm=r["arm"], screen=r["screen"],
                                 pub=pm.get(r["screen"], "?" + r["screen"]),
                                 hits=int(r["hits"]), picks=picks, n_measured=n_meas,
                                 n_genes=int(r["n_genes"]), total_hits=int(r["total_hits"]),
                                 ef=r["hits"] / (eff * prev),
                                 ef_nvr=float(r.get("ef", np.nan))))
        if not recs:
            missing.append((sset, pool, proto, files))
            continue
        D = pd.DataFrame(recs)
        base = D[D.arm == BASE].set_index("screen")
        for arm in ARMS:
            g = D[D.arm == arm]
            if g.empty:
                continue
            _, _, _ = None, None, None
            ef_pub, _, _ = boot({p: v.ef.tolist() for p, v in g.groupby("pub")}, a.n)
            d_by_pub, wins, tot, dh = defaultdict(list), 0, 0, 0
            if arm != BASE and not base.empty:
                for r in g.itertuples():
                    if r.screen not in base.index:
                        continue
                    d = r.ef - base.at[r.screen, "ef"]
                    d_by_pub[r.pub].append(d)
                    wins += d > 0
                    tot += 1
                    dh += r.hits - base.at[r.screen, "hits"]
            dm, dlo, dhi = boot(d_by_pub, a.n) if d_by_pub else (0.0, np.nan, np.nan)
            rows.append(dict(
                set=sset, pool=pool, protocol=proto, llm_budget_share=share,
                arm=SHORT[arm], arm_recorded=arm,
                screens=g.screen.nunique(), publications=g.pub.nunique(),
                ef_screen=g.ef.mean(), ef_pub=ef_pub, ef_nvr=g.ef_nvr.mean(),
                hits=int(g.hits.sum()),
                hits_available=int(g.groupby("screen").total_hits.first().sum()),
                hits_per_1000_slots=g.hits.sum() / (g.screen.nunique() * BUDGET) * 1000,
                mean_measured_picks=g.n_measured.mean(),
                mean_off_library_picks=(g.picks - g.n_measured).mean(),
                d_ef_vs_af_depmap=dm, d_ef_lo=dlo, d_ef_hi=dhi,
                d_hits_vs_af_depmap=dh, win_rate=(wins / tot) if tot else np.nan,
                significant="*" if (dlo == dlo and dhi == dhi
                                    and (dlo > 0 or dhi < 0)) else "",
                source_file=src.get(arm, "")))

    T = pd.DataFrame(rows)
    T["pct_of_available_hits"] = 100 * T.hits / T.hits_available.replace(0, np.nan)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    T.round(4).to_csv(a.out, index=False)
    if missing:
        print("MISSING (row omitted):")
        for m in missing:
            print("  ", m)
    print(f"wrote {a.out}  ({len(T)} rows, {T['set'].nunique()} sets)")
    print(T[["set", "pool", "protocol", "arm", "screens", "ef_screen", "hits",
             "d_ef_vs_af_depmap", "d_hits_vs_af_depmap", "significant"]]
          .round(3).to_string(index=False))


if __name__ == "__main__":
    main()
