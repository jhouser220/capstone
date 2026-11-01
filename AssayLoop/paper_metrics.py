"""Compute the paper's Table-1 columns for sweeps you have run, using the authors' own metric code.

    # which sweeps do I have?
    uv run python paper_metrics.py --list

    # score one or more sweeps (10-round handoff sweeps are the comparable ones)
    uv run python paper_metrics.py --screen-set paper_test \\
        --sweeps sweep-aaaa1111 sweep-bbbb2222 --labels "CLI handoff" "fullgenome handoff"

How each column is computed (all means over the screens that have a run; values taken from the
saved result.json files, exactly as scripts/full_genome_table.py does):

  EF      adjusted enrichment factor. Picks in the screen's library count; real genes outside the
          library but inside the shared universe are forgiven; out-of-universe picks AND unfilled
          slots are charged. Imported from the repo (_adj_ef_from_run).
  nAUC %  adjusted normalised area under the cumulative-hits curve (_adj_nauc_from_run).
  FH %    fraction of the screen's hits recovered (final_metrics.frac_hits).
  SF %    the paper's "shortfall": share of acquired genes OUTSIDE the screen's library.
          NOT the same as `shortfall_frac` in result.json (that is the share of requested slots
          left unfilled; it is reported here as "unfilled %").
  %ess    share of recovered hits that are DepMap common-essential genes (needs --essentials and
          the DepMap file that assaybench asks you to supply).
  EP-B/S/D  effective number of Reactome pathway groups per batch / screen / dataset.

EF and nAUC charge unfilled slots against a budget (default 1000 = 10 rounds x 100). A sweep that
ran fewer rounds (e.g. the 3-round Gemini warm start alone) is therefore NOT comparable to the paper
rows unless you pass --budget to match; the script warns when it spots this.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from assayloop import config
from assayloop.scripts.full_genome_table import _adj_ef_from_run, _adj_nauc_from_run
from assayloop.tasks import gene_universe, load_screens

# Paper Table 1 (20 held-out TEST screens). Only meaningful for --screen-set paper_test.
PAPER_ROWS = [
    ("[paper] Random",                          1.04,  3.4,  4.9, 10.0, 18.6, 21.7, 56.9, 82.6),
    ("[paper] Gemini-3.1-Pro (10 rounds)",      4.71, 20.6, 24.6,  3.6, 30.6, 14.6, 34.6, 66.2),
    ("[paper] AssayFormer alone",               4.83, 17.2, 23.2, 10.0, 40.7, 19.3, 46.7, 65.0),
    ("[paper] AssayLoop: Gemini-3.1-Pro -> AF", 5.67, 21.7, 27.7,  7.6, 30.1, 18.4, 46.6, 70.4),
    ("[paper] AssayLoop: GPT-5.6 Sol -> AF",    5.67, 21.7, 27.6,  7.8, 29.4, 18.1, 45.9, 69.8),
]
COLS = ["label", "sweep", "screens", "EF", "nAUC %", "FH %", "SF %", "unfilled %", "%ess",
        "EP-B", "EP-S", "EP-D", "harness EF"]


def run_dir_name(sweep_id: str, idx: int, dataset_name: str) -> str:
    """Same naming rule as the runner: {sweep}-{NN}-{name[:40] with '/' -> '_'}."""
    return f"{sweep_id}-{idx:02d}-{dataset_name[:40].replace('/', '_')}"


def find_result(runs_dir: Path, sweep_id: str, idx: int, dataset_name: str):
    exact = runs_dir / run_dir_name(sweep_id, idx, dataset_name) / "result.json"
    if exact.is_file():
        return exact
    safe = dataset_name[:40].replace("/", "_")
    # fall back to a recursive search so run folders you moved into subfolders (e.g. old/) are still found
    hits = list(runs_dir.rglob(f"{sweep_id}-[0-9][0-9]-{safe}/result.json"))
    return hits[0] if len(hits) == 1 else None


def list_sweeps(sweeps_dir: Path):
    rows = []
    for sj in sorted(sweeps_dir.glob("*/sweep.json"), key=lambda p: p.stat().st_mtime, reverse=True):
        try:
            d = json.loads(sj.read_text())
        except (OSError, ValueError):
            continue
        cfg = d.get("config") or {}
        rows.append(dict(
            sweep=d.get("sweep_id", sj.parent.name), tag=cfg.get("tag"), screen_set=cfg.get("screen_set"),
            model=cfg.get("model"), acq=cfg.get("acq"), n_steps=cfg.get("n_steps"),
            screens=len(d.get("per_screen") or []),
            mean_harness_EF=(d.get("aggregate") or {}).get("mean_n_hits_vs_random"),
        ))
    return pd.DataFrame(rows)


def score_sweep(sweep_id, screens, universe_set, runs_dir, budget, essentials=None, with_ep=True):
    efs, naucs, fhs, sfs, unf, hef, pess, screen_batches, missing, n_rounds = [], [], [], [], [], [], [], [], [], 0
    for i, s in enumerate(screens):
        fp = find_result(runs_dir, sweep_id, i, s.dataset_name)
        if fp is None:
            missing.append(s.dataset_name)
            continue
        lib, H = set(s.genes), int(sum(s.hits))
        rd = json.loads(fp.read_text())
        steps = rd.get("steps") or []
        n_rounds = max(n_rounds, len(steps))
        efs.append(_adj_ef_from_run(fp, lib, universe_set, H, len(lib), budget))
        naucs.append(_adj_nauc_from_run(fp, lib, universe_set, H, len(lib), budget))
        fm = rd.get("final_metrics") or {}
        picks = [(g, h) for st in steps for g, h in zip(st.get("acquired_batch", []), st.get("hits", []))]
        fhs.append(fm["frac_hits"] if "frac_hits" in fm else (sum(1 for _, h in picks if h) / H if H else 0.0))
        sfs.append(sum(1 for g, _ in picks if g not in lib) / len(picks) if picks else np.nan)
        unf.append(1 - len(picks) / budget)
        if fm.get("n_hits_vs_random") is not None:
            hef.append(fm["n_hits_vs_random"])
        if essentials is not None:
            ess_hits = sum(1 for g, h in picks if h and g in lib and g in essentials)
            all_hits = sum(1 for g, h in picks if h and g in lib)
            pess.append(ess_hits / all_hits if all_hits else 0.0)
        batches = [st["acquired_batch"] for st in steps if len(st.get("acquired_batch") or []) >= 2]
        if batches:
            screen_batches.append(batches)

    ep = {}
    if with_ep and screen_batches:
        try:
            from assayloop.metrics.effective_pathways import effective_pathways
            ep = effective_pathways(screen_batches)
        except Exception as e:                                   # keep the other columns
            print(f"  (EP skipped for {sweep_id}: {type(e).__name__}: {e})")
    m = lambda xs, k=1.0: float(np.nanmean(xs)) * k if len(xs) else np.nan
    return dict(
        sweep=sweep_id, screens=f"{len(efs)}/{len(screens)}",
        EF=m(efs), **{"nAUC %": m(naucs, 100), "FH %": m(fhs, 100), "SF %": m(sfs, 100),
                      "unfilled %": m(unf, 100), "%ess": m(pess, 100) if pess else np.nan,
                      "EP-B": ep.get("ep_batch", np.nan), "EP-S": ep.get("ep_screen", np.nan),
                      "EP-D": ep.get("ep_dataset", np.nan), "harness EF": m(hef)},
    ), missing, n_rounds


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--list", action="store_true", help="list sweeps in output/sweeps and exit")
    ap.add_argument("--sweeps", nargs="+", help="sweep ids, e.g. sweep-1a2b3c4d")
    ap.add_argument("--labels", nargs="+", help="display names, same order as --sweeps")
    ap.add_argument("--screen-set", default="paper_test")
    ap.add_argument("--budget", type=int, default=1000, help="requested genes per screen (10 x 100)")
    ap.add_argument("--min-screen-freq", type=int, default=2, help="2 = the paper's f2 universe")
    ap.add_argument("--runs-dir", type=Path, default=config.OUTPUT_PATH / "runs")
    ap.add_argument("--essentials", action="store_true", help="also compute %%ess (needs the DepMap file)")
    ap.add_argument("--no-ep", action="store_true", help="skip the pathway-diversity columns")
    ap.add_argument("--no-paper", action="store_true", help="do not append the paper's reference rows")
    ap.add_argument("--csv", type=Path, default=None)
    args = ap.parse_args(argv)

    if args.list:
        df = list_sweeps(config.OUTPUT_PATH / "sweeps")
        print(df.to_string(index=False) if len(df) else f"no sweeps under {config.OUTPUT_PATH / 'sweeps'}")
        return df
    if not args.sweeps:
        ap.error("pass --sweeps (or --list to see what exists)")
    labels = args.labels or args.sweeps
    if len(labels) != len(args.sweeps):
        ap.error("--labels must have one entry per --sweeps")

    screens = load_screens(target_set=args.screen_set)
    universe = gene_universe(screens, min_screen_freq=args.min_screen_freq)
    universe_set = set(universe)
    print(f"{len(screens)} screens in {args.screen_set!r}; shared universe = {len(universe)} genes; "
          f"budget = {args.budget}")

    essentials = None
    if args.essentials:
        from assayloop.scripts.full_genome_table import _load_essentials
        essentials = _load_essentials()                          # raises with setup steps if absent

    rows = []
    for label, sweep_id in zip(labels, args.sweeps):
        row, missing, n_rounds = score_sweep(sweep_id, screens, universe_set, args.runs_dir,
                                             args.budget, essentials, with_ep=not args.no_ep)
        if missing:
            print(f"WARNING {sweep_id}: no run found for {len(missing)} screen(s): {missing}")
        if n_rounds and n_rounds * 100 < args.budget:
            print(f"WARNING {sweep_id}: runs have {n_rounds} rounds (~{n_rounds * 100} genes) but "
                  f"--budget is {args.budget}; EF/nAUC charge the missing slots. Not comparable to "
                  "the paper's 10-round rows (pass --budget to change it).")
        rows.append({"label": label, **row})

    df = pd.DataFrame(rows)
    if args.screen_set == "paper_test" and not args.no_paper:
        ref = pd.DataFrame(PAPER_ROWS, columns=["label", "EF", "nAUC %", "FH %", "SF %", "%ess",
                                                "EP-B", "EP-S", "EP-D"])
        df = pd.concat([df, ref], ignore_index=True)
    df = df.reindex(columns=COLS)
    with pd.option_context("display.width", 250, "display.max_columns", 30):
        print("\n" + df.to_string(index=False, float_format=lambda x: f"{x:.2f}"))
    if args.csv:
        df.to_csv(args.csv, index=False)
        print(f"\nwrote {args.csv}")
    return df


if __name__ == "__main__":
    main()
