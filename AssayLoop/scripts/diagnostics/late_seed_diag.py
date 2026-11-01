"""Late-seed pathway-propagation diagnostic.

Question: after 7 rounds of a campaign, do the confirmed hits still carry enough
pathway signal to rank the REMAINING untested pool better than base rate -- and
better than AssayFormer's own rounds 8-10?

Seeds from rounds 1..SEED_ROUNDS hits; evaluates on the untested pool at that point.
Ground-truth labels come from the screen records, not the run.
"""
import argparse, json, sys
from pathlib import Path

import numpy as np

REPO = Path("/home/jeff/Documents/Textbooks/DATASCI_210/capstone/AssayLoop")
GMT = REPO / "src/assayloop/data/gene_sets/ReactomePathways.gmt"
MIN_P, MAX_P = 2, 150          # module size window (the cap is load-bearing)
MIN_SEEDS = 2                  # hot module needs >=2 confirmed hits


def load_modules():
    mods = []
    for line in GMT.read_text().splitlines():
        parts = line.rstrip("\n").split("\t")
        genes = {g for g in parts[2:] if g}
        if MIN_P <= len(genes) <= MAX_P:
            mods.append((parts[0], genes))
    return mods


def score_pool(mods, seeds, pool):
    """score(g) = sum over hot modules containing g of c(p)/|p|."""
    seeds, pool = set(seeds), set(pool)
    scores = {}
    n_hot = 0
    for _name, genes in mods:
        c = len(seeds & genes)
        if c < MIN_SEEDS:
            continue
        n_hot += 1
        w = c / len(genes)
        for g in genes & pool:
            scores[g] = scores.get(g, 0.0) + w
    return scores, n_hot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sweep-dir", required=True)
    ap.add_argument("--sweep", help="only run dirs whose name starts with this")
    ap.add_argument("--screen-set", default="paper_validation")
    ap.add_argument("--seed-rounds", type=int, default=7)
    ap.add_argument("--eval-rounds", type=int, default=3)
    ap.add_argument("--batch", type=int, default=100)
    ap.add_argument("--out")
    a = ap.parse_args()

    sys.path.insert(0, str(REPO / "src"))
    from assayloop.tasks import load_screens

    screens = {s.dataset_name: s for s in load_screens(target_set=a.screen_set)}
    mods = load_modules()
    print(f"modules {MIN_P}-{MAX_P} genes: {len(mods)}   screens: {len(screens)}", file=sys.stderr)

    rows = []
    for run in sorted(Path(a.sweep_dir).iterdir()):
        rj, cj = run / "result.json", run / "config.json"
        if not rj.is_file():
            continue
        if a.sweep and not run.name.startswith(a.sweep):
            continue
        res = json.loads(rj.read_text())
        name = json.loads(cj.read_text()).get("screen_name") or run.name.split("-", 2)[-1]
        sc = screens.get(name)
        if sc is None:
            print(f"  skip {name}: not in {a.screen_set}", file=sys.stderr)
            continue

        universe = json.loads(cj.read_text())["config"].get("universe_genes") or sc.genes
        steps = res.get("steps") or []
        if len(steps) < a.seed_rounds + a.eval_rounds:
            print(f"  skip {name}: only {len(steps)} rounds", file=sys.stderr)
            continue

        hit_of = dict(zip(sc.genes, sc.hits))
        rel_of = dict(zip(sc.genes, sc.relevance_scores))
        # strong hits = top quartile of relevance among this screen's hits
        hit_rels = [rel_of[g] for g in sc.genes if hit_of.get(g)]
        thr = float(np.percentile(hit_rels, 75)) if hit_rels else np.inf
        strong = {g for g in sc.genes if hit_of.get(g) and rel_of[g] >= thr}

        seed_batches = [steps[i].get("acquired_batch") or [] for i in range(a.seed_rounds)]
        eval_batches = [steps[i].get("acquired_batch") or []
                        for i in range(a.seed_rounds, a.seed_rounds + a.eval_rounds)]
        tested = {g for b in seed_batches for g in b}
        seeds = [g for g in tested if hit_of.get(g)]

        pool = [g for g in universe if g not in tested]
        pool_hits = {g for g in pool if hit_of.get(g)}
        if not pool or not pool_hits or len(seeds) < MIN_SEEDS:
            print(f"  skip {name}: seeds={len(seeds)} pool_hits={len(pool_hits)}", file=sys.stderr)
            continue
        base = len(pool_hits) / len(pool)
        base_strong = len(strong & set(pool)) / len(pool)

        scores, n_hot = score_pool(mods, seeds, pool)
        cands = [g for g, v in scores.items() if v > 0]
        ranked = sorted(cands, key=lambda g: (-scores[g], g))
        k = a.batch * a.eval_rounds

        def prec(genes, target):
            genes = list(genes)
            return (sum(1 for g in genes if g in target) / len(genes)) if genes else float("nan")

        af_eval = [g for b in eval_batches for g in b]
        af_hits = sum(1 for g in af_eval if hit_of.get(g))

        rows.append(dict(
            screen=name, n_seeds=len(seeds), n_hot=n_hot, pool=len(pool),
            pool_hits=len(pool_hits), base_pct=100 * base,
            n_cands=len(cands),
            cand_pct=100 * prec(cands, pool_hits),
            top_k_pct=100 * prec(ranked[:k], pool_hits),
            af_late_pct=100 * af_hits / max(1, len(af_eval)),
            af_late_hits=af_hits,
            top_k_hits=sum(1 for g in ranked[:k] if g in pool_hits),
            base_strong_pct=100 * base_strong,
            cand_strong_pct=100 * prec(cands, strong),
            overlap_af=len(set(ranked[:k]) & set(af_eval)),
        ))
        print(f"  {name:18s} seeds={len(seeds):4d} hot={n_hot:4d} cands={len(cands):5d} "
              f"lift={rows[-1]['cand_pct']/rows[-1]['base_pct']:5.2f}x "
              f"topK={rows[-1]['top_k_pct']:5.2f}% vs AF {rows[-1]['af_late_pct']:5.2f}%",
              file=sys.stderr)

    if a.out:
        Path(a.out).write_text(json.dumps(rows, indent=1))
    print(json.dumps(rows))


if __name__ == "__main__":
    main()
