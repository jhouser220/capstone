"""Turn the Abt & Wang 2026 PARCB CRISPR dependency screen into a ScreenRecord.

A second EXTERNAL held-out test: PNAS, published 2026-03-20, so like Schukken it
post-dates every AssayBench year fold.

Source: Abt, Wang, Varuzhanyan, Freeland et al., "Synthetic lethality between
RB-loss and E2F3 inhibition in small cell cancers targeted by pyrimidine
synthesis blockade", PNAS 2026, DOI 10.1073/pnas.2532814123 (PMC13012052).
Data: ``Zhang_count/merged_spacer_counts_mageck_test.gene_summary_FILT.txt``
from github.com/jfreeland01/Abt_Wang_E2F3_CRISPR -- the authors' own final
gene-level MAGeCK output, the one the manuscript's figures are built on.

**This dataset is ONE screen, not a panel.** The authors ran ``mageck test`` in
paired mode across all nine replicates at once (PARCB1, PARCB3 and PARCB9 x
three replicates, 2 weeks vs starting material), so the published analysis is a
single pooled dependency call rather than three per-line screens. Deriving
per-line screens would mean re-running MAGeCK in a configuration the paper does
not report, and re-deciding its >27-good-sgRNA filter, which is keyed to the
nine pooled replicates. That is a different experiment, so it is not done here.
The consequence is that this test has n=1 and therefore no error bars; what it
does give is the ten budgets of one screen, and hit calls that need no judgment
from me.

**The hit definition is the authors' own, exactly.** ``neg.p.value < 0.05``
reproduces their reported count of 1398 essential genes to the gene, and they
chose a nominal p-value over FDR deliberately (FDR < 0.05 recovered only
68/1398 of common-essential genes). It needs no AssayBench-matching either: the
resulting hit fraction, 0.0967 of 14459 genes, already sits on the AssayBench
genome-scale knockout-fitness median of 0.094.

Gene symbols are remapped through HGNC by the same collision-free routine as
the Schukken build -- GeCKO v2 is a 2015 library, so it has the same old-symbol
problem.
"""

from __future__ import annotations

import argparse
import json
import logging
import pickle
from pathlib import Path

import pandas as pd

from assaybench.data.screens import ScreenRecord
from build_schukken import hgnc_synonyms, remap_library

log = logging.getLogger("build_abt")

HIT_P = 0.05


def build(summary: Path, hgnc: Path, vocab: Path, hit_p: float) -> list[ScreenRecord]:
    d = pd.read_csv(summary, sep="\t")
    syn = hgnc_synonyms(hgnc)
    known = set(json.loads(vocab.read_text())["itos"])

    raw = [str(g).strip() for g in d["id"]]
    genes = remap_library(raw, known, syn)
    n_oov = sum(1 for g in genes if g not in known)
    log.info("library %d genes, %d unique; %d (%.1f%%) outside the model vocabulary",
             len(genes), len(set(genes)), n_oov, 100 * n_oov / len(genes))
    if len(set(genes)) != len(genes):
        raise RuntimeError(f"gene list not unique ({len(set(genes))} of {len(genes)})")

    p = pd.to_numeric(d["neg.p.value"], errors="coerce").fillna(1.0)
    lfc = pd.to_numeric(d["neg.lfc"], errors="coerce").fillna(0.0)
    hits = [bool(v < hit_p) for v in p]
    log.info("hits at neg.p.value < %g: %d (frac %.4f)",
             hit_p, sum(hits), sum(hits) / len(hits))

    condition = (
        "Genome-wide CRISPR-Cas9 knockout dependency screen in PARCB small cell "
        "neuroendocrine prostate carcinoma models (PARCB1, PARCB3, PARCB9), "
        "GeCKO v2 library, sgRNA abundance after fourteen days of culture "
        "versus the starting material"
    )
    return [ScreenRecord(
        dataset_name="ABT_PARCB_dependency",
        split="external",
        organism="Homo sapiens",
        gene_symbol_convention="HGNC",
        phenotype=("Genes whose knockout reduces proliferation and is depleted "
                   "from the population over fourteen days (MAGeCK)"),
        cleaned_phenotype="Fitness / Proliferation / Viability",
        cell_line="PARCB1 / PARCB3 / PARCB9",
        cell_type="small cell neuroendocrine prostate carcinoma",
        library_type="GeCKO v2 whole-genome CRISPR knockout library",
        library_methodology="Knockout",
        direction_str="decrease",
        condition_clause=condition,
        num_genes=len(genes),
        genes=list(genes),
        # Dropout magnitude for hits, masked to 0 for non-hits (AssayBench
        # convention); the unmasked column keeps every gene's LFC.
        relevance_scores=[float(max(0.0, -v)) if h else 0.0
                          for v, h in zip(lfc, hits)],
        hits=hits,
        description=condition,
        contrast_label="2 weeks vs starting material",
        notes=("Abt & Wang et al. 2026 (PNAS, PMC13012052). Hits are the "
               f"authors' own dependency calls, MAGeCK neg.p.value < {hit_p}, "
               "pooled over nine replicates from three PARCB lines."),
        duration="fourteen days",
        unmasked_relevance_scores=[float(-v) for v in lfc],
    )]


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", required=True, type=Path)
    ap.add_argument("--hgnc", required=True, type=Path)
    ap.add_argument("--vocab", required=True, type=Path)
    ap.add_argument("--hit-p", type=float, default=HIT_P)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()

    screens = build(a.summary, a.hgnc, a.vocab, a.hit_p)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "wb") as fh:
        pickle.dump(screens, fh)
    s = screens[0]
    print(f"\n{s.dataset_name}: {s.num_genes} genes, {s.total_hits} hits "
          f"(frac {s.total_hits / s.num_genes:.4f}), "
          f"EF ceiling {s.num_genes / max(s.total_hits, 1000):.2f}")
    print(f"-> {a.out}")


if __name__ == "__main__":
    main()
