"""Turn the Wilson et al. 2024 MERCS FACS screen into ScreenRecords.

A third EXTERNAL test, and the first one that is NOT a fitness screen. Where
Schukken and Abt are knockout *dropout* screens on proliferation, this is a
FACS-sorted *enrichment* screen on a reporter phenotype -- the opposite corner
of AssayBench's type space, and the corner where the DepMap gain was measured
strongest (enrichment +0.423 [+0.247, +0.634], molecular-output phenotypes).
It is also SPARSE: ~2% of the library are hits, against ~10% for the two
fitness screens, which puts it in the same regime as AssayBench's own test
screens (median 220 hits).

Source: Wilson, Yu, Leal, Woodward, Patikas et al., "Genome-wide CRISPR/Cas9
screen shows that loss of GET4 increases mitochondria-endoplasmic reticulum
contact sites and is neuroprotective", Cell Death & Disease 15:203 (2024),
PMC10928201. Data: ``Supplementary Table 2_Whole Genome Pooled Screening
Results.csv`` from github.com/m1gus/MitoER-CRISPR-Screen -- the authors' own
MAGeCK RRA gene-level output.

**Leakage check.** This one is 2024, inside AssayBench's year range, so it had
to be checked rather than assumed: no screen in any fold mentions GET4, BAG6,
MERCS, mVenus or this study's authors, and the only HeLa reporter screens in
the train fold are mitophagy (mt-mKeima) screens from a different study.

**Three screens, following AssayBench's own convention** for sorted screens
(it splits them into ``_inc`` / ``_dec`` / ``_merged``):

- ``_inc``  -- genes whose loss INCREASES contact sites (MAGeCK ``pos``), 394
- ``_dec``  -- genes whose loss DECREASES contact sites (MAGeCK ``neg``), 229
- ``_merged`` -- the union, 623

``_inc`` and ``_dec`` are disjoint (overlap 0) and are the two independent
screens; ``_merged`` is a derived view of the same data, so it must not be
counted as a third independent unit in any interval.

**Hit threshold: FDR < 1e-4, which is the paper's stated "FDR = 0.01%".** That
gives 394 up / 229 down against the 410 / 230 the paper reports. The gap is not
a disagreement about the cut: MAGeCK RRA's permutation FDRs are heavily TIED
(357 genes share the single value 1.4e-05), so no threshold yields 410 at all --
the achievable counts step 394 -> 412 as the cut passes 1.08e-4. 230 is
reachable, at FDR <= 1.29e-4. The stated threshold is used as stated, the
residual is under 4%, and ``--hit-fdr`` moves it.

**Symbols are resolved from Ensembl IDs first.** This file carries an ENSG per
row, so the authoritative route is ENSG -> current HGNC symbol, with the
old-symbol remap only as a fallback. That also repairs the four duplicated
symbols, two of which are Excel date corruption (``Mar-01``, ``Mar-02`` for
MARCHF genes) and would otherwise be unmappable; the other two (RABGEF1,
STPG4) carry two Ensembl IDs each. Any symbol still duplicated after
resolution keeps its most significant row, because a duplicated gene list makes
``reveal()`` collapse entries and corrupt the hit labels.
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import pickle
from pathlib import Path

import pandas as pd

from assaybench.data.screens import ScreenRecord
from build_schukken import hgnc_synonyms, remap_library

log = logging.getLogger("build_wilson")

HIT_FDR = 1e-4

COLS = ["ensg", "gene", "num",
        "neg_score", "neg_p", "neg_fdr", "neg_rank", "neg_good", "neg_lfc",
        "pos_score", "pos_p", "pos_fdr", "pos_rank", "pos_good", "pos_lfc"]

ARMS = {
    # suffix: (fdr column, lfc column, phenotype, direction)
    "inc": ("pos_fdr", "pos_lfc",
            "Genes whose knockout increases mitochondria-endoplasmic reticulum "
            "contact sites, enriched in the high-reporter-fluorescence sorted "
            "population", "increase"),
    "dec": ("neg_fdr", "neg_lfc",
            "Genes whose knockout decreases mitochondria-endoplasmic reticulum "
            "contact sites, enriched in the low-reporter-fluorescence sorted "
            "population", "decrease"),
}

CONDITION = (
    "Genome-wide CRISPR-Cas9 knockout screen in HeLa-Cas9 cells carrying a "
    "doxycycline-inducible split-mVenus mitochondria-endoplasmic reticulum "
    "contact site reporter, sorted by median reporter fluorescence intensity "
    "and compared high versus low population"
)


def ensg_to_symbol(path: Path) -> dict[str, str]:
    """``{ensembl_gene_id: approved symbol}`` from the HGNC complete set."""
    out: dict[str, str] = {}
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            e, s, st = (r.get("ensembl_gene_id") or "").strip(), r["symbol"].strip(), r["status"].strip()
            if e and s and st == "Approved":
                out[e] = s
    return out


def build(csv_path: Path, hgnc: Path, vocab: Path, hit_fdr: float) -> list[ScreenRecord]:
    d = pd.read_csv(csv_path, sep=r"\s+", engine="python", skiprows=1, names=COLS)
    for c in COLS[2:]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d["pos_fdr"] = d["pos_fdr"].fillna(1.0)
    d["neg_fdr"] = d["neg_fdr"].fillna(1.0)

    by_ensg = ensg_to_symbol(hgnc)
    syn = hgnc_synonyms(hgnc)
    known = set(json.loads(vocab.read_text())["itos"])

    # Ensembl id is authoritative; the given symbol is the fallback.
    n_by_ensg = 0
    sym = []
    for e, g in zip(d["ensg"], d["gene"]):
        s = by_ensg.get(str(e).strip())
        if s:
            n_by_ensg += 1
            sym.append(s)
        else:
            sym.append(str(g).strip())
    d["sym"] = sym
    log.info("resolved %d/%d symbols from Ensembl ids", n_by_ensg, len(d))

    # Drop duplicates, keeping each symbol's most significant row in either
    # direction, so the gene list is unique and the hit labels stay intact.
    d["best"] = d[["pos_fdr", "neg_fdr"]].min(axis=1)
    before = len(d)
    d = d.sort_values("best").drop_duplicates("sym", keep="first").reset_index(drop=True)
    if before != len(d):
        log.info("dropped %d duplicate-symbol rows (kept most significant)", before - len(d))

    genes = remap_library([str(g) for g in d["sym"]], known, syn)
    n_oov = sum(1 for g in genes if g not in known)
    log.info("library %d genes, %d unique; %d (%.1f%%) outside the model vocabulary",
             len(genes), len(set(genes)), n_oov, 100 * n_oov / len(genes))
    if len(set(genes)) != len(genes):
        raise RuntimeError(f"gene list not unique ({len(set(genes))} of {len(genes)})")

    masks = {}
    out: list[ScreenRecord] = []
    for suffix, (fdr_col, lfc_col, phenotype, direction) in ARMS.items():
        masks[suffix] = [bool(v < hit_fdr) for v in d[fdr_col]]
    masks["merged"] = [a or b for a, b in zip(masks["inc"], masks["dec"])]
    overlap = sum(1 for a, b in zip(masks["inc"], masks["dec"]) if a and b)
    log.info("inc %d, dec %d, overlap %d, merged %d",
             sum(masks["inc"]), sum(masks["dec"]), overlap, sum(masks["merged"]))

    for suffix in ("inc", "dec", "merged"):
        hits = masks[suffix]
        if suffix == "merged":
            phenotype = ("Genes whose knockout changes mitochondria-endoplasmic "
                         "reticulum contact sites in either direction")
            direction = "both"
            # magnitude of whichever direction called the gene
            mag = [max(abs(p), abs(n)) for p, n in zip(d["pos_lfc"].fillna(0.0),
                                                       d["neg_lfc"].fillna(0.0))]
        else:
            _, lfc_col, phenotype, direction = ARMS[suffix]
            mag = [abs(v) for v in d[lfc_col].fillna(0.0)]

        out.append(ScreenRecord(
            dataset_name=f"WILSON_MERCS_{suffix}",
            split="external",
            organism="Homo sapiens",
            gene_symbol_convention="HGNC",
            phenotype=phenotype,
            cleaned_phenotype="Molecular Output / Reporter / Pathway Activity",
            cell_line="HeLa",
            cell_type="cervical adenocarcinoma, split-mVenus MERCS reporter",
            library_type="whole-genome CRISPR knockout library (91,536 gRNAs)",
            library_methodology="Knockout",
            direction_str=direction,
            condition_clause=CONDITION,
            num_genes=len(genes),
            genes=list(genes),
            relevance_scores=[float(m) if h else 0.0 for m, h in zip(mag, hits)],
            hits=hits,
            description=CONDITION,
            contrast_label="high versus low reporter fluorescence",
            notes=("Wilson et al. 2024 (Cell Death & Disease 15:203, PMC10928201). "
                   f"Hits are the authors' MAGeCK RRA calls at FDR < {hit_fdr:g} "
                   "in the stated direction."),
            unmasked_relevance_scores=[float(m) for m in mag],
        ))
    return out


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True, type=Path)
    ap.add_argument("--hgnc", required=True, type=Path)
    ap.add_argument("--vocab", required=True, type=Path)
    ap.add_argument("--hit-fdr", type=float, default=HIT_FDR)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()

    screens = build(a.csv, a.hgnc, a.vocab, a.hit_fdr)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "wb") as fh:
        pickle.dump(screens, fh)
    print(f"\n{len(screens)} screens -> {a.out}")
    for s in screens:
        print(f"  {s.dataset_name:24s} {s.num_genes:6d} genes {s.total_hits:5d} hits "
              f"frac {s.total_hits/s.num_genes:.4f}  EF ceiling "
              f"{s.num_genes/max(s.total_hits,1000):5.2f}")


if __name__ == "__main__":
    main()
