"""Turn the Sanger organoid biobank dependency screens into ScreenRecords.

A sixth EXTERNAL test, and by far the largest: 178 genome-wide knockout
dependency screens over 162 patient-derived tumour organoids, against 1-14
screens for every earlier external dataset. Enough to test the method on
external data with something close to the benchmark's own power.

Source: Herranz-Ors, Bhosle, Beck, Gilbert, Picco ... Garnett, "A tumour-derived
organoid biobank maps cancer gene dependencies", Nature (2026), doi
10.1038/s41586-026-10830-y, PMC13581617. Data: Supplementary Table 5
(``sample_ID gene ensembl_id HGNC_ID library is_depleted``, 2 948 036 rows),
with tissue labels from the Figshare item's ``data_covariates_ALL-organoids.csv``.

Libraries: MinLibCas9 (primary) and Human CRISPR Library Yusa v1.1 for a subset.
CRISPRcleanR-corrected LFCs, fitness genes called by BAGEL2.

**A TRAP WORTH RECORDING.** The obvious-looking file on the paper's Figshare
item, ``07_EssMatrix_bDepletionsB2.Rdata``, is NOT this paper's organoids. It is
17 647 genes x 930 columns whose names are ``ACH-001421``, ``ACH-000888`` and so
on -- DepMap/Achilles CELL LINE identifiers, the cell-line reference matrix from
Pacini et al. that the paper compares against. Building screens from it would
have scored the DepMap co-essentiality atlas against the very cell lines the
atlas is computed from: total leakage, and it would have manufactured a large
fake gain. The organoid calls live only in Supplementary Table 5, whose
``sample_ID`` values are ``HCM-SANG-*``.

**The unit is (sample_ID, library), giving 178 screens.** 16 of the 162
organoids were screened with both libraries, so grouping by sample alone
double-counts their genes (33 124 rows where a single-library model has 16 562)
and inflates their hit counts past the paper's stated range. Split by library,
hits per screen run 297 to 2160, which reproduces the paper's quoted range
exactly. The authors do pick one library per model for their own figures, using
the composite ``order`` metric in Supplementary Table 4, but that metric
disagrees with AUC_ROC on some models and its direction is not documented in
the sheet, so rather than guess at their rule this build keeps both screens and
leaves the duplication to the bootstrap: the resampling cluster is the ORGANOID
(162 clusters over 178 screens), so a dual-library model cannot count twice.

**Hits are the authors' own BAGEL2 calls** (``is_depleted``), not re-derived.

**Leakage.** Nature, August 2026, so it post-dates every AssayBench year fold.
These are organoids, not Achilles cell lines, so they are absent from the
DepMap atlas the method reads -- but note that this is a FITNESS/dependency
dataset, the one phenotype where a co-essentiality atlas is measuring the same
kind of quantity as the screen itself.
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

log = logging.getLogger("build_organoid")

TISSUE = {
    "COLO": ("colorectal", "colorectal carcinoma"),
    "OESO": ("oesophageal", "oesophageal carcinoma"),
    "ESCA": ("oesophageal", "oesophageal carcinoma"),
    "OV":   ("ovarian", "ovarian carcinoma"),
    "PAAD": ("pancreatic", "pancreatic adenocarcinoma"),
    "STAD": ("gastric", "gastric adenocarcinoma"),
    "GAST": ("gastric", "gastric adenocarcinoma"),
}
LIBNAME = {"minLib": "MinLibCas9 minimal genome-wide CRISPR knockout library",
           "v1_1": "Human CRISPR Library Yusa v1.1"}


def build(table5: Path, cov: Path, hgnc: Path, vocab: Path,
          limit: int = 0) -> list[ScreenRecord]:
    syn = hgnc_synonyms(hgnc)
    known = set(json.loads(vocab.read_text())["itos"])

    tis = {}
    c = pd.read_csv(cov)
    for sid, t in zip(c["sample_ID"].astype(str), c["tissue"].astype(str)):
        tis[sid.strip()] = t.strip()

    log.info("reading %s ...", table5.name)
    d = pd.read_csv(table5, sep=r"\s+", engine="c",
                    usecols=["sample_ID", "gene", "library", "is_depleted"])
    log.info("  %d rows, %d samples, %d sample x library groups",
             len(d), d["sample_ID"].nunique(),
             d.groupby(["sample_ID", "library"]).ngroups)

    out: list[ScreenRecord] = []
    for (sid, lib), grp in d.groupby(["sample_ID", "library"], sort=True):
        raw = [str(g).strip() for g in grp["gene"]]
        dep = [int(v) == 1 for v in grp["is_depleted"]]
        # one row per gene within a screen; drop any repeat, keeping a call
        seen: dict[str, bool] = {}
        for g, h in zip(raw, dep):
            seen[g] = seen.get(g, False) or h
        gene_list = list(seen)
        hits_by_raw = seen

        genes = remap_library(gene_list, known, syn)
        if len(set(genes)) != len(genes):
            raise RuntimeError(f"{sid}/{lib}: gene list not unique")
        hits = [hits_by_raw[g] for g in gene_list]

        tcode = tis.get(str(sid).strip(), "")
        short, ctype = TISSUE.get(tcode, ("tumour", "patient-derived tumour organoid"))
        condition = (
            f"Genome-wide CRISPR-Cas9 knockout dependency screen in {sid}, a "
            f"patient-derived {short} cancer organoid, using the "
            f"{LIBNAME.get(lib, lib)}; sgRNA depletion between the early "
            "time point and the end of screening, CRISPRcleanR-corrected and "
            "called with BAGEL2"
        )
        out.append(ScreenRecord(
            dataset_name=f"ORGANOID_{sid}_{lib}",
            split="external",
            organism="Homo sapiens",
            gene_symbol_convention="HGNC",
            phenotype=("Genes whose knockout reduces organoid fitness and is "
                       "depleted over the course of the screen (BAGEL2 fitness "
                       "genes)"),
            cleaned_phenotype="Fitness / Proliferation / Viability",
            cell_line=str(sid),
            cell_type=f"{ctype}, patient-derived organoid",
            library_type=LIBNAME.get(lib, lib),
            library_methodology="Knockout",
            direction_str="decrease",
            condition_clause=condition,
            num_genes=len(genes),
            genes=list(genes),
            # BAGEL2 gives a binary call here, so a hit's relevance is flat; the
            # score-aware metrics are not used in this comparison.
            relevance_scores=[1.0 if h else 0.0 for h in hits],
            hits=hits,
            description=condition,
            contrast_label="end of screen versus early time point",
            notes=("Herranz-Ors et al. 2026 (Nature, PMC13581617), Sanger "
                   "organoid biobank. Hits are the authors' own BAGEL2 "
                   f"fitness-gene calls. Library {lib}."),
            unmasked_relevance_scores=[1.0 if h else 0.0 for h in hits],
        ))
        if limit and len(out) >= limit:
            break
    return out


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    ap = argparse.ArgumentParser()
    ap.add_argument("--table5", required=True, type=Path)
    ap.add_argument("--cov", required=True, type=Path)
    ap.add_argument("--hgnc", required=True, type=Path)
    ap.add_argument("--vocab", required=True, type=Path)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()

    screens = build(a.table5, a.cov, a.hgnc, a.vocab, a.limit)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "wb") as fh:
        pickle.dump(screens, fh)

    import numpy as np
    H = np.array([s.total_hits for s in screens])
    G = np.array([s.num_genes for s in screens])
    print(f"\n{len(screens)} screens -> {a.out}")
    print(f"  genes per screen : min {G.min()}  median {int(np.median(G))}  max {G.max()}")
    print(f"  hits per screen  : min {H.min()}  median {int(np.median(H))}  max {H.max()}"
          f"   (paper: median 1440, range 297-2160)")
    print(f"  hit fraction     : min {(H/G).min():.4f}  median {np.median(H/G):.4f}  "
          f"max {(H/G).max():.4f}")
    print(f"  sparse (H<500)   : {int((H < 500).sum())} screens")


if __name__ == "__main__":
    main()
