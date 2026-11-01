"""Turn the Shan et al. 2025 cholesterol-homeostasis screens into ScreenRecords.

A fifth EXTERNAL test. Three genome-scale knockout screens for synthetic
lethality with the cholesterol-synthesis axis, all in the SPARSE regime
(1.2-1.9% of the library are hits), which is where both the benchmark
stratification (+0.391 for H<500) and the external screens (+0.595) say
AF + DepMap should gain.

Source: Shan, Fan, Li, Liang, Chen, Wang, Wang, Li, Chen, Yu, Fei, "Systematic
interrogation of functional genes underlying cholesterol and lipid
homeostasis", Genome Biology 26 (2025), doi 10.1186/s13059-025-03531-8,
PMC11912599. Data: Additional file 2 (Table S1, "CRISPR screen data").

Human H3 knockout library (Addgene 133914), 117 606 sgRNAs over >18 000
protein-coding genes, MAGeCK beta scores. Three screens:

    1  HepG2, lovastatin vs vehicle                          368 hits
    2  HeLa,  lovastatin vs vehicle                          359 hits
    3  HeLa,  LDLR-ablated clones vs AAVS1-KO control        227 hits

**Hits are the authors' own, read from their ``_hits`` sheets, not re-derived.**
Their rules are stated in the workbook's own legend -- for screens 1 and 2,
``beta(lovastatin) - beta(vehicle) <= -0.2 and beta(lovastatin) <= -0.3``; for
screen 3, ``score(LDLR KO clone) - score(AAVS1 KO) <= -2 and score(LDLR KO
clone) <= -1``. Screen 1's 368 reproduces the count quoted in the paper's text.
Taking the sheets directly means the thresholds need no reconstruction, which
matters here because screen 3's rule references a "tRRA" score whose column is
not in the workbook under that name.

**Leakage: this study is not in AssayBench**, but unlike the earlier external
sets its BIOLOGY is partly represented in the train fold -- statin-response
screens in other lines (``U_1182`` K-562, ``U_1321`` hTERT-RPE1, ``U_1503``
HAP-1) and a cholesterol/LDLR screen (``1434``, SV589J). No HepG2 or HeLa
lovastatin screen, and nothing authored by this group, so the screens
themselves are held out; the model has nonetheless seen the pathway. That makes
this a weaker novelty test than Wilson's reporter screen and it is reported as
such.

Screens 2 and 3 share the HeLa background, so the resampling cluster is the
cell line (n=2), not the screen (n=3).
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

log = logging.getLogger("build_shan")

# sheet stem -> (name, cell line, cell type, condition, contrast, phenotype tail)
SCREENS = {
    "Screen 1": dict(
        name="HepG2_lovastatin", cell_line="Hep G2",
        cell_type="hepatocellular carcinoma",
        condition=("Genome-scale CRISPR-Cas9 knockout screen in Hep G2 "
                   "hepatocellular carcinoma cells treated with lovastatin, an "
                   "HMG-CoA reductase inhibitor, and compared with vehicle"),
        contrast="lovastatin versus vehicle",
        phen=("Genes whose knockout is synthetically lethal with lovastatin, "
              "depleting the cells that carry it"),
        condition_name="lovastatin",
        diff="Lovastatin - Vehicle"),
    "Screen 2": dict(
        name="HeLa_lovastatin", cell_line="HeLa",
        cell_type="cervical adenocarcinoma",
        condition=("Genome-scale CRISPR-Cas9 knockout screen in HeLa cells "
                   "treated with lovastatin, an HMG-CoA reductase inhibitor, "
                   "and compared with vehicle"),
        contrast="lovastatin versus vehicle",
        phen=("Genes whose knockout is synthetically lethal with lovastatin, "
              "depleting the cells that carry it"),
        condition_name="lovastatin",
        diff="Lovastatin - Vehicle"),
    "Screen 3": dict(
        name="HeLa_LDLRko", cell_line="HeLa",
        cell_type="cervical adenocarcinoma, LDLR-ablated clones",
        condition=("Genome-scale CRISPR-Cas9 knockout screen in LDLR-ablated "
                   "HeLa single-cell clones compared with an AAVS1-knockout "
                   "control, probing synthetic lethality with loss of the LDL "
                   "receptor"),
        contrast="LDLR-knockout clones versus AAVS1-knockout control",
        phen=("Genes whose knockout is synthetically lethal with loss of the "
              "LDL receptor, depleting the cells that carry it"),
        condition_name="LDLR ablation",
        # two clones; magnitude takes whichever contrast is stronger
        diff=["LDLR KO clone1 - AAVS1 KO", "LDLR KO clone2 - AAVS1 KO"]),
}


def build(xlsx: Path, hgnc: Path, vocab: Path) -> list[ScreenRecord]:
    syn = hgnc_synonyms(hgnc)
    known = set(json.loads(vocab.read_text())["itos"])
    out: list[ScreenRecord] = []

    for stem, spec in SCREENS.items():
        alls = pd.read_excel(xlsx, sheet_name=f"{stem}_all")
        hits_sheet = pd.read_excel(xlsx, sheet_name=f"{stem}_hits")
        gcol = alls.columns[0]
        alls[gcol] = alls[gcol].astype(str).str.strip()
        hit_names = set(hits_sheet[hits_sheet.columns[0]].astype(str).str.strip())

        before = len(alls)
        alls = alls.drop_duplicates(gcol, keep="first").reset_index(drop=True)
        if before != len(alls):
            log.info("%s: dropped %d duplicate-gene rows", stem, before - len(alls))

        genes = remap_library(list(alls[gcol]), known, syn)
        n_oov = sum(1 for g in genes if g not in known)
        if len(set(genes)) != len(genes):
            raise RuntimeError(f"{stem}: gene list not unique")

        # Hit membership is decided on the ORIGINAL symbols, because both sheets
        # of the workbook use the same spellings; the HGNC remap applies after.
        hits = [g in hit_names for g in alls[gcol]]
        missing = len(hit_names) - sum(hits)
        if missing:
            log.warning("%s: %d of %d called genes not found in the _all sheet",
                        stem, missing, len(hit_names))

        diff = spec["diff"]
        cols = [diff] if isinstance(diff, str) else diff
        mag = (pd.concat([pd.to_numeric(alls[c], errors="coerce").abs()
                          for c in cols], axis=1).max(axis=1).fillna(0.0))

        log.info("%s (%s): %d genes, %d hits (frac %.4f), %d (%.1f%%) OOV",
                 stem, spec["name"], len(genes), sum(hits),
                 sum(hits) / len(genes), n_oov, 100 * n_oov / len(genes))

        out.append(ScreenRecord(
            dataset_name=f"SHAN_{spec['name']}",
            split="external",
            organism="Homo sapiens",
            gene_symbol_convention="HGNC",
            phenotype=spec["phen"],
            cleaned_phenotype="Drug / Chemical / Environmental Response",
            cell_line=spec["cell_line"],
            cell_type=spec["cell_type"],
            library_type="Human H3 whole-genome CRISPR knockout library",
            library_methodology="Knockout",
            direction_str="decrease",
            condition_clause=spec["condition"],
            condition_name=spec["condition_name"],
            num_genes=len(genes),
            genes=list(genes),
            relevance_scores=[float(m) if h else 0.0 for m, h in zip(mag, hits)],
            hits=hits,
            description=spec["condition"],
            contrast_label=spec["contrast"],
            notes=("Shan et al. 2025 (Genome Biol 26, PMC11912599). Hits are "
                   "the authors' own synthetic-lethal calls, taken from the "
                   f"workbook's {stem}_hits sheet."),
            unmasked_relevance_scores=[float(m) for m in mag],
        ))
    return out


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    ap = argparse.ArgumentParser()
    ap.add_argument("--xlsx", required=True, type=Path)
    ap.add_argument("--hgnc", required=True, type=Path)
    ap.add_argument("--vocab", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()

    screens = build(a.xlsx, a.hgnc, a.vocab)
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
