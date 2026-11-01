"""Turn the Nguyen et al. 2025 ADE chemotherapy screens into ScreenRecords.

A fourth EXTERNAL test, and the first DRUG-RESPONSE one. Schukken and Abt are
knockout fitness screens; Wilson is a sorted reporter screen; this is
"Drug / Chemical / Environmental Response", the other phenotype cut where the
DepMap gain was measured strongest (+0.452 [+0.228, +0.707]).

Source: Nguyen, Rafiee, Parcha, Tagmount, Rubnitz, Ribeiro, Cao, Pounds, Vulpe,
Lamba, "Genome-wide CRISPR/Cas9 screen identifies AraC-daunorubicin-etoposide
response modulators associated with outcomes in pediatric AML", Blood Advances
9(5):1078 (2025), doi 10.1182/bloodadvances.2024014157, PMC11914169. Data:
supplemental tables ``mmc2.xlsx``, sheets S1 (genome-wide diff-beta) and S2
(the significant calls).

Brunello knockout library, 19 114 genes, K-562, three drugs each against an
untreated control, day 0 to day 4, MAGeCK-MLE beta scores. ``diff-beta`` is
treatment beta minus control beta.

**Nine screens from three drugs.** Each drug gives two directional screens plus
their union, following AssayBench's own ``_inc``/``_dec``/``_merged`` convention
for two-tailed screens. The paper's own words fix the labels, and they are kept
rather than reinterpreted: positive diff-beta is "sensitive" (positive
selection), negative is "resistant" (negative selection).

    AraC          resistant  240   sensitive 1800
    Daunorubicin  resistant 1229   sensitive  258
    Etoposide     resistant  296   sensitive 1410

The six directional screens are the units; the three ``_merged`` ones are
derived unions of the same data and must not be counted as independent. The
natural resampling cluster is the DRUG (n=3), since all nine come from one cell
line, one library and one experiment, and the two tails of a drug come from a
single diff-beta vector, so they are anti-correlated rather than independent.

**Hits are the authors' own calls, taken from sheet S2, not re-derived.** Their
stated rule is "diff-beta beyond 2 SD", but no plain reading of it reproduces
their counts from sheet S1: mean +/- 2 SD of each drug's diff-beta gives
8/976, 731/221 and 120/771 against the published 240/1800, 1229/258 and
296/1410, and centring on zero instead gives 4/1066, 899/165 and 81/889. Their
threshold is evidently computed some other way (MAGeCK-Flute offers several).
Rather than guess at it, this build reads the per-drug ``Response`` labels out
of sheet S2, which reproduces all six published counts EXACTLY. Sheet S1
supplies the 19 114-gene library and the diff-beta magnitudes.

**Leakage check.** 2025, inside AssayBench's year range, so it was checked: the
AraC screens already in the benchmark are hTERT-RPE1 (test fold, U_2334-U_2342)
and the one etoposide screen is NCI-H69AR; none of the 93 K-562 screens in the
train fold mentions any ADE drug, and no fold holds a screen authored by Lamba.
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

log = logging.getLogger("build_nguyen")

DRUGS = {
    # key: (full name, condition phrase)
    "AraC":  ("cytarabine (AraC)", "cytarabine"),
    "Dauno": ("daunorubicin", "daunorubicin"),
    "Etop":  ("etoposide", "etoposide"),
}

# The paper's own direction labels, and what each selects for.
DIRECTIONS = {
    "sensitive": ("Sensitive", "increase",
                  "whose knockout increases sensitivity to {drug} (positive "
                  "selection, sgRNA enrichment)"),
    "resistant": ("Resistant", "decrease",
                  "whose knockout confers resistance to {drug} (negative "
                  "selection, sgRNA depletion)"),
}


def build(xlsx: Path, hgnc: Path, vocab: Path) -> list[ScreenRecord]:
    lib = pd.read_excel(xlsx, sheet_name="Supplementary Table S1", skiprows=1)
    lib.columns = [str(c).strip() for c in lib.columns]
    sig = pd.read_excel(xlsx, sheet_name="Supplementary Table S2", header=1)
    sig.columns = [str(c).strip() for c in sig.columns]

    # Library: dedupe on gene, keeping the row with the largest |diff-beta| in
    # any drug, so a duplicated symbol cannot collapse in reveal().
    dcols = [f"{d} Diff-Beta" for d in DRUGS]
    for c in dcols:
        lib[c] = pd.to_numeric(lib[c], errors="coerce")
    lib["Gene"] = lib["Gene"].astype(str).str.strip()
    lib["_mag"] = lib[dcols].abs().max(axis=1).fillna(0.0)
    before = len(lib)
    lib = (lib.sort_values("_mag", ascending=False)
              .drop_duplicates("Gene", keep="first")
              .reset_index(drop=True))
    if before != len(lib):
        log.info("dropped %d duplicate-gene rows from the library", before - len(lib))

    syn = hgnc_synonyms(hgnc)
    known = set(json.loads(vocab.read_text())["itos"])
    genes = remap_library(list(lib["Gene"]), known, syn)
    n_oov = sum(1 for g in genes if g not in known)
    log.info("library %d genes, %d unique; %d (%.1f%%) outside the model vocabulary",
             len(genes), len(set(genes)), n_oov, 100 * n_oov / len(genes))
    if len(set(genes)) != len(genes):
        raise RuntimeError(f"gene list not unique ({len(set(genes))} of {len(genes)})")

    # Map the authors' calls onto the deduped library order, by ORIGINAL symbol
    # (both tables use the same spellings; the remap happens after).
    pos_of = {g: i for i, g in enumerate(lib["Gene"])}
    sig["Gene"] = sig["Gene"].astype(str).str.strip()

    # Sheet S2 went through an R step that mangled '-' to '__' in composite and
    # hyphenated symbols -- NKX3__1 for NKX3-1, HLA__C for HLA-C,
    # ATP5J2__PTCD1 for ATP5J2-PTCD1. Sheet S1 kept the hyphens, so 34 of the
    # authors' 3748 calls would otherwise fail to find their library row and
    # silently drop out of the hit sets. Undoing the mangling recovers all 34
    # and makes the per-arm counts reproduce the paper exactly.
    def unmangle(g: str) -> str:
        return g if g in pos_of else g.replace("__", "-")

    sig["Gene"] = [unmangle(g) for g in sig["Gene"]]

    out: list[ScreenRecord] = []
    for dkey, (dname, dphrase) in DRUGS.items():
        resp_col, beta_col = f"{dkey} Response", f"{dkey} Diff-Beta"
        masks: dict[str, list[bool]] = {}
        for dirkey, (label, _, _) in DIRECTIONS.items():
            m = [False] * len(genes)
            called = sig.loc[sig[resp_col] == label, "Gene"]
            missing = 0
            for g in called:
                i = pos_of.get(g)
                if i is None:
                    missing += 1
                else:
                    m[i] = True
            if missing:
                log.warning("%s/%s: %d called genes absent from the library",
                            dkey, dirkey, missing)
            masks[dirkey] = m
            log.info("%s %s: %d hits", dkey, dirkey, sum(m))
        masks["merged"] = [a or b for a, b in zip(masks["sensitive"], masks["resistant"])]

        mag = lib[beta_col].abs().fillna(0.0).tolist()
        for dirkey in ("sensitive", "resistant", "merged"):
            hits = masks[dirkey]
            if dirkey == "merged":
                phen = (f"Genes whose knockout changes {dname} response in "
                        "either direction")
                direction = "both"
            else:
                _, direction, tmpl = DIRECTIONS[dirkey]
                phen = "Genes " + tmpl.format(drug=dname)
            condition = (
                f"Genome-wide CRISPR-Cas9 knockout screen in K-562 myeloid "
                f"leukaemia cells treated with {dphrase} for four days and "
                f"compared with an untreated control, Brunello library, "
                f"MAGeCK-MLE differential beta score"
            )
            out.append(ScreenRecord(
                dataset_name=f"NGUYEN_{dkey}_{dirkey}",
                split="external",
                organism="Homo sapiens",
                gene_symbol_convention="HGNC",
                phenotype=phen,
                cleaned_phenotype="Drug / Chemical / Environmental Response",
                cell_line="K-562",
                cell_type="chronic myeloid leukaemia, AML model",
                library_type="Brunello whole-genome CRISPR knockout library",
                library_methodology="Knockout",
                direction_str=direction,
                condition_clause=condition,
                condition_name=dphrase,
                num_genes=len(genes),
                genes=list(genes),
                relevance_scores=[float(m) if h else 0.0 for m, h in zip(mag, hits)],
                hits=hits,
                description=condition,
                contrast_label=f"{dphrase} versus untreated control, day 4",
                notes=("Nguyen et al. 2025 (Blood Adv 9:1078, PMC11914169). Hits "
                       "are the authors' own 2 SD diff-beta calls, read from "
                       "supplemental Table S2."),
                duration="four days",
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
        print(f"  {s.dataset_name:26s} {s.num_genes:6d} genes {s.total_hits:5d} hits "
              f"frac {s.total_hits/s.num_genes:.4f}  EF ceiling "
              f"{s.num_genes/max(s.total_hits,1000):5.2f}")


if __name__ == "__main__":
    main()
