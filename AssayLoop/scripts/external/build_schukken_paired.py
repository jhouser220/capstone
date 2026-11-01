"""The Schukken PAIRED aneuploid-vs-near-euploid contrast -- the paper's own analysis.

`build_schukken.py` reads each of the 14 WGL columns as what it is on its own: a
knockout dropout screen, hit = `beta < -0.75`. That is NOT the paper's analysis.
The paper's result is the PAIRED contrast: for each gene, the mean difference in
beta score between an aneuploid clone and its own isogenic near-euploid parent,
over all seven pairs. One screen, not fourteen, and it asks a different question
-- aneuploid-SPECIFIC dependency rather than essentiality.

Source: Schukken, Aitken, Sheltzer et al., DOI 10.64898/2026.04.26.720636
(PMC13142536), `SupData/SD3- CRISPR Screen Beta Scores.xlsx`, sheet **`WGL`**
(quantile-normalised, 19 112 genes) -- not `Raw WGL` (20 111), which has no
paired columns. Columns used: `ANvEU.Beta.Paired` (aneuploid minus euploid),
`ANvEU.P.Paired` (paired t-test), `AneuploidRank` (the authors' own ranking;
"low values are aneuploid dependencies", and it is exactly the ascending rank of
`ANvEU.Beta.Paired`, Spearman 1.000).

**Correction to `build_schukken.py`'s docstring.** It says the paired cut
"leaves 15 hits genome-wide -- far too few to measure EF over 1000 picks". That
number was computed on the `DG + WGL` sheet, which is 3 377 druggable-genome
genes, not genome-wide. On the genome-wide `WGL` sheet the paper's own rule
leaves **188**, which is a perfectly usable hit count. The paired contrast was
viable all along.

**Two hit rules are built, because neither can be assumed.**

  `paper`   `ANvEU.P.Paired < 0.05 and ANvEU.Beta.Paired < -0.2`  -> 188 hits.
            The rule as the paper states it.
  `rank194` the top 194 genes by `AneuploidRank`  -> 194 hits, equivalently
            `ANvEU.Beta.Paired < -0.2725`.

194 is the count in the team's `external_sources.json` for
`E224_SCHUKKEN_PAIRED_WGL`, and no round threshold reproduces it: not p < 0.05
with any 2-dp beta cut (188, or 197 if beta is rounded first), not BH-FDR at any
level (0 hits at q < 0.3 -- the paired t-test over seven pairs has no gene
surviving multiplicity correction), not a z-score cut on beta (151 at -3 SD, 296
at -2.5 SD), and not p < 0.01 (206). The `rank194` arm matches their count
exactly using the authors' own ranking column; the `paper` arm matches the
authors' stated rule. They overlap heavily, so if the two arms agree the
ambiguity does not matter.

**Leakage: none.** Posted April 2026, after every AssayBench year fold.
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

log = logging.getLogger("build_schukken_paired")

P_CUT = 0.05
BETA_CUT = -0.2
RANK_N = 194

CONDITION = (
    "Paired whole-genome CRISPR-Cas9 knockout screen across seven isogenic "
    "pairs of engineered aneuploid clones and their own near-euploid parental "
    "controls (DLD-1, HCT-116, Vaco432, A2780, A2058, AGS, MCF10A), comparing "
    "the MAGeCK beta score of each gene in the aneuploid clones against the "
    "matched near-euploid controls"
)
PHENOTYPE = (
    "Genes whose knockout reduces fitness specifically in aneuploid cells, "
    "more depleted in the aneuploid clones than in their matched near-euploid "
    "controls (aneuploid-specific dependencies)"
)


def build(xlsx: Path, hgnc: Path, vocab: Path) -> list[ScreenRecord]:
    d = pd.read_excel(xlsx, sheet_name="WGL")
    d["Gene"] = d["Gene"].astype(str).str.strip()
    for c in ("ANvEU.Beta.Paired", "ANvEU.P.Paired", "AneuploidRank"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    log.info("WGL sheet: %d genes, %d unique", len(d), d["Gene"].nunique())

    before = len(d)
    d = (d.sort_values("ANvEU.Beta.Paired")
          .drop_duplicates("Gene", keep="first")
          .reset_index(drop=True))
    if before != len(d):
        log.info("dropped %d duplicate-gene rows (kept the strongest)", before - len(d))

    syn = hgnc_synonyms(hgnc)
    known = set(json.loads(vocab.read_text())["itos"])
    genes = remap_library(list(d["Gene"]), known, syn)
    n_oov = sum(1 for g in genes if g not in known)
    log.info("library %d genes, %d unique; %d (%.1f%%) outside the model vocabulary",
             len(genes), len(set(genes)), n_oov, 100 * n_oov / len(genes))
    if len(set(genes)) != len(genes):
        raise RuntimeError(f"gene list not unique ({len(set(genes))} of {len(genes)})")

    beta = d["ANvEU.Beta.Paired"].fillna(0.0)
    pval = d["ANvEU.P.Paired"].fillna(1.0)
    rank = d["AneuploidRank"]

    rules = {
        "paper": [bool(p < P_CUT and b < BETA_CUT) for p, b in zip(pval, beta)],
        "rank194": [bool(r <= RANK_N) for r in rank],
    }
    agree = sum(1 for a, b in zip(rules["paper"], rules["rank194"]) if a and b)
    log.info("paper rule %d hits, rank194 %d hits, %d in both",
             sum(rules["paper"]), sum(rules["rank194"]), agree)

    out: list[ScreenRecord] = []
    for name, hits in rules.items():
        rule_text = (f"paired t-test p < {P_CUT} and aneuploid-minus-euploid "
                     f"beta < {BETA_CUT} (the paper's stated rule)"
                     if name == "paper" else
                     f"the top {RANK_N} genes by the authors' AneuploidRank")
        out.append(ScreenRecord(
            dataset_name=f"SCHUKKEN_PAIRED_WGL_{name}",
            split="external",
            organism="Homo sapiens",
            gene_symbol_convention="HGNC",
            phenotype=PHENOTYPE,
            cleaned_phenotype="Fitness / Proliferation / Viability",
            cell_line="isogenic aneuploid / near-euploid clone panel",
            cell_type=("engineered aneuploid human cancer and immortalised "
                       "epithelial clones with matched near-euploid parents"),
            library_type="Brunello whole-genome CRISPR knockout library",
            library_methodology="Knockout",
            direction_str="decrease",
            condition_clause=CONDITION,
            num_genes=len(genes),
            genes=list(genes),
            relevance_scores=[float(-b) if h else 0.0 for b, h in zip(beta, hits)],
            hits=hits,
            description=CONDITION,
            contrast_label="aneuploid clones versus matched near-euploid controls",
            notes=("Schukken et al. 2026 (PMC13142536), the paper's PAIRED "
                   f"aneuploid-vs-near-euploid contrast. Hits are {rule_text}. "
                   "Quantile-normalised WGL sheet."),
            duration="ten population doublings",
            unmasked_relevance_scores=[float(-b) for b in beta],
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
        print(f"  {s.dataset_name:30s} {s.num_genes:6d} genes {s.total_hits:5d} hits "
              f"frac {s.total_hits / s.num_genes:.4f}  EF ceiling "
              f"{s.num_genes / max(s.total_hits, 1000):5.2f}")


if __name__ == "__main__":
    main()
