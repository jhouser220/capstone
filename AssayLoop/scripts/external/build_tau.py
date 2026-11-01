"""Turn the Samelson et al. 2026 tau-proteostasis CRISPRi screen into ScreenRecords.

A seventh EXTERNAL test, and the only CRISPRi one. Every earlier external
dataset is CRISPR knockout; AssayBench's train fold holds 1303 Knockout screens
against 28 Inhibition (CRISPRi) screens, so this is the under-represented
modality as well as a sorted reporter readout -- the same class as Wilson 2024,
the one external dataset where AF + DepMap clearly helped.

Source: Samelson, Ariqat, McKetney et al., "CRISPR screens in iPSC-derived
neurons reveal principles of tau proteostasis", Cell 189:1517-1534.e19 (2026),
doi 10.1016/j.cell.2025.12.038, PMID 41610849 (Kampmann lab). The article is
not open access and its PMC record is still inside the NIH embargo, so the two
supplementary workbooks were supplied directly: ``mmc2.xlsx`` (all genes) and
``mmc3.xlsx`` (the called hits), each with one sheet per screen.

**Only the genome-wide screen is usable.** The six other sheets are the
follow-up screens on a focused ~1000-gene sublibrary (984-986 genes). The
budget in this comparison is 1000 picks, so a run on any of them would simply
acquire the entire library and both arms would score identically by
construction; they are read and reported, not scored.

**The library is TSS-level, as Kampmann CRISPRi libraries are.** 20 436 entries
for 18 825 genes: 1567 genes carry two or more transcription-start-site entries
(``HMBS_P1``, ``PRAF2_P2``, ``NDUFA9_ENST00000266544.5``). They are collapsed to
the base symbol, keeping the strongest entry and OR-ing the hit calls, because a
duplicated gene list makes ``reveal()`` collapse entries and corrupt the hit
labels. A further 20 436 rows are ``NTC_*`` negative-control quasi-genes, five
to seven copies each for the library's sublibraries; they are not genes and are
dropped from both the library and the hit list.

**Hits are the authors' own calls, read from mmc3, not re-derived.** They are
NOT reconstructible from mmc2's columns: ``gene score`` is
``epsilon * -log10(pvalue)``, but no threshold on it reproduces the lists --
1017 genes clear the weakest negative hit's score while only 729 are called,
and genes at |score| 1.2 are left out while genes at 0.51 are kept. That is the
signature of the Kampmann hit rule, whose discriminant threshold is calibrated
against the NTC distribution SEPARATELY PER SUBLIBRARY, and the table does not
say which sublibrary a gene sits in. The authors' lists are therefore taken as
given, which is this project's policy for every other external dataset too.
45 of mmc3's 1142 genome-wide rows are NTC quasi-genes -- the authors' own
false-discovery background -- and are dropped with the rest of the controls.
Four hit rows have their Ensembl transcript suffix truncated
(``NDUFA9_ENST00000266544`` for ``...544.5``) and are matched back by prefix.

**Direction is fixed by an internal control, not by guesswork.** MAPT itself is
the single most negative gene in the screen (epsilon -2.454, score -9.51):
knocking down tau must lower the tau signal, so negative epsilon is a DECREASE
in the sorted readout and positive is an INCREASE. The abstract's own finding
that loss of CRL5/SOCS4 raises tau agrees -- CUL5 is the second-strongest
positive gene.

**Leakage: none.** Published January 2026, after every AssayBench year fold, and
no screen in any fold mentions tau, MAPT, tauopathy, oligomers, iPSC-derived
neurons or this group. The two iPSC screens in the benchmark are a neural
progenitor viral-infection screen (train, 1172) and a tunicamycin-resistance
screen (test, 2454).
"""

from __future__ import annotations

import argparse
import json
import logging
import pickle
import re
from pathlib import Path

import pandas as pd

from assaybench.data.screens import ScreenRecord
from build_schukken import hgnc_synonyms, remap_library

log = logging.getLogger("build_tau")

SHEET = "GenomeWide"
COLS = ["gene", "epsilon", "pvalue", "score"]

CONDITION = (
    "Genome-wide CRISPR interference (CRISPRi) knockdown screen in "
    "iPSC-derived human neurons carrying a doxycycline-inducible CRISPRi "
    "machinery, stained for tau and sorted by antibody fluorescence, "
    "comparing the high- and low-signal sorted populations"
)

ARMS = {
    "inc": ("Genes whose knockdown increases the neuronal tau level, enriched "
            "in the high-tau sorted population", "increase"),
    "dec": ("Genes whose knockdown decreases the neuronal tau level, enriched "
            "in the low-tau sorted population", "decrease"),
}


def tss_base(g: str) -> str:
    """``HMBS_P1`` -> ``HMBS``, ``NDUFA9_ENST00000266544.5`` -> ``NDUFA9``.

    HGNC symbols never contain an underscore, so the suffix is unambiguous.
    """
    g = re.sub(r"_ENST.*$", "", g)
    return re.sub(r"_P\d.*$", "", g)


def read_sheet(path: Path, sheet: str) -> pd.DataFrame:
    d = pd.read_excel(path, sheet_name=sheet)
    d.columns = COLS
    d["gene"] = d["gene"].astype(str).str.strip()
    for c in COLS[1:]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    return d


def build(all_xlsx: Path, hits_xlsx: Path, hgnc: Path, vocab: Path) -> list[ScreenRecord]:
    d = read_sheet(all_xlsx, SHEET)
    h = read_sheet(hits_xlsx, SHEET)

    ntc = d["gene"].str.startswith("NTC_")
    log.info("mmc2: %d rows, %d NTC control rows (%d distinct quasi-genes)",
             len(d), int(ntc.sum()), d.loc[ntc, "gene"].nunique())
    d = d[~ntc].reset_index(drop=True)

    h_ntc = h["gene"].str.startswith("NTC_")
    log.info("mmc3: %d hit rows, %d of them NTC quasi-genes (the authors' own "
             "false-discovery background)", len(h), int(h_ntc.sum()))
    called = set(h.loc[~h_ntc, "gene"])

    # Four hit rows lost the decimal part of their Ensembl transcript id.
    unknown = called - set(d["gene"])
    if unknown:
        by_prefix = {}
        for g in d["gene"]:
            by_prefix.setdefault(g, g)
        for u in sorted(unknown):
            cand = [g for g in d["gene"] if g.startswith(u)]
            if len(cand) == 1:
                called.discard(u)
                called.add(cand[0])
                log.info("  repaired truncated hit name %s -> %s", u, cand[0])
            else:
                log.warning("  hit %s matches %d library entries; dropped", u, len(cand))
                called.discard(u)

    d["hit"] = d["gene"].isin(called)
    log.info("genome-wide: %d library entries, %d called", len(d), int(d["hit"].sum()))

    # Collapse the TSS-level library to one row per gene: strongest entry wins
    # the score, and the gene is a hit if ANY of its entries was called.
    d["base"] = [tss_base(g) for g in d["gene"]]
    d["mag"] = d["score"].abs()
    hit_any = d.groupby("base")["hit"].any()
    pos_any = d.groupby("base").apply(
        lambda g: bool((g["hit"] & (g["epsilon"] > 0)).any()), include_groups=False)
    neg_any = d.groupby("base").apply(
        lambda g: bool((g["hit"] & (g["epsilon"] < 0)).any()), include_groups=False)
    n_multi = int((d.groupby("base").size() > 1).sum())
    g = (d.sort_values("mag", ascending=False)
          .drop_duplicates("base", keep="first")
          .sort_values("score", ascending=False)
          .reset_index(drop=True))
    log.info("collapsed %d TSS entries to %d genes (%d had >1 start site)",
             len(d), len(g), n_multi)

    syn = hgnc_synonyms(hgnc)
    known = set(json.loads(vocab.read_text())["itos"])
    genes = remap_library(list(g["base"]), known, syn)
    n_oov = sum(1 for x in genes if x not in known)
    log.info("library %d genes, %d unique; %d (%.1f%%) outside the model vocabulary",
             len(genes), len(set(genes)), n_oov, 100 * n_oov / len(genes))
    if len(set(genes)) != len(genes):
        raise RuntimeError(f"gene list not unique ({len(set(genes))} of {len(genes)})")

    masks = {
        "inc": [bool(pos_any[b]) for b in g["base"]],
        "dec": [bool(neg_any[b]) for b in g["base"]],
    }
    masks["merged"] = [bool(hit_any[b]) for b in g["base"]]
    both = sum(1 for a, b in zip(masks["inc"], masks["dec"]) if a and b)
    log.info("inc %d, dec %d, overlap %d, merged %d", sum(masks["inc"]),
             sum(masks["dec"]), both, sum(masks["merged"]))

    mag = [float(v) for v in g["mag"].fillna(0.0)]
    out: list[ScreenRecord] = []
    for suffix in ("inc", "dec", "merged"):
        hits = masks[suffix]
        if suffix == "merged":
            phenotype = ("Genes whose knockdown changes the neuronal tau level "
                         "in either direction")
            direction = "both"
        else:
            phenotype, direction = ARMS[suffix]
        out.append(ScreenRecord(
            dataset_name=f"TAU_CRISPRi_{suffix}",
            split="external",
            organism="Homo sapiens",
            gene_symbol_convention="HGNC",
            phenotype=phenotype,
            cleaned_phenotype="Molecular Output / Reporter / Pathway Activity",
            cell_line="iPSC-derived neuron",
            cell_type="human iPSC-derived glutamatergic neuron, CRISPRi",
            library_type="hCRISPRi-v2 genome-wide CRISPRi sgRNA library",
            library_methodology="Inhibition",
            direction_str=direction,
            condition_clause=CONDITION,
            num_genes=len(genes),
            genes=list(genes),
            relevance_scores=[m if hh else 0.0 for m, hh in zip(mag, hits)],
            hits=hits,
            description=CONDITION,
            contrast_label="high versus low tau antibody signal",
            notes=("Samelson et al. 2026 (Cell 189:1517, doi "
                   "10.1016/j.cell.2025.12.038). Hits are the authors' own "
                   "calls from supplementary table mmc3; library and "
                   "phenotype scores from mmc2."),
            unmasked_relevance_scores=mag,
        ))
    return out


def secondary_report(all_xlsx: Path, hits_xlsx: Path) -> None:
    """Print, but do not build, the focused follow-up screens."""
    a = pd.read_excel(all_xlsx, sheet_name=None)
    b = pd.read_excel(hits_xlsx, sheet_name=None)
    print("\nfollow-up screens on the focused sublibrary (NOT scored: the "
          "library is smaller than the 1000-pick budget)")
    for k in a:
        if k == SHEET:
            continue
        ga = {tss_base(str(x).strip()) for x in a[k].iloc[:, 0]}
        gb = {tss_base(str(x).strip()) for x in b[k].iloc[:, 0]}
        ga = {x for x in ga if not x.startswith("NTC_")}
        gb = {x for x in gb if not x.startswith("NTC_")}
        print(f"  {k:24s} {len(ga):5d} genes {len(gb):5d} hits  "
              f"frac {len(gb & ga) / max(len(ga), 1):.3f}")


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    ap = argparse.ArgumentParser()
    ap.add_argument("--all", required=True, type=Path, help="mmc2.xlsx")
    ap.add_argument("--hits", required=True, type=Path, help="mmc3.xlsx")
    ap.add_argument("--hgnc", required=True, type=Path)
    ap.add_argument("--vocab", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()

    screens = build(a.all, a.hits, a.hgnc, a.vocab)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "wb") as fh:
        pickle.dump(screens, fh)
    print(f"\n{len(screens)} screens -> {a.out}")
    for s in screens:
        print(f"  {s.dataset_name:22s} {s.num_genes:6d} genes {s.total_hits:5d} hits "
              f"frac {s.total_hits / s.num_genes:.4f}  EF ceiling "
              f"{s.num_genes / max(s.total_hits, 1000):5.2f}")
    secondary_report(a.all, a.hits)


if __name__ == "__main__":
    main()
