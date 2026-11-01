"""Turn Schukken et al. 2026 whole-genome CRISPR screens into ScreenRecords.

An EXTERNAL held-out test set: posted April 2026, so it post-dates every
AssayBench year fold and cannot be in AssayFormer's training data.

Source: Schukken, Aitken, Sheltzer et al., "Paired CRISPR screens identify
mitochondrial metabolism and UBE2H ...", DOI 10.64898/2026.04.26.720636
(PMC13142536). Data: ``SupData/SD3- CRISPR Screen Beta Scores.xlsx``, sheet
``Raw WGL``, from github.com/kschukken/Aneuploid_CRISPR_Screens.

Three decisions are baked in here, each one checkable:

**Only the whole-genome library.** The workbook also holds 26 Druggable-Genome
screens at 3377 genes. A 1000-pick budget against a 3377-gene library caps EF
at ``D/budget = 3.4``, which compresses the very quantity the comparison is
about; the WGL sheet's 20111 genes sit in the same regime as AssayBench's own
genome-scale screens. This is §4.29's genome-scale filter, applied up front.

**A hit is ``beta < HIT_BETA`` (-0.75).** The paper's own -0.2/p<0.05 cut is for
its *paired* aneuploid-vs-euploid contrast, which leaves 15 hits genome-wide --
far too few to measure EF over 1000 picks. Each column is instead read as what
it is on its own: a knockout dropout screen. -0.75 is not tuned on any model;
it is the threshold that puts these 14 screens on top of the AssayBench
genome-scale knockout *Fitness/Proliferation/Viability* hit-fraction
distribution (median 0.094, p10 0.051, p90 0.121), which is what makes EF
comparable to the numbers this repo already quotes. ``--hit-beta`` re-runs it
anywhere else; -0.5 is the obvious sensitivity check.

**Old gene symbols are remapped through HGNC.** Brunello is a 2016 library, so
11.2% of its symbols are not in AssayFormer's 24217-gene vocabulary -- `AARS`
(now AARS1), `ADCK3` (COQ8A), `ADRBK1` (GRK2). Every one of those collapses onto
``<unk>`` (id 0), where the model cannot rank it at all, so leaving them alone
would handicap AssayFormer on a ninth of the library for a bookkeeping reason.
Remapping each symbol to its current approved form -- only when exactly one
approved symbol claims it, never on an ambiguous alias -- recovers 1130 of 2246
and takes vocabulary coverage to 94.4%. The residual 5.5% stays in the library
and stays unrankable: that is a real cost of deploying on an older library, it
is reported, and both arms carry it identically.
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

log = logging.getLogger("build_schukken")

HIT_BETA = -0.75

# Hand-written, because the column names carry the karyotype and the sheet has
# no metadata rows. "Cont" = the near-euploid parental control; everything else
# is an engineered aneuploid clone of that same parent. Ts = trisomy,
# Di = disomy (a single-copy loss from a near-triploid background), p/q = arm.
CELL_LINES = {
    "DLD1": ("DLD-1", "colorectal adenocarcinoma"),
    "HCT116": ("HCT 116", "colorectal carcinoma"),
    "Vaco432": ("VACO432", "colorectal adenocarcinoma"),
    "A2780": ("A2780", "ovarian carcinoma"),
    "A2058": ("A2058", "melanoma"),
    "AGS": ("AGS", "gastric adenocarcinoma"),
    "MCF10A": ("MCF 10A", "non-tumorigenic mammary epithelium"),
}

KARYOTYPE = {
    "Cont": "the near-euploid parental control",
    "Ts13": "an engineered clone trisomic for chromosome 13",
    "Di1q": "an engineered clone with a single-copy loss of 1q",
    "7pDi": "an engineered clone with a single-copy loss of 7p",
}


def hgnc_synonyms(path: Path) -> dict[str, set[str]]:
    """``{old_or_alias_symbol: {approved symbols claiming it}}``."""
    syn: dict[str, set[str]] = {}
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            sym, status = r["symbol"].strip(), r["status"].strip()
            if not sym or status != "Approved":
                continue
            for col in ("prev_symbol", "alias_symbol"):
                for alias in (r.get(col) or "").split("|"):
                    alias = alias.strip()
                    if alias and alias != sym:
                        syn.setdefault(alias, set()).add(sym)
    return syn


def remap_library(raw: list[str], target: set[str],
                  syn: dict[str, set[str]]) -> list[str]:
    """Remap a whole library's symbols to current HGNC form, collision-free.

    Per gene: keep it if the model already knows the symbol, else move it to its
    current approved symbol. Ambiguity is left alone on purpose -- if two
    approved symbols both claim an old alias, picking one invents a gene
    identity.

    The collision guard is the part that matters. A remap is skipped whenever
    its target symbol is ALREADY in this library, either as an original symbol
    or because an earlier row was remapped onto it. Without it, 22 rows of this
    library land on a symbol it already carries (the library measured one gene
    twice, under an old alias and a current one), the gene list stops being
    unique, and ``reveal()`` silently collapses the pair to the last index --
    which would quietly corrupt the hit labels rather than fail. Skipping
    instead leaves the duplicate row under its old symbol, where it is OOV and
    unrankable; that is the honest reading, since the gene is still represented
    once, by the row that already had the right name.
    """
    taken = set(raw)
    out: list[str] = []
    for g in raw:
        if g in target:
            out.append(g)
            continue
        cands = [c for c in syn.get(g, ()) if c in target]
        if len(cands) == 1 and cands[0] not in taken:
            taken.add(cands[0])
            out.append(cands[0])
        else:
            out.append(g)
    return out


def describe(line: str, cell_type: str, karyo: str) -> tuple[str, str]:
    """``(phenotype, condition_clause)`` for one column."""
    phenotype = (
        "Genes whose knockout reduces cell proliferation and is depleted from "
        "the population over ten doublings (MAGeCK beta score)"
    )
    condition = (
        f"Genome-wide CRISPR-Cas9 knockout proliferation screen in {line} "
        f"({cell_type}), {karyo}, harvested after ten population doublings"
    )
    return phenotype, condition


def build(xlsx: Path, hgnc: Path, vocab: Path, hit_beta: float) -> list[ScreenRecord]:
    df = pd.read_excel(xlsx, sheet_name="Raw WGL")
    syn = hgnc_synonyms(hgnc)
    known = set(json.loads(vocab.read_text())["itos"])

    raw = [str(g).strip() for g in df["Gene"]]
    genes = remap_library(raw, known, syn)
    n_oov = sum(1 for g in genes if g not in known)
    log.info(
        "library %d genes, %d unique; %d (%.1f%%) still outside the model vocabulary",
        len(genes), len(set(genes)), n_oov, 100 * n_oov / len(genes),
    )
    if len(set(genes)) != len(genes):
        raise RuntimeError(
            f"gene list is not unique ({len(set(genes))} of {len(genes)}); "
            "reveal() would collapse duplicates and corrupt the hit labels"
        )

    out: list[ScreenRecord] = []
    for col in [c for c in df.columns if c != "Gene"]:
        stem = col[: -len(".Beta")] if col.endswith(".Beta") else col
        prefix, _, karyo_key = stem.partition(".")
        line, cell_type = CELL_LINES[prefix]
        karyo = KARYOTYPE[karyo_key]
        beta = pd.to_numeric(df[col], errors="coerce")

        hits = [bool(b < hit_beta) for b in beta.fillna(0.0)]
        # Hits carry dropout magnitude so a score-aware metric sees a gradient;
        # non-hits are masked to 0, which is the AssayBench convention.
        rel = [float(max(0.0, -b)) if h else 0.0 for b, h in zip(beta.fillna(0.0), hits)]
        phenotype, condition = describe(line, cell_type, karyo)

        out.append(ScreenRecord(
            dataset_name=f"SCHUKKEN_{stem.replace('.', '_')}",
            split="external",
            organism="Homo sapiens",
            gene_symbol_convention="HGNC",
            phenotype=phenotype,
            cleaned_phenotype="Fitness / Proliferation / Viability",
            cell_line=line,
            cell_type=cell_type,
            library_type="Brunello whole-genome CRISPR knockout library",
            library_methodology="Knockout",
            direction_str="decrease",
            condition_clause=condition,
            num_genes=len(genes),
            genes=list(genes),
            relevance_scores=rel,
            hits=hits,
            description=condition,
            contrast_label="day 10 vs day 0",
            notes=(
                "Schukken et al. 2026 (PMC13142536), aneuploid vs near-euploid "
                "isogenic clone panel. Hits are genes depleted in this screen on "
                f"its own (MAGeCK beta < {hit_beta}), not the paper's paired "
                "aneuploid-specific dependencies."
            ),
            duration="ten population doublings",
            unmasked_relevance_scores=[float(-b) for b in beta.fillna(0.0)],
        ))
    return out


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    ap = argparse.ArgumentParser()
    ap.add_argument("--xlsx", required=True, type=Path)
    ap.add_argument("--hgnc", required=True, type=Path)
    ap.add_argument("--vocab", required=True, type=Path)
    ap.add_argument("--hit-beta", type=float, default=HIT_BETA)
    ap.add_argument("--out", required=True, type=Path)
    a = ap.parse_args()

    screens = build(a.xlsx, a.hgnc, a.vocab, a.hit_beta)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "wb") as fh:
        pickle.dump(screens, fh)

    print(f"\n{len(screens)} screens -> {a.out}")
    print(f"{'screen':34s} {'genes':>7} {'hits':>6} {'frac':>7}")
    for s in screens:
        print(f"{s.dataset_name:34s} {s.num_genes:7d} {s.total_hits:6d} "
              f"{s.total_hits / s.num_genes:7.4f}")


if __name__ == "__main__":
    main()
