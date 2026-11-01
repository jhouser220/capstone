"""Build the positive gene sets for GATE 0, including ones appropriate to a
GAIN-OF-FUNCTION atlas.

Why this exists. §4.39 failed JUMP ORF on gate 0 because co-complex genes are no more
similar than random pairs, and then explained the failure by stoichiometry: a
gain-of-function phenotype depends on dosage and balance, so overexpressing ONE subunit
of a complex creates an IMBALANCE. That mechanism PREDICTS the failure, which
disqualifies the control. Co-complex membership tests whether an atlas carries
complex-membership signal; on an overexpression atlas it is the wrong ground truth, and
failing it says nothing about whether the atlas carries FUNCTIONAL signal.

So this builds three groups, which behave differently under loss and gain of function:

  complex  Physical complexes. Under KNOCKOUT, subunits converge -- removing any one
           disables the whole machine, which is why DepMap scores +4 to +11 sd here.
           Under OVEREXPRESSION, the stoichiometry argument predicts divergence.
  paralog  Sequence/domain families: substitutes, not partners. Overexpressing either
           of two paralogues should do a SIMILAR thing, because they act on the same
           substrates -- no imbalance is created. The prediction is the reverse of the
           complex case, and it is reversed for DEPENDENCY too: paralogues BUFFER each
           other, so knocking out one of a redundant pair often has no phenotype, and
           dependency profiles of paralogues are routinely uncorrelated. This group is
           therefore the discriminating one: a double dissociation (paralogues beat
           complexes in JUMP while complexes beat paralogues in DepMap) would show the
           overexpression atlas carries real functional signal of the kind it SHOULD.
  pathway  Same-pathway membership, the weakest of the three: pathway sets mix
           complexes, substrates and regulators, so neither prediction is clean. Kept
           as a middle reference, not as a decision gate.

Sources, both free:

  HGNC `hgnc_complete_set.txt` (17 MB, `storage.googleapis.com/public-download-files/
  hgnc/tsv/tsv/`), column `gene_group` -- curated families given as gene SYMBOLS, so no
  id mapping. The column mixes families with complexes ("DIP2 family" next to "MTOR
  complex 2"), which is why COMPLEX_RE splits them by name rather than trusting it.

  MSigDB C2 CP (already in `src/assayloop/data/gene_sets/`) for the pathway group.

TWO EXCLUSIONS, both decided on grounds independent of any result, and both stated here
because they are the only hand choices in the file:

  - Families not expressed in a cultured line (olfactory and taste receptors,
    immunoglobulin and T-cell receptor segments, keratin-associated proteins, late
    cornified envelope). Overexpressing them produces no phenotype in either atlas, so
    they are noise in the pooled statistic rather than evidence about either.
  - Families above MAX_SIZE genes. A 300-member class is a naming convention, not a
    claim that its members do the same thing.

    uv run python scripts/diagnostics/build_gate0_sets.py \
        /home/jeff/hgnc/hgnc_complete_set.txt output/analysis/gate0_sets.json
"""
from __future__ import annotations

import collections
import csv
import json
import re
import sys
from pathlib import Path

MIN_SIZE, MAX_SIZE = 2, 12

# A family whose NAME says "complex" is a complex whatever column it arrived in.
COMPLEX_RE = re.compile(
    r"complex|subunits|components|holoenzyme|ribosom|proteasom|spliceosom|condensin|"
    r"cohesin|nucleoporin|chaperonin|exosome|signalosome|inflammasome|"
    r"\bcomplexes\b|particle|anaphase-promoting",
    re.I)
# Not expressed in a cultured cell line -> no ORF or knockout phenotype in any atlas.
NOT_EXPRESSED_RE = re.compile(
    r"olfactory receptor|taste receptor|immunoglobulin|T cell receptor|"
    r"keratin associated|late cornified|vomeronasal|pheromone",
    re.I)
# Printed individually: paralogue families whose members have a textbook convergent
# phenotype WHEN OVEREXPRESSED, which is exactly the regime being tested.
SPOTLIGHT = {
    "AKT kinase family", "E2F transcription factors", "Src family tyrosine kinases",
    "Jak family tyrosine kinases", "Aurora kinase family", "Polo like kinase family",
    "Sirtuins", "SMAD family", "Caspases", "Tubulin beta family",
    "Signal transducer and activator of transcription family",
    "Histone deacetylases, class I",
}


def main() -> None:
    hgnc, out = Path(sys.argv[1]), Path(sys.argv[2])
    gmt = Path(sys.argv[3]) if len(sys.argv) > 3 else Path(
        "src/assayloop/data/gene_sets/c2.cp.v2023.2.Hs.symbols.gmt")

    fams: dict[str, list[str]] = collections.defaultdict(list)
    with open(hgnc, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["locus_group"] != "protein-coding gene":
                continue
            for g in (r["gene_group"] or "").split("|"):
                g = g.strip()
                if g:
                    fams[g].append(r["symbol"])
    print(f"HGNC: {len(fams)} gene groups over protein-coding genes")

    groups: dict[str, dict[str, list[str]]] = {"complex": {}, "paralog": {}, "pathway": {}}
    n_big = n_unexpressed = 0
    for name, genes in fams.items():
        genes = sorted(set(genes))
        if NOT_EXPRESSED_RE.search(name):
            n_unexpressed += 1
            continue
        if not (MIN_SIZE <= len(genes) <= MAX_SIZE):
            n_big += 1
            continue
        groups["complex" if COMPLEX_RE.search(name) else "paralog"][name] = genes

    # The six complexes §4.39 used, so the new table is comparable to the old one.
    groups["complex"].setdefault("ribosome (S4.39)", ["RPL3", "RPL4", "RPL5", "RPL6", "RPL7", "RPL8"])
    groups["complex"].setdefault("RNA pol II (S4.39)", ["POLR2A", "POLR2B", "POLR2C", "POLR2D", "POLR2E"])
    groups["complex"].setdefault("eIF3 (S4.39)", ["EIF3A", "EIF3B", "EIF3C", "EIF3D", "EIF3E"])
    groups["complex"].setdefault("proteasome (S4.39)", ["PSMA1", "PSMA2", "PSMA3", "PSMB1", "PSMB2"])
    groups["complex"].setdefault("SF3B (S4.39)", ["SF3B1", "SF3B2", "SF3B3", "SF3B5"])
    groups["complex"].setdefault("mediator (S4.39)", ["MED11", "MED17", "MED18", "MED20", "MED27"])

    if gmt.is_file():
        for line in gmt.read_text().splitlines():
            p = line.split("\t")
            genes = sorted(set(p[2:]))
            # Larger than the family cap: a pathway is not claimed to be a family, and
            # too few small CP sets exist to pool at MAX_SIZE.
            if 5 <= len(genes) <= 30 and not COMPLEX_RE.search(p[0]):
                groups["pathway"][p[0]] = genes
    else:
        print(f"WARNING: {gmt} missing, pathway group empty", file=sys.stderr)

    spot = sorted(SPOTLIGHT & set(groups["paralog"]))
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"groups": groups, "spotlight": spot,
                               "min_size": MIN_SIZE, "max_size": MAX_SIZE}, indent=1))
    for k, v in groups.items():
        npair = sum(len(g) * (len(g) - 1) // 2 for g in v.values())
        print(f"  {k:8s} {len(v):5d} sets, {npair:7d} within-set pairs, "
              f"median size {sorted(len(g) for g in v.values())[len(v) // 2]}")
    print(f"  dropped {n_big} groups outside size {MIN_SIZE}-{MAX_SIZE}, "
          f"{n_unexpressed} not expressed in culture")
    print(f"  spotlight families present: {len(spot)} -> {spot}")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
