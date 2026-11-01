"""Reading a screen's selection direction from its phenotype sentence, including the inverse-quantity cases."""
from __future__ import annotations

import pytest

from assayloop.screen_direction import BOTH, DROPOUT, ENRICHMENT, UNKNOWN, classify_selection, should_use_llm


@pytest.mark.parametrize("phenotype,expected", [
    # readout goes UP on knockout: the condition needed the gene
    ("increases resistance to HCoV-229E coronavirus infection.", ENRICHMENT),
    ("increase resistance to SARS-CoV-2 infection", ENRICHMENT),
    ("increases accumulation of NMD-reporter protein as measured by fluorescence.", ENRICHMENT),
    ("increases rAAV (AAV2-GFP) transgene expression as measured by GFP protein accumulation.", ENRICHMENT),
    ("increases T-cell proliferation as measured by enrichment of highly dividing T cells", ENRICHMENT),
    ("decreases sensitivity to erlotinib", ENRICHMENT),                 # inverse quantity, down verb
    # readout goes DOWN on knockout: the cells needed the gene
    ("decrease cell proliferation in a KRASG12V mutant background", DROPOUT),
    ("decreases mTORC1 activity as measured by reduced rpS6 phosphorylation.", DROPOUT),
    ("decrease TNF-alpha production by CD4+ T cells", DROPOUT),
    ("increases sensitivity to Cpd-391 as measured by decreased cell proliferation", DROPOUT),
    ("decrease resistance to hydroxychloroquine", DROPOUT),
    ("increases susceptibility to infection", DROPOUT),
    # not decidable from the sentence
    ("either increases or decreases IL-2 production by CD4+ T cells", BOTH),
    ("either increase or decrease cell proliferation", BOTH),
    ("impacts response to chemicals", UNKNOWN),
    ("impacts cell proliferation", UNKNOWN),
    ("", UNKNOWN),
])
def test_classify_selection(phenotype, expected):
    assert classify_selection(phenotype)[0] == expected


def test_why_string_quotes_the_evidence():
    d, why = classify_selection("increases sensitivity to Cpd-391 as measured by decreased proliferation")
    assert d == DROPOUT and "inverse quantity" in why and "increases" in why


def test_only_the_leading_clause_decides():
    # the trailing "decreased cell proliferation" is how the gain in sensitivity was measured, not a second call
    assert classify_selection("increases resistance to RSL3, as measured by decreased cell death")[0] == ENRICHMENT


@pytest.mark.parametrize("direction,mode,amb,expected", [
    (ENRICHMENT, "all", False, True),
    (DROPOUT, "all", False, True),
    (BOTH, "all", False, True),
    (ENRICHMENT, ENRICHMENT, False, True),
    (DROPOUT, ENRICHMENT, False, False),
    (BOTH, ENRICHMENT, False, False),
    (BOTH, ENRICHMENT, True, True),
    (UNKNOWN, ENRICHMENT, True, True),
    (DROPOUT, ENRICHMENT, True, False),            # a known opposite direction is never admitted
    (DROPOUT, DROPOUT, False, True),
])
def test_should_use_llm(direction, mode, amb, expected):
    assert should_use_llm(direction, mode, amb) is expected


def test_matches_the_curated_screen_sets():
    pytest.importorskip("assaybench")
    from assayloop.tasks import load_screens

    got = {s.dataset_name: classify_selection(s.phenotype)[0] for s in load_screens(target_set="paper_validation")}
    # the 6 host-pathogen screens and the 2 reporter-style screens select for gain of readout
    for n in ("1561", "1581", "1584", "1815", "1817", "1818", "2048", "2104", "1644", "TR_1975_1976"):
        assert got[n] == ENRICHMENT, (n, got[n])
    for n in ("1725", "U_1793_dec", "U_1794_dec", "U_1795_dec", "U_1796_dec", "U_1797_dec"):
        assert got[n] == DROPOUT, (n, got[n])
    assert got["U_1793_merged"] == BOTH and got["1990"] == UNKNOWN
    assert sum(v == ENRICHMENT for v in got.values()) == 10
