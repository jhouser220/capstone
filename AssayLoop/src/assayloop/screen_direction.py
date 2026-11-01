"""Which way a screen selects: does a hit's knockout turn the readout UP or DOWN?

The benchmark's labels are binary, so this is not recorded anywhere; it is read from the screen's own
phenotype sentence (``ScreenRecord.phenotype``, the text the LLM is shown too). It matters because the two
directions are different biological questions, and on the validation screens AssayFormer is far better at one
of them:

  dropout     "decrease cell proliferation", "decreases mTORC1 activity" — a hit is a gene the cells NEED.
              Largely the essential-gene set, shared across screens; AssayFormer ran at 66-95% hit rate.
  enrichment  "increases resistance to MERS-CoV", "increases accumulation of NMD-reporter" — a hit is a gene
              the CONDITION needs, and losing it rescues the cell. Condition-specific, barely shared across
              screens; AssayFormer ran at ~5%, close to each screen's base rate.

Polarity comes from the leading verb, flipped when the measured quantity is itself inverse (sensitivity,
susceptibility): "increases sensitivity" means growth goes DOWN, so it is a dropout screen, as is "decrease
resistance". Screens whose sentence states both directions ("either increases or decreases") or none
("impacts cell proliferation") are reported as such rather than guessed at.
"""

from __future__ import annotations

import re

ENRICHMENT = "enrichment"
DROPOUT = "dropout"
BOTH = "both"
UNKNOWN = "unknown"

_RX_BOTH = re.compile(r"\beither\s+(increase|decrease)s?\s+or\s+(increase|decrease)s?\b", re.I)
_RX_VERB = re.compile(r"\b(increas|decreas|reduc|enhanc|impair|suppress|elevat|diminish)(?:e|es|ed|ing)?\b", re.I)
# the measured quantity, when it runs opposite to the phenotype (knockout raising "sensitivity" lowers growth)
_RX_INVERSE = re.compile(r"\b(sensitivit|sensitis|sensitiz|susceptib|depletion|dropout|death|toxicit|lethalit)", re.I)
_RX_CLAUSE_END = re.compile(r",|\bas\s+(?:measured|indicated|determined|assessed|shown)\b", re.I)
_DOWN_VERBS = {"decreas", "reduc", "impair", "suppress", "diminish"}


def classify_selection(phenotype: str) -> tuple[str, str]:
    """``(direction, why)`` for a screen's phenotype sentence; direction is one of the four constants above."""
    text = (phenotype or "").strip()
    if not text:
        return UNKNOWN, "no phenotype text"
    if _RX_BOTH.search(text):
        return BOTH, "phenotype states both directions"
    m = _RX_VERB.search(text)
    if not m:
        return UNKNOWN, f"no direction verb in {text[:60]!r}"
    stem = m.group(1).lower()
    down = stem in _DOWN_VERBS
    # Only the clause the verb governs decides polarity. The trailing "..., as measured by decreased cell death"
    # names how the change was observed and would otherwise flip the reading, so cut it off first.
    tail = _RX_CLAUSE_END.split(text[m.end():], maxsplit=1)[0][:60]
    inverse = bool(_RX_INVERSE.search(tail))
    up_readout = down == inverse            # "decreases sensitivity" -> up; "increases sensitivity" -> down
    word = m.group(0)
    obj = " ".join(tail.split()[:4])
    why = f"{word!r} + {obj!r}" + (" (inverse quantity)" if inverse else "")
    return (ENRICHMENT if up_readout else DROPOUT), why


def should_use_llm(direction: str, mode: str, include_ambiguous: bool = False) -> bool:
    """Whether the LLM gets slots on a screen with this ``direction``.

    ``mode`` is ``"all"`` (every screen) or one of the direction constants (only screens of that direction;
    ``include_ambiguous`` also admits ``both``/``unknown`` screens).
    """
    if mode == "all":
        return True
    if direction == mode:
        return True
    return include_ambiguous and direction in (BOTH, UNKNOWN)


__all__ = ["classify_selection", "should_use_llm", "ENRICHMENT", "DROPOUT", "BOTH", "UNKNOWN"]
