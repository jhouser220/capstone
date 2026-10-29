# Failure-Aware Reflection: IL-2 Pilot

## Research question
Can an LLM agent that critically examines unsuccessful
experiments improve subsequent gene selection compared
with standard feedback and general-purpose reflection?

## Experimental setup
- Model: gemini-3.5-flash-lite
- Dataset: IL-2
- Five rounds of 128 genes each
- Identical first batch across all three approaches
- One pilot run per approach
- Previously measured outcomes only

## Results
- Standard full feedback: 68 hits (10.40% recall)
- General reflection: 56 hits (8.56% recall)
- Failure-aware reflection: 109 hits (16.67% recall)
- Failure-aware advantage: 41 hits over standard

## Biological sensitivity analysis
- All benchmark hits: 109 vs. 68
- Excluding RPL/RPS families: 75 vs. 68
- Excluding broader RNA/translation/transport families:
  57 vs. 68

Failure-aware discovered 34 RPL/RPS hits and 18
additional hits in the broader excluded families.

Its final round discovered 46 hits, including many
ribosomal and RNA-processing-related genes.

## Interpretation
The failure-aware agent shifted gene selection toward
RNA processing, translation, and cellular transport.

This shift was associated with higher benchmark recall.
However, most of the advantage was concentrated in
broad cellular-process gene families.

The exclusions are post-hoc sensitivity analyses based
on gene-name prefixes, not formal essentiality or
pathway annotations. Filtered selection budgets differ.

One matched-start run is insufficient to establish a
general performance advantage or causal superiority.

No new wet-lab experiments were performed.

## Reproducibility
- research_assistant.py: failure-aware mode
- tools.py: failure-analysis agent and matched start
- failure_aware_IL2/: gene selections and agent logs
- failure_aware_comparison.csv: headline comparison
- biological_sensitivity_results.csv: filtered results
- hit_family_composition.csv: family composition
- family_hits_by_round.csv: round-level results
