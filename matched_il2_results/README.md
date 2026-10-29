# Matched-Start IL-2 Reflection Experiment

## Research question
Does an additional LLM reflection/synthesis agent improve
iterative gene discovery compared with standard full feedback?

## Configuration
- Model: gemini-3.5-flash-lite
- Dataset: IL-2
- Candidate pool: 18,939 genes
- Known hits: 654
- Rounds: 5
- Genes per round: 128
- Replications: 1 matched pair
- Both approaches reuse the same initial 128 genes.

## Results
- Standard full feedback: 68 hits
- Full feedback + reflection: 56 hits
- Difference: -12 hits
- Both approaches selected 640 unique genes.
- Both started with 21 hits in round 1.

## Interpretation
Reflection did not improve hit discovery in this pilot.
It changed the selected genes but resulted in fewer hits.

This is one matched-start comparison, not definitive
evidence that reflection generally reduces performance.

## Reproducibility
- research_assistant.py and tools.py contain the
  matched-start implementation and feedback-mode fix.
- shared_seed_IL2/common_seed_1 contains the initial batch.
- standard_full_IL2 contains baseline results.
- with_reflection_IL2 contains reflection results.
- round_comparison.csv contains the round-level analysis.

The first batch is reused; later model responses are stochastic.
No new wet-lab experiments were performed.
