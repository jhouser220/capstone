# GPT-OSS 120B matched handoff experiment

## Purpose

Extend the team's existing LLM-to-numerical handoff analysis to
GPT-OSS 120B.

Protocol:

1. GPT-OSS 120B selects genes for rounds 1-3.
2. The exact same three-round warm start is replayed into:
   - AssayFormer
   - HSA
3. The numerical method controls rounds 4-10.

This matched design holds the LLM starting condition fixed so the
downstream difference is attributable to the continuation method.

## Current status

13 of the planned 15 paper-test screens are complete.

Still missing because of Groq quota limits:

- U_2423_merged
- U_2427_dec

## Preliminary 13-screen result

- AssayFormer total hits: 1,301
- HSA total hits: 1,441
- HSA delta: +140 hits
- HSA wins on total hits: 9 / 13 screens
- AssayFormer measurable selections: 11,602
- HSA measurable selections: 12,407
- AssayFormer pooled measured hit rate: 11.21%
- HSA pooled measured hit rate: 11.61%

A substantial portion of HSA's total-hit advantage comes from keeping
more proposed genes within the measurable assay universe. Ranking
performance remains heterogeneous across screens.

These results are preliminary until the remaining two screens are
completed.

## Files

- `run_matched_methods.py`
  - matched AssayFormer/HSA continuation runner

- `warm_start_rounds_13_of_15.json`
  - saved GPT-OSS rounds 1-3 for the 13 completed screens

- `gptoss_13screen_handoff_results.csv`
  - current matched handoff results

- `collect-gpt-oss-120b-groq.yaml`
  - GPT-OSS/Groq configuration
  - API key is referenced through the `GROQ_API_KEY` environment
    variable and is not stored in this repository
