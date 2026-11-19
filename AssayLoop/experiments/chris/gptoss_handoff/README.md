# GPT-OSS 120B Matched Handoff Experiment

## Purpose

Extend the team's existing LLM-to-numerical handoff evaluation to
GPT-OSS 120B as a new LLM condition.

## Protocol

For each AssayLoop paper-test screen:

1. GPT-OSS 120B performs rounds 1-3.
2. The same initial GPT-OSS rounds are used for both downstream arms.
3. AssayFormer or HSA controls rounds 4-10.

## Status

Primary paper-test experiment complete:

- 15 / 15 paper-test screens
- 3 GPT-OSS rounds + 7 numerical rounds
- AssayFormer vs HSA

## Final results

| Metric | AssayFormer | HSA |
|---|---:|---:|
| Total hits | 1,373 | 1,509 |
| Mean hits / screen | 91.53 | 100.60 |
| Mean recall | 0.2503 | 0.2567 |
| Mean enrichment factor | 5.2931 | 5.4276 |
| Pooled measured hit rate | 10.23% | 10.50% |

HSA found 136 more total hits and produced more hits on 10 of 15 screens.

HSA also produced 950 more measurable selections
(14,365 vs 13,415), so part of the total-hit advantage reflects better
coverage of the measurable gene universe. HSA had the higher
measured-gene hit rate on 5 of 15 screens.

## Files

- `run_matched_methods.py`
- `collect-gpt-oss-120b-groq.yaml`
- `warm_start_rounds_15_of_15.json`
- `gptoss_15screen_handoff_results.csv`

## Next experiments

1. AssayLoop validation set
2. External reconstructed datasets
3. Secondary hybrid condition
4. AssayLLM confirmation
