> Historical research record, copied 2026-10-01. Some linked raw artifacts remain in the local archive. See the package README for included records and runnable entry points.

# Cross-dataset numerical learning: protocol before execution

2026-09-29. User explicitly requested broader testing now. New API cost: $0. No LLM calls. Use existing fixed numerical methods, released outcomes and feature cache; do not retune to each screen's unobserved labels.

## Question

Do the IFNG/IL2 differences recur on other single-gene screens? Does their direction depend on representation, prediction target, or initial observation composition? This extends method analysis rather than launching another rich-critic method.

## Sources and conceptual moves

- [BioDiscoveryAgent, ICLR2025](https://arxiv.org/abs/2405.17631): multi-screen sequential gene selection → extend our small numerical campaigns across its available single-gene tasks. This is not a reproduction of its full budgets/agent.
- [LLMNN / LLMs for Bayesian Optimization, EMNLP2025](https://arxiv.org/abs/2509.21403): separate priors from numerical acquisition and test dependence on feedback → retain geometry, random and learned controls; follow up with frozen histories where informative. This is an experimental-design transfer, not implementation of LLMNN.
- [GOLLuM, NMI2026](https://www.nature.com/articles/s42256-026-01283-z): task-relevant representation learning with GP objective → compare existing compressed/full representations across screens; no learned encoder reproduction or new-method claim.

## Frozen phase1 design

Seven released single-gene tables: IFNG, IL2, Carnevale, Sanchez day21, Sanchez day30, Scharenberg, CAR-T. The two Sanchez tables are related measurements, not independent studies or just alternative relabelings. Horlbeck is excluded because candidates are gene pairs, whereas this representation and model interface are for single genes. A valid pairwise feature construction and equal-budget experiment is separate work.

Three seeds8101–8103 per table. Matched256 starting measurements (64 geometry-prior +192 random), three new batches of32. Arms: geometry-centrality, random, isotropic GP on absolute effect, median-effect tree on PCA16, same tree on full778 features, full-feature binary classifier. In each campaign use its own newly revealed results, keep unselected genes available, and remove tested genes. Fit only on observed measurements. All dataset labels and feature coverage are held constant across arms within a screen. Random and geometry are controls, not BDA LLM baselines.

No target direction is silently redefined: exact released hit membership is the evaluator. Effect models retain the previously tested absolute-score objective. This may mismatch some screens, especially Sanchez day21's two tails on a statistical score. That is an explicit transfer diagnostic, not a claim that absolute scores universally encode biological utility. Binary classification supplies a label-aligned control. Interpret score semantics using the existing verified [EDA source record](../../EDA/SOURCE_RECORD.md), not the task prose alone; CAR-T headers0/1 must not lose the first row.

Report per-screen hits/288 new choices per arm, per-seed differences, hit rate, remaining-hit recall, mean endpoint discovery curve and time. Do not pool raw counts across prevalence/pool-size differences to declare a winner. Scharenberg's256 history is a larger fraction of its pool; report the resulting initial hit depletion. Candidate feature coverage and distributions are evaluator diagnostics, not agent inputs.

After phase1, choose a limited confirmation based on the results: repeat key arms with random-only histories and/or compare updated/frozen histories. Record this decision before the continuation. Keep all negative results. These tests consume formerly reserved screens for exploration; none is thereafter a pristine final holdout. Carnevale already appears in earlier project experiments.

Validity limits: Achilles features have known cached lineage but unverified independent rebuild/assay overlap; seeds share genes and biological datasets. No claim of independent seven-study replication, causal biological explanation, or publication novelty. Results may identify better next tests, not guarantee a new method.
