> Historical research record, copied 2026-10-01. Some linked raw artifacts remain in the local archive. See the package README for included records and runnable entry points.

# Objective, representation and useful feedback

2026-09-29. Follow-on to the rich numerical study and GPT-6 Astra spot check. Method analysis first. Prior conservative project spend$4.001451; full cap$50. No charge for local numerical rounds. Refresh provider prices and reconcile all ledgers before paid work.

## Source mechanisms and transfers

| Source | Original mechanism | Our transfer / fidelity limit |
|---|---|---|
|[GOLLuM, NMI2026](https://www.nature.com/articles/s42256-026-01283-z)|GP objective adapts representations to outcome prediction|Projection: isolate representation versus target/loss in matched numerical models. Fixed raw biological profiles versus16 PCA components are a diagnostic, not GOLLuM or encoder adaptation.|
|[LLMNN, Findings EMNLP2025](https://arxiv.org/abs/2509.21403)|Separates scientific prior from numerical acquisition; tests response to shuffled outcomes|Design transfer: matched no-report/report and real/shuffled-new-result controls; preserve available untested genes.|
|[Learning to Rank for In-Context Example Retrieval, NeurIPS2025](https://proceedings.nips.cc/paper_files/paper/2025/hash/cf68c8b43e692c63deec7b1242bbb3a0-Abstract-Conference.html)|Relative usefulness ranking rather than useful/useless classification|Distant projection: compare an established gene-ranking objective with classification. Different items and labels; no reproduction of learned example retriever. [XGBoost LambdaMART implementation](https://xgboost.readthedocs.io/en/stable/tutorials/learning_to_rank.html) is the actual algorithm.|
|[MESA,2026 preprint](https://arxiv.org/abs/2608.10108), [JitMem,2026 preprint](https://arxiv.org/abs/2609.27334)|Select task-dependent evidence views; curate memory at use time|Projection: separate score-only feedback from score-plus-measured-examples on the current candidate menu. No trained memory router.|
|[LGBO,2026](https://arxiv.org/abs/2605.17976)|LLM preferences modify GP mean|Inspiration: test incremental LLM judgment over numerical selection. Do not rebrand generic mixing as its mechanism.|

A search result labelled arXiv2026, *Metric-agnostic Learning-to-Rank via Boosting and Rank Approximation*, is actually ICDM2023 according to its abstract page; excluded from the2025–26 source shortlist. Standard numerical algorithms remain legitimate baselines.

## Round1 fixed comparison

IFNG/IL2, seeds7101–7102,256 revealed measurements,32 new tests. Same100 depth3 trees,learning rate.05,L2=5. Factor2 feature views(full778-dimensional normalized Achilles versus16 standardized PCA coordinates) ×5 objectives: binary logistic; binary squared error; absolute-effect squared error; absolute-effect median quantile; binary LambdaMART NDCG ranking (one query containing the observed campaign genes, top32 pair sampling). Geometry and isotropic GP controls. The two squared-error arms isolate target choice; the binary logistic/squared/ranking arms isolate objective at fixed binary information. Rank scores and effect scores are NOT probabilities. All feature-covered untested genes remain eligible.

Expected usefulness: task-aligned target/ranking and less destructive compression might improve discovery. Weakening result: gains change sign across seeds/horizons or fail against fixed geometry. Rounds2–5 selected sequentially from actual results and logged in DECISIONS.md before execution. No exhaustive combined system.

## Validity

Model fit sees revealed labels/effects only. Unlabeled feature PCA is transductive. Predictions/selected genes saved before evaluator scoring. Fresh starting seeds reuse biological datasets and many genes; not independent held-out biology. Feature source lineage known but no independent rebuild/assay-overlap audit;32-gene short campaigns are not full BDA reproduction. Numerical baselines here are newly executed; archived inference is never described as a fresh intervention.
