> Historical research record, copied 2026-10-01. Some linked raw artifacts remain in the local archive. See the package README for included records and runnable entry points.

# Rich numerical feedback: sequential pilot

2026-09-29. Method analysis first; all novelty claims provisional.

## Question and fidelity

Do richer predictors make better discoveries, and does an LLM add useful judgment beyond their scores? The prior restricted-menu improvements are motivation, not controls for this new experiment.

| Source | Original mechanism | Our move | Expected benefit / weakening result |
|---|---|---|---|
| [LLMNN, Findings EMNLP 2025](https://arxiv.org/abs/2509.21403) | LLM prior coupled to nearest-neighbor acquisition; tests shuffled outcomes | Experimental-design transfer: numerical baselines, genuine-versus-no-report LLM decisions, same candidates | LLM may supply complementary biological prior; no incremental gains weakens that role |
| [GOLLuM, Nature Machine Intelligence 2026](https://www.nature.com/articles/s42256-026-01283-z), [author code](https://github.com/schwallergroup/gollum) | GP marginal likelihood trains a projection and/or language encoder | Limited projection: optimize GP automatic relevance determination (ARD) on 16 fixed PCA components of Achilles profiles; compare isotropic GP | Task-dependent distances may improve ranking. This is NOT GOLLuM reproduction or a trained LLM embedding. Failure to beat isotropic GP weakens this inexpensive adaptation |
| [XGBoost quantile regression](https://xgboost.readthedocs.io/en/stable/parameter.html) | Gradient boosted tree classification and conditional quantiles | Established baselines: hit classification versus absolute effect quantiles (10/50/90%) | Continuous non-hit measurements may help; classification winning weakens this particular continuous objective, not all continuous feedback |
| [LGBO, 2026](https://arxiv.org/abs/2605.17976) | LLM preferences modify GP mean while retaining covariance | Inspiration: compare numerical probabilities, LLM probabilities, fixed and observed-error-weighted blends | Complementarity must beat each component; our blending is not the original GP-mean mechanism |
| [BATCHIE, Nature Communications 2025](https://www.nature.com/articles/s41467-024-55287-7) | Bayesian information-based batch experimental design for drug combinations | Conditional projection: uncertainty/diversity acquisition if useful numerical predictions exist | Better later discoveries, possibly lower immediate yield; variance alone must not be called expected information gain |
| [Confidence vs. Critique, ACL 2025](https://aclanthology.org/2025.acl-long.203/) | Separates correction from preservation | Diagnostic transfer: hits added and discarded by each revision | Improved factual reporting alone does not establish better choices |

## Frozen first comparison

Six matched states: IFNG and IL2, seeds6101–6103.256 observed genes (64 geometry-prior +192 random).32 new selections per model from ALL feature-covered untested genes, retaining every unselected candidate. Models: geometry-only prior, k20, XGBoost hit classifier, XGBoost median absolute-effect predictor plus quantile intervals, isotropic GP and ARD GP. All models except raw-feature kNN/geometry share16 standardized unsupervised PCA coordinates. PCA uses unlabeled candidate features; supervised fitting only receives revealed measurements. No target-specific hyperparameter search in this round. GP Gaussian tail probabilities use a threshold estimated from observed hit/effect pairs; this can be misspecified, especially with asymmetric screen cutoffs. Report this separately from effect ranking.

Round2 will be chosen after round1, emphasizing own-history campaigns and retention of available genes. Round3–4 will test fresh LLM inference with numerical reports and component/blend controls; exact variants frozen before each round. Round5 will confirm or diagnose using fresh seeds. No attempt to implement every proposed tool at once.

## Accounting and validity

$1 new-batch cap; conservative project stop at$40, leaving reserve below$50. Refreshed OpenRouter prices and prior ledgers in PRIOR_ACCOUNTING.json. Unresolved old calls remain reserved. Predictions and selections saved before evaluator outcomes are consulted. Hit sums across seeds count discoveries per run, NOT globally distinct genes. Experimental budget excludes the shared256-gene initial history; total is256+32×rounds per campaign. Small retrospective campaigns are not published BDA benchmark reproduction. Feature cache hash pinned; source lineage known but independent rebuild and assay-overlap audit remain unresolved. No causal biological explanations from SHAP/PCA importance.
