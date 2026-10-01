> Historical research record, copied 2026-10-01. Some linked raw artifacts remain in the local archive. See the package README for included records and runnable entry points.

# Rich numerical feedback: methods, controls and complete results

## Mechanisms actually implemented

- **Geometry baseline:** normalized Achilles gene-profile centrality, no outcome fitting. This is not an LLM prior. kNN uses20 nearest observed raw profiles and Laplace-smoothed hit frequency.
- **XGBoost:**100 depth3 trees, learning rate.05, L2 penalty5; binary classifier versus10/50/90% absolute-effect quantile regressor. Selection uses classifier probability or effect median respectively. This changes both target and loss, so it does not isolate the value of continuous labels by itself.
- **GP:** standardized16-dimensional PCA of unlabeled Achilles profiles; constant×RBF+noise; observed absolute effects. Isotropic versus16 ARD length scales, both optimized through marginal likelihood, max35 iterations. ARD learns a diagonal distance metric; it is not GOLLuM's learned linear projection or LLM fine-tuning. No assertion of optimizer convergence. Hit probability is a Gaussian tail above a threshold inferred from observed hit/non-hit effects, a misspecified approximation when cutoffs are asymmetric.
- **Reports:** GP probability, kNN probability, tree probability, tree effect quantiles, nearest observed hit AND non-hit (gene, raw signed effect, label, similarity), and top2 TreeSHAP contributions in PCA-coordinate log odds. Attributions describe the fitted classifier, not biological causes. All candidates receive the same types of information. No hidden candidate result is supplied.
- **Mixtures:** fixed50/50 or softmax of negative Brier errors (temperature coefficient10). R3 calibrates using32 masked old observations. R5 starts50/50 and updates only on outcomes obtained after saved predictions; these selected-sample errors can be biased for future menus. This is an inspiration from LGBO, not its GP-mean intervention.
- **Acquisition:**24 highest GP probabilities plus8 maximum latent posterior variances, versus8 greedily maximizing0.5 log(1+conditional latent variance/noise). The latter conditions covariance on earlier batch choices without knowing their future outcomes. Exact under the fitted Gaussian model, not a guarantee about biology or information quality. Compare32 greedy at the same test budget.
- **Other executed tools:** Spearman feature/effect associations with30 bootstrap resamples; feature/effect mutual information; grouped permutation importance on64 masked already-observed genes with20 permutations/group. Results in OBSERVED_HISTORY_TOOL_DIAGNOSTICS.json. These diagnostics did not choose the tested arms. Tree contributions were included in LLM reports; MI/permutation results were not.

## Metrics and limitations

Hit rate=hits/new gene tests (precision at the chosen batch size). Cumulative unique hits counts each selected gene once per campaign. Discovery AUC is the mean cumulative-hit count at the3 checkpoints. Brier/log loss measure probability quality; Spearman measures rank agreement with absolute effects. Nominal80% interval coverage/width diagnose uncertainty. Added/removed hits compare each LLM choice with the GP proposal in that same state; cumulative arm differences also include changed later histories.

First-round mean pool diagnostics: XGBoost classifier Brier0.0406, GP0.0428, kNN0.0687. The classifier's better global probability score did not make it the best selector. Effect XGB80% intervals cover69.1%; GP94.8%; ARD91.5%. These are evaluator-only whole-unrevealed-pool diagnostics, not validation supplied to agents. Probability errors on mostly non-hits need not predict top32 utility.

All numerical acquisition sees the entire feature-covered pool and retains unselected genes. LLM ranking sees a shuffled union of GPtop48 and geometrytop48; this is still a shortlist experiment, not unrestricted BDA gene generation. The no-report LLM has no explicit observed outcomes, but its shortlist already depends on history through the GP. Its comparison isolates extra report value, not history-free system behavior.

Six/four small matched states, one completion per model-state, overlapping gene pools,256 initial observations and32×3 campaigns limit generalization. Fresh seeds are not new biological datasets. No new Carnevale result. Achilles source lineage is documented but independent rebuild and assay overlap remain unverified. No wet-lab outcomes, no training of LLM weights, no full BDA reproduction. R5 Flash uses the exact R3 cases intentionally for a capacity check, not an independent biological replication.

## Recovery and output validity

An initial list-format response supplied94 values for96 candidates and was rejected; the old calibration and failed response remain archived and charged. Switched to explicit candidate-ID maps before scoring candidate outcomes. One later response had a single stray quote after an unquoted number; this was repaired mechanically without changing any numeric value, with a flag. Numeric strings in another response were converted to floats with a count recorded. Exact IDs, ranges and finite values are validated. A local NumPy-int JSON serialization error was fixed; identical already-paid responses were reused on restart, not counted as new experiments. API calls are unique ledger records. Protocols are execution records, not formal preregistration; source files/configs/logs and outputs are retained.

## Source fidelity

See [protocol](PROTOCOL.md). GOLLuM author README and GP implementation were inspected and archived at commit c418d7ed3c17e5995f503f6df7c26e9f7c58d09f. Its BoTorch/GPyTorch training code was not run. This pilot uses sklearn/XGBoost implementations to isolate inexpensive numerical mechanisms. Source algorithms are not refuted by a failed projection.

## r1_full_pool_models

| Method | Hits / selections | Hit rate | IFNG / IL2 | Paired Δ vs prior | Win/tie/loss |
|---|---:|---:|---|---|---|
|gp_ard|41/192|21.35%|{'IFNG': 13, 'IL2': 28}|[-1, -3, 2, 4, 4, 2]|[4, 0, 2]|
|gp_iso|46/192|23.96%|{'IFNG': 22, 'IL2': 24}|[-3, 5, 5, 1, 3, 2]|[5, 0, 1]|
|knn|13/192|6.77%|{'IFNG': 7, 'IL2': 6}|[-2, -2, -4, -4, -3, -5]|[0, 0, 6]|
|prior|33/192|17.19%|{'IFNG': 15, 'IL2': 18}|[0, 0, 0, 0, 0, 0]|[0, 6, 0]|
|xgb_effect|47/192|24.48%|{'IFNG': 15, 'IL2': 32}|[-2, -1, 3, 7, 2, 5]|[4, 0, 2]|
|xgb_hit|29/192|15.10%|{'IFNG': 15, 'IL2': 14}|[2, -3, 1, -2, 0, -2]|[2, 1, 3]|

## r2_own_history_campaigns

| Method | Hits / selections | Hit rate | IFNG / IL2 | Paired Δ vs prior | Win/tie/loss |
|---|---:|---:|---|---|---|
|gp_ard|62/384|16.15%|{'IFNG': 29, 'IL2': 33}|[-13, -2, -22, -17]|[0, 0, 4]|
|gp_iso|115/384|29.95%|{'IFNG': 51, 'IL2': 64}|[4, 3, -1, -7]|[2, 0, 2]|
|knn|42/384|10.94%|{'IFNG': 14, 'IL2': 28}|[-17, -13, -24, -20]|[0, 0, 4]|
|prior|116/384|30.21%|{'IFNG': 44, 'IL2': 72}|[0, 0, 0, 0]|[0, 4, 0]|
|xgb_effect|68/384|17.71%|{'IFNG': 18, 'IL2': 50}|[-9, -17, -12, -10]|[0, 0, 4]|
|xgb_hit|75/384|19.53%|{'IFNG': 20, 'IL2': 55}|[-13, -11, -5, -12]|[0, 0, 4]|

## r3_llm_numeric_reports

| Method | Hits / selections | Hit rate | IFNG / IL2 | Paired Δ vs numeric | Win/tie/loss |
|---|---:|---:|---|---|---|
|calibrated_blend|25/128|19.53%|{'IFNG': 16, 'IL2': 9}|[4, 0, 0, 0]|[1, 3, 0]|
|fixed_blend|25/128|19.53%|{'IFNG': 16, 'IL2': 9}|[4, 0, 0, 0]|[1, 3, 0]|
|llm_prior|27/128|21.09%|{'IFNG': 12, 'IL2': 15}|[4, -4, 0, 6]|[2, 1, 1]|
|llm_report|29/128|22.66%|{'IFNG': 16, 'IL2': 13}|[4, 0, 0, 4]|[2, 2, 0]|
|numeric|21/128|16.41%|{'IFNG': 12, 'IL2': 9}|[0, 0, 0, 0]|[0, 4, 0]|

## r4_information_acquisition

| Method | Hits / selections | Hit rate | IFNG / IL2 | Paired Δ vs gp_iso | Win/tie/loss |
|---|---:|---:|---|---|---|
|gp_eig8|74/384|19.27%|{'IFNG': 30, 'IL2': 44}|[-19, -11, -8, -6]|[0, 0, 4]|
|gp_iso|118/384|30.73%|{'IFNG': 60, 'IL2': 58}|[0, 0, 0, 0]|[0, 4, 0]|
|gp_uncertainty8|74/384|19.27%|{'IFNG': 30, 'IL2': 44}|[-19, -11, -8, -6]|[0, 0, 4]|

## r5_llm_campaigns

| Method | Hits / selections | Hit rate | IFNG / IL2 | Paired Δ vs numeric | Win/tie/loss |
|---|---:|---:|---|---|---|
|llm_prior|102/384|26.56%|{'IFNG': 47, 'IL2': 55}|[-6, -4, -1, 4]|[1, 0, 3]|
|llm_report|103/384|26.82%|{'IFNG': 49, 'IL2': 54}|[-3, -5, 0, 2]|[1, 1, 2]|
|numeric|109/384|28.39%|{'IFNG': 57, 'IL2': 52}|[0, 0, 0, 0]|[0, 4, 0]|
|online_blend|101/384|26.30%|{'IFNG': 49, 'IL2': 52}|[-4, -4, 0, 0]|[0, 2, 2]|
|prior|116/384|30.21%|{'IFNG': 43, 'IL2': 73}|[-9, -5, 10, 11]|[2, 0, 2]|

## r5_flash_spot

| Method | Hits / selections | Hit rate | IFNG / IL2 | Paired Δ vs numeric | Win/tie/loss |
|---|---:|---:|---|---|---|
|calibrated_blend|25/128|19.53%|{'IFNG': 16, 'IL2': 9}|[4, 0, 0, 0]|[1, 3, 0]|
|fixed_blend|25/128|19.53%|{'IFNG': 16, 'IL2': 9}|[4, 0, 0, 0]|[1, 3, 0]|
|llm_prior|30/128|23.44%|{'IFNG': 15, 'IL2': 15}|[4, -1, 0, 6]|[2, 1, 1]|
|llm_report|31/128|24.22%|{'IFNG': 16, 'IL2': 15}|[4, 0, 0, 6]|[2, 2, 0]|
|numeric|21/128|16.41%|{'IFNG': 12, 'IL2': 9}|[0, 0, 0, 0]|[0, 4, 0]|

Full accounting: FINAL_ACCOUNTING.json. All per-case predictions, prompts, responses, histories and decisions are retained in stage folders and calls/.

## Final validity and acquisition diagnostics

144 saved numerical selection records and60 own-history decision records passed eligibility/uniqueness checks. Ten hidden-outcome mutation checks left model predictions and selections identical.6314 supplied neighbouring measurements were checked against the corresponding revealed histories and labels/effects. This tests leakage boundaries in these implementations, not upstream feature provenance.

Each exploration arm found0 hits among its96 exploratory selections, while its288 greedy selections found74. The GP assigned those exploratory candidates about8.9% mean hit probability versus21.4% for the greedy portion. This suggests a poorly targeted/miscalibrated exploration region in this setup; overlapping genes and adaptive sampling preclude treating96 as independent Bernoulli trials. No evidence that these exploratory observations paid back within three rounds.
