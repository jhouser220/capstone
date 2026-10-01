> Historical research record, copied 2026-10-01. Some linked raw artifacts remain in the local archive. See the package README for included records and runnable entry points.

# Numerical objectives and feedback views: methods and results

2026-09-29. See [five-round summary](FIVE_ROUND_SUMMARY.md) for the main conclusions. All numerical fits and30 unique paid requests in this directory were newly executed. Resuming interrupted inference reused identical saved responses; these are not new independent interventions.

## 1. Design and paper-to-method mapping

| Source | Original mechanism | What we actually tested | Boundary |
|---|---|---|---|
| [GOLLuM, Nature Machine Intelligence2026](https://www.nature.com/articles/s42256-026-01283-z) | GP marginal likelihood optimizes representations for the downstream objective | Cross feature representation (778 normalized Achilles features vs16 PCA components) with target and loss; compare existing GP | Diagnostic projection. No trained GOLLuM encoder, no reproduction. |
| [Learning to Rank for In-Context Example Retrieval / SeDPO, NeurIPS2025](https://proceedings.nips.cc/paper_files/paper/2025/hash/cf68c8b43e692c63deec7b1242bbb3a0-Abstract-Conference.html) | Relative example usefulness trains a retrieval ranker | Rank candidate genes instead of estimating their class labels | Distant inspiration across components. Implemented standard [XGBoost LambdaMART](https://xgboost.readthedocs.io/en/stable/tutorials/learning_to_rank.html), not SeDPO. |
| [LLMNN, Findings EMNLP2025](https://arxiv.org/abs/2509.21403) | Separates LLM priors and numerical acquisition; uses feedback controls | Frozen/genuine/shuffled numerical updates; no-report/genuine/shuffled LLM reports | Experimental-control transfer. Our feature/model system is different. |
| [MESA,2026 preprint](https://arxiv.org/abs/2608.10108) | Learns selection/composition of memory structures | Observed-history expert selection and scores-versus-examples views | Projections; neither is the learned MESA policy. Generic cross-validation is not novel. |
| [Just-in-Time Memory,2026 preprint](https://arxiv.org/abs/2609.27334) | Builds task-specific memory at retrieval time from retained experience | For each current candidate, retrieve nearest already-observed hit and non-hit | Read-time curation projection; fixed similarity rather than trained retrieval. |
| [LGBO,2026 preprint](https://arxiv.org/abs/2605.17976) | LLM preferences guide the GP mean/acquisition process | Ask LLM to select from numerical advice and observed examples | Combination inspiration; direction of advice differs. No claim to reproduce LGBO. |

Sources were checked against primary pages. These experiments do not refute their original findings. We deliberately did not cite arXiv2604.15101 as a new2026 advance: its abstract identifies ICDM2023.

## 2. Shared setup and outcome separation

- IFNG and IL2, feature-covered candidate pool. Published screen hit labels and signed effects remain separate from model fitting.
- Initial256 observations: geometry-top64 plus192 random other genes, shuffled. Seeds change the random observations, but histories share top64 and the underlying biological dataset. They are not independent biological replications.
- One batch=32 distinct untested genes. Numerical campaigns run3 batches; unchosen candidates remain available. Initial256 observations do not count toward new discoveries. New hits are unique within a campaign; totals across campaigns can count the same gene more than once.
- Raw features are778 normalized Achilles profiles. PCA16 is fit without target labels on the full feature-covered pool, so this is transductive feature preprocessing. Cache provenance and hash are known; independent feature rebuild/assay overlap remain unverified. No feature causal interpretation is made.
- Trees:100 estimators, depth3, learning rate0.05, L2=5, no row/column subsampling, random state0. Objectives: binary logistic; squared error on binary labels; squared error on absolute effects; median quantile on absolute effects; rank:ndcg on binary labels in one query group. Ranking pair setting32 does **not** make this direct precision@32 optimization.
- Existing isotropic GP and geometry-centrality baselines are reused. The geometry baseline uses features but no measured target outcomes.
- Only revealed outcomes enter fitting. Predictions and chosen IDs are saved before evaluator scoring. Evaluation-only full-pool diagnostics never enter agent prompts. R4 validation outcomes are drawn from already-observed initial measurements.

## 3. Numerical results

Each R1 arm has128 choices; R2–R4 each384. Rows are matched within a round, not across rounds.

| R1 model | Total hits | IFNG | IL2 |
|---|---:|---:|---:|
| Geometry |22|10|12|
| GP |28|15|13|
| PCA binary logistic |23|12|11|
| PCA binary squared error |24|6|18|
| PCA effect squared error |29|10|19|
| PCA median effect |22|2|20|
| PCA binary ranking |24|11|13|
| Full binary logistic |27|11|16|
| Full binary squared error |17|5|12|
| Full effect squared error |21|5|16|
| Full median effect |40|10|30|
| Full binary ranking |29|15|14|

At fixed squared-error loss, continuous targets helped slightly in this one step (PCA29 vs24; full21 vs17). This is not evidence that every continuous objective is better. Full features help median strongly, but worsen effect squared error: representation and loss interact.

| R2 three-batch policy | Total | IFNG | IL2 |
|---|---:|---:|---:|
| Geometry |117|44|73|
| GP |105|54|51|
| PCA median |78|22|56|
| Full median |103|25|78|
| Full ranker |74|35|39|
| Full binary classifier |96|32|64|
| PCA effect squared error |63|17|46|

Full median's gain over compressed median persists, but its advantage over geometry does not. GP helps IFNG and loses IL2; full median shows the opposite aggregate relation in this round. We did not construct a post-hoc dataset-name oracle and call it an achievable method.

| R3 policy/history | Total | IFNG | IL2 | Four case totals: IFNG1,IFNG2,IL2_1,IL2_2 |
|---|---:|---:|---:|---|
| GP updated |108|48|60|24,24,28,32|
| GP frozen |88|25|63|22,3,31,32|
| GP shuffled new |102|37|65|22,15,32,33|
| Full median updated |90|29|61|10,19,35,26|
| Full median frozen |76|18|58|7,11,34,24|
| Full median shuffled new |71|22|49|11,11,30,19|

Frozen policies ignore new measurements but remove tested genes. Shuffled policies retain initial256 correct observations and jointly permute only new hit/effect pairs. The GP updated–frozen gain is concentrated in one IFNG case; its +6 versus shuffled is weaker than the tree's +19. Updated trajectories diverge, so these estimate whole-policy effects rather than same-state marginal treatment effects at every later step. Geometry was not rerun in R3; do not compare to R2 geometry as if paired.

R4 uses first192 observed genes for training and last64 for validation. Choose between GP and full median using top16 validation hits or absolute-effect MAE; freeze the choice before new experiments and refit that expert after every batch. Hit-based choice picked GP all4 states:113/384. MAE picked median all4:76/384. Geometry118. Routers reuse the corresponding deterministic expert trajectories, explicitly labeled in records; these are not extra independent campaigns. Validation MAE can favor a predictor with poor discovery utility. The small biased validation panel may also be unrepresentative; this does not establish that a better utility-based router is impossible.

## 4. LLM report views and stronger-model check

Four states: IFNG/IL2 ×7501/7502. A shared shuffled menu is the union of GPtop48, median-treeTop48 and geometryTop16 (at most112 genes). Both numerical top32 sets are included. No-report still receives this history-dependent shortlist, so it is not a completely feedback-free system.

Views:

1. No report: task and candidate names.
2. Scores: GP hit probability and predicted effect magnitude, tree median magnitude, both rank percentiles; observed validation hits/MAE/Spearman. Probability is uncalibrated; ranks/effect predictions are not probabilities.
3. Examples: same scores plus each candidate's nearest observed hit and non-hit, their signed effects and similarities. Both successes and failures are represented.
4. Shuffled scores: jointly permute per-candidate score fields, retaining real validation summary. This is a **partial** information control, not fully meaningless feedback.

This batch tests selection from supplied analyses, not autonomous tool invocation, live information acquisition, or adaptive probability blending. Content and prompt length vary together; do not label this a pure representation/format effect. Extra examples are not guaranteed to be counterexamples to any stated hypothesis.

| Model/view | IFNG7501 | IFNG7502 | IL2_7501 | IL2_7502 | Total hits/tests |
|---|---:|---:|---:|---:|---:|
| Lite no report |5|5|12|9*|31/128 (24.22%)|
| Lite scores |10|10*|8*|7*|35/128 (27.34%)|
| Lite examples |8|5|7|10|30/128 (23.44%)|
| Lite shuffled |8*|10*|9*|invalid|27/96; incomplete arm|
| Astra no report |9|—|6|—|15/64 (23.44%)|
| Astra scores |8|—|5|—|13/64 (20.31%)|
| Astra examples |9|—|7|—|16/64 (25.00%)|
| GP alone |5|11|4|9|29/128 (22.66%)|
| Tree alone |5|5|13|10|33/128 (25.78%)|

*One count-correction request used, without new outcomes. Eight of16 initial Lite responses had incorrect counts/duplicates;7 repairs succeeded and1 remained invalid. Astra required no repairs. Original outputs, requests, corrections and their charges are retained. Do not impute the invalid arm, count it as zero biological hits, or compare its96 choices with other arms'128 choices.

**Matched controls and sensitivity:**

- On the three complete shuffled-control states, genuine scores28 vs shuffled27 vs no report22 vs examples20 per96. The apparent scores-versus-no-report benefit therefore does not establish meaningful evidence use.
- Removing every state with a repair in either Lite no-report or scores leaves only IFNG7501: scores10 vs no report5/32. This single case is too weak to confirm a benefit. All Lite shuffled states required repair, so there is no repair-free genuine-versus-shuffled comparison.
- For examples versus no report, excluding the repaired no-report IL2_7502 state gives examples20 vs no report22/96.
- Same two states as Astra: Lite no-report17, scores18, examples15/64; numerical GP9, tree18. Astra15/13/16. This tiny sample does not rank general model capability. Astra used medium reasoning, Lite disabled reasoning: model and reasoning budget vary together.
- Relative to each state's GP draft, Lite scores added23 hits and discarded17 (net+6); examples added19/discarded18 (net+1); no-report25/23 (net+2). These are set differences, not observed interactive revisions. Astra examples13/6 (net+7), but no-report15/9 (net+6): only+1 incremental report benefit.

## 5. Metrics, audit and spending

- Primary: unique new hits within a campaign, hit rate=hits/new tests. `METRICS.json` includes mean cumulative hits at batch endpoints (discrete discovery AUC; not normalized or a continuous integral). Rankings/effect magnitudes are not treated as calibrated probabilities.
- Saved diagnostics include Spearman versus absolute effect and full-pool evaluation Brier/log loss for probability-like outputs. Binary squared-error scores are clipped for those diagnostic metrics; they are not calibrated. Full-pool outcomes are evaluator-only. These are not selective-sampling-corrected estimates of deployment calibration.
- Audit passed240 numerical selection checks,120 campaign-record checks (including shared deterministic router records),16 shuffled-pair preservation checks,22 hidden-outcome-mutation checks, and the additional valid LLM/menu and observed-example checks in `AUDIT.json`. Mutation tests cover the new tree/selection code; this is not an independent provenance audit of Achilles features.
- API prices checked before inference: Lite $0.10/$0.40 per million input/output tokens; Astra $10/$50, with any actual caching reflected in charges. Hard new-batch reservation cap$2.50, conservative project working stop$40 to preserve reserve below$50. No unresolved new calls.
- Lite24 calls (16 original+8 correction):$0.01622037. Astra6 calls:$0.68814250. **Total$0.70436287**, provider account usage difference matches. Entire recorded project actual$4.69536828; conservative$4.70581388 including two older unresolved reservations. API spend is not local compute cost.
- A postprocessing import stalled; finalization was changed to use saved arrays and a lightweight read-only audit. This made no additional inference calls and did not change outputs.

## 6. What this supports and the next focused tests

**Method analysis first:** the strongest candidate is a controlled account of how feature compression and prediction/validation objectives alter the value of numerical feedback across tasks and horizons. R3 is positive evidence, not another blanket negative. R2 and R4 show why lower prediction error or an attractive one-step result can still fail to yield more discoveries.

**Expected effectiveness versus demonstrated effectiveness:** full-feature median and numerical updating have plausible expected utility and local measured gains. No broadly superior new selection method, learned router, or useful rich LLM report has been demonstrated. All components are established methods or simple projections; novelty would require identifying a general mechanism beyond this small dataset-specific case study.

**Next, without a large combined system:**

1. Freeze the GP/full-median/geometry comparison; repeat using random-only initial histories as well as the current geometry-enriched histories. Test whether task differences are caused by observation composition or persist beyond it. Source: LLMNN prior/feedback separation + GOLLuM task-dependent representation projection. Weakening result: gains vanish outside the enriched histories.
2. Compare same-sized observed validation panels and global-error versus top-batch discovery criteria, recording uncertainty in their expert ordering. Use a disjoint retrospective panel acquired at equal experimental cost, and retain a no-switch baseline. Do not use hidden full-pool labels to choose the expert. Source: SeDPO relative-utility inspiration + MESA conditional selection projection. Generic model selection remains a control; novelty is unestablished.
3. Only if the score-report signal warrants confirmation: enforce exactly32 unique IDs with validated structured output, compare genuine and shuffled scores with equal evidence length, then run matched continuations. Source: LLMNN evidence controls + JitMem/MESA evidence-view projections. The current shuffled-control result is a reason to test information dependence, not to add more critics.

No new independent dataset, full128×5 BDA campaign, or publication-ready claim was produced in this batch. Fresh random seeds do not erase adaptivity from repeatedly studying IFNG/IL2. Reserve a genuinely held-out setting for a pre-specified final comparison.
