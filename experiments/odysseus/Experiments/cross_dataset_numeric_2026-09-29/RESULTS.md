> Historical research record, copied 2026-10-01. Some linked raw artifacts remain in the local archive. See the package README for included records and runnable entry points.

# Cross-dataset numerical learning: results and interpretation

Completed 2026-09-29T23:04:01.236680+00:00. [Summary](FIVE_ROUND_SUMMARY.md) · [Protocol](PROTOCOL.md) · [Adaptive decisions](DECISIONS.md).

## Scope and paper connections

We extended the existing GP/tree experiments to every available compatible single-gene table: seven tables from six screen families, with related Sanchez day21/day30 measurements. Horlbeck's gene-pair task is excluded because a single-gene feature/predictor interface is not a valid pairwise baseline. Nothing is claimed about Horlbeck or general biological-agent performance.

[BioDiscoveryAgent, ICLR2025](https://arxiv.org/abs/2405.17631) supplies the released task/outcome setting and multi-screen motivation. [LLMNN, EMNLP2025](https://arxiv.org/abs/2509.21403) separates priors from numerical acquisition and uses feedback controls; our transfer is to compare fixed priors, learned selections and matched information states. [GOLLuM, NMI2026](https://www.nature.com/articles/s42256-026-01283-z) trains representations through a GP objective; our much smaller projection tests fixed feature compression, target definitions and acquisition across tasks. We did not reproduce either learned hybrid, train an LLM encoder, or establish method novelty.

## Definitions and shared design

- Geometry means ranking by centrality in normalized Achilles gene profiles. It is a feature-based, label-free baseline, not the published BDA agent.
- GP fits absolute observed scores using an isotropic GP on16 standardized PCA components, then ranks by probability of exceeding a threshold estimated only from observed hit/non-hit data. This is a model-based score, not a calibrated hit probability. GP mean ranks by predicted absolute score using the same fit.
- Median trees predict absolute observed scores; PCA and full versions use16 versus778 feature dimensions. Binary XGBoost predicts released hit labels. Hyperparameters are fixed as in the preceding batch. A one-class observed history gets constant classifier scores rather than an invalid fit; tied candidates follow fixed feature-file order. Some later classifier predictions are constant too. This is a limitation of the tested implementation, not evidence that every classifier fails with scarce hits.
- Phase1 has256 observations (64 geometry-prior+192 random), phase2 has256 uniformly random observations. Three seeds8101–8103 per table; three new32-gene batches with own-history updates. Full feature-covered untested pools remain available. Methods within a condition start from identical histories. Reported totals below are **hits among288 new choices per method per dataset**, summing three campaigns. Genes can repeat across seeds; they cannot repeat within a campaign.
- Initial observed hits are excluded from discovery totals. `SUMMARY.json` also records per-seed gains, hit rates, recall among initially remaining hits, and mean cumulative hits at batch endpoints (discrete discovery AUC). Initial history acquisition contributes to total experimental budget:352 measurements per main campaign,112 for phase4. We do not compare the two history-size conditions as equal-total-budget competitors.
- Preserve original released hit labels. Task prose alone does not define the oracle: Carnevale hits are both tails despite a directional task description; Sanchez scores are statistical RSA_Down values, and its two tables differ in day and hit definition. CAR-T numeric column headers are handled explicitly without dropping its first record. See [verified local EDA](../../EDA/SOURCE_RECORD.md).

## Round1: transfer with the previous starting-history design

| Dataset | Geometry baseline | Random | GP probability | PCA median | Full median | Binary classifier |
|---|---:|---:|---:|---:|---:|---:|
| IFNG | 64 | 19 | 54 | 27 | 53 | 65 |
| IL-2 | 108 | 10 | 80 | 89 | 100 | 100 |
| Carnevale | 45 | 14 | 45 | 15 | 19 | 24 |
| Sanchez day21 | 21 | 12 | 18 | 10 | 16 | 15 |
| Sanchez day30 | 36 | 17 | 16 | 9 | 23 | 16 |
| Scharenberg | 19 | 10 | 32 | 32 | 36 | 52 |
| CAR-T | 3 | 4 | 1 | 1 | 2 | 3 |


The previous IFNG/IL2 generalization is weakened. IFNG GP cases26,0,28 versus geometry21,21,22 show that two improvements coexist with a catastrophic failure. IL2 GP loses in all three mixed-history cases. Scharenberg classifier wins all three (19,16,17 vs6,6,7). Carnevale GP's45 total ties geometry but consists of17,10,18 vs15,15,15. The Sanchez models generally fail to beat geometry in this condition. CAR-T's tiny counts cannot establish useful relative superiority.

The full-feature median improves over PCA median on most screens here, but is still worse than geometry on many. Representation gain relative to a weaker learner does not imply a competitive discovery method.

## Round2: random-only starting histories

| Dataset | Geometry baseline | GP probability | Full median | Binary classifier |
|---|---:|---:|---:|---:|
| IFNG | 47 | 53 | 44 | 28 |
| IL-2 | 76 | 64 | 96 | 76 |
| Carnevale | 26 | 42 | 16 | 19 |
| Sanchez day21 | 20 | 26 | 29 | 14 |
| Sanchez day30 | 18 | 31 | 27 | 21 |
| Scharenberg | 29 | 48 | 43 | 49 |
| CAR-T | 0 | 4 | 2 | 3 |


Carnevale GP wins all three paired cases (15,17,10 vs9,8,9); Sanchez day30 GP also wins all three (10,8,13 vs6,6,6). IL2 full median wins all three (38,27,31 vs25,25,26). Scharenberg classifier49vs29 remains positive. These are real local improvements, not a story in which nothing worked.

However, changing initial observations also removes different genes from the remaining action pool. Raw differences between phase1 and phase2 therefore combine information effects, candidate depletion and stochastic history variation. Fresh seeds are not independent biological datasets. No dataset-name routing rule was validated.

## Round3: does a different GP selection rule repair the failures?

The following are totals per288 choices. GP probability controls are reused from phase1; GP-mean campaigns are newly run with their own subsequent histories.

| Dataset | GP probability | GP mean |
|---|---:|---:|
| IFNG | 54 | 58 |
| IL-2 | 80 | 81 |
| Carnevale | 45 | 34 |
| Sanchez day21 | 18 | 18 |
| Sanchez day30 | 16 | 14 |
| Scharenberg | 32 | 31 |
| CAR-T | 1 | 1 |


IFNG8102's failed probability policy had a first-step score range only0.0599164–0.0606911. Scores were not all equal or underflowed to zero. Mean ranking changes the complete-campaign yield from0 to1; it does not rescue the failure. It also loses on Carnevale. This narrows the diagnosis beyond the acquisition rule: fitting, representation, target and observed-state support need examination. We have not isolated which causes the near-flat fit.

## Round4: was the small screen's advantage due to much more observed coverage?

| Scharenberg random starting history | Geometry | GP | Binary classifier | Initial covered-pool fraction |
|---|---:|---:|---:|---:|
|256 observations|29|48|49|25.42%|
|16 observations|31|24|29|1.59%|

Each row counts288 new choices per method. With16 observations the initial hit counts were1,0,2; the classifier no longer beats geometry. This supports a dependence on available training support. It does **not** identify coverage fraction as the sole cause: observation count, positive examples, composition and remaining candidates change together. Total measured budgets differ intentionally. More support and a smaller candidate pool are plausible explanations for part of the cross-screen pattern; this elementary learning-curve phenomenon is not new by itself.

## Round5: hold the remaining action pool identical

Adaptive one-step replay on three screens: each model sees either the mixed256 or random256 history; both arms select from the same pool excluding the union of both histories. Each cell is hits/96 (three32-choice batches). Geometry is identically reused across the two history labels, not independently rerun.

| Dataset | Model | Mixed history | Random history |
|---|---|---:|---:|
| IL-2 | Geometry | 19 | 19 |
| IL-2 | GP | 18 | 18 |
| IL-2 | Full median | 36 | 17 |
| Carnevale | Geometry | 5 | 5 |
| Carnevale | GP | 18 | 16 |
| Carnevale | Full median | 8 | 3 |
| Sanchez day30 | Geometry | 3 | 3 |
| Sanchez day30 | GP | 3 | 12 |
| Sanchez day30 | Full median | 5 | 13 |


History composition matters even with identical available actions. For IL2, the median model does better with the mixed history (36vs17, paired differences+1,+6,+12); for Sanchez day30 the GP does better with random history (12vs3, +6,+1,+2). The strongest effect is not always in the same direction or the same model. This is a concrete next method-analysis question: which already-observed evidence is useful for this particular candidate pool and predictor? It is not yet a learned selection rule or a full-campaign gain.

## Dataset properties that help interpret the differences

| Dataset | Covered unique genes | Covered hits | Hit prevalence | Initial256 / pool | Absolute score versus hit AUC |
|---|---:|---:|---:|---:|---:|
| IFNG | 16794 | 857 | 5.10% | 1.52% | 1.000 |
| IL-2 | 17457 | 628 | 3.60% | 1.47% | 1.000 |
| Carnevale | 16921 | 849 | 5.02% | 1.51% | 1.000 |
| Sanchez day21 | 17073 | 849 | 4.97% | 1.50% | 0.497 |
| Sanchez day30 | 17073 | 882 | 5.17% | 1.50% | 1.000 |
| Scharenberg | 1007 | 48 | 4.77% | 25.42% | 1.000 |
| CAR-T | 17309 | 139 | 0.80% | 1.48% | 1.000 |


**The final column is an evaluator-only description of score/label alignment, not predictive model performance.** It uses known measured scores and released labels over the full covered screen. Near1 means larger absolute observed scores align with hit membership. Sanchez day21 is near0.5 because its labels cover both tails of a nonpositive statistical score: large absolute magnitude is not an adequate monotonic utility proxy. This identifies a target-definition problem in transferring our absolute-score regression, not an intrinsically unlearnable biological task. Binary classification is an aligned control but was not sufficient to win here.

CAR-T has about0.80% hits; in mixed phase1 initial histories it had1,0,0 known hits, compared with16–17 on Scharenberg. This makes its feedback regime qualitatively different. The two additional CAR-T hits in one arm are not persuasive method evidence.

Feature-neighbour diagnostics (512 uniformly sampled query genes, nearest10 other genes) are saved in `DATASET_PROPERTIES.json`. They are evaluator-only and were never used for selection. They suggest clustering for some screens, but only2 CAR-T hits appeared in the probe; do not infer a universal biological geometry law or causality from this small descriptive sample.

## Integrity, deviations and limits

-240 newly executed three-batch campaigns;54 additional one-step comparison records, including9 duplicated deterministic geometry records. No LLM/API requests. Resume reused completed campaign files rather than relabeling them fresh runs.
-Input validation stopped before Scharenberg because six gene IDs had conflicting duplicate scores: PSMA1, TM9SF1, PI4K2A, TLR9, CLN3, IDS. All six were excluded from candidate/history pools in every arm. The feature-covered pool is1007 unique unambiguous genes, with48 hits. Original data/hit files are unchanged; neither score was chosen opportunistically. Published full-screen results are not directly comparable.
-Numerical loader/runtime issues were resolved in a temporary environment with the same pinned package versions; installation records and execution snapshots are retained. Original stopped logs and first protocol snapshots are preserved. This changes no outcomes or hyperparameters.
-Audit passed774 selected-batch checks,8 hidden-outcome-mutation checks, common-query/identical-geometry checks for phase5 and unchanged input-file hashes. Training receives only revealed histories. Candidate effects/full-screen diagnostics remain evaluator-only.
-Achilles feature cache provenance is known but independent rebuild and assay-overlap validation remain open. The feature-covered subsets exclude some original genes/hits. PCA uses unlabeled full-pool features, so preprocessing is transductive. No test label was used to learn the feature projection.
-Three seeds per condition are exploratory. Follow-ups were selected adaptively; no independent biological replication, multiplicity-corrected statistical claim or novelty claim is made. These screens, including Carnevale, are now development/exploration data; choose an audited genuinely unseen study or protocol for final confirmation. The prior blanket suggestion to treat Carnevale as untouched is stale.
-No paid API cost: $0 new, conservative project total$4.70581388/$50. Local compute is not priced as API spending. No stronger-model check was needed to answer this numerical cross-screen question; previous frontier checks remain in their own records.

## Research decision

**Most useful analysis lead:** feedback usefulness depends on observed support, objective alignment and the candidate pool—not merely dataset name or whether a model receives more history. The fixed-pool replay is a sharper test than comparing different campaigns without controlling the available genes.

**Best next controlled tests:**

1. On identical remaining pools, vary history size and hit/non-hit support independently where feasible; preserve total acquisition budgets for policy comparisons. Learn whether an observable support signal predicts when to trust numerical learning. Inspired by LLMNN's prior/acquisition separation, this would need a new, validated mechanism to become a method paper.
2. Correct the Sanchez day21 utility representation: model the two hit tails rather than assume absolute statistical score equals desired biological effect. Estimate any boundaries from revealed data; retain binary classification and geometry controls. This is objective alignment motivated by GOLLuM's task-relevant representation principle, not an implementation of it or an intrinsically novel method.
3. Diagnose GP fit stability on failed histories before adding a complex LLM mixture. Posterior-mean acquisition alone failed. Check alternative fit initializations and calibrated/binary surrogates on fixed histories, then campaigns; do not optimize on hidden outcomes.

A support-aware combination of feature priors and numerical learners is plausible, but **not demonstrated** here. Better performance on both IFNG and IL2, or on every screen, remains open. No new named framework is warranted yet.
