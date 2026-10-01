> Historical research record, copied 2026-10-01. Some linked raw artifacts remain in the local archive. See the package README for included records and runnable entry points.

# Five sequential rounds: numerical objectives, representations and feedback views

Completed 2026-09-29. These are new experiments following `rich_numeric_feedback_2026-09-29` and `frontier_spot_2026-09-29`; their old results are not counted again. Method analysis remains the first goal.

**Main finding:** numerical feedback can help, but improving average prediction error, adding examples, or using a stronger LLM does not automatically improve discovery. No general new-method advantage is established.

| Round/question | Source papers and research move | Comparison | Actual result | API cost | Evidence strength | Decision |
|---|---|---|---|---|---|---|
| 1. Does the prediction objective interact with feature compression? | [GOLLuM](https://www.nature.com/articles/s42256-026-01283-z): task-dependent representation projection; [SeDPO](https://proceedings.nips.cc/paper_files/paper/2025/hash/cf68c8b43e692c63deec7b1242bbb3a0-Abstract-Conference.html): distant ranking-utility inspiration | Same XGBoost capacity; five objectives × full778 features/PCA16, plus GP and geometry | Full-feature median40/128 vs compressed22/128; geometry22, GP28. Gain mostly IL2. | $0 | 1/5 exploratory | Advance full/compressed median and controls to campaigns. |
| 2. Does the gain persist over three batches? | GOLLuM projection; [LLMNN](https://arxiv.org/abs/2509.21403) prior-versus-acquisition separation | Own-history campaigns, fresh starting seeds, same budget | Full median103/384 vs compressed78, but geometry117 and GP105. Full median IL2:78 vs geometry73; IFNG:25 vs44. Direct ranker74. | $0 | 2/5 preliminary | Study task dependence rather than advertise one superior predictor. |
| 3. Is new experimental feedback actually useful? | LLMNN feedback-dependence control transferred to numerical learners | Updated vs frozen vs jointly shuffled new outcomes | GP108 vs88 vs102/384; full median90 vs76 vs71/384. **Real positive numerical feedback signal**, uneven across cases. | $0 | 2/5 preliminary | Check whether observed validation can identify the useful expert. |
| 4. Can observed-history metrics choose the right predictor? | [MESA](https://arxiv.org/abs/2608.10108): distant adaptive-view-selection projection to numerical experts | Select GP/tree by validation hits or effect MAE; geometry control | Hit criterion113/384; MAE criterion76; geometry118. Criteria picked the same expert in all four states. **No successful adaptive routing.** | $0 | 2/5 diagnostic | Separate prediction accuracy from discovery utility; do not deploy this router. |
| 5. Do scores or candidate-specific examples improve LLM decisions? | MESA + [JitMem](https://arxiv.org/abs/2609.27334): view/read-time evidence projection; LLMNN controls; [LGBO](https://arxiv.org/abs/2605.17976) numerical/LLM inspiration | Matched menus: no report, scores, scores+observed examples; shuffled scores for Lite; Astra spot check | Lite31/35/30 hits per128. On three valid shuffled-control states: genuine28 vs shuffled27/96. Astra15/13/16 per64. GP29/tree33 per128. **No reliable added report benefit.** | $0.704363 | 1/5, repairs and tiny sample | Fix output schema and test genuine versus shuffled evidence before more rich-report variants. |

Strength scale: 1 = exploratory/confounded; 2 = controlled but small and overlapping; 3 = repeated independent confirmation; 4 = broad transfer; 5 = strong external replication. It rates evidence, not novelty or effect size.

## Plain-language verdict

- **Worked locally:** retaining the full features helped this median-effect model; using new measurements improved both numerical models over frozen history in R3.
- **Did not work reliably:** direct ranking, choosing an expert by average effect error, and adding richer reports to the LLM. A larger model did not rescue the report intervention on the two tested states.
- **Most useful analysis direction:** how the prediction target, representation and validation objective interact with discovery—not merely whether an agent can restate feedback accurately.
- **Next experiment:** repeat the feature/objective contrast with genuinely different initial-history composition and separate validation panels, then compare global-error and top-batch utility criteria. Avoid a dataset-name routing rule. If an incremental LLM effect remains, confirm it with schema-enforced32 selections and a length-matched shuffled-information control before campaigns.

**Spending:** 30 new paid calls including8 formatting corrections; $0.70436287. Whole-project recorded charges $4.69536828; conservative total including two older unresolved reservations **$4.70581388/$50**. Provider usage change matches new charges. Numerical local runtime is recorded, not priced as API spending.

See [full methods/results](RESULTS.md), [sequential decisions](DECISIONS.md), [metrics](METRICS.json), [audit](AUDIT.json), [accounting](FINAL_ACCOUNTING.json). All gains are within this restricted retrospective setting; no top-conference novelty or full-BDA superiority claim follows.
