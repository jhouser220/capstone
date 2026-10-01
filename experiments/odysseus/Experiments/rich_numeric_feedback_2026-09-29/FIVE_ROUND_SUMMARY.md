> Historical research record, copied 2026-10-01. Some linked raw artifacts remain in the local archive. See the package README for included records and runnable entry points.

# Five-round summary: richer numerical feedback

2026-09-29T22:14:11.394859+00:00

**Completed new experiments: six numerical models, multi-round numerical campaigns, new Lite/Flash inference, and own-history LLM campaigns.** Cost this batch **$0.135683**, 58 unique charged requests (including format failures); conservative project total **$3.320381/$50**. Shared/cached identical calls during recovery or common initial states are not independent experiments.

| Round / question | Specific source and move | Exact comparison | Actual result | API cost | Evidence strength / next decision |
|---|---|---|---|---:|---|
|1. Do richer numerical models help? |[GOLLuM](https://www.nature.com/articles/s42256-026-01283-z) limited metric-adaptation projection; [LLMNN](https://arxiv.org/abs/2509.21403) numerical/prior separation |Full available pool; kNN, XGBoost hit/effect, isotropic/ARD GP, geometry;6 states |Effect XGB47, GP46, ARD41, geometry33, classifier29, kNN13 per192 |$0|Early positive for continuous models;2/5, test campaigns |
|2. Does the gain last? |LLMNN sequential numerical acquisition; GOLLuM effect-model projection |Four new3-round cases;32 genes/round; own updated history |GP115 vs geometry116; effect XGB68, classifier75, ARD62, kNN42 per384 |$0|Early XGB gain failed to persist;2/5 descriptive, GP strongest learned comparator |
|3. Does the LLM use the numerical report usefully? |[LGBO](https://arxiv.org/abs/2605.17976) preference/numerical combination inspiration; LLMNN no-report control |GP vs Lite without/with report vs fixed/history-calibrated mixtures;4 matched states |GP21, no-report27, report29, mixtures25 per128 |$0.017440|Report+2 over no-report is weak;1/5 method benefit. Check later value and stronger model |
|4. Does planned exploration help? |[BATCHIE](https://www.nature.com/articles/s41467-024-55287-7) information-based acquisition projection |32 greedy vs24 greedy+8 marginal-uncertainty vs24+8 conditional information gain;4×3 rounds |Greedy118 vs74 for each exploration rule per384 |$0|Both lose all4 paired campaigns;2/5 local negative. Defer this25% exploration allocation |
|5. Does feedback/blending survive actual campaigns and capacity check? |LGBO/LLMNN combination and incremental-value analysis |Four new Lite3-round campaigns; plus Flash on the4 fixed R3 one-step states |Campaign: geometry116, GP109, no-report102, report103, online blend101 per384. Separate Flash: report31 vs no-report30, GP21; mixtures25 per128 |$0.118243|Small pilot; inspect dataset-specific reversals. No established new method |

## What this means

- **Repeated dataset split:** GP beats geometry on IFNG51vs44 in round2 and57vs43 in round5; it loses on IL264vs72 then52vs73. All four campaign seeds per dataset have the corresponding direction, but overlapping genes limit independence.
- Richer numerical tools changed results substantially. Their initial advantages did not automatically become campaign advantages.
- GP versus XGBoost and binary versus continuous targets are useful analysis axes. Isotropic GP was stronger than the more flexible ARD adaptation in the campaign development test. This is not evidence against learned representations generally.
- Reports beat GP locally, but their increment over the LLM without a report was only2 hits with Lite and1 with Flash. Much of the apparent advantage can therefore come from LLM prior selection or the shortlist, not demonstrated feedback use.
- Expected information gain is distinct from feature–outcome mutual information. Our Gaussian EIG implementation reduced model variance, but did not improve short-horizon discoveries. EIG and marginal variance had identical hit totals, not identical selections (6/12 batches identical).
- **Priority remains method analysis.** Established GP/tree/blending tools do not constitute a novel method. The useful question is which prediction objective, representation and uncertainty signal improves discovery at different campaign horizons. Novelty of that specific claim still needs comparison against related work; no publication-ready superiority is established.

**Best next experiment:** compare continuous-magnitude, binary-hit and target-tail/ranking objectives with the same tree features and training budget, including full Achilles coordinates versus16 PCA components. Keep the GP and geometry baselines; then test only a small score/report ablation where the report beats a no-report LLM on matched states. This addresses the actual objective/representation bottleneck instead of adding more critics or trusting uncalibrated confidence.

Counts are unique within each campaign, summed across runs; they are not globally distinct genes. Initial256 observed genes are shared and excluded from the reported new-hit totals. Different rounds have different cases and budgets; do not pool their totals. No BDA benchmark-superiority claim.
