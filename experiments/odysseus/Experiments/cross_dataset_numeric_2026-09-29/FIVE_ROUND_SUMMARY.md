> Historical research record, copied 2026-10-01. Some linked raw artifacts remain in the local archive. See the package README for included records and runnable entry points.

# Cross-dataset experiments: five-round summary

Completed 2026-09-29T23:04:01.236680+00:00; project local date2026-09-30. Folder retains its2026-09-29 start date.

**Yes, broader testing changed the conclusion.** There is no stable rule that GP helps IFNG and hurts IL2. Initial history, the remaining action pool, training-data support and score meaning materially change the comparisons. A classifier has a substantial local advantage on Scharenberg with256 observations; that advantage disappears with16. No universal new method is demonstrated.

| Round | Question and source move | Exact comparison | Result | Cost | Evidence / next decision |
|---|---|---|---|---|---|
|1|Cross-screen transfer; [BDA](https://arxiv.org/abs/2405.17631) task breadth + [GOLLuM](https://www.nature.com/articles/s42256-026-01283-z) representation projection|Seven tables,3 seeds,6 methods,256 mixed-history observations +3×32 choices|Scharenberg classifier52 vs geometry19/288; GP Carnevale45vs45; Sanchez day30 GP16vs36; CAR-T GP1vs3. IFNG GP54vs64 includes a zero-hit campaign.|$0|Preliminary,2/5. Test starting-history dependence.|
|2|Prior versus learned evidence; [LLMNN](https://arxiv.org/abs/2509.21403) experimental-design transfer|Same seven tables/3 seeds; four core methods with random-only256 histories|Carnevale GP42vs26; Sanchez day30 GP31vs18; IL2 full median96vs76; Scharenberg classifier49vs29. Different histories also leave different eligible genes.|$0|2/5. Positive effects exist, but some signs depend on initialization.|
|3|Prediction versus acquisition; GOLLuM/LLMNN diagnostic projection|Same GP fit, posterior mean versus threshold-exceedance score, mixed histories; reuse phase1 controls|IFNG58vs54, but failed case1vs0; Carnevale34vs45. No general rescue.|$0|2/5. Do not treat acquisition change alone as a new method.|
|4|Observation fraction; BDA screen-size and GOLLuM sample-efficiency inspiration|Scharenberg16-random-history versus phase2's256; same3×32 new choices within each condition|Classifier29vsgeometry31 with16, compared with49vs29 with256. Different total experimental budgets are intentional.|$0|2/5. More observed support plausibly contributes; biology alone is insufficient.|
|5|History information versus action availability; LLMNN control projection|Two256-history views, identical query pool excluding their union; one-step3×32 per arm on3 selected screens|IL2 median mixed36vsrandom17; Sanchez day30 GP mixed3vsrandom12. Equal available actions do not remove history dependence.|$0|2/5 adaptive replay. Prioritize information composition and target alignment; confirm with fresh states.|

Evidence2/5 means controlled but exploratory: three starting seeds on overlapping biological data, adaptive follow-ups, and no independent study replication. It is not a rating of novelty or effect size.

**Executed:**240 new three-batch numerical campaigns (126+84+21+9), plus54 one-step comparison records in phase5;9 of those54 records are the identical geometry selection repeated under the other history label. Controls reused between phases are not new independent runs. These are retrospective simulations against recorded screen outcomes, not new laboratory experiments or LLM campaigns.

**Spending:** no new API calls, $0 additional API cost. Conservative recorded project total$4.70581388/$50; two older unresolved reservations remain included. Local runtime is recorded separately and is not billed as API usage.

**Best next tests:** (a) learning curves varying observed hit/non-hit support and representative coverage while holding available candidates fixed; (b) task-aligned two-tail modeling for Sanchez day21, where absolute score is a poor proxy for hit membership; (c) diagnose GP fit instability before adding an LLM router. These are proposed follow-ups, not demonstrated methods. The elementary controls are not novel by themselves.

See [full results](RESULTS.md), [all metrics](SUMMARY.json), [dataset properties](DATASET_PROPERTIES.json), [decisions](DECISIONS.md), [audit](AUDIT.json).
