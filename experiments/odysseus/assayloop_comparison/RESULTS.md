# Results history

Keep a short dated entry here for each result-affecting commit: what changed, the comparison and settings, results, and the linked evidence. Commit code and its entry together; Git history provides the commit ID. For changes without a new run, say “results unchanged” or “evaluation pending.” Preserve earlier entries.

## All benchmark and external results together

**Mean unique measured hits after 10 rounds × 100 selected candidates.** HSA means **Historical Subspace Adaptation**. Both systems use AWS Gemma 4 31B and Qwen3 32B, two runs per model per task, frozen AssayFormer, and three language-assisted rounds followed by seven numerical rounds. Original AssayLoop proposal shortfalls are filled to match the selection budget. These are our runs of the original policy with substituted actors, not numbers copied from the original paper.

| Dataset / cohort | Evidence | Original AssayLoop + budget fill | AssayLoop + HSA | Gain in hits | Relative gain |
|---|---|---:|---:|---:|---:|
| Original benchmark: initial test (3 publications) | Historical E160 | 86.58 | 98.33 | +11.75 | +13.6% |
| Original benchmark: additional test (6 publications) | Historical E160 | 147.50 | 173.83 | +26.33 | +17.9% |
| Original benchmark: validation (6 publications) | Historical E160 | 145.63 | 160.83 | +15.21 | +10.4% |
| Original benchmark: training (3 publications) | Historical E160 | 43.17 | 59.58 | +16.42 | +38.0% |
| Original benchmark: initial test (3 publications) | Fresh rerun (same exposed cohorts) | 93.08 | 96.75 | +3.67 | +3.9% |
| Original benchmark: additional test (6 publications) | Fresh rerun (same exposed cohorts) | 148.38 | 175.96 | +27.58 | +18.6% |
| Original benchmark: validation (6 publications) | Fresh rerun (same exposed cohorts) | 147.08 | 161.17 | +14.08 | +9.6% |
| Original benchmark: training (3 publications) | Fresh rerun (same exposed cohorts) | 44.67 | 57.50 | +12.83 | +28.7% |
| External: Schukken (1 assay) | Historical external | 47.75 | 50.00 | +2.25 | +4.7% |
| External: Abt–Wang (1 assay) | Historical external | 322.00 | 365.50 | +43.50 | +13.5% |
| External: Wilson (1 assay) | Historical external | 97.25 | 107.75 | +10.50 | +10.8% |
| External: Nguyen (1 assay) | Historical external | 13.50 | 15.00 | +1.50 | +11.1% |
| External: Shan (3 assays) | Historical external | 76.42 | 86.83 | +10.42 | +13.6% |
| External: Herranz-Ors (3 assays) | Historical external | 657.08 | 712.00 | +54.92 | +8.4% |
| External: Samelson (1 assay) | Historical external | 113.00 | 155.00 | +42.00 | +37.2% |

The historical E160 and external rows are the results used in the presentation. Fresh-rerun rows are a separate completed replication on the same exposed original cohorts, not additional datasets or independent confirmation. No new inference was performed for this table. **Training rows are explicitly training-cohort evaluations, not held-out generalization evidence and not LLM weight training.** The original panel contains 18 selected publications (9 test, 6 validation, 3 training), not the full original paper's splits. The external panel contains 7 publications / 11 selected assays, not all conditions available in their source papers.

Average the four actor/repetition campaigns within each task, then weight publications equally within original cohorts and assays equally within each external publication. Each original panel has 144 campaigns; the external panel has 88. Never sum or pool these table rows into one grand score. Historical E160 and external language-feedback conventions differ; see [dataset contracts](DATASETS.md). Gains use unrounded means. The 1,000 selections are not 1,000 measured outcomes; availability differs between policies.

All four original cohort averages improve, but historical publication means have **11 wins / 7 losses**, and fresh reruns **12 wins / 6 losses**. External publication means all improve, but individual paths can lose, and strong numerical controls often match or beat language + HSA. These whole-system comparisons do not isolate HSA's subspace component, certify absence of pretraining exposure, or establish benefit from language alone.

[Machine-readable combined table](recorded_results/combined_results.csv) · [historical campaign records](recorded_results/campaign_metrics.csv) · [fresh campaign records](recorded_results/fresh_replication_campaigns.csv) · [external campaign records](recorded_results/external_campaigns.json) · [reproduction report](REPRODUCTION_REPORT.md).

## 2026-10-27 — consolidated results navigation (no new experiments)

Original historical test, validation and training results, the completed fresh rerun, and all seven external sources are now in one table above and in the package README. Absolute means and relative gains were recomputed from 376 paired-arm campaign records. Results and implementation are unchanged. Detailed historical entries follow.

## 2026-10-22 — Portable package verified

**Code:** snapshot `c8eca0f81f6f` ([full source hashes](SOURCE_PROVENANCE.json)); release `odysseus-assayloop-reproduction-2026-10-22`. The public and private copies share this snapshot. Use `git log -- RESULTS.md` from this directory to find the commit containing an entry.

**Change:** extracted the existing method into a portable package and added complete replay and fresh-run commands. The method is unchanged. A repeated-worker transport-wrapper bug found during verification was fixed before the successful runs.

**Comparison:** our complete extension versus native AssayLoop with proposal shortfalls filled by AssayFormer. Both select 1,000 candidates. Gemma4 31B and Qwen3 32B; two repetitions; three language plus seven numerical decisions; 18 publications. Gains are additional unique measured hits, averaged equally across publications after averaging models/repetitions within publication.

| Cohort | Historical gain, reproduced exactly | Fresh rerun gain |
|---|---:|---:|
| Original exposed test (3 publications) | +11.75 | +3.67 |
| Additional exposed test (6) | +26.33 | +27.58 |
| Validation (6) | +15.21 | +14.08 |
| Training (3) | +16.42 | +12.83 |

**Checks and cost:** 342 historical campaigns replayed exactly in 4.61 minutes; 144 fresh campaigns and 432 calls in 31.27 minutes, costing approximately $0.2043. Zero failures in the successful full runs; 11 tests passed.

**Conclusion:** the whole-system gain persists in all four cohort averages; the original-test gain is smaller in the fresh run. Fresh publication averages: 12 wins, 6 losses. These are exposed development cohorts. Selected budgets match; measured availability varies. Attribution to individual components remains an ablation question.

**Evidence:** [full report](REPRODUCTION_REPORT.md), [historical campaign results](recorded_results/campaign_metrics.csv), [fresh campaign results](recorded_results/fresh_replication_campaigns.csv), [reproduction commands and asset requirements](REPRODUCING.md).

## 2026-10-27 — external transfer packaging (no new experiments)

Added seven external publication sources / eleven selected assays, paired HSA and budget-completed Original AssayLoop records, an explicit external measured-only feedback adapter, and a saved-proposal full-loop replay. Current code is based on c72a034 and preserves its original six-arm implementation and fresh-reproduction evidence. All 88 external campaigns matched at every selection/outcome; 15 tests passed. [Dataset contracts and limitations](DATASETS.md), [validation](VALIDATION.md), and [source/code hashes](SOURCE_PROVENANCE.json) identify this version. HSA128 remains the incumbent; recent RL refinements are not promoted.
