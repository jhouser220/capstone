# Results history

Keep a short dated entry here for each result-affecting commit: what changed, the comparison and settings, results, and the linked evidence. Commit code and its entry together; Git history provides the commit ID. For changes without a new run, say “results unchanged” or “evaluation pending.” Preserve earlier entries.

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
