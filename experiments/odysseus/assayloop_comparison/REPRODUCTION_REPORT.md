# 2026-10-22T15:35:09+01:00 — Public-package reproduction verification

The portable package **reproduces the presentation's complete historical results exactly**. A separate fresh AWS replication completed **144/144 campaigns** with **432 new provider calls**, across the same 18 screens, two actors, two repetitions and ten decisions per campaign.

## Exact reproduction of the presentation

All 342 historical campaigns, 3,420 decisions and 864 recorded language requests passed. The package reconstructed each prompt from its actual current history, passed the saved response through the original parser, recomputed the numerical continuation and matched every selected candidate and revealed outcome. All final hits, selected and measured counts and official EF values matched the public campaign CSV. The presentation's equal-selected-budget margins are therefore reproduced: **+11.75, +26.33, +15.21 and +16.42 hits**, respectively for original exposed test, additional exposed test, validation and training.

Replay runtime: **4.61 minutes**; fresh provider calls during replay: **zero**. Using genuine recorded responses makes this an exact implementation check. It also avoids treating random changes in model completions as a code discrepancy.

## Fresh primary comparison

| Cohort | Our extension: hits | Native with budget completion: hits | Fresh difference | Historical difference |
|---|---:|---:|---:|---:|
| Original exposed test: 3 publications | 96.75 | 93.08 | +3.67 | +11.75 |
| Additional exposed test: 6 publications | 175.96 | 148.38 | +27.58 | +26.33 |
| Validation: 6 publications | 161.17 | 147.08 | +14.08 | +15.21 |
| Training: 3 publications | 57.50 | 44.67 | +12.83 | +16.42 |

Values are final unique measured hits averaged across actors and repetitions within each publication, then equally across publications. Both systems selected 1,000 distinct candidates in every campaign. The number with available measured outcomes differs by policy; campaign records include measured counts, measured-hit precision and EF. Across the 18 publication means, the fresh extension wins **12**, loses **6** and ties **0** against native with budget completion.

Fresh runtime for the full matrix, including both actors and repetitions: **31.27 minutes**. Usage-based inference cost: **$0.204319**; uncertain billing reserve: **$0.000000**; failures: **0**. The incremental cap was $1 and concurrency was two. Cost is calculated from returned token usage and the verified AWS price catalogue, rather than a settled invoice. Published rates checked on 22 October 2026 were $0.14/$0.40 per million input/output tokens for Gemma4 31B and $0.15/$0.60 for Qwen3 32B. [AWS pricing catalogue](https://pricing.us-east-1.amazonaws.com/offers/v1.0/aws/AmazonBedrock/current/us-east-1/index.json).

Each language call was newly requested. Once trajectories diverged, each branch used its own feedback. The original three-language/seven-numerical schedule, 100 candidates per decision, source checkpoint, historical banks, candidates, prompts and parser remained fixed. This rerun used the complete extension and the equal-budget native comparator. The four other historical controls were fully replayed from their genuine recorded responses.

## What this establishes

The prepared public implementation can regenerate the presentation results from the frozen companion assets. Fresh results provide an additional stochastic replication on these same exposed cohorts. Method selection has already used these cohorts, so independent generalization confirmation remains a separate experiment. A whole-system comparison includes both our proposal handling and our numerical scoring; isolating the subspace's contribution requires matched ablations.

The native comparator uses the authors' implementation with Gemma/Qwen actor substitutions and explicit completion of proposal shortfalls. Equal selected budget is the stated comparison; availability and precision remain visible. These numbers should be compared with the declared controls and cohorts, rather than the paper's different actor/data headline.

## Checks, assets and handoff

- 11 unit/integration tests passed. Three altered-receipt checks rejected changed history, request and parser evidence.
- Six critical upstream files match the authors' pinned GitHub source byte for byte. All 547 input/upstream source files in the companion archive passed hash verification.
- An initial replay attempt found a worker-wrapper installation bug in the new harness. The wrapper was corrected and the entire replay then passed. The scientific method and frozen data were unchanged.
- Code and detailed instructions: [REPRODUCING.md](REPRODUCING.md); historical results: [campaign metrics](recorded_results/campaign_metrics.csv); fresh results: [campaign metrics](recorded_results/fresh_replication_campaigns.csv) and [publication-balanced summary](recorded_results/fresh_replication_summary.json).
- Companion archive: `assayloop-presentation-reproduction-2026-10-22.tar.gz`, 186,856,583 bytes, SHA-256 `19d409b2ab048ff1b8f402183d33bd0073efb0995a5ffa8049b6e08dc0e93903`. The archive is distributed with the [versioned GitHub Release](https://github.com/jhouser220/capstone/releases/tag/odysseus-assayloop-reproduction-2026-10-22). The Git package alone needs this companion data to execute the replay.
- Worktree branch: `reproduce/presentation-2026-10-22`. The original public checkout's staged patch was preserved byte for byte. The active research checkout was read for frozen input export and was left unchanged by this work. Verification was completed before publication.

Raw requests, responses, per-decision trajectories, costs, failed-attempt evidence and runtime manifests are retained under the isolated worktree's `_reproduction/` directory. The public package contains compact summaries and executable reproduction instructions; the companion archive supplies the larger frozen inputs.
