# Validation of the public package — 22 October 2026

## Complete historical execution replay

All **342 campaigns passed**, covering **3,420 decisions and 864 genuine recorded language requests**. The portable package regenerated requests from the current state, used the original parser, recomputed numerical predictions and selections, revealed outcomes after selection, and matched each historical trajectory exactly. Final unique hits, selected/measured counts and enrichment factors also matched the compact presentation CSV. This took **276.62 seconds** with two workers and made zero provider calls.

The tested matrix contains 72 campaigns each for the complete extension, original native, equal-budget native and same-proposal ordinary system; 36 historical-scorer-only campaigns; and 18 AssayFormer-only campaigns. Each campaign has ten decisions. Language-led systems have three language decisions followed by seven numerical decisions.

The first full attempt exposed repeated installation of the verification transport in reused workers. The wrapper was made idempotent, a regression test was added, and the complete matrix was rerun successfully. This was a harness defect; the experiment's method and historical data were unchanged.

## Additional checks

- **11 tests passed** with pinned upstream dependencies. They cover numerical fitting, the original frozen recipient fixture to absolute tolerance 1e-12, donor exclusion, duplicated-record rejection, cohort aggregation, full ten-decision integration for all six arms, proposal preservation, explicit live-call opt-in and repeat-safe transport setup.
- Three altered-evidence checks independently changed a saved history, HTTP request and parsed proposal. All were rejected by the replay before a result could be accepted.
- Six critical upstream files were fetched from the declared GitHub commit and matched byte for byte: language acquisition, greedy acquisition, AssayFormer wrapper, task, language client and ranker network.
- Every input and upstream Python file in the companion archive was hash-verified: 547 files. The bundle was scanned for credential patterns and private account/path markers.
- The original public checkout's staged patch remained byte-identical. Work was performed in the separate reproduction worktree.

## Fresh replication

The fresh AWS comparison completed 144/144 campaigns and 432 calls in 31.27 minutes, at a usage-based cost of $0.204319. The extension gained hits in all four cohort means and won 12 of 18 publication means against equal-budget native. Full results are recorded in [REPRODUCTION_REPORT.md](REPRODUCTION_REPORT.md). Historical replay and fresh replication answer different questions: replay establishes exact implementation reproduction; fresh calls measure another stochastic realization on the same exposed cohorts.

See [REPRODUCING.md](REPRODUCING.md) for commands, asset requirements and runtime versions. The companion bundle is prepared separately from Git. The [GitHub Release](https://github.com/jhouser220/capstone/releases/tag/odysseus-assayloop-reproduction-2026-10-22) ties the bundle to the code version. Verification was completed before publication.

## External extension — 27 October 2026

Preparation started with independent raw recounts of 88 external and 144 original benchmark HSA/budget-matched campaigns, with zero recount errors. The portable adapter then reproduced all 88 external campaigns (880 decisions across eleven assays / seven publications): every selected identity and outcome matched, as did HSA posterior modes within absolute tolerance 2e-5. The replay consumes saved parsed proposals only after exact preceding-history verification; it does not regenerate HTTP requests or rerun the source parser. Earlier original-benchmark transport-level validation above remains a separate stronger check of that layer.

No provider calls, new scientific campaigns or paid compute were launched during this preparation. Replay does not establish independence from pretraining exposure, absence of adaptive dataset selection, or future stochastic language-model equivalence. The external results are retrospective development evidence, not untouched confirmation.

The package is reconciled onto team remote c72a034, preserving its six-arm adapter, original reproduction harness and released assets. On this reconciled checkout, all 15 unit/integration tests passed and the full 88-campaign replay passed again. All 44 external source assets matched the declared SHA-256 inventory. See [machine-readable replay verification](recorded_results/external_replay_verification.json). The original benchmark compact recount and external compact recount both passed.

The exported measured-only feedback adapter also produced identical prompt tuples to the original E209 adapter on all 264 saved external language histories. This is prompt-construction equivalence, not transport/parser replay or a fresh inference run.
