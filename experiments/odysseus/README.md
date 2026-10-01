# Numerical learning and LLM feedback experiments

Odysseus's focused research package, 1 October 2026. These are retrospective, exploratory experiments using released BioDiscoveryAgent screens. They are **not a reproduction of the complete BDA agent** or a demonstrated new method. The original repository's agent and the other collaborators' code are unchanged.

## What is included

| Folder under `Experiments/` | Main comparisons | Main entry points |
|---|---|---|
| `rich_numeric_feedback_2026-09-29` | GP, trees, kNN, numerical reports, LLM revisions/blending, information-based acquisition | `study.py`, `llm_study.py`, `llm_campaign.py`, `acquisition.py`, `diagnostics.py` |
| `objective_representation_2026-09-29` | Prediction targets, full/compressed features, updated/frozen history, validation-based routing, report views | `study2.py`, `router.py`, `views.py` |
| `cross_dataset_numeric_2026-09-29` | Seven screen tables, random/mixed starting histories, smaller initial history, fixed candidate-pool replay | `run.py`, `matched_pool.py` |

Keep these three folders together: the later studies import the earlier learners. `recorded_results/` contains bundled per-case results, protocols and summaries from the completed experiments. These are **historical outputs**, not freshly rerun results. Original descriptions and paper links are in each study's `PROTOCOL.md`, `RESULTS.md` and `FIVE_ROUND_SUMMARY.md`. Links to full raw logs, predictions and accounting artifacts in those archived notes may refer to the original local archive; those large artifacts and account records are not included here.

## Setup

From the repository root, use a separate environment (the original BDA environment has different dependencies):

```sh
python3.11 -m venv experiments/odysseus/.venv
source experiments/odysseus/.venv/bin/activate
python -m pip install -r experiments/odysseus/requirements.txt
mkdir -p experiments/odysseus/data
curl -L --fail https://figshare.com/ndownloader/files/49843176 -o experiments/odysseus/data/achilles.csv
python experiments/odysseus/prepare_features.py --csv experiments/odysseus/data/achilles.csv
python experiments/odysseus/check.py
```

The CSV source is the same public Figshare file referenced by this repository's `achilles.py`. The builder preserves gene-column order, drops incomplete cell-line rows, transposes to gene profiles and applies float32 L2 normalization. It records source/output checksums. This reconstructs the feature transformation; it does **not** establish absence of assay overlap or biological leakage. Historical feature-file SHA256: `ff7787507c4b79e98087415c71a2221941c13e6f6d979164a1c0f962e80acc16`. ZIP serialization can change the file hash; compare arrays as well when auditing a rebuild.

By default, inputs come from the repository's existing `datasets/` and the generated `experiments/odysseus/data/achilles_normalized.npz`. Override with `BDA_DATA_DIR` and `BDA_FEATURE_FILE`; optionally require an exact checksum with `BDA_FEATURE_SHA256`. No private source project or absolute user path is needed. PCA features and outputs are generated locally and ignored by Git.

## Run the numerical comparisons — no paid calls

```sh
cd experiments/odysseus/Experiments/cross_dataset_numeric_2026-09-29
python run.py phase1.json phase2.json phase3.json phase4.json
python matched_pool.py
```

The four configs run 240 three-batch campaigns. `matched_pool.py` runs 54 one-step comparison records, including nine repeated identical baseline selections. For a small four-method IFNG check, run `python run.py smoke_check.json` (one seed, three batches). Always use a new output `name` when changing a config: the historical runner resumes existing case files and does not validate all changed settings. Do not count resumed outputs as new experiments. Runs may take substantial local CPU time; they do not call LLMs.

Earlier numerical entry points accept their corresponding JSON configurations:

```sh
# From rich_numeric_feedback_2026-09-29:
python study.py r1.json
python study.py r2.json
python acquisition.py r4.json
python diagnostics.py
# From objective_representation_2026-09-29:
python study2.py r1.json
python study2.py r2.json
python study2.py r3.json
python router.py r4.json
```

## Optional LLM experiments

**New paid requests are disabled by default.** This packaging work makes no paid requests. LLM entry points are `llm_study.py r3.json`, `llm_study.py r5_flash.json`, `llm_campaign.py r5.json` (rich study) and `views.py` (objective study). The model names, provider restrictions and price ceilings preserve the historical protocols; verify current availability/prices before rerunning them.

To opt in, provide `OPENROUTER_API_KEY`, set `BDA_ENABLE_PAID_CALLS=1`, and set `BDA_BUDGET_ROOT` to the **shared, reconciled experiment-ledger directory** containing one study directory per ledger (`*/budget.sqlite`). Do not use an empty directory as a substitute for previous spending. Reservations use that shared location; unresolved calls remain reserved, provider prices are capped, and the working allocation stops below $40, preserving room within the project's total $50 cap. The ledger cannot account for other accounts, external charges or concurrent reservations in another study database; serialize paid runs and reconcile before each batch. Individual historical caps are $1 and $2.50. Never commit keys, responses containing account metadata, or ledger databases.

## Interpretation and paper roots

- [BDA](https://arxiv.org/abs/2405.17631): released screen outcomes and sequential gene selection.
- [LLMNN](https://arxiv.org/abs/2509.21403): inspiration for separating prior knowledge, numerical learning and incremental feedback value.
- [GOLLuM](https://www.nature.com/articles/s42256-026-01283-z): projection to feature compression, GP distance adaptation and target comparisons; we do not reproduce its learned encoder.
- [BATCHIE](https://www.nature.com/articles/s41467-024-55287-7): projection of information-based experimental design to Gaussian batch acquisition.
- [LGBO](https://arxiv.org/abs/2605.17976): inspiration for combining LLM preferences and numerical scores; our probability mixtures differ from its mechanism.
- [MESA](https://arxiv.org/abs/2608.10108) and [Just-in-Time Memory](https://arxiv.org/abs/2609.27334): inspiration for report-view and candidate-specific evidence comparisons.

Findings: numerical learning sometimes improves discoveries, but direction depends on dataset and starting history; richer LLM reports have not shown reliable incremental gains. Scores, binary hits and ranking objectives are not interchangeable. Sanchez day21 includes both score tails, and absolute magnitude can target the wrong region. Seven tables include related Sanchez screens, not seven independent studies. The fixed geometry baseline is not full BDA.

Hit totals exclude initial observations, count a gene once per campaign, and sum across overlapping runs. No outcome of an unselected gene enters the fitting interface. Whole-pool diagnostics are evaluator-only. LLM menus are restricted numerical shortlists; even the no-report arm therefore inherits some history dependence. Small/adaptively selected cases, feature provenance and incomplete independent replication limit claims.

`SOURCE_MANIFEST.json` maps the copied execution files to their original hashes. Packaging changes cover paths, missing-feature errors, optional checksum validation, paid-call opt-in and shared ledger configuration. Learners, prompts, experiment seeds and scoring rules are preserved. See `VALIDATION.md` for checks run before pushing.
