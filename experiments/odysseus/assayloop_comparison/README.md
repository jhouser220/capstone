# Original AssayLoop versus Historical Subspace Adaptation

This package contains Historical Subspace Adaptation (HSA), adapters for the authors' full sequential loop, and compact historical comparison data. The comparison uses **Qwen3 32B and Gemma4 31B**, two repetitions, three language decisions followed by seven numerical decisions, and 100 candidate slots per decision. Recorded evidence covers the original E160 panel and seven external publications (eleven selected assays). These are existing development results. The original benchmark package also has a previously completed fresh 144-campaign replication; see [REPRODUCTION_REPORT.md](REPRODUCTION_REPORT.md) and [REPRODUCING.md](REPRODUCING.md). This update adds offline external replay, not another fresh language-model trial. See [datasets and contracts](DATASETS.md) and [validation](VALIDATION.md).

## Experimental setup behind the presentation results

| Setting | Recorded protocol |
|---|---|
| Language models | **Gemma 4 31B** (`google.gemma-4-31b`) and **Qwen3 32B** (`qwen.qwen3-32b-v1:0`), served through AWS |
| Numerical foundation | The same **frozen AssayFormer checkpoint** in both systems; no LLM fine-tuning. HSA fits task-specific numerical weights using selected measured outcomes. |
| Repetitions | **2 runs per model, per assay, per method**: 4 campaigns per assay per method |
| Rounds | **10 sequential rounds per campaign** |
| Selection budget | **100 genes per round × 10 rounds = 1,000 gene selections per campaign**, matched between the two presented systems |
| Language phase | Rounds **1–3**: language proposals plus numerical completion to 100 selections |
| Numerical phase | Rounds **4–10**: numerical selection only |
| Proposal handling | Original AssayLoop retains valid language proposals and fills short batches using pre-batch AssayFormer. HSA retains up to 25 ranked language proposals and numerically fills the remainder. |
| Sampling configuration | Temperature 1.0; 8,192-token output cap; preserve recorded model identities, prompts and parser for reruns |
| Main metric | **Unique measured hits after the full 1,000-selection campaign** |
| Reported means | Average the 4 model/run campaigns within each assay; then average assays equally within each publication (3 assays each for Shan and Herranz-Ors). |
| External comparison size | **7 publications, 11 selected assays, 88 campaigns** = 11 assays × 2 methods × 2 models × 2 runs; 880 decision rounds |

A “run” or “campaign” means the complete ten-round trajectory; it is not one language-model call. The 1,000-selection budget is neither 1,000 measured outcomes nor 1,000 hits: some selected genes lack a recorded measurement. It is also not a dollar budget or an API-call count. Historical inference costs and provider attempts require the original receipts; this package's offline recount/replay makes no new LLM calls. Each arm receives its own selected-history feedback after divergence.

These external results use measured-only language feedback. The original E160 benchmark panel uses the same two actors and two repetitions, but an earlier declared feedback convention; do not silently pool the protocols or present the local 18-publication panel as the complete original paper test split. See [DATASETS.md](DATASETS.md).

## Current result and scope

**Original AssayLoop + matched-budget completion** versus **AssayLoop + HSA** is the main comparison. “Original” identifies the source model/policy; our language actors, measured-only feedback and budget completion are declared substitutions/extensions, not exact reproduction of the paper's original actors. HSA remains the incumbent after E252. This is a lean method package, not the full research archive.

| External source | Original + budget fill | AssayLoop + HSA | Relative gain |
|---|---:|---:|---:|
| Schukken | 47.75 | 50.00 | 4.7% |
| Abt–Wang | 322.00 | 365.50 | 13.5% |
| Wilson | 97.25 | 107.75 | 10.8% |
| Nguyen | 13.50 | 15.00 | 11.1% |
| Shan (3 assays) | 76.42 | 86.83 | 13.6% |
| Herranz-Ors (3 cultures) | 657.08 | 712.00 | 8.4% |
| Samelson | 113.00 | 155.00 | 37.2% |

Means are unique measured hits per 1,000 requested slots: four actor/draw paths per assay, then equal assay means within each publication. Seven publications are not eleven independent sources. Strong cold HSA/MAP/frequency controls often match or exceed language+HSA; Nguyen loses on two of four paired paths. These gains do not establish a benefit from language, isolated subspace causality, no training exposure, or independent confirmation. [Full comparison boundaries](DATASETS.md).

## What differs between the systems?

| Arm | Language decisions1–3 | Numerical completion / decisions4–10 | Question |
|---|---|---|---|
| `native` | Original parsed proposals retained | Original AssayFormer continuation; proposal shortfalls preserved | Behavior of native implementation with our actors |
| `native_filled` | Same native proposals, all valid proposals retained | Fill missing slots using pre-batch AssayFormer; original continuation | Full-system comparison at equal selected-candidate budget |
| `incumbent` (ours) | Source/MAP ranking within proposal pool, retain up to25 | Historical Subspace Adaptation fills to100 and continues | Benefit of our complete extension |
| `ordinary_full` | Same proposal handling as ours | Source/MAP rank fusion; no historical32D, q or posterior-probability integration | Strong simpler system; differs in several scoring components |
| `assayformer_only` | Numerical from the first decision | AssayFormer for all ten decisions | Numerical-only reference |
| `historical_only` | Numerical from the first decision | HSA for all ten decisions | Historical-scorer-only reference |

Our changes include **both proposal handling and numerical scoring**. A whole-system gain is supported by comparing ours with native_filled. The historical-subspace component alone needs the separately planned matched ablation. The paper's original actors are different; this is a comparison using our substituted actors and declared cohorts.

## Historical results: additional final unique measured hits

| Cohort | Ours − original native | Ours − native filled | Ours − same-proposal ordinary |
|---|---:|---:|---:|
| Original test3 | +21.50 | +11.75 | −2.75 |
| Additional exposed test6 | +47.33 | +26.33 | +8.92 |
| Validation6 | +31.21 | +15.21 | +3.79 |
| Training3 | +23.00 | +16.42 | +15.33 |

At equal1,000 selected candidates, ours wins11 and loses7 publication means against native_filled. Cohort means improve, while enrichment factor and measured-hit precision remain heterogeneous. **Selected candidates and measured outcomes are different counts.** Availability differs between policies. Keep requested, selected and measured columns visible. All test cohorts here have been exposed during the research programme; independent confirmation is outstanding.

[Campaign records](recorded_results/campaign_metrics.csv) include all methods, actors, repeats and scientific artifact hashes. [Cohort table](recorded_results/cohort_means.csv) includes AssayFormer-alone, EF, counts and measured-hit precision. Rows for numerical-only controls have different repetition structure; average within publication before comparing cohort means. The raw trajectory archives and model weights are external assets, rather than public Git payloads.

## Inspect and verify without inference

From the repository root:

```sh
python -m pip install -r experiments/odysseus/assayloop_comparison/requirements.txt
python experiments/odysseus/assayloop_comparison/audit_results.py
python experiments/odysseus/assayloop_comparison/audit_external.py
python -m unittest discover -s experiments/odysseus/assayloop_comparison/tests -v
```

The two recount commands use compact saved numbers; unit tests use synthetic fixtures and pinned-loop integration. The real-data replay below separately verifies full selections and outcomes. They require no credentials, GPU, network, biological experiments or LLM calls. [VALIDATION.md](VALIDATION.md) records the checks actually performed for this packaging change.

## Code map

- [numerics.py](numerics.py): preserved damped-Newton fit; current42D Laplace scoring with128 antithetic samples and one availability multiplier.
- [assets.py](assets.py): explicit portable bank schema and32D construction from approved, aligned historical donor arrays. Query-publication exclusion is enforced before fitting.
- [adapter.py](adapter.py): method recipient and acquisition adapters using the authors' `SequentialLoop`, `GreedyFromModel`, supplied AssayFormer and genuine actor acquisition.
- [audit_results.py](audit_results.py): independent cohort recount and paired budget checks.
- [SOURCE_PROVENANCE.json](SOURCE_PROVENANCE.json): original component paths relative to the research archive, hashes and extraction changes.

The current method appends32 historical directions to the source's10 features, uses a history-conditioned source center and a fixed historical mean, fits on measured selected labels, and ranks by availability times integrated hit probability. The source encoder retains the recorded native unavailable-outcome convention; the likelihood filters unavailable labels. Later research variants are documented separately; this package retains the frozen incumbent rather than promoting them.

## Running the full system

Install the authors' AssayLoop at commit `e7251afacda88c803864896e0d1444196f482471` and its AssayBench dependency, following their installation instructions. [Upstream repository](https://github.com/Genentech/AssayLoop). Preserve the checkpoint `model_last.pt`, vocabulary, cached description embeddings, approved candidate universe and per-screen metadata. Source and donor assets must follow the frozen manifests used for the intended comparison. A rebuilt bank can differ numerically across SVD/library versions; use frozen factors for exact historical replay.

The array-only bank constructor needs historical mask/hit matrices, donor publication IDs, aligned source gene weights and cached donor/query description embeddings. It has no input for query outcomes. `HistoricalBank.load()` accepts a pickle-free NPZ with `genes`, `residual` (unscaled32D), `mean` (already0.1-scaled), `availability`, `scalar_c`, `query_publication`, and `donor_publications`. The scorer applies0.25 to residual directions exactly once. Verify donor provenance and embedding normalization before generating or sharing these assets.

Use a **fresh task, source instance and actor per campaign**. The following integration pattern assumes those objects have been constructed using the pinned upstream classes and approved data:

```python
from experiments.odysseus.assayloop_comparison.adapter import SourceState, run_campaign
from experiments.odysseus.assayloop_comparison.assets import HistoricalBank

state = SourceState(source_model, cached_description_vector)
bank = HistoricalBank.load(bank_path)
result, audit = run_campaign(
    arm="incumbent", task=task, source_model=source_model,
    screen_id=screen_id, draw=0, source_state=state, bank=bank,
    actor_acquisition=original_llm_acquisition, allow_live_calls=True,
)
```

For the newer external comparison, construct `feedback.MeasuredOnlyAcquisition`, a subclass of the original `LLMSingleAcquisition`, with an explicitly configured authorized provider. It removes unavailable labels from language feedback, lists those candidates neutrally, and removes the query total-hit count from context. Retain the original numerical history. Historical E160 uses its own declared earlier feedback contract; do not silently pool these versions. The recorded comparison used temperature1.0, an8,192-token output limit and the exact served Qwen/Gemma identities in the campaign CSV; preserve these settings and the source parser/prompt for a matched rerun. Caller-owned transport must enforce the existing budget/quota policy and preserve requests, responses and failures. The package never reads private account files. Credentials stay in the caller's environment. Full campaigns are opt-in and should use the same publication/actor/repetition matrix across arms. No changed-history response reuse; feedback is revealed only after selection. Save returned histories, audit, exact assets/configuration and cost records under ignored `outputs/`.

## Replay real saved campaigns

The packaged adapter can run the authors' full ten-decision loop against exact saved language proposals. It refuses changed preceding selections or feedback, verifies every subsequent selected identity and outcome, and refits HSA numerically. This is stronger than a summary recount but **is not fresh LLM sampling**. Its source arrays/checkpoints and recorded proposal inputs must be the frozen originals.

```sh
python experiments/odysseus/assayloop_comparison/replay_saved.py \
  --research-root /path/to/restored/Capstone \
  --checkpoint /path/to/trusted/assayformer/checkpoint \
  --all-paths --output /tmp/hsa-replay.json
```

Install the pinned upstream dependencies and `torch`/`threadpoolctl` first. The private research snapshot supplies the relative source paths in `recorded_results/external_campaigns.json`; it is not downloaded automatically and no private account configuration is read. Only load trusted checkpoint/factor files. Original datasets, source weights and raw provider receipts are not redistributed in this small public package; consult their original licenses and [external asset inventory](EXTERNAL_REPRODUCTION_ASSETS.json). A fresh external genuine-LLM equivalence campaign remains a separate, unperformed validation level. The original benchmark replay/reproduction bundle remains available via [REPRODUCTION_ASSETS.json](REPRODUCTION_ASSETS.json).

## Recent refinements, not promoted replacements

E249 immediate-profile learning was a small validation development lead. E250 fixed-direction attribution was inconclusive; E251 showed actor transfer of unchanged HSA; E252's longer training did not improve the lead. These remain in the private research archive. They do not replace this package's HSA128 method or justify publishing only favorable variants.

## Relationship to the team repository

The root team project retains its BioDiscoveryAgent implementation. This AssayLoop extension is scoped to Odysseus's experiment directory and does not replace Jeff or Chris's workflows. The existing older numerical/RL package is retained separately. No private sessions, account ledgers, credentials or raw provider responses are part of this addition.
