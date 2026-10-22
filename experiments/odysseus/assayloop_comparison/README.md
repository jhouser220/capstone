# Native AssayLoop versus Historical Response Subspace Adaptation

This package contains our numerical method, adapters for the authors' full sequential loop, and compact historical comparison data. The comparison uses **Qwen3 32B and Gemma4 31B**, two repetitions, three language decisions followed by seven numerical decisions, and 100 candidate slots per decision. The historical data are from completed E160 and its verified reused references. The portable extraction now reproduces all 342 complete historical campaigns exactly. A fresh 144-campaign replication also completed successfully; see the [verification report](REPRODUCTION_REPORT.md) and [reproduction instructions](REPRODUCING.md).

[Results history for this code](RESULTS.md) records each evaluated change and links to its evidence. Download the [frozen reproduction inputs](https://github.com/jhouser220/capstone/releases/download/odysseus-assayloop-reproduction-2026-10-22/assayloop-presentation-reproduction-2026-10-22.tar.gz) to rerun it.

## What differs between the systems?

| Arm | Language decisions1–3 | Numerical completion / decisions4–10 | Question |
|---|---|---|---|
| `native` | Original parsed proposals retained | Original AssayFormer continuation; proposal shortfalls preserved | Behavior of native implementation with our actors |
| `native_filled` | Same native proposals, all valid proposals retained | Fill missing slots using pre-batch AssayFormer; original continuation | Full-system comparison at equal selected-candidate budget |
| `incumbent` (ours) | Source/MAP ranking within proposal pool, retain up to25 | Historical Response Subspace Adaptation fills to100 and continues | Benefit of our complete extension |
| `ordinary_full` | Same proposal handling as ours | Source/MAP rank fusion; no historical32D, q or posterior-probability integration | Strong simpler system; differs in several scoring components |
| `assayformer_only` | Numerical from the first decision | AssayFormer for all ten decisions | Numerical-only reference |

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
python -m unittest discover -s experiments/odysseus/assayloop_comparison/tests -v
```

The recount and tests use saved numbers or synthetic fixtures. They require no credentials, GPU, network, biological experiments or LLM calls. [VALIDATION.md](VALIDATION.md) records the checks actually performed for this packaging change.

## Code map

- [numerics.py](numerics.py): preserved damped-Newton fit; current42D Laplace scoring with128 antithetic samples and one availability multiplier.
- [assets.py](assets.py): explicit portable bank schema and32D construction from approved, aligned historical donor arrays. Query-publication exclusion is enforced before fitting.
- [adapter.py](adapter.py): method recipient and acquisition adapters using the authors' `SequentialLoop`, `GreedyFromModel`, supplied AssayFormer and genuine actor acquisition.
- [audit_results.py](audit_results.py): independent cohort recount and paired budget checks.
- [SOURCE_PROVENANCE.json](SOURCE_PROVENANCE.json): original component paths relative to the research archive, hashes and extraction changes.

The current method appends32 historical directions to the source's10 features, uses a history-conditioned source center and a fixed historical mean, fits on measured selected labels, and ranks by availability times integrated hit probability. The source encoder retains the recorded native unavailable-outcome convention; the likelihood filters unavailable labels. New R160 proposals are future experiments and are absent from the implementation.

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

The actor is the original `LLMSingleAcquisition` with an explicitly configured authorized provider. The recorded comparison used temperature1.0, an8,192-token output limit and the exact served Qwen/Gemma identities in the campaign CSV; preserve these settings and the source parser/prompt for a matched rerun. Caller-owned transport must enforce the existing budget/quota policy and preserve requests, responses and failures. The package never reads private account files. Credentials stay in the caller's environment. Full campaigns are opt-in and should use the same publication/actor/repetition matrix across arms. No changed-history response reuse; feedback is revealed only after selection. Save returned histories, audit, exact assets/configuration and cost records under ignored `outputs/`.

Use [REPRODUCING.md](REPRODUCING.md) for executable historical replay and fresh AWS replication commands. Exact replay checks every generated request, parsed proposal, selection, outcome and final metric against the frozen evidence. Fresh model responses are stochastic, so their results are reported separately.

## Relationship to the team repository

The root team project retains its BioDiscoveryAgent implementation. This AssayLoop extension is scoped to Odysseus's experiment directory and does not replace Jeff or Chris's workflows. The existing older numerical/RL package is retained separately. No private sessions, account ledgers, credentials or raw provider responses are part of this addition.
