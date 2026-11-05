# DepMap / HSA comparison — results and how to reproduce them

`main_results.csv` is the one table behind every DepMap and HSA comparison in this work:
50 rows, one per **(screen set, candidate pool, LLM protocol, arm)**. It is regenerated,
never hand-edited:

```bash
uv run python scripts/diagnostics/build_main_table.py        # writes results/main_results.csv
```

## The three axes — a row is only comparable to another sharing all three

| axis | values | what differs |
|---|---|---|
| **set** | `paper_test`, `paper_validation` (curated 20-screen panels), `validation fold`, `TEST fold` (whole AssayBench folds, genome-scale subset), `team11 external` | which screens |
| **pool** | `f2`, `library` | what the policy could pick from |
| **protocol** | `cold`, `retain25`, `warm` | what rounds 1–3 did |

**pool.** `f2` is the shared 21,147-gene universe. The paper uses it *"for fair LLM
comparison"* — classical models normally pick from the screen's own library, LLMs name
genes from the whole genome, so the models are handicapped up to the LLM's pool. Picks that
land outside a screen's library have no label and are **refunded** from the EF denominator.
That refund is **per-arm**: `AF + DepMap` sends ~85 picks per 1,000 off-library, `dm.50`
~38, so it discounts `AF + DepMap` twice as hard and **reorders the arms**. `library` gives
each screen its own genes, nothing can land outside, and every EF convention agrees.
With no LLM in a comparison, `library` is the setting with nothing to argue about.

**protocol.** `cold` = no LLM. `retain25` = keep GLM-5's 25 best-scoring proposals per warm
round and let the arm fill the other 75 (7.5% of the budget). `warm` = the whole recorded
round, shortfall lost (29.1%). Rounds 4–10 are identical in all three. Both LLM protocols
are **replayed from recorded proposals** — no API call is made.

## The columns that need explaining

| column | meaning |
|---|---|
| `ef_screen` | **canonical adjusted EF**, plain screen mean. `hits / (eff · H/D)`, `D` = the screen's own library, `eff = 1000 − n2` (off-library-but-in-pool picks refunded, unfilled slots charged). This is the paper's metric: our AssayFormer reads 4.828 against their published 4.83. |
| `ef_pub` | the same metric averaged over publications. Differs from `ef_screen` by up to 1.4 where one publication dominates a set. |
| `ef_nvr` | `n_hits_vs_random` = `hits / (1000 · H/pool)`. What every `result.json` stores. **Not** the paper's metric — never put it beside 4.83 / 5.27 / 5.67. |
| `d_ef_vs_af_depmap`, `d_ef_lo/hi` | paired contrast against `AF + DepMap` within the same (set, pool, protocol) block. **Publication-clustered**: the point estimate and the interval are both means of publication means. |
| `significant` | `*` when the interval excludes zero |
| `mean_off_library_picks` | per 1,000. Zero for every `library` row by construction; it is the size of the refund on `f2` rows. |
| `source_file` | which JSON supplied that arm, so mean=0 and real-mean bank variants can never be silently mixed |

**Clustering matters.** Screens from one publication are near-duplicates. `paper_validation`
is 20 screens over 10 publications with one supplying 8; the test fold is 146 over 34.
Screen-level intervals are too narrow, so everything here clusters by publication.

## What the table says

- **`library` pool, cold:** `DepMap + HSA dm.50` beats `AF + DepMap` on the two curated
  panels (+0.13, +0.15) and on team11 (+0.45\*), loses on the validation fold (−0.20) and
  the held-out test fold (−0.23). Only team11's interval excludes zero. **HSA finds more
  hits on every set.**
- **Protocol, paper_test:** every arm peaks at `retain25` on both EF and hits, and falls at
  the full warm start. The LLM is worth having at 7.5% of the budget and costs EF at 29%.
- **Pool:** moving f2 → library costs AF-based arms 0.19–0.28 EF and HSA arms ≤0.11, because
  the refund was discounting them unequally.

## Reproducing from scratch

Everything below is **free** — no API calls. The LLM protocols replay recorded proposals and
the runners assert that no `llm_calls.jsonl` appears.

**External assets, none of them in git:**

| asset | where it comes from |
|---|---|
| AssayFormer checkpoint (`model_last.pt`, vocab, config) | `Genentech/assayformer` on HuggingFace — use `model_last.pt`, not `model.pt` |
| AssayBench screens | `Genentech/assaybench` on HuggingFace, pulled by `load_screens` |
| reduced DepMap atlas `achilles_noeval.npz` | built by `scripts/diagnostics/build_noeval_atlas.py` |
| HSA banks | rebuilt locally — see below |
| GLM-5 recorded rounds | `output/runs/baseline/glm_51/` and `output/runs/.../sweep-glm5tu-pt-af-*` |

```bash
# 1. donor matrices, then the banks each target set needs (~5 s per target)
uv run python scripts/external/build_donor_matrices.py \
    --checkpoint _reproduction/bundle/checkpoint \
    --reference-bank _reproduction/bundle/banks/1566.npz \
    --out _reproduction/donor/train_1349.npz
uv run python scripts/external/rebuild_banks.py \
    --donor _reproduction/donor/train_1349.npz --verify 1566        # must print 1.0000
uv run python scripts/external/rebuild_banks.py \
    --donor _reproduction/donor/train_1349.npz \
    --manifest _reproduction/donor/test_fold_334.yaml --target-set test \
    --skip-existing --out _reproduction/banks_real_test_fold

# 2. the own-library runs (3 min for 20 screens x 5 arms; ~20 min for 146)
ARMS="AF alone,AF + DepMap,HSA only (ours),DepMap + HSA dm.25,DepMap + HSA dm.50"
PYTHONPATH=_reproduction:scripts/external uv run python scripts/external/run_paper_test_hsa.py \
    --banks _reproduction/banks_real_test_fold --checkpoint _reproduction/bundle/checkpoint \
    --screen-set test --per-screen-pool --restrict-to-bank --min-genes 18000 \
    --arms "$ARMS" --out output/runs/external/test_fold_hsa_libpool.json

# 3. the LLM protocols (paper_test; replayed from disk, nothing billed)
PYTHONPATH=_reproduction:scripts/external uv run python scripts/external/run_glm5_handoff_hsa.py \
    --banks _reproduction/banks_real_paper_test --checkpoint _reproduction/bundle/checkpoint \
    --retain 25 --out output/runs/external/glm5_retain25_hsa.json

# 4. rebuild this table
uv run python scripts/diagnostics/build_main_table.py
```

`--min-genes 18000` is required on a whole fold: the 218-screen validation fold has median
library 2,318 genes and 98 screens at or under 1,001, where a 1,000-pick budget exhausts the
library and EF is meaningless. It keeps 69 of 218 and 146 of 334.

## The code

| file | role |
|---|---|
| `scripts/external/run_paper_test_hsa.py` | every own-library and f2 run; `--per-screen-pool`, `--restrict-to-bank`, `--min-genes` |
| `scripts/external/run_glm5_handoff_hsa.py` | the GLM-5 warm start and the `--retain` proposal cap |
| `scripts/external/rebuild_banks.py`, `build_donor_matrices.py` | HSA banks from public data |
| `scripts/external/run_bundle18.py` | the `ARMS` table — the weight triple behind each arm name |
| `src/assayloop/acquisitions/glm_handoff_acq.py` | the `retain`/`rank` proposal cap |
| `scripts/diagnostics/build_main_table.py` | this table |
| `scripts/diagnostics/score_libpool.py` | one own-library run, with the one-metric identity asserted per row |
| `scripts/diagnostics/libpool_three_sets.py` | all four sets as levels, plus each arm's f2→library move |
| `scripts/diagnostics/score_retain25.py`, `score_glm5_handoff_hsa.py` | the LLM protocols |
| `notebooks/compare_hsa_depmap.ipynb` | figures; section 7 reads these runs directly |
| `tests/test_glm_handoff_retain.py` | pins the proposal cap |

## Two traps that already cost time

1. **`n_hits_vs_random` is not adjusted EF.** It is what every `result.json` stores and what
   notebook cells 0–37 plot. The two rank arms differently — `dm.50` reads 5.863 vs
   `AF + DepMap`'s 5.670 in that unit, and 5.472 vs 5.532 in the paper's. Use `ef_screen`
   or `ef_pub` from this table.
2. **Bank variants.** HSA banks exist as `mean=0` and real-mean versions. Every own-library
   run uses `banks_real_*`; comparing one against a mean=0 file mixes variants and moved one
   arm by 0.10 EF. `source_file` in this table is the audit trail.
