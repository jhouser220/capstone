# Packaging validation — 1 October 2026

No paid API requests were made. This validates the portable package, not a new research hypothesis.

- Parsed all 14 Python files successfully.
- Existing repository screen CSV/NPY files match their local research counterparts byte for byte.
- Rebuilt Achilles features from the existing public-source CSV using `prepare_features.py`: 18,119 genes × 778 features, float32. Gene order and all feature values are exactly equal to the original normalized cache. The compressed file checksum differs because archive serialization differs.
- Source CSV SHA256: `71f32e3515dfa24d3bff559ea99fd5a0fa102c264064a6377e13d8eecd54298f`.
- Rebuilt NPZ SHA256: `9492ddd06b89c4f9da90c3e2937b54e0c12b396829760dc694b5bb9337507820`.
- `check.py` passed with the existing research feature file and with the rebuilt default file: hidden-outcome mutation leaves predictions unchanged for four learners; selections contain 32 distinct eligible genes; report neighbours belong to revealed history; both API clients reject fresh paid requests before network access when disabled.
- Ran four three-batch IFNG campaigns, random initial history, seed 8101, with the rebuilt features and default repository datasets. Initial histories and all selected batches exactly match archived phase2 results:

| Policy | Newly selected hits / 96 | Historical match |
|---|---:|---|
| Geometry baseline | 16 | Exact choices in all 3 batches |
| GP | 15 | Exact choices in all 3 batches |
| Full-feature effect median | 9 | Exact choices in all 3 batches |
| Full-feature hit classifier | 4 | Exact choices in all 3 batches |

These are reproducibility checks of existing cases, not independent confirmation of effectiveness. Runtime outputs and features remain locally ignored; historical scientific results are bundled separately in `recorded_results/`.

The feature builder ran with pandas 2.3.3 / numpy 2.4.6. Learner checks ran in the existing research environment (Python 3.11.1, numpy 2.4.6, scipy 1.17.1, scikit-learn 1.9.1, XGBoost 3.2.0). A clean-environment install was attempted but ran out of disk space; that incomplete environment was removed. Therefore clean installation and other operating systems are not yet verified. No LLM response reproduction is claimed, and current API pricing/model availability was not checked because no paid calls were made.

Independent biological feature provenance and assay-overlap checks remain open despite the successful numerical rebuild. Full historical response logs, raw prediction matrices and old accounting ledgers are intentionally not in this focused package.
