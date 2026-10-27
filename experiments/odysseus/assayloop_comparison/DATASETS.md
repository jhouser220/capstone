# Dataset provenance and comparison boundaries

The common evaluation system is AssayLoop with a fixed 20,893-gene candidate pool. Each external source supplies measured outcomes for selected assays; no claim is made to evaluate every assay in each publication. Candidate-pool coverage is not the fraction experimentally queried: each campaign requests only 1,000 unique candidates. Missing outcomes remain unknown and are charged to the request budget.

| Source / experiment | Biological endpoint | Selected assays | Primary source | Contract limitation |
|---|---|---:|---|---|
| Schukken / E224 | Aneuploid cancer-cell dependencies | 1 | [Primary article](https://pmc.ncbi.nlm.nih.gov/articles/PMC13142536/) | Declared source-derived RPE1-WGD/ESL endpoint, not every author screen |
| Abt–Wang / E225 | PARCB cell-fitness dependencies | 1 | [PNAS](https://doi.org/10.1073/pnas.2532814123) | Endpoint/sign and high-prevalence qualifications retained |
| Wilson / E226 | Mitochondria–ER contacts | 1 | [Cell Death & Disease](https://doi.org/10.1038/s41419-024-06568-y) | Source table versus paper hit-count mismatch; not exact biological label parity |
| Nguyen / E232 | AraC sensitization | 1 | [Primary article](https://pmc.ncbi.nlm.nih.gov/articles/PMC11914169/) | Modest mean gain (+1.5 hits); two of four matched paths lose |
| Shan / E236 | Lovastatin response and LDLR knockout | 3 | [Genome Biology](https://doi.org/10.1186/s13059-025-03531-8) | Three assays from one publication; strong-control promotion gate failed |
| Herranz-Ors / E239 | Colorectal, oesophageal and ovarian tumour organoids | 3 | [Nature](https://doi.org/10.1038/s41586-026-10830-y) | Three cultures from one publication; strong-control promotion gate failed |
| Samelson / E245 | Neuronal tau lowering | 1 | [Primary article](https://pmc.ncbi.nlm.nih.gov/articles/PMC12978015/) | CRISPRbrain T22 release used; final-paper annotation/version parity unresolved |

The JSON source inventory records exact screen IDs, publication identities and source-table counts; frozen target inputs and transformation scripts live in the corresponding experiment directories of the private snapshot. Rebuild from the declared source release and ID/threshold contract before substituting inputs. This public package does not grant redistribution rights for source datasets or weights.

## Strong-control context

| Source | Language + HSA | HSA without language | MAP in the same 42D representation |
|---|---:|---:|---:|
| Schukken | 50.00 | 52.50 | 44.00 |
| Abt–Wang | 365.50 | 369.00 | 363.00 |
| Wilson | 107.75 | 115.00 | 111.00 |
| Nguyen | 15.00 | 18.00 | 20.00 |
| Shan | 86.83 | 90.67 | 90.33 |
| Herranz-Ors | 712.00 | 717.00 | 713.33 |
| Samelson | 155.00 | 152.00 | 153.00 |

These are descriptive publication means with different repetition structures for cold controls. Do not treat them as independent paired LLM samples. Frequency-based proposals are another strong comparator, especially on PARCB and lovastatin; no universal HSA superiority is established.

## Original AssayLoop panel

E160 covers 18 selected, development-exposed publications (3 original test, 6 additional test, 6 validation, 3 training), not the original paper's full official 20-test-screen benchmark. HSA beats the budget-filled comparator on 11/18 publication means; positive cohort averages hide heterogeneous and negative cases. Training is diagnostic, and extensively inspected test data are no longer independent confirmation.

## Research basis

- [AssayLoop: Biology-in-the-loop](https://arxiv.org/abs/2609.11877): original adaptive loop, AssayFormer and source acquisition machinery.
- Empirical Gaussian Processes (ICML 2026): historical response structure as a prior; inspiration, not an implementation-equivalence claim.
- Residual Feature Integration is Sufficient to Prevent Negative Transfer (ICLR 2026): source-preserving additional features; a comparator/inspiration, not evidence that our normalized basis is strictly orthogonal.
- Many Needles in a Haystack (ICML 2026): hit-discovery objectives and evaluation.

HSA appends historical response directions, fits task coefficients to selected measured outcomes and integrates hit probability with availability. Whole-package comparisons also change proposal retention and numerical scoring. A matched intervention is needed to assign the gain specifically to one component.

## Audit limits

Fresh packaging checks reproduce recorded arithmetic and selected actions. They do not certify unknown checkpoint/LLM training ancestry, paper-level source-label equivalence on every source, future-publication generalization or absence of every possible form of leakage. Requested-slot AUC and source measured-axis nAUC are different metrics. Never use one label for both.
