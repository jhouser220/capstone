# Where every external screen came from

Seven papers, 206 directional screens. Each row below is the exact file that
`build_*.py` reads, the exact column or sheet the hit calls come from, and the
route that worked when the publisher's own site blocked.

Every file is now collected under **`output/runs/external/sources/`** (gitignored,
so nothing here bloats the repo). The builders take explicit `--` paths, so point
them anywhere.

---

## 1. Schukken 2026 — aneuploid vs near-euploid knockout fitness

| | |
|---|---|
| paper | Schukken, Aitken, Sheltzer et al., *Paired CRISPR screens...*, PMC13142536 |
| data | `SupData/SD3- CRISPR Screen Beta Scores.xlsx`, sheet **`Raw WGL`** |
| from | https://github.com/kschukken/Aneuploid_CRISPR_Screens |
| local | `sources/schukken_SD3_beta_scores.xlsx` |
| screens | 14 (7 cell-line pairs) |
| hit rule | **mine, not theirs**: MAGeCK `beta < −0.75`. The paper reports no gene-level call list. −0.75 was chosen to match the AssayBench genome-scale knockout-fitness median hit fraction (0.094; this gives 0.086). `schukken_b05_results.json` is the same thing at −0.5, as a sensitivity check. |

This is the **only** dataset whose hit labels I defined. Every other one uses the
authors' own calls.

## 2. Abt 2026 — PARCB small-cell prostate, E2F3 synthetic lethality

| | |
|---|---|
| paper | Abt, Wang, Varuzhanyan, Freeland et al., PNAS 2026, doi 10.1073/pnas.2532814123, PMC13012052 |
| data | `Zhang_count/merged_spacer_counts_mageck_test.gene_summary_FILT.txt` |
| from | https://github.com/jfreeland01/Abt_Wang_E2F3_CRISPR |
| local | `sources/abt_gene_summary_FILT.txt` |
| screens | 1 |
| hit rule | authors' own: `neg.p.value < 0.05` → exactly 1,398 hits |

pnas.org returns 403 to scripted requests; the authors' GitHub has the final
gene-level MAGeCK output, already filtered to ≥3 high-quality sgRNAs, which is
the file their own figures use. Its `README.md` documents the whole pipeline.
Only one screen because the published table is the 9-replicate pooled MAGeCK
test — a per-line re-run would give n=3 but needs the raw counts.

## 3. Wilson 2024 — MERCS split-mVenus reporter, FACS-sorted

| | |
|---|---|
| paper | Wilson, Yu, Leal, Woodward, Patikas et al., *Cell Death & Disease* 15:203 (2024), PMC10928201 |
| data | `Supplementary Table 2_Whole Genome Pooled Screening Results.csv` |
| from | https://github.com/m1gus/MitoER-CRISPR-Screen |
| local | `sources/wilson_ST2_screen_results.csv` |
| screens | 3 (`_inc` 394, `_dec` 229, `_merged` 623 — 2 independent) |
| hit rule | authors' stated **FDR < 1e-4** ("FDR = 0.01%") on MAGeCK RRA `pos_fdr` / `neg_fdr` |

Their quoted 410/230 is **not reachable at any threshold**: MAGeCK RRA
permutation FDRs are heavily tied (357 genes share 1.4e-05) and the achievable
counts step 394 → 412. The stated cut is used as stated; residual under 4%;
`--hit-fdr` moves it. Symbols are resolved **ENSG → HGNC first**, which also
repairs Excel date corruption (`Mar-01`/`Mar-02` for MARCHF genes).

## 4. Nguyen 2025 — AraC / daunorubicin / etoposide drug response

| | |
|---|---|
| paper | Nguyen, Rafiee, Parcha ... Lamba, *Blood Advances* 9(5):1078 (2025), doi 10.1182/bloodadvances.2024014157, PMC11914169 |
| data | supplemental `mmc2.xlsx` — sheet **S1** (19,114-gene diff-beta) + sheet **S2** (the calls) |
| from | Europe PMC `supplementaryFiles` for PMC11914169 (ashpublications.org → 403) |
| local | `sources/nguyen_mmc2_supplemental_tables.xlsx` |
| screens | 9 (3 drugs × `sensitive` / `resistant` / `merged` — 6 independent) |
| hit rule | authors' own per-drug **`Response`** labels in sheet S2 |

Their stated "beyond 2 SD" rule reproduces **none** of their six counts (any
plain reading gives 8/976, 731/221, 120/771 against the published 240/1800,
1229/258, 296/1410). Sheet S2's explicit labels reproduce all six exactly.
Sheet S2 also went through an R step that mangled `-` to `__` (`NKX3__1`,
`HLA__C`, `ATP5J2__PTCD1`); 34 of 3,748 calls silently drop without the
unmangle step.

## 5. Shan 2025 — cholesterol / lipid homeostasis synthetic lethality

| | |
|---|---|
| paper | Shan, Fan, Li ... Fei, *Genome Biology* 26 (2025), doi 10.1186/s13059-025-03531-8, PMC11912599 |
| data | **Additional file 2** (Table S1, "CRISPR screen data") — `Screen N_all` + `Screen N_hits` sheets |
| from | Europe PMC `supplementaryFiles` for PMC11912599 |
| local | `sources/shan_additional_file2_TableS1.xlsx` |
| screens | 3 (HepG2 + lovastatin 368, HeLa + lovastatin 359, HeLa LDLR-KO 227) |
| hit rule | authors' own, read straight off the `_hits` sheets |

Taking the sheets avoids reconstructing their thresholds, which matters because
screen 3's stated rule references a "tRRA" score whose column is not in the
workbook under that name. Screen 1's 368 reproduces the paper's text.
**Weakest novelty test of the seven**: statin-response screens in other lines
(`U_1182` K-562, `U_1321` hTERT-RPE1, `U_1503` HAP-1) and a cholesterol/LDLR
screen (`1434`) are in the train fold, so the model has seen the pathway even
though it has not seen these screens.

## 6. Herranz-Ors 2026 — Sanger tumour organoid biobank (the big one)

| | |
|---|---|
| paper | Herranz-Ors, Bhosle, Beck ... Garnett, *Nature* (2026), doi 10.1038/s41586-026-10830-y, PMC13581617 |
| data | **Supplementary Table 5** (`sample_ID gene ensembl_id HGNC_ID library is_depleted`, 2,948,036 rows) |
| from | the article's Supplementary Information (Springer `static-content`); covariates from the paper's Figshare item, doi **10.6084/m9.figshare.28339340** |
| local | `sources/organoid_suppl_table5.csv` (164 MB), `sources/organoid_data_covariates_ALL-organoids.csv`, plus tables 1 and 4 |
| screens | 178 = `(sample_ID, library)` over 162 organoids; 162 bootstrap clusters |
| hit rule | authors' own **BAGEL2** fitness calls (`is_depleted`) |

**The trap.** The obvious-looking file on that Figshare item,
`07_EssMatrix_bDepletionsB2.Rdata` (ndownloader `63130879`), is **not** this
paper's organoids: it is 17,647 × 930 with column names like `ACH-001421` —
DepMap/Achilles **cell-line** identifiers, the Pacini et al. reference matrix the
paper compares against. Scoring against it would pit the co-essentiality atlas
against the very cell lines it is computed from and manufacture a large fake
gain. Organoid calls live **only** in Supplementary Table 5, whose `sample_ID`
values are `HCM-SANG-*`.

The unit must be `(sample_ID, library)`: 16 of the 162 organoids were screened
with both MinLibCas9 and Yusa v1.1, and grouping by sample alone doubles their
genes (33,124 rows where a single-library model has 16,562). Split by library,
hits run 297–2160 (median 1440), which is the paper's stated range exactly.

## 7. Samelson 2026 — tau proteostasis, genome-wide CRISPRi *(the new one)*

| | |
|---|---|
| paper | Samelson, Ariqat, McKetney et al., *CRISPR screens in iPSC-derived neurons reveal principles of tau proteostasis*, **Cell** 189:1517–1534.e19 (2026), doi 10.1016/j.cell.2025.12.038, PMID 41610849 (Kampmann lab) |
| data | `mmc2.xlsx` (all genes, 7 sheets) + `mmc3.xlsx` (called hits, same 7 sheets) |
| from | **supplied by hand** — not open access, and PMC12978015 is still inside the NIH embargo |
| local | `sources/tau_mmc2_all_genes.xlsx`, `sources/tau_mmc3_hits.xlsx` (also `mmc2.xlsx` / `mmc3.xlsx` in the repo root) |
| screens | 3 from the `GenomeWide` sheet (`_inc` 366, `_dec` 729, `_merged` 1094 — 2 independent) |
| hit rule | authors' own, read from `mmc3` |

Routes that **failed**, for the record: cell.com and sciencedirect 403;
Europe PMC `fullTextXML` HTTP 500; Europe PMC `supplementaryFiles` → *"not open
access"*; PMC web page → browser challenge; GEO → only a different tauopathy
paper (GSE255902/255904); CRISPRbrain `/api/screens` → 401 and `/api/v1/screens`
→ an SPA shell. **crisprbrain.org** is the authors' stated portal and would add
sgRNA-level data, but the gene-level phenotypes and calls needed here are
entirely in mmc2/mmc3.

Only the `GenomeWide` sheet is scorable — the other six sheets are the follow-up
screens on a 984–986-gene focused sublibrary, and the budget is 1,000 picks, so
both arms would acquire the whole library and tie by construction.

The hit calls are **not reconstructible** from mmc2's columns: `gene score` is
`epsilon · −log10(pvalue)`, but 1,017 genes clear the weakest called negative
while only 729 are called, and genes at |score| 1.2 are excluded while genes at
0.51 are kept. That is the Kampmann discriminant rule, calibrated against the
NTC distribution **separately per sublibrary**, and the table does not say which
sublibrary a gene is in.

---

## The benchmark comparison row

`benchlib_results.json` is the 348 held-out AssayBench screens re-run through the
same driver with `universe_genes` set to each screen's own library, so it sits in
the same metric regime as the seven papers above rather than on its native
fold-union pool (where 562 of every 1,000 picks are forgiven and `b_eff ≈ 437`).
Screens come from `dump_benchmark_screens.py`; families from
`benchmark_stratified.json`.

## Shared reference files

| file | path |
|---|---|
| HGNC complete set (`symbol`, `status`, `prev_symbol`, `alias_symbol`, `ensembl_gene_id`) | `/home/jeff/hgnc/hgnc_complete_set.txt` |
| AssayFormer vocabulary (24,217 tokens) | `<assayformer snapshot>/vocab.json` |
| DepMap reduced atlas — **the only atlas used, every run** | `/home/jeff/depmap/achilles_noeval.npz` |

## Rebuilding any of them

```bash
cd AssayLoop
export PYTHONPATH=src:scripts/external
SRC=output/runs/external/sources
V=/home/jeff/.cache/huggingface/hub/models--Genentech--assayformer/snapshots/4927703e5903984ae0254a9fe8e93ec0b29c27e4/vocab.json

.venv/bin/python scripts/external/build_tau.py \
    --all $SRC/tau_mmc2_all_genes.xlsx --hits $SRC/tau_mmc3_hits.xlsx \
    --hgnc /home/jeff/hgnc/hgnc_complete_set.txt --vocab $V --out /tmp/tau.pkl

.venv/bin/python scripts/external/run_external.py \
    --screens /tmp/tau.pkl --n-steps 10 \
    --arms "AF alone,AF alone (greedy),AF + DepMap" \
    --tag tau --out output/runs/external/tau_results.json
```

Each builder prints its gene counts, hit counts, hit fraction and EF ceiling, and
raises if its gene list is not unique — duplicate symbols make `reveal()` collapse
entries and corrupt the hit labels.

---

## The team's `external_sources.json` (11 frozen screens)

`scripts/external/team_subset.py` maps every entry to a screen built here and
prints the AF vs AF + DepMap comparison for just those 11. Ten matched on
`(source_genes, source_hits)` with no ambiguity. The eleventh,
`E224_SCHUKKEN_PAIRED_WGL`, is the Schukken paper's **paired**
aneuploid-vs-near-euploid contrast — a screen this project had not built, because
`build_schukken.py` reads the 14 WGL columns as independent dropout screens
instead. `build_schukken_paired.py` builds it from the quantile-normalised `WGL`
sheet (19,112 genes).

Its 194-hit count is not reproducible from any round threshold, so that builder
carries two readings — the paper's stated `p<0.05 & beta<-0.2` (188 hits) and the
top 194 by the authors' `AneuploidRank`. They share only 84 genes and they
disagree in sign, so E224 needs the teammate's actual hit list before it is quoted.
