# AssayLoop warm-start investigation: handoff

Last updated **2026-10-27**. Companion to `CLAUDE.md`, which covers the original reproduction setup
and is still accurate for environment and patch details. This file covers the warm-start improvement
work: what was built, what was measured, and what the numbers mean.

## Current state

**THE TWO RESULTS.**

- **The paper's LLM warm start does not reproduce.** Thirteen variants, two models, nine prompts,
  all at or below AssayFormer alone (§1, §4). Three competing explanations were refuted
  experimentally: model capability (GLM-5 ~ Flash-Lite), unfilled slots (top-up fixed the shortfall,
  EF unchanged), prompt engineering (nine variants).
- **A non-LLM intervention works: `--acq codep_seed` with
  `query=diffmeans weight=0.5 knn_weight=0.5`.** AssayFormer (frozen) fused by reciprocal rank with
  (a) DepMap co-essentiality of the hits it has already confirmed and (b) hit frequency in the most
  similar training screens. **CONFIRMED ON TWO INDEPENDENT HELD-OUT FOLDS (§4.46):** validation
  198 screens **+0.175 [+0.077, +0.290]** per screen / **+0.331** per family; test 314 screens
  **+0.257 [+0.050, +0.476]** / **+0.398**. Genome-scale only (the honest filter, §4.29) the two
  folds land on top of each other: **+0.485 [+0.281, +0.730] and +0.484 [+0.130, +0.917] per
  family.** Every interval excludes zero. On the curated `paper_test` 20 it reads +0.562, but those
  screens chose the weight. **Quote the two-fold figures, or +0.485/+0.484 with the filter stated.**
- **BUT THE GAIN IS NOT ROBUST ACROSS BUDGETS, and that is new as of 2026-10-23 (§4.46).** The
  endpoint replicates; the per-round curve does not. On the test fold the method is **BEHIND from
  300 to 800 picks, by up to −1.44 EF at round 7**, overtaking only in rounds 9-10 — while on
  validation it leads at every budget. dnAUC over the ten budgets: validation **+0.144
  [+0.053, +0.241]**, test **−0.397 [−0.855, −0.008]**, both excluding zero with OPPOSITE signs.
  **The headline holds at the paper's metric — ten rounds, EF at 1000 — and any claim must say so.
  Someone stopping early would be worse off.** Report EF-at-1000 and dnAUC together; a
  single-endpoint metric hid this for the whole project.

**WHY THE NEGATIVE HOLDS.** AssayFormer compresses a screen to 10 numbers, and the best possible
such compression — fitted to the full answer key — is worth only +0.41/+0.84 EF, so every warm
start, re-ranker and context trick competed for a prize that barely exists (§4.7). The LLM's
gene-level knowledge is already inside the BPMF gene factors: AUROC 0.55-0.58 standalone, **zero**
incremental (§4.11). One mechanism explains every null.

**WHY THE POSITIVE HOLDS.** Nested baseline (`weight=0` reproduces AssayFormer pick for pick, all
1000 picks verified); coherent dose-response with the same optimum on both screen sets; measured out
of sample; and three leakage audits passed — training labels (§4.22: gain is strongest on the
*least* contaminated screens), DepMap (§4.35: deleting each screen's own cell-line column changes
the margin by +0.010 — and since §4.56 the **default atlas has those 40 lines removed outright**, so
no screen's own cell line is in it by construction rather than by audit), text embedding (no label
contact by construction).

**THE TASK HAS HUGE ROOM, AND IT IS EVERYWHERE.** On the 100 genome-scale held-out screens
AssayFormer averages EF 4.67 against perfect play's 14.81 — a gap of **~10 EF per screen, about 3x
the current score**. Very-hard screens carry 37% of it, middling 27%, and **screens where
AssayFormer already does well still carry 28%** (§4.12 correction). Caveat for any ceiling claim:
**98 of the fold's 198 screens are metric-degenerate**, their library smaller than the 1000-pick
budget, so perfect-play EF computes below 1 — exclude them.

**ONE HEADER CLAIM WAS WITHDRAWN, 2026-10-23.** §4.39 said gain-of-function "may not support
co-functional clustering at all". §4.41 repaired that control — it was both mis-specified (complexes
are the wrong ground truth for an overexpression atlas, by §4.39's own stated mechanism) and
underpowered (six complexes) — and the clustering is present and significant. JUMP is still rejected,
on a measured **4.5x** signal deficit against DepMap rather than on a failed control.

**WHY DEPMAP WORKS AND NOTHING ELSE DOES** — five conditions, each measured (§4 synthesis):
perturbation phenotype not baseline state (expression, same machinery, is **harmful** at −0.109,
§4.24); readout matching the screen (fitness screens +0.766 vs non-fitness +0.302, §4.35);
perturbation DIRECTION matching (knockout +0.546, CRISPRi +0.397, **CRISPRa +0.022**, §4.30); many
contexts (808 lines; cutting to 34 costs −0.40, §4.37); and it supplies a SIMILARITY SPACE pointed
by the screen's own revealed labels, not an answer — PRISM uses the same DepMap data pointed by the
drug name instead and scores AUROC 0.500 (§4.27). **DepMap is the closest available proxy for the
experiment being predicted, and the hard screens are precisely those furthest from it.**

**ELEVEN FEATURE CANDIDATES SCREENED, ONE WINNER.** Expression harmful (§4.24); Perturb-seq
(§4.17), co-hit (§4.21), BioPlex PPI (§4.28), PRISM (§4.27), condition-matched transfer (§4.34),
PPMI/SVD co-hit (§4.32), direction-filtered co-hit (§4.31), lineage-restricted DepMap (§4.37), JUMP
ORF (§4.39) and LINCS OE (§4.36) all fail, each on a named condition. **"More training screens" is
not available: AssayBench already IS BioGRID ORCS** — 1901 of 1952 human screens (§4.33).

**EVERY TUNING KNOB IS ALREADY OPTIMAL.** Both fusion weights peak at exactly 0.5 on 198 held-out
screens (§4.19); `min_seeds = 2` is best in both directions (§4.20, §4.23); the saturated-screen
gate is a no-op, so **deploy ungated** (§4.13 correction); a rate-adaptive weight has no basis
(Spearman −0.045, §4.23). Weights are **absolute doses, not shares**: three scorers at 1/3 carry the
deployed total pressure and are significantly worse (§4.21).

**SIX RULES EARNED THE HARD WAY.**
1. **Offline margins overstate deployed gains — six times, once with a sign flip.** co-hit +0.198
   offline (CI excluding zero) -> +0.006 deployed; expression +0.87 -> −0.109. An offline CI
   excluding zero still predicts nothing. Use offline as a filter, never an estimate.
2. **Filtering data to the "matched" subset always loses — three for three.** Enrichment-only
   training screens (§4.31), rank-truncated co-hit (§4.32), pancreatic-only cell lines (§4.37).
   **Use everything.** Sample size beats specificity by 6-13x here.
3. **A matching DIAGNOSIS does not license a matching PRESCRIPTION** (§4.31). Hard screens ask
   enrichment questions, but enrichment-direction training data is worse than the pooled corpus.
4. **Check gate 2 (does the atlas have values for the hits you need to rank?) BEFORE recommending a
   feature.** It is structural, costs minutes, and killed LINCS after three recommendations (§4.36).
5. **Check that a data source holds something you do not already have** before arguing for it —
   three sections went into ORCS before anyone compared screen ids (§4.33).
6. **Use `prescreen_atlas.py` for any new atlas**, and read GATE 0 first: a positive control on
   known complexes. It is what diagnosed JUMP (co-complex genes no more similar than random pairs)
   and what would have caught §4.17's 6832 infinities immediately.

**WE MATCH THE PAPER'S HEADLINE WITHOUT AN LLM, AND ITS OWN DATA SAYS WHY (§4.52, 2026-10-27).**
The user supplied `paper_per_round_data.xlsx`, the paper's 59 methods round by round on its 20 test
screens — which are exactly `paper_test`, and on which **our AssayFormer baseline reproduces the
paper's own baseline row to the hit on all 20** (mean |diff| 0.000). On those screens the paper's
Gemini-3.1-Pro + AssayLoop scores **5.48 [5.265, 5.699]** and **our method scores 5.389 with no
language model at all.** And **Gemini 3.1 Pro ALONE scores 4.774, BELOW AssayFormer's 4.828** — the
paper's gain comes from the handoff structure, not from the LLM being a good picker. **Two traps
when reading that file: its row named `AssayFormer` is NOT the released model (EF ~2.4; the
released one is `AssayFormer (BPMF + GRPO)`), and its EF can only be bracketed, not computed
exactly.**

**THE LLM AND OUR METHOD ARE SUBSTITUTES, CONFIRMED FROM BOTH DIRECTIONS (§4.48-§4.51).** Give the
campaign PERFECT early hits and our method's benefit vanishes (+0.331 -> −0.074, §4.50). Give it
WORSE early hits by letting an LLM pick rounds 1-3 and our method contributes the same as ever
(+0.549 from an LLM start vs +0.562 from a cold start, §4.51). The method is insensitive to its
seeds, so the LLM's displacement cost passes straight through: **LLM + our method beats the paper's
handoff clearly, and still loses to not using the LLM at all.**

**THE INPUT-SIDE LEVER IS SPENT, 2026-10-23.** §4.26(a) called the unread description
metadata "the largest untried lever in the project". §4.45 found `direction_str` is empty on every
screen, and §4.47 measured the rest: `notes` is flat with a ±0.045 interval, the full field set is
significantly HARMFUL, and in round 1 — where the description is all the screen latent has — the
variants cost **−0.52 to −0.86 EF**. A richer description produces a worse `u`, because the
checkpoint was trained on the renderer as it stands. **Nothing on the input side remains that does
not require retraining.**

**WHAT IS LEFT — REOPENED 2026-10-23 by external review (§12).** The claim that "the candidate list
is closed" was too strong. It was closed against *new data sources*, which is where ten of eleven
candidates were spent. Two axes were never touched at all: **the geometry of the atlas** (cosine is
computed on raw 808-dim profiles; no PC removal, no variance weighting) and **the statistical
instrument** (EF at one budget, unadjusted, on a fold that has now served as both screening and
confirmation set). A third, §4.39's gain-of-function null, rests on a positive control that may be
wrong *by its own stated mechanism*. Six live items, ranked in §12; the `test`-fold confirmation
(§10 item 1) is still the only one that strengthens what already exists.

**DO NOT REBUILD** the long list at the end of §10. It includes four failures of the "hits span
multiple modules" idea, five of "look the condition up in a database", and three of "filter the data
to the relevant subset".

Read `CLAUDE.md` first for repo background, then this.

---

## 1. The question and the answer

**Question.** The paper reports that an LLM warm start (rounds 1-3) followed by AssayFormer
(rounds 4-10) beats AssayFormer alone: EF 5.67 vs 4.83 on `paper_test`, i.e. +0.84. Does that
reproduce, and can the warm start be improved by aiming it at genes AssayFormer misses?

**Answer: the LLM warm start does not reproduce, and thirteen interventions did not fix it — but a
non-LLM one does work, and it has been confirmed out of sample.** Headline, on 198 screens that
played no part in choosing any setting (§4.15): AssayFormer fused with **DepMap co-essentiality of
its own confirmed hits** plus **hit frequency in the most similar training screens** gives
**+0.175 EF per screen [+0.075, +0.293], +0.331 per screen family** — and **+0.485 per family on the
100 genome-scale screens** once the metric-degenerate ones are excluded (§4.12 correction, §4.29).
On the 20 curated `paper_test` screens it reads +0.562, but those screens chose the weight, so the
held-out figure is the one to quote. No LLM, no retraining, and a nested baseline verified pick for
pick — all 1000 picks. Every LLM warm start tested is at best level with AssayFormer alone and
usually worse; model capability and unfilled slots were both ruled out by direct experiment. For
scale the paper claims +0.84 for its LLM warm start, and the oracle-û ceiling is +0.84 — so this
works only because it adds a dimension outside the model's own gene representation.

**Why it works, and why nothing else does** (§4 synthesis, five conditions each measured): DepMap is
a *perturbation* phenotype rather than a baseline state, its *readout* (fitness) matches what the
benefitting screens measure, its perturbation *direction* (loss of function) matches theirs, it has
*many contexts* so the gene-gene structure is estimable, and it supplies a *similarity space pointed
by the screen's own revealed labels* rather than an answer. Eleven candidate atlases were screened
against those conditions and each failure violates at least one. **DepMap is the closest available
proxy for the experiment being predicted, and the hard screens are precisely those furthest from
it** — wrong perturbation direction (every CRISPRa screen is hard, kernel worth +0.022), wrong
readout (surface marker, protein binding), or a mutant background nothing else shares.

**The decisive result is a ceiling argument** (§4 diagnostic 7). AssayFormer scores
`a(g) = û·v_g + b_g` with `d_gene = 10`, so a screen's reachable rankings are indexed by a 10-vector.
Fitting û to a screen's *complete* label vector — an oracle with the answer key — scores **6.510 on
validation and 5.668 on test**, against the released model's 6.100 and 4.828. So **every intervention
that works by giving AssayFormer a better û is bounded by +0.41 / +0.84 EF, and that bound requires
knowing all the answers.** Warm starts, rerankers, exploration schemes and context surgery are all in
that family. **But the bound is on û, not on the architecture** — §4 diagnostic 9 shows that appending
one informative additive feature raises the ceiling by +9.6 / +13.1 EF, because only 21% / 3.5% of the
hit vector is expressible in span(V, b, 1) at all. The constraint is feature quality, not capacity. Two supporting ablations: a screen-agnostic `b_g`-only ranking already scores 3.017 on
test, the description adds ~0.91, and all ten rounds of labels add ~0.90.

**This bound is on û, not on the task, and that distinction matters** (§4 diagnostic 12). Ranking
every true hit first would give EF **16.77** on validation and **18.92** on test, so AssayFormer
reaches only **42.7%** and **25.8%** of what the benchmark allows — and on the hard screens just
9-12%, finding 4.7-9.9% of the hits. The task has roughly **+10.7 EF** of headroom on average and
**+17.7** on the hard screens. The ~10 EF between the oracle-û ceiling and perfect play is the part
of the task this 10-number-per-gene representation **cannot express**; it is reachable only by
changing the representation (added features, or retraining `v_g`), never by a better latent. Note
the budget binds on the easy fitness screens, where AssayFormer already reaches 73-87% of ceiling,
so essentially all the available room is on the hard, low-density screens.

Two further consequences. The paper's reported AssayLoop figure, **5.67 on `paper_test`, equals the
oracle-û ceiling of 5.668** — an observation, not an accusation, but a label-free three-round warm
start matching an oracle is implausible, and it explains twelve null interventions with one fact
instead of twelve. And the only *input* with headroom that needs no retraining is the **screen
description text**, worth ~0.9 EF and currently rendered from 9 fields that omit `direction_str`.

**Headline numbers** (adjusted EF via the authors' own `_adj_ef_from_run`, budget 1000; "change"
is a mean over independent screen families with a 95% bootstrap CI — see §6 on families):

`paper_validation` (AssayFormer alone = **6.10**):

| Policy | EF | Change vs AF alone | Total hits | Better on |
|---|---|---|---|---|
| GLM-5 + top-up -> AF (n=3) | 5.22 | −0.52 [−0.94, −0.13] | −1076 | 5/20 |
| Flash-Lite -> AF (n=3) | 5.17 | −0.52 [−0.92, −0.18] | −1177 | 4/20 |
| GLM-5 plain -> AF (n=3) | 5.15 | −0.54 [−0.91, −0.21] | −1197 | 4/20 |

`paper_test` (AssayFormer alone = **4.83**, exactly the paper's reported 4.83):

| Policy | EF | Change vs AF alone |
|---|---|---|
| Decline-triggered handoff (r=0.5, max 3) | 4.90 | +0.08 [−0.51, +0.56] |
| Flash-Lite top-up -> AF | 4.70 | −0.13 [−0.79, +0.43] |
| Flash-Lite categories(8) -> AF | 4.68 | −0.15 [−0.74, +0.39] |
| GLM-5 -> AF (n=3) | 4.65 | −0.18 [−0.77, +0.30] |
| Flash-Lite categories(12) -> AF | 4.62 | −0.20 [−0.68, +0.31] |
| Flash-Lite -> AF (n=3), paper method | 4.57 | −0.26 [−0.85, +0.25] |
| Gemini 3.5 Flash-Lite alone, 10 rounds | 2.19 | −2.64 [−3.79, −1.50] |

---

## 2. Why: the diagnosis

**(a) The bar is AssayFormer's marginal pick, which is very high on validation.**

| Rounds 1-3, 20 validation screens | Hit rate |
|---|---|
| AssayFormer alone | **39.6%** |
| GLM-5 warm start | 13.4% |
| Flash-Lite warm start | 15.9% |

Each slot the LLM takes displaces a pick that hits ~40% of the time. On `paper_test` AssayFormer's
early-round rate is only 8.3%, which is why the losses there are much smaller.

**(b) Model capability is NOT the constraint.** GLM-5 (Bedrock `zai.glm-5`) and Gemini 3.5
Flash-Lite at minimal reasoning are statistically indistinguishable on both sets (5.15 vs 5.17
validation; 4.65 vs 4.57 test). High reasoning effort also did not help: Flash-Lite's hit rate went
9.1% -> 10.4% overall but 7.3% -> 5.8% on the enrichment screens that matter.

**(c) Unfilled slots are NOT the constraint.** `llm_single_topup` cut GLM-5's shortfall from 20.0%
to 2.1% (the paper reports 3.6% for Gemini-3.1-Pro, so this matches the paper's regime). EF moved
5.15 -> 5.22, still −0.52 vs AssayFormer alone. Hit rate stayed flat at 14.9%.

**(d) Most LLM hits are redundant with AssayFormer.** On validation, 62% of GLM-5's warm-start hits
were also found by AssayFormer alone **in the same screen** (Flash-Lite: 67%). Earlier cross-screen
analysis gave 65%, but that was inflated by near-duplicate screens; the within-screen figure from a
real AssayFormer-alone run is the honest one.

**(e) The genuinely new contribution looks model-independent.**

| Warm start, validation | Hits AF alone never found in that screen |
|---|---|
| GLM-5 | 245 |
| Flash-Lite | 247 |

Two models of very different capability converge on ~245 new hits over 20 screens (~12/screen)
against ~1730 hits given up. This suggests the ceiling is set by how much signal the **screen
description** carries, not by the model. Consistent with the benchmark's own eligibility rule
(screens admitted only if an LLM cleared AnDCG@100 >= 0.05).

**(f) Both models fail the same structural way.** Repeats/invalid symbols in rounds 1-3 on
validation: GLM-5 1137, Flash-Lite 1138. Accepted-per-round collapses (GLM-5: 1844 -> 1671 -> 1286).
The models re-propose genes they were explicitly shown. This is what top-up papers over.

---

## 3. The most interesting finding: selection direction

A screen's **selection direction**, readable from its own phenotype sentence, predicts almost
everything. Implemented in `src/assayloop/screen_direction.py` (`classify_selection`), verified
correct on all 40 curated screens.

- **dropout** — "decrease cell proliferation": a hit is a gene the *cells* need. Largely the
  essential-gene set, recurs across screens. AssayFormer: 66-95% hit rate.
- **enrichment** — "increases resistance to MERS-CoV": a hit is a gene the *condition* needs.
  Condition-specific, barely recurs. AssayFormer: ~5%, near each screen's base rate.

| Validation, 10 screens each | AF EF | AF hit rate | LLM warm-start hit rate |
|---|---|---|---|
| dropout / both | 9.5-9.7 | 66-67% | 25-31% |
| enrichment | **2.56** | **5.3%** | 4.3% |

Counts: validation 10 enrichment / 6 dropout / 3 both / 1 unknown; test 9 / 5 / 4 / 2.

**The LLM only ever helps on enrichment screens.** Restricting mixing to them (`--selection-mode
enrichment`) flips total hits from −145/−224 to +13/+9, and EF from −0.01/+0.05 to +0.02/+0.05 —
i.e. it removes the damage but does not create a gain. On `paper_test` the handoff is positive on
host-pathogen screens (+0.32 mean) and negative on drug (−0.73) and molecular-output (−0.60).

This is a cleaner version of E4 in the user's improvement brief (selection direction rather than
phenotype class), and it needs no knowledge of the hit count.

---

## 4. Interventions tried, and what each did

| Intervention | Where | Result |
|---|---|---|
| Paper method, n=3 handoff | `handoff_fullgenome.py` | −0.52 val, −0.26 test |
| Top-up (refill short rounds) | `llm_single_topup` (pre-existing) | shortfall 20% -> 2%, EF unchanged |
| Category planning (8 and 12 cats) | `llm_category_acq.py` (new) | ≈0 |
| Decline-triggered handoff | `decline_handoff_acq.py` (new) | +0.08 test, best single result |
| Mixed batches, 10-40 LLM genes/batch | `mixed_batch_acq.py` (new) | ≈0; 10/batch best (+0.03) |
| AF-coverage filter, horizon 100-1000 | `MixedBatchAcquisition.af_horizon` | hit rate 18% -> 9.7%, new-hit share 39% -> 84%, EF ≈0 |
| Category-coverage reweighting | `LLMCategoryComplementAcquisition` | **backwards** — see below |
| Enrichment prompt (direction + no hedging) | `ENRICHMENT_RULES` | wrong-sign picks 30 -> 0, hit rate 7.3% -> 5.7% |
| High reasoning effort | `-high.yaml` config | no gain where it matters |
| Enrichment-only gating | `--selection-mode` | stops losses, no gain |
| Stronger model (GLM-5) | Bedrock | indistinguishable from Flash-Lite |
| Pathway propagation (offline scoring) | diagnostic only | loses to AssayFormer on 16/20 screens — never built |

**Diagnostics run on 2026-10-21 (free, no LLM calls) on the hard screens:**

First, the framing fact: on the 7 `paper_test` screens where AssayFormer is near-random (EF < 2.5),
**Gemini alone is also near-random** — EF 1.18 vs 1.28, beating it on 2/7. The hard screens are hard
for both, so reallocating budget between them cannot help. Gene sets were fetched for these
(`scripts/fetch_gene_sets.sh`, public MSigDB + Reactome, now in `src/assayloop/data/gene_sets/`).

1. **Pathway propagation has real signal on enrichment screens, but only for real hits.** Seeding
   from the hits AssayFormer finds in rounds 1-3, taking small Reactome modules (<=150 genes) that
   contain >=2 of those hits, and scoring the untested members (~465 candidates/screen):

   | Target | Mean lift over base rate | Median | Screens >1.5x |
   |---|---|---|---|
   | all hits | 1.78x | 1.40x | 8/16 |
   | strong hits only (top quartile by relevance) | **3.53x** | **2.31x** | **11/16** |

   Best: 2466 10.1x, 1830 8.7x, 2048 6.8x, 2076 6.2x, TR_1975_1976 5.8x, 1561 4.2x. No lift on
   1815, 1818, 2054, 2104. Dropout screens get 2.1-3.7x, so the method is *also* good where
   AssayFormer is already strong — the value is specifically that it is complementary (Appendix F.3:
   BPMF embeddings organize by phenotype, not pathway). The PPI half of the test never ran:
   `reactome_interactions.txt` has UniProt/Ensembl/Entrez ids but no gene symbols. The user's own
   `../candidate_retrieval.py` does symbol-level PPI neighbours with 187 cached entries.

2. **Essentiality gating is BACKWARDS - do not build it.** Hypothesis was that a common-essential
   gene cannot produce a resistance phenotype. Using "hit in >=3 dropout screens" as the proxy
   (1560 genes), essential-proxy genes hit at **2.81% vs a 2.29% base rate on enrichment screens,
   ratio 1.56** (11/19 screens above 1.0; TR_1975_1976 6.05x). The LLM already spends 21-24% of its
   enrichment-screen picks there and that is mildly favourable. Gating would delete hits.

3. **The weak tail of the hit lists is largely noise**, which caps achievable EF on these screens:

   | | Top quartile of hits | Bottom quartile | Random genes |
   |---|---|---|---|
   | in any Reactome pathway | 63.9% | 53.4% | 50.9% |
   | shares a pathway with another hit | 41.9% | 24.4% | - |

   The bottom quartile is indistinguishable from random on annotation coverage. 1817 labels 888 of
   20,671 genes as hits (p <= 0.05 AND log2FC >= 2) including an olfactory receptor, a
   keratin-associated protein and a microRNA. EF/FH score prediction of that tail too, so part of
   the tier-3 gap is an unreachable target rather than a method failure. Only inversion: 1644, the
   one CRISPRa screen.

4. **Late seeding has no headroom, and this RETIRES the pathway scorer (diagnostic 1 used the
   wrong baseline).** Prompted by the idea of reversing the handoff — AssayFormer rounds 1-7, LLM
   rounds 8-10, fed the confirmed hits as seeds. Reproduced AssayFormer alone on validation
   (`sweep-afalone-val`, adjusted EF **6.100**, matching the recorded 6.10 exactly), then seeded the
   pathway scorer from its rounds 1-7 hits and scored the untested pool, comparing the top 300
   against AssayFormer's actual rounds 8-10 picks on the same 300 slots:

   | Seeds from | Subset | base | pathway candidate pool | pathway top-300 | **AssayFormer 8-10** |
   |---|---|---|---|---|---|
   | rounds 1-3 | all 20 | 3.80% | 8.29% (2.18x) | 21.30% (5.61x) | **37.42% (9.85x)** |
   | rounds 1-3 | 12 non-family | 2.32% | 3.40% (1.47x) | 4.31% (1.86x) | **6.86% (2.96x)** |
   | rounds 1-7 | all 20 | 3.14% | 5.80% (1.84x) | 13.93% (4.43x) | **30.73% (9.77x)** |
   | rounds 1-7 | 12 non-family | 2.22% | 3.01% (1.35x) | 4.53% (2.04x) | **5.89% (2.65x)** |

   Three conclusions, in order of importance:

   - **Diagnostic 1's 1.78x/3.53x was measured against random, not against AssayFormer.** Against
     random the scorer wins; against the policy it would replace or reweight, it loses everywhere
     (pathway top-300 beats AssayFormer on 4/20 screens, 836 hits vs 1844). Lift over base rate is
     the wrong yardstick when the incumbent already clears base by 2.65-9.85x. This is the single
     most important correction in this file.
   - **More seeds do not beat a more depleted pool.** Late seeding is *weaker* than early
     (mean lift 1.83x at round 7 vs 2.07x at round 3) even though seed counts rise 3-10x. The easy
     pathway neighbours are already taken by round 7.
   - **No deployable gate exists.** Spearman(AssayFormer's own observed rounds 5-7 hit rate,
     pathway's advantage in 8-10) = **-0.85**, so the gate direction is right — but every threshold
     nets zero or negative: `<2%` 2 screens -4 hits, `<3%` 4 screens -2 hits, `<5%` 6 screens
     -14 hits, `<10%` 10 screens -13 hits. The screens where AssayFormer collapses are screens where
     the pathway scorer is also near-random (2048: 3 hits vs 1 in 300 picks).

   The four screens where pathway wins are 1817 (+6 hits), 1644 (+6), 1818 (+3), 2048 (+2) — all
   tiny, and 1644 is the lone CRISPRa screen that inverted on diagnostic 3 as well.

   **Implication for the reversed handoff specifically.** The LLM at round 8 faces the same depleted
   pool *plus* a penalty the pathway scorer does not pay: its favourite genes are the famous ones
   AssayFormer already took, so they come back "already tested" and the slot is spent deeper in its
   ranking. The `af_horizon` dose-response measures exactly this — excluding AssayFormer's top 800-1000
   dropped Gemini from 17.7% to 9.7% — so the prediction for a round-8 handoff is ~11-12% against
   AssayFormer's 30.7%. The strong-hit lift does survive (2.52x on the 12 non-family screens), but EF
   counts a weak hit the same as a strong one, so that is not a usable win.

   Scripts are now in the repo at `scripts/diagnostics/` (they were nearly lost to a scratchpad
   wipe once already): `late_seed_diag.py` (scorer + per-screen table, `--seed-rounds` picks the
   horizon), `gate.py` (the gate test), `report.py` (aggregates a `--out` json), `pooled.py`
   (the pooled table above), `ef.py` (adjusted-EF check of an AssayFormer-alone sweep).

**Diagnostics run on 2026-10-21, part 2 — test-time interventions on AssayFormer itself
(§5 documents the code).** This is the first work here that modified the model rather than the warm
start, prompted by the observation that `relevance_score` is revealed per acquisition
(`tasks/gene_batch.py:116`) and thrown away by `predict`. Protocol: nested baseline (defaults are
bit-identical to the released model and score the recorded 6.100 / 4.828), paired delta, 95%
bootstrap CI over screen **families**, dose-response, and replication across both screen sets.

5. **Context ensembling does nothing. Context denoising does not survive the dose-response.**

   | Validation (11 families) | dEF | 95% CI | | Test (20 families) | dEF | 95% CI |
   |---|---|---|---|---|---|---|
   | ensemble K=8, sub 0.8, mean | -0.058 | [-0.164, +0.016] | | drop weakest 10% of hits | +0.066 | [-0.005, +0.169] |
   | ensemble K=8, sub 0.5, mean | -0.006 | [-0.148, +0.122] | | drop weakest 25% | **+0.122** | **[+0.022, +0.243]** |
   | ensemble K=8, sub 0.8, rrf | -0.013 | [-0.047, +0.019] | | drop weakest 35% | -0.025 | [-0.110, +0.051] |
   | ensemble K=16, sub 0.5, mean | -0.008 | [-0.084, +0.075] | | drop weakest 50% | +0.021 | [-0.079, +0.145] |
   | drop weakest 10% of hits | **-0.046** | **[-0.107, -0.002]** | | drop weakest 75% | -0.096 | [-0.389, +0.162] |
   | drop weakest 25% | +0.030 | [-0.120, +0.197] | | | | |
   | drop weakest 35% | +0.041 | [-0.107, +0.214] | | | | |
   | drop weakest 50% | +0.067 | [-0.077, +0.235] | | | | |
   | drop weakest 75% | -0.050 | [-0.313, +0.167] | | | | |
   | weak hits -> negative, 25% | +0.030 | [-0.110, +0.184] | | | | |

   **Read it as a null.** Across the 10 denoising cells, 6 are positive and 4 negative; the two whose
   CI excludes zero point in **opposite directions** (validation q=0.10 is significantly *negative*,
   test q=0.25 significantly positive); and there is no dose-response — validation's best dose (0.50)
   is test's near-zero, and test's best (0.25) is validation's near-zero. The runs are deterministic,
   so this spread is screen heterogeneity, not seed noise, and the SE on dEF is ~0.08-0.10. Do not
   report the +0.122 as a result: it is one cell of ten, chosen after looking.

   **Where it helps does not replicate either.** On test the gain concentrates on strong screens
   (Spearman(baseline EF, delta) = +0.67; weaker half -0.016 vs stronger half +0.261); on validation
   that correlation is -0.30. Poor-screen delta is ~0 on both (+0.007 test, +0.046 validation), so
   **this lever does not support per-screen routing** — the hoped-for use case. Direction gating adds
   nothing on test either (enrichment +0.116 vs overall +0.122).

   **The useful finding is the insensitivity itself.** Dropping three-quarters of the revealed hits
   from the context moves EF by ~0.1, and bagging half the observed set moves it less than 0.06.
   AssayFormer's û is saturated with respect to large perturbations of its observed-gene set. That is
   an independent explanation of the project's central negative result: a front end that only changes
   *which genes are revealed early* cannot move the outcome much, whatever it nominates. It also
   kills the whole family of "feed AssayFormer better context" ideas, not just these two.

6. **How much is AssayFormer's context conditioning worth at all? About as much as the paper claims
   for the warm start.** `--model-param ignore_context=true` (an ablation hook that already existed)
   makes the model score from the screen description alone, so the ten rounds of revealed labels
   contribute nothing and the policy is a static ranking walked 1000 deep:

   | Set | AssayFormer | context removed | value of ALL context | paper's claimed warm-start gain |
   |---|---|---|---|---|
   | `paper_validation` | 6.100 | 5.542 | **-0.478** [-0.852, -0.149] | — |
   | `paper_test` | 4.828 | 3.930 | **-0.897** [-1.720, -0.229] | **+0.84** |

   On the paper's own test set, in the paper's own metric, **everything AssayFormer gains from ten
   rounds of label feedback is worth ~0.9 EF — and the paper asks an LLM warm start, acting in three
   rounds with no labels available to it, to add +0.84 on top of that.** The claimed gain is
   comparable to the entire adaptive contribution of the model it is bolted onto. The CI is wide, so
   the claim is "same order", not "larger"; but it is an independent reason to doubt +0.84 that does
   not rest on any of the nine interventions, and it explains why none of them moved anything.

   It also bounds every remaining idea in this direction: no scheme that changes *which genes are
   revealed* can be worth more than ~0.9 EF on test, because that is the whole value of the labels.


7. **THE CEILING RESULT. AssayFormer is within 0.41-0.84 EF of the best possible member of its own
   hypothesis class, and that bound needs the answer key.** The score is `a(g) = û·v_g + b_g` with
   `d_gene = 10`, so for any screen the reachable rankings are indexed by a 10-vector û. Fitting û to
   a screen's **complete** label vector by logistic MLE (`scripts/diagnostics/oracle_u.py`) gives the
   best member of that family — an oracle, not a policy — and therefore bounds *every* intervention
   that works by handing AssayFormer a better û: warm starts, rerankers, exploration schemes,
   context surgery, all of it.

   | | `paper_validation` | `paper_test` |
   |---|---|---|
   | `b_g` only (ignores the screen entirely) | 3.685 | 3.017 |
   | + description, no AL context (`ignore_context`) | 5.542 | 3.930 |
   | **AssayFormer as released** | **6.100** | **4.828** |
   | û fitted to ALL labels (oracle ceiling) | 6.510 | 5.668 |
   | headroom inside the family | **+0.410** | **+0.840** |
   | û fitted to only rounds 1-3 labels (a real policy) | 3.948 | 3.457 |

   Four things follow, and they close out the project's question.

   - **The family is nearly exhausted.** A cheating oracle with every label gains +0.41 on validation.
     Ceilings are soft in the regulariser but not much: on test, `l2` 0.1 / 1.0 / 10.0 give
     5.572 / 5.668 / 5.468.
   - **The paper's claimed AssayLoop number, 5.67 on `paper_test`, equals the oracle ceiling (5.668).**
     Stated as an observation, not an accusation — their universe, metric variant and
     Gemini-3.1-Pro may all differ. But a three-round, label-free warm start reaching what an oracle
     with the full answer key reaches is implausible on its face, and it is a far more economical
     explanation of eleven null interventions than eleven separate failures.
   - **The amortized encoder is not underusing the labels — it beats direct fitting.** û fitted to 300
     revealed labels scores 3.948 / 3.457, *far worse* than the released model's 6.100 / 4.828. The
     transformer's learned prior is worth more than a well-posed 10-parameter MLE at that sample size.
     This refutes the hypothesis that motivated diagnostic 5.
   - **`b_g` is decorative.** sd 0.037 against the bilinear term's 3-7, and deleting it changes nothing
     (+0.013 validation, +0.039 test). Any plan to recalibrate the global gene prior is dead; the
     bias-only row is low *because* it ignores the screen, not because `b_g` matters.

   **Where headroom does exist:** the description channel. It is worth ~0.9 EF on test
   (3.017 -> 3.930), it is the one input you can change without retraining, and
   `amortized/text_embed.py:81` renders only 9 ScreenRecord fields — `direction_str` is **not** among
   them, although selection direction is this project's single most predictive screen variable (§3).
   Caveats before anyone spends on it: the shipped embedding cache is keyed on `sha1(model, text)`, so
   changing the text **misses the cache and requires fresh `text-embedding-3-small` calls** (Azure
   credentials, small but nonzero cost — hand the command over, do not run it), and the checkpoint was
   trained on the current renderer, so new fields are an off-distribution risk.

   Beyond that, raising the ceiling itself means better gene factors `v_g` — a different or larger
   representation, i.e. retraining, i.e. out of this project's compute budget.


8. **Hit strength carries no ranking information beyond the binary label — measured at the oracle, so
   no plumbing can rescue it.** The task reveals `relevance_score` per acquisition, and the ordering
   is real: across screens 91% (validation) / 82% (test) of hits have distinct relevance values, and
   several screens span a wide range (1644 0.35-1.00, 2104 0.50-1.00, 1817 0.57-1.00). It is *not*
   degenerate — 1561 (sd 0.0028, 14% tied) is the exception, not the rule.

   But feeding the oracle that ordering changes nothing. Fitting û to three targets with the
   estimator held fixed (`oracle_u.py --target {binary,revealed,rank} --estimator {logistic,ridge}`):

   | target | validation (logistic) | test (logistic) | validation (ridge) | test (ridge) |
   |---|---|---|---|---|
   | binary hit/no-hit | **6.510** | **5.668** | 5.640 | 3.990 |
   | revealed relevance (hits; non-hits 0) | 6.503 | 5.665 | 5.641 | 3.989 |
   | rank among hits, scale-free | 6.463 | 5.374 | 5.524 | 3.750 |

   Strength is worth **-0.007 EF on validation and -0.003 on test** — zero to three decimals — and the
   rank transform actively hurts. The reason is structural: with only 10 parameters, spending them on
   "how strong a hit is" instead of "whether it is a hit" is a bad trade, and adjusted EF counts a
   weak hit exactly like a strong one.

   **Hold the estimator fixed when reading this.** The first pass compared logistic-on-binary (6.510)
   against ridge-on-relevance (5.641) and looked like strength hurting by 0.87. That gap is entirely
   the loss function: ridge-on-binary is also 5.640. Target and estimator must be varied separately.

   **This explains diagnostic 5 mechanistically.** `ctx_weak_hits` was null not because the dropping
   was clumsy but because there was no strength signal to exploit in the first place. It also tempers
   diagnostic 3: the relevance ordering distinguishing "strong" from "weak" hits is real, but it is
   not information this model can rank on.

   **For the record, how strength *could* be fed in, since the mechanisms are not obvious and someone
   will ask again.** `_encode_hidden` does `self.hit_emb(ctx_hit.clamp(0, 1))` over an
   `nn.Embedding` with exactly 2 rows. You would not add rows (that is the retired graded-labels
   idea, which asks the model to consume embeddings it never learned); you would pass a convex
   combination `(1-w)·hit_emb[0] + w·hit_emb[1]` with `w` the normalised strength, which stays inside
   the trained span. Cruder alternative: duplicate strong hits in the context to upweight them. Both
   are a few lines. Neither is worth writing, because the oracle above bounds what any of them could
   deliver at ~0.00 EF.


9. **IMPORTANT SCOPE CORRECTION to diagnostic 7. The ceiling binds û, not the architecture.**
   Diagnostic 7 showed the best û is worth only +0.41/+0.84 over the released model, and §1 drew from
   that the conclusion that "the hypothesis class is nearly exhausted". **That is true of û-based
   interventions and false of the architecture.** Appending features,

       a(g) = û·v_g + b_g + Σ_k w_k·φ_k(g),

   lives in a (10+k)-dimensional subspace *strictly containing* the original, so `w = 0` recovers the
   released model exactly and the oracle ceiling does not apply. `scripts/diagnostics/augment_u.py`
   fits (û, w) jointly to all labels and reports the augmented ceiling:

   | | `paper_validation` | `paper_test` |
   |---|---|---|
   | R² of the hit vector on span(V, b, 1) | 0.214 (median 0.051) | **0.035 (median 0.008)** |
   | plain oracle ceiling | 6.502 | 5.666 |
   | + one *oracle* feature (unmasked relevance) | **16.089 (+9.59)** | **18.762 (+13.10)** |
   | + pathway propagation | 6.557 (+0.054) | 5.783 (+0.117) |
   | + pathway propagation × hit strength | 6.558 (+0.056) | 5.784 (+0.119) |
   | + `b_g` (sanity: already in the model) | 6.515 (+0.013) | 5.707 (+0.042) |
   | + Gaussian noise (negative control) | 6.499 (-0.003) | 5.619 (-0.047) |

   - **Only 3.5% of "which genes are hits" is expressible in AssayFormer's basis on test** (21% on
     validation, and the medians are far lower still — the mean is carried by the fitness family).
     The remaining structure is orthogonal to span(V, b, 1) and simply unreachable for any û.
   - **A single extra additive dimension is worth up to +13 EF** when it is informative. So capacity
     is emphatically *not* the binding constraint, and "nothing can improve AssayFormer" is the wrong
     conclusion to draw from diagnostic 7.
   - **The binding constraint is feature quality.** The one biologically motivated feature available
     locally, pathway propagation, is worth **+0.05 to +0.12** — consistent with its failure as a
     competing ranker in diagnostic 4, and now explained: it is nearly inside the existing span.
   - **Hit strength adds +0.002 on top of unweighted pathway propagation**, closing out §4 diagnostic 8
     from the additive-feature side as well. Strength is not useful as a target, as a context signal,
     or as a feature weight.

   **What the existing 11 numbers per gene are.** `V` is 20,893 x 10 and its variance is spread
   almost evenly over nine directions (0.152 ... 0.077; the tenth is degenerate at 0.001, so the
   usable rank is 9). The directions are not simply "how often is this gene a hit": the strongest
   Spearman between any singular direction and cross-screen hit frequency is -0.289, and for `||v_g||`
   it is +0.068. No gene has an all-zero factor. The decisive structural fact for choosing features is
   that **`v_g` is screen-independent** — only û varies by screen — so a gene-only feature competes
   with `v_g` for the same job and is likely partly redundant, whereas a gene x screen feature
   (expression in *this* cell line, an LLM score for *this* phenotype) is something the architecture
   can represent only through the 10-number bottleneck. Necessary but not sufficient: pathway
   propagation is gene x screen and still scored +0.05, so screen the feature rather than trusting
   the argument.

   **`augment_u.py` is the instrument to use before building anything else.** Propose a gene feature,
   measure its augmented ceiling, and only implement if it clears, say, +0.5 EF. It costs seconds per
   feature and it would have killed the pathway scorer before I wrote 500 lines for it. Candidates not
   yet screened, all available locally or free: DepMap/Achilles essentiality (the user's own
   `../achilles.py` already downloads it), cell-line expression, PPI degree and shortest-path distance
   to seeds (`../candidate_retrieval.py`), GO-BP membership as a block of indicators, paralog family,
   and an LLM-assigned per-gene score. The negative control and the oracle row are both in the script,
   so a new feature is always read against a floor and a ceiling.


10. **DepMap essentiality and expression are both null as additive features.** Screened with
    `scripts/diagnostics/screen_depmap.py` (data: DepMap 24Q4 Public, figshare article 27993248,
    free and public — `achilles.csv` 297 MB via the user's own `../achilles.py` file id, plus
    `Model.csv` and `OmicsExpressionProteinCodingGenesTPMLogp1.csv` via
    `https://api.figshare.com/v2/file/download/<file_id>`; the plain `figshare.com/ndownloader`
    route returns HTTP 202 forever and the DepMap portal API is bot-gated).

    | added feature | validation | test | screens (val/test) |
    |---|---|---|---|
    | essentiality, pan-cancer mean | -0.016 | +0.046 | 20 / 20 |
    | essentiality, fraction of lines essential | +0.005 | +0.040 | 20 / 20 |
    | essentiality in **this screen's cell line** | -0.005 | +0.058 | 20 / 20 |
    | expression in **this screen's cell line** | +0.057 | -0.126 | 18 / 12 |
    | "not expressed" flag (log2 TPM+1 < 1) | +0.020 | -0.002 | 18 / 12 |
    | *oracle relevance (positive control)* | *+9.586* | *+13.097* | 20 / 20 |
    | *Gaussian noise (negative control)* | *-0.003* | *-0.047* | 20 / 20 |

    Every candidate sits inside the noise band and **flips sign between the two screen sets**. The
    positive control works, so these are true nulls rather than a broken instrument. Nothing comes
    near the +0.5 threshold set in diagnostic 9.

    - Cell-line matching: 17/20 validation and 10/20 test screens resolve to a DepMap model via
      `Model.csv`; of those, 10 and 7 appear in this gene-effect release and 18 and 12 in the
      expression release. The misses are real — HEK293T and primary T-cells have no DepMap model.
      `ALIAS` in the script maps engineered derivatives to their parent line (Huh-7.5.1 -> HuH-7,
      HeLa RelA-mNeonGreen -> HeLa, Calu1-ACE2 -> Calu-1).
    - **Bookkeeping trap worth not repeating.** The first pass averaged a feature over the 18 screens
      it covered and compared that to a baseline averaged over all 20, which read as expression being
      worth **-0.131** on validation. Compared on the same 18 screens it is **+0.057** — the sign
      flipped. Any partially-covered feature must be scored against the baseline restricted to its
      own screens; the script now does this and labels the row.
    - Why they fail, most likely: screens only ever test genes in their own library, so non-expressed
      and pan-essential genes are already largely absent from the places AssayFormer spends its
      budget. The information is not new relative to span(V).

    Running tally of additive features screened: oracle relevance +9.6/+13.1, pathway propagation
    +0.05/+0.12, pathway x strength +0.06/+0.12, `b_g` +0.01/+0.04, essentiality ~0, expression ~0.
    **Every real, locally available feature tried so far is nearly inside span(V).** The one candidate
    class not yet screened is a per-gene score conditioned on the screen's *phenotype text* — i.e. an
    LLM score as a feature rather than as a picker, the one role it has never been given. That costs
    API calls, so it must be handed over, not run.

    Reusable tooling: `screen_depmap.py` (prefers a `.npz` cache), `depmap_to_npz.py` (build it once;
    the 18k-column CSV parse was the entire runtime), `run_depmap_screens.sh` (both sets).


11. **The LLM as an additive feature: real signal, zero incremental value. This closes the project's
    central question with a mechanism.** Flash-Lite scored every gene 0-10 for the phenotype over
    AssayFormer's top 3000, in two pilots (§5 for the tooling). Data quality was clean — 3000/3000
    genes returned a score on nearly every screen, 30 batches of ~100, no parse failures, and a
    sensible right-skewed distribution (~60% zeros, a tail out to 10).

    | pilot | screens | AUROC of the raw score | oracle dEF | deployable dEF | noise floor |
    |---|---|---|---|---|---|
    | mid-tier enrichment | 1561, 1581, 1584 | 0.584 | **-0.012** | -0.584 | +0.031 |
    | hard enrichment | 1817, 1818, 1815, 1644, 2104, 2048 | 0.553 | **-0.077** | -0.408 | -0.122 |

    Per screen on the hard set: 1818 +0.116 (the only positive), 2104 -0.008, 1644 -0.017,
    1815 -0.030, 1817 -0.030, 2048 -0.494.

    **The two numbers disagree, and that disagreement is the finding.** AUROC 0.553-0.584 says the
    LLM's scores genuinely separate hits from non-hits with no fitting at all — Flash-Lite is not
    guessing. The oracle delta says that information adds nothing once `û·v_g` is already in the
    model, even with the weight fitted to the full answer key. AUROC measures standalone
    discriminative power; the augmented oracle measures *incremental* power. So **the LLM's
    gene-level knowledge is already inside the BPMF gene factors.**

    This is the strongest form of the redundancy result. The original measurement (62-67% of
    warm-start hits are genes AssayFormer finds anyway) could be dismissed as an artefact of the
    picking mechanism — both nominate famous genes, so of course they overlap. This demonstrates it at
    the *representation* level, which means no re-plumbing of the LLM's output can help: not as a
    picker, not as a warm start, not as a mixed batch, not as a feature. One mechanism explains every
    null in this project.

    **Do not spend more on model capability here.** AUROC > 0.5 with complete coverage shows the model
    engaging sensibly, so the failure mode is redundancy, not effort or capability — the same
    conclusion the GLM-5 ≈ Flash-Lite comparison reached, now with a cause. A stronger model or
    `reasoning_effort: high` cannot fix redundancy.


12. **SCOPE CORRECTION, and the most important number in this file. There is enormous room for
    improvement in the TASK; there is almost none reachable by improving û.** I had been summarising
    diagnostic 7 as "little room for improvement", which is wrong. Two different ceilings were being
    conflated. `scripts/diagnostics/task_ceiling.py` computes the second one:

    - **oracle-û ceiling** — the best ranking expressible as `a(g) = û·v_g + b_g` for *any* û.
      Bounds every intervention that works through the screen latent. **+0.41 / +0.84.**
    - **perfect-play ceiling** — rank every true hit first. This is what the benchmark allows, and it
      says nothing about the architecture. With budget 1000 over a D-gene library holding H hits,
      `EF_max = min(1000, H) / (1000H/D)`, which is `D/1000` (about 19) whenever `H <= 1000`.

    | | `paper_validation` | `paper_test` |
    |---|---|---|
    | AssayFormer | 6.10 | 4.83 |
    | oracle û (best possible latent) | 6.51 | 5.67 |
    | **perfect play** | **16.77** | **18.92** |
    | fraction of the task ceiling reached | **42.7%** | **25.8%** |
    | hard screens (EF < 5) | 9.5% of ceiling, 9.0% of hits | 12.4% of ceiling, 11.2% of hits |
    | easy screens (EF >= 5) | 69.9% of ceiling, 49.3% of hits | 45.8% of ceiling, 41.1% of hits |

    On the nine hard validation screens AssayFormer finds **4.7% to 9.9% of the hits** (1817: 42 of
    888; 1818: 20 of 371; 2048: 11 of 128). The task headroom there is roughly **+17 EF**, forty times
    the oracle-û headroom. So the correct statement is not "there is little room" but:

    > **Everything we tried was confined to a 10-dimensional family whose ceiling is +0.41, while the
    > task has +10.7 of headroom on average and +17.7 on the hard screens. The ~10 EF gap between the
    > two is the part of the task this gene representation cannot express.**

    Three measurements agree on that gap, which is why it is believable:
    - only **3.5%** (test) / 21% (validation) of the hit vector lies in `span(V, b, 1)` (§4.9);
    - oracle û reaches 6.51 on validation, but oracle û **plus one informative feature** reaches
      **16.09**, and perfect play is **16.77** — i.e. a single genuinely new feature would very nearly
      close the whole gap (§4.9);
    - every real feature screened so far (pathway, DepMap essentiality, DepMap expression, LLM score)
      landed inside span(V) at ~0, so the gap is real but we have not found a key to it.

    **The budget binds on the easy screens, and only there.** The U_179x fitness screens hold
    1287-1632 hits, more than the 1000-pick budget, so their perfect-play EF is only 11.7-14.8 and
    AssayFormer already reaches 73-87% of it. Those screens genuinely are near-saturated. **All of the
    available headroom is on the hard, low-hit-density screens** — which is also where the user's
    interest has been throughout, and it means a targeted method does not have to beat AssayFormer
    everywhere to be worth having.

    **CORRECTION 2026-10-22, measured on the 198-screen fold instead of the curated 20.** Two claims
    above do not generalise, and the second changes where to look for room.

    - **Half the fold is metric-DEGENERATE.** 98 of 198 screens have a library SMALLER than the
      1000-pick budget (many are 681-gene focused libraries), so `min(1000, H)/(1000H/D)` puts
      perfect play below 1 and the metric is meaningless there. The curated sets never exposed this
      because their manifests filter to 18-22k-gene libraries (§7). **Exclude them from any ceiling
      or fraction-of-ceiling statement.**
    - **The headroom is BROADLY DISTRIBUTED, not concentrated on the hard screens.** On the 100
      non-degenerate screens:

      | group | screens | hits held | mean EF | perfect play | mean gap | share of headroom |
      |---|---|---|---|---|---|---|
      | zero hits found | 4 | 9 | 0.00 | 18.93 | 18.93 | 7.5% |
      | very hard (EF < 2) | 32 | 22573 | 1.23 | 12.94 | 11.71 | **37.0%** |
      | middling | 24 | 4664 | 3.35 | 14.87 | 11.52 | 27.3% |
      | does well (EF >= 5) | 40 | 28218 | 8.68 | 15.85 | 7.18 | **28.3%** |
      | **ALL** | **100** | **55464** | **4.67** | **14.81** | **10.14** | 100% |

      Only **23 of 100** screens hold more hits than the budget, and even they reach just **54%** of
      their ceiling; the other 77 reach **40%**. So "the easy screens are saturated" described
      `paper_validation`'s eight-screen PC-9 fitness family, not the population. **Even where
      AssayFormer does well it leaves ~7 EF unclaimed**, and that group carries 28% of the total
      headroom — so a method does not have to specialise in hard screens to find room. The headroom
      is about **10 EF per screen, roughly 3x the current score**, essentially everywhere.

    **What this does and does not change.** It does not rescue any of the twelve null interventions:
    they were either û-based (hard-bounded at +0.41) or used features already inside span(V), and
    those measurements stand. It does change the conclusion to draw. "AssayFormer cannot be improved"
    is wrong; "AssayFormer cannot be improved through its screen latent, and the one lever with real
    headroom is the gene representation — new features, or retraining `v_g`" is right, and is a
    sharper finding than a flat negative.


13. **A POSITIVE RESULT: AssayFormer fused with DepMap co-essentiality gains +0.540 EF on the
    paper's test set, CI excluding zero, with no LLM and no retraining.** This is the first
    intervention in the project to clear zero with a confidence interval, and it came from correcting
    my own mistake in diagnostic 10 — essentiality as a *scalar mean* is null, essentiality as a
    *profile across cell lines* is not.

    The feature. Each round, average the row-centred unit-normalised DepMap knockout profiles of
    every gene revealed as a hit so far into a query vector `q`, and score each untested gene by
    `corr(profile(g), q)` — genes that co-vary across 808 cell lines with what already worked in
    *this* screen. Fused with AssayFormer by reciprocal rank fusion. It is pathway propagation with a
    data-derived graph instead of a curated one, and unlike `v_g` it is conditioned on the screen.

    | config | `paper_validation` dEF | `paper_test` dEF |
    |---|---|---|
    | weight 0 (nested baseline) | 6.100 exactly, pick for pick | 4.828 exactly |
    | **weight 0.5, RRF** | +0.202 [-0.055, +0.475], 8/20 | **+0.540 [+0.183, +0.941], 13/20** |
    | weight 1.0, RRF | -0.009 [-0.373, +0.365] | +0.406 [-0.055, +0.930] |
    | weight 1.0, zsum | -0.078 [-0.588, +0.426] | +0.446 [-0.175, +1.072] |
    | weight 2.0, RRF | **-0.715 [-1.335, -0.187]** | -0.029 [-0.685, +0.632] |

    **What the weight actually controls, since "0.5 vs 2.0" is not a matter of degree.** Under
    `fused = 1/(60 + r_AF) + weight/(60 + r_codep)`, the weight sets how deep in AssayFormer's
    ranking a gene may sit and still be promoted. A gene AssayFormer ranks 5000th needs a codep rank
    of ~#20 at weight 0.5, ~#97 at 1.0 and ~#238 at 2.0. More importantly there is a threshold
    between 1.0 and 2.0: below it the codep scorer **cannot** seat its own favourite at the head of a
    batch (its best possible contribution, 1/60 x weight, stays under AssayFormer's top pick), and at
    2.0 it can (0.0334 against 0.0173). So 0.5 is a surgical re-ranker and 2.0 is the primary ranker
    with AssayFormer demoted to a tiebreaker. Since the codep direction is estimated from ~20-30
    hits, its top ~20 are worth trusting and its top ~240 are mostly noise, and every promotion
    displaces an AssayFormer pick.

    That also explains the asymmetry in the harm: weight 2.0 costs **-0.715 on validation** but only
    **-0.029 on test**, and AssayFormer's early-round hit rate is 39.6% on validation against 8.3% on
    test. **The damage from over-weighting scales with the quality of the picks being displaced** --
    which is the §2(a) mechanism behind the LLM warm-start failure, recurring in a method with no LLM
    in it. One diagnosis covers both, and it is the reason to prefer a narrow promotion window over a
    wide one whenever the incumbent policy is strong.

    **Why I believe it, in order of weight:**
    - **A coherent dose-response with the same optimum on both sets.** 0.5 best, 1.0 smaller, 2.0
      gone (significantly negative on validation). This is exactly what the context-denoising family
      lacked, and it is the single strongest discriminator between signal and noise here.
    - **The baseline is nested and verified.** `weight=0` reproduces the AssayFormer-alone campaign
      pick for pick — all 1000 picks, 47 hits on screen 1561 — so the delta cannot be a harness
      artefact.
    - **Replication in sign on the held-out set.** Validation's CI includes zero, but it has 11
      families against test's 20, so less power; the point estimate and the optimum agree.
    - **On test: +120 hits** (1611 -> 1731), 15.3% of picks changed, feature active in 8.2/10 rounds.

    **Where it helps, and the mechanism.** It helps where the 1000-pick budget leaves headroom and
    hurts where it does not:

    | | screens | hits before -> after | better on |
    |---|---|---|---|
    | validation, budget-saturated (H > 1000) | 8 | 6361 -> 6203 (**-158**) | **0/8** |
    | validation, has headroom (H <= 1000) | 12 | 800 -> 865 (+65) | 8/12 |
    | test, has headroom | 19 | 1150 -> 1238 (+88) | 13/19 |

    The 8 losers are the near-duplicate PC-9 fitness family, where AssayFormer already reaches 73-87%
    of a ceiling the budget caps (§4.12) and any perturbation can only cost hits. Validation's raw hit
    count therefore FALLS by 93 overall even as family-weighted EF rises — report both numbers.

    **CORRECTION 2026-10-22: there is nothing to gate, and "do not deploy ungated" below is wrong.**
    Re-measured for the DEPLOYED configuration (`diffmeans`, weight 0.5, knn 0.5) on the 198 fresh
    screens of the `validation` fold:

    | group | n | families | dEF/screen | dEF/family | hits |
    |---|---|---|---|---|---|
    | budget-saturated (total hits > 1000) | 23 | 8 | **+0.012** | **+0.023** | −5 |
    | has headroom (total hits <= 1000) | 175 | 102 | +0.197 | +0.355 | +291 |
    | what the gate actually catches (rounds 1-3 rate > 40%) | 18 | 2 | +0.004 | −0.003 | −15 |

    Dropping all 23 saturated screens moves the family-weighted result from +0.331 to +0.355, so they
    cost **0.024 EF**, and applying the gate changes the reported dEF by **0.000 at every threshold
    from 20% to 70%** (`codep_gated.py --screen-set validation --off sweep-dm-validation-w0 --on
    sweep-dm-validation-dm05knn --exclude-set paper_validation`).

    **Why the original result did not hold.** It was measured on `paper_validation`, where 8 of 20
    screens are ONE near-duplicate PC-9 fitness family — 40% of the set and most of its raw hit count.
    On the fold those screens are 23 of 198 but only 8 of 110 families, so family weighting already
    neutralises them and their per-screen delta is mildly POSITIVE. **A per-screen-group effect
    measured on a 20-screen set with an 8-screen family in it is a statement about that family.**

    The gate is still sound as a classifier and still leak-free — early hit rate separates the groups
    at 49.4% against 2.4%, using only rounds 1-3 and never `|H|` — there is simply nothing worth
    classifying. **The remaining harm is diffuse:** 75 screens improve and 75 worsen, and the ones
    that worsen are not the saturated ones, so they are not identifiable from saturation or early hit
    rate. Three gating attempts have now failed this way (this one, direction gating in §4.5, early-rate
    gating for the pathway scorer in §4.4); prefer reducing the feature's variance over switching it
    off selectively.

    The original, superseded reasoning follows.

    **A deployable gate exists, unlike every previous attempt.** Total hits is not observable mid-campaign
    but the early hit *rate* is, and it separates the groups cleanly: AssayFormer hits at 71-95% in
    rounds 1-3 on the saturated screens and 2-23% on the rest. `codep_gated.py` scores the policy
    "disable the feature if rounds 1-3 hit above 40%", deciding from round-1-3 information only:
    validation **+0.229 [-0.014, +0.488]** (12 of 20 screens gated on), test **+0.540 [+0.183, +0.941]**
    (all 20 gated on, so unchanged). The gate removes most of validation's damage without touching test.

    **Scale.** The paper claims +0.84 on test for its LLM warm start; this gets **+0.540 with no LLM
    at all**, from a public dataset, on a frozen checkpoint. It also sits inside the room diagnostic 12
    identified: the oracle-û ceiling is +0.84, so a feature beating what *any* screen latent can do is
    only possible because it adds a dimension outside `span(V)`.

    **What is not established.** Validation does not reach significance on its own. Eight configs
    across two sets is a modest multiplicity. Weight 0.5 was picked from this sweep, so a fresh
    confirmation at a pre-specified weight 0.5 on new screens would be the clean test. And the
    saturated-screen damage is real: do not deploy ungated. *(Superseded — see the correction at the
    top of this diagnostic. On held-out data the damage is +0.012 per screen, i.e. absent, and the
    gate is a no-op. Deploy ungated.)*

    Code: `src/assayloop/acquisitions/codep_seed_acq.py` (registered as `--acq codep_seed`, 11 tests in
    `tests/test_codep_seed.py`), `scripts/diagnostics/run_codep_campaigns.sh`, `codep_gated.py`.
    Needs `--acq-param profile_npz=<achilles.npz>` from `depmap_to_npz.py`. A new `--acq-param` flag
    was added to `assayloop run`; acquisition overrides were previously unreachable from the CLI.

    **One design trap the unit tests caught.** Under RRF, *possessing* a co-dependency rank is itself
    worth up to `weight/60`, so leaving the 3.4k genes missing from DepMap unranked would have
    systematically outranked them — the coverage confound re-entering through the fusion rule, right
    after §4.10 measured it. Uncovered genes now take the **mean** covered contribution, so coverage
    conveys no advantage either way.


14. **Held-out confirmation, plus a second feature. Best result AT THE TIME: +0.145 EF per screen
    [+0.053, +0.253] on 198 screens that chose nothing** — superseded by §4.15's +0.175 once the
    query direction changed from `centroid` to `diffmeans`; this section is the `centroid`-era record.** Two ideas were screened offline first; one
    was refuted before any campaign was built, the other works and is additive.

    **The confirmation (§4.13's result, on held-out data).** `paper_test`'s +0.540 was measured on the
    20 screens that chose weight 0.5. On the 198 fresh screens of the `validation` fold:

    | arm | per-screen dEF | per-family dEF | up / down / same |
    |---|---|---|---|
    | codep only (the §4.13 method) | +0.106 [+0.015, +0.203] | +0.245 [+0.118, +0.383] | 69 / 80 / 49 |
    | screen-kNN only | +0.033 [-0.021, +0.091] | +0.081 [+0.005, +0.165] | 58 / 84 / 56 |
    | **both, 0.5 / 0.5** | **+0.145 [+0.053, +0.253]** | **+0.278 [+0.141, +0.445]** | 75 / 77 / 46 |
    | both, 0.5 / 0.25 | +0.121 [+0.029, +0.221] | +0.261 [+0.129, +0.404] | 70 / 81 / 47 |

    It replicates, at **about a fifth of the curated-set figure** — the usual shrinkage when you stop
    measuring on the data that selected the setting. Report +0.106 per screen / +0.245 per family for
    the single-feature method and +0.145 / +0.278 for the pair. Note the per-family number is roughly
    double the per-screen one because the three largest families (17, 16 and 14 screens) all have
    negative means; family weighting is correct given near-duplicates but it is doing a lot of work,
    so quote both.

    **Idea A — better seed aggregation — REFUTED, and do not retry.** The deployed feature scores a
    gene by correlation with the MEAN of the confirmed hits' profiles. I predicted that centroid
    destroys multi-module structure and that max-over-seeds or mean-of-top-3 would beat it. Both are
    *worse*, on both screen sets:

    | aggregation | validation | test |
    |---|---|---|
    | centroid (deployed) | +0.148 [-0.087, +0.435] | +0.435 [+0.093, +0.772] |
    | max over seeds | +0.136 | +0.082 |
    | mean of top 3 | +0.118 | +0.251 |

    With ~30 noisy profiles the mean is a far better estimate of the shared signal than any single
    one, and max-over-seeds latches onto whichever seed is noisiest in a gene's direction. Averaging
    denoises more than it blurs.

    **Idea C — screen-kNN — WORKS and is additive.** `v_g` is a 10-number compression of 1349
    training screens. This scores a gene by how often it is a hit among the training screens whose hit
    sets most resemble what the current screen has revealed, read off the raw 1349 x 20874 matrix
    (`build_train_hits.py`, 1.5M hit entries, 5.32% dense). Using training labels at inference is not
    leakage: they are training data and the folds are disjoint. On held-out data
    0.106 + 0.033 = 0.139 predicted against 0.145 measured — **the offline additivity survived
    deployment**, which is the first time in this project an offline margin transferred. The two
    features help in different places: screen-kNN contributes +0.032 on the poor screens against
    codep's +0.43, so they are complements, not substitutes.

    Caveat worth keeping: on `paper_test` the pair (+0.515 to +0.566) was indistinguishable from codep
    alone (+0.540), i.e. the additivity did NOT show there. With n=20 that is the less reliable
    measurement, and the 198-screen fold is what the conclusion rests on.

    **A methodology fix that changes how to read §4.9 and §4.13.** The width-matched Gaussian control
    was a SINGLE draw per feature. Changing only the order the controls were drawn moved margins by
    ~0.1 EF — a single draw is itself a noisy estimate of the capacity baseline. `screen_seed_aggs.py`
    now averages **5 draws**, with a generator seeded per (screen, feature) so results cannot depend
    on how many features ran earlier. **The 50-column numbers in §4.9 (+0.60 / +0.84) used single
    draws and carry more uncertainty than stated there.**

    **Where the remaining headroom is.** The 50-column profile block scored +0.596 offline against the
    1-column centroid's +0.148, so most of the profile's information still never reaches the campaign.
    Getting several profile dimensions into the acquisition — not a better scalar — is the larger
    remaining prize, and the harder one (online weight fitting on ~300 labels was measured noisy
    enough that a pure confound control scored +0.58).

    Code: `codep_seed_acq.py` now carries two scorers fused as
    `1/(k+r_AF) + w_codep*phi_codep + w_knn*phi_knn`, with `--acq-param train_hits=<prefix>
    knn_weight=0.5 knn_k=50`. Both scorers apply the same neutrality rule (a gene missing from either
    source takes the mean contribution of those that have one), refactored into shared `_z` and
    `_rrf_term` helpers. 17 tests in `tests/test_codep_seed.py` pin that `knn_weight=0` is bit-identical
    to codep-only, both weights 0 is exactly greedy, and either scorer fires independently.
    Runners: `run_codep_twofeat.sh`, `screen_seed_aggs.py`, `build_train_hits.py`.


15. **Direction estimation: the single biggest improvement to the method, and it is one line.**
    §4.14 left a gap — a 50-column profile block scored +0.596 offline against the deployed
    1-column centroid's +0.148 — and the obvious reading was "get more dimensions into the campaign
    with per-screen weights". That reading was wrong, and the right one is simpler.

    **The reframing.** A *direction* in profile space collapses to ONE column:
    `corr(profile(g), w)`. The 50-column oracle margin does not say the campaign needs capacity; it
    says **a better single direction exists**. That turns an unstable 50-parameter online fit (which
    §4.9 measured as hopeless — a pure confound control scored +0.58 EF under it) into a closed-form
    direction estimate, and it sidesteps the capacity/noise-control problem entirely because
    everything stays 1-column and directly comparable to the deployed feature.

    **What was screened** (`screen_seed_aggs.py`, 1-column features, 5-draw width-matched noise
    controls, both screen sets, two independent runs each):

    | query direction | validation | test |
    |---|---|---|
    | hit centroid (the old default) | +0.126 … +0.158 | +0.444 … +0.467 |
    | **hits mean MINUS non-hits mean** | **+0.342 … +0.365** | **+0.559 … +0.656** |
    | diagonal variance shrinkage of the above | +0.346 … +0.350 | +0.627 … +0.647 |
    | first principal direction of the seeds | -0.022 … +0.006 | **-0.167 … -0.110** |
    | best 2- and 3-direction blocks | +0.298 … +0.418 | +0.577 … +0.807 |

    - **Difference of means roughly DOUBLES the centroid's margin**, on both sets, with intervals
      excluding zero where the centroid's does not on validation. Two reasons it works: it uses all
      ~300 labels per three rounds instead of only the ~20 hits, and subtracting the non-hit mean
      cancels whatever is merely *typical of the region AssayFormer already picks from* — which the
      hit mean otherwise bakes in, promoting genes that are typical rather than distinguishing.
      It is the two-class centroid difference, i.e. the Fisher direction under isotropic covariance.
    - **Diagonal shrinkage adds nothing** over it. Expected in hindsight: same numerator.
    - **The seeds' principal direction is actively HARMFUL** on test. This is the THIRD failure of the
      multi-module hypothesis (after max-over-seeds and mean-of-top-3 in §4.14). One target direction
      is the right model; it was simply being computed badly.
    - **Multi-direction blocks buy capacity, not signal** once run-to-run variation is accounted for.
      That variation is still +-0.03 to +-0.10 *with* 5 noise draws, which is the right error bar for
      comparing any two rows in this table.

    **Deployed, on `paper_test`:**

    | arm | dEF | wins | poor-screen dEF |
    |---|---|---|---|
    | centroid + kNN (the §4.14 best) | +0.566 [+0.226, +0.943] | 15/20 | +0.410 |
    | diffmeans, no kNN | +0.497 [+0.153, +0.901] | 14/20 | +0.495 |
    | **diffmeans + kNN** | **+0.562 [+0.274, +0.878]** | **18/20** | +0.462 |

    The mean is unchanged but **consistency improves sharply — 18 of 20 screens, the best in the
    project — and the interval tightens** ([+0.274, +0.878] against [+0.226, +0.943]). Same expected
    gain, much less variance. The held-out fold is the number to quote; see §1.

   **DEPLOYED, HELD OUT — the final numbers.** All five configurations on the same 198 fresh screens
   of the `validation` fold, each against its own weight-0 baseline:

   | configuration | per-screen dEF | per-family dEF | up / down / same |
   |---|---|---|---|
   | centroid only (§4.13) | +0.106 [+0.015, +0.203] | +0.245 [+0.118, +0.383] | 69 / 80 / 49 |
   | centroid + kNN (§4.14) | +0.145 [+0.053, +0.253] | +0.278 [+0.141, +0.445] | 75 / 77 / 46 |
   | diffmeans only | +0.113 [+0.016, +0.219] | +0.258 [+0.124, +0.402] | 68 / 81 / 49 |
   | **diffmeans + kNN (BEST)** | **+0.175 [+0.075, +0.293]** | **+0.331 [+0.178, +0.514]** | 75 / 75 / 48 |
   | diffmeans at weight 1.0 | **-0.109 [-0.257, +0.043]** | +0.068 [-0.136, +0.279] | 62 / 90 / 46 |

   **The offline doubling did not transfer, and that is the headline caveat.** Offline, diffmeans beat
   the centroid 2.3x (+0.35 vs +0.15). Deployed and alone it beats it by 1.07x (+0.113 vs +0.106) —
   essentially nothing. The real gain comes from the COMBINATION: +0.175 against +0.145 for
   centroid + kNN, so diffmeans contributes about +0.03 on top of the second scorer rather than the
   +0.2 the offline screen implied.

   **Weight 1.0 is significantly harmful** (-0.109 per screen), against +0.175 at 0.5. The dose
   response is sharp, so do not raise the weight on the strength of an offline margin.

   **The pattern across three attempts: offline oracle margins overstate deployed gains, and the
   ORDERING of candidates does not transfer reliably either.** The 50-column blocks, the
   paper_test additivity, and now diffmeans all behaved this way. Offline screening remains the right
   filter — it is seconds per candidate and it killed max-over-seeds, seed-PC1 and diagonal shrinkage
   correctly — but treat a positive offline margin as a licence to run the campaign, never as an
   estimate of what the campaign will give.

   Net from the whole line of work: AssayFormer alone 6.100 on `paper_validation`, and
   **+0.175 EF per screen / +0.331 per family on 198 screens that chose nothing**, with 75 screens
   improved against 75 worsened — the gain is in the size of the wins, not their frequency.

    Code: `codep_seed_acq.py` takes `query=diffmeans|centroid`, defaulting to `diffmeans` from
    2026-10-22. `centroid` is kept because every result before that date used it and must stay
    reproducible. Two guards matter: `min_nonhits` (default 10) falls back to the centroid when too
    few non-hits are available, and it counts non-hits **present in the DepMap atlas**, not labelled
    non-hits — a failing test caught that those are different numbers. `max_labels` caps the means.
    34 tests in `tests/test_codep_seed.py`.


16. **Hierarchical shrinkage toward similar training screens — REFUTED offline, in both of its
    forms, on both screen sets, with controls. No campaign was run.** This was §10 item 1, the
    identified next lever. The reasoning was sound and the measurement is clean, so the record of
    *why* it fails is the useful part.

    **The idea.** The deployed query direction is a difference of means estimated from the ~20-30
    hits revealed so far, in 808-dimensional space, and §4.14/§4.15 established that the binding
    constraint is sample size rather than the estimator (three attempts to extract more structure
    from those few profiles all lost). The 1349 training screens have complete labels, so each one's
    own optimal direction can be computed offline from hundreds of hits. Average the directions of
    the most similar training screens — reusing the similarity the screen-kNN scorer already
    computes — into a prior `p`, and shrink the live estimate toward it:
    `q(alpha) = normalise((1 - alpha) * q_own + alpha * p)`, with `alpha = 0` the deployed feature.

    `scripts/diagnostics/build_dir_bank.py` built the bank: **1289 of 1349 training screens have a
    usable direction, median 1294 atlas-covered hits each** against the ~20-30 being shrunk, so there
    was real strength to borrow. The directions also genuinely differ — pairwise cosine mean +0.433,
    sd 0.379, 5th/95th percentile −0.092/+0.949 — so *which* screens you borrow from could matter,
    and the selection was not decoration.

    | query direction | `paper_validation` | `paper_test` |
    |---|---|---|
    | **alpha 0.00 (deployed)** | **+0.343 [+0.039, +0.762]** | **+0.610 [+0.290, +0.967]** |
    | alpha 0.25 | +0.324 [+0.043, +0.701] | +0.522 [+0.196, +0.856] |
    | alpha 0.50 | +0.155 [−0.105, +0.486] | +0.387 [+0.125, +0.659] |
    | alpha 0.75 | +0.129 [−0.091, +0.399] | +0.339 [+0.074, +0.651] |
    | alpha 1.00 (prior only) | +0.019 [−0.078, +0.133] | +0.321 [+0.099, +0.598] |
    | *[ctl] shuffled bank, alpha 0.50* | *+0.055* | *+0.054* |
    | *[ctl] shuffled bank, alpha 1.00* | *−0.044* | *−0.012* |
    | *[ctl] drop nearest 5, alpha 0.50* | *+0.193* | *+0.340* |
    | *[ctl] drop nearest 5, alpha 1.00* | *+0.029* | *+0.157* |
    | 2 columns: own + prior | +0.390 [+0.101, +0.775] | +0.659 [+0.313, +1.036] |
    | *[ctl] 2 columns: own + shuffled* | *+0.399 [+0.090, +0.781]* | *+0.705 [+0.301, +1.138]* |

    **Four readings, in order of how much they settle.**

    - **As a blend it is refuted: the margin declines monotonically in alpha on both sets, and
      alpha = 0 is the best cell of five.** Replicated across two screen sets with the same shape.
      This is the cleanest negative dose-response in the file — compare §4.5, where the denoising
      family had no dose-response at all and was therefore only arguably null.
    - **As a separate column it is refuted by its own control, which is the decisive measurement.**
      Two columns beat one (+0.659 against +0.610 on test), which looks like additivity — but the
      *shuffled* bank scores the same or higher (+0.705 test, +0.399 validation). The second column's
      apparent lift is width, not information. Without that control this would have been written up
      as a small positive and taken to a campaign.
    - **The prior does carry standalone signal on test and it does not replicate.** alpha = 1 uses
      none of the live screen's own profile labels and still scores +0.321 on test against a shuffled
      control of −0.012 — real, screen-specific, borrowed information. On validation the same arm is
      +0.019. One set out of two, and `drop nearest 5` halves the test figure (+0.321 → +0.157), so
      part of what little there is comes from the closest neighbours, which is where near-duplicate
      contamination would sit (§6: hit-set Jaccard reaches 1.00 within a curated set).
    - **Why it fails, and it is not the sample size.** `cos(q_own, p)` is **+0.237 validation /
      +0.246 test** — the borrowed direction points somewhere substantially different. Borrowing
      strength only pays when the borrowed quantity estimates the *same* thing; here the training
      screens' optimal directions genuinely differ (cosine sd 0.379), so a neighbour's direction is a
      biased estimate of this screen's, and the bias swamps the variance it removes. Shrinking a noisy
      estimate toward a precise estimate of something else is not shrinkage, it is substitution.

    **Note the regime this was measured in favours shrinkage.** The seeds come from rounds 1-3, i.e.
    the ~20-30-hit regime where `q_own` is at its noisiest and a prior has the most to offer. By round
    8 the live screen has hundreds of labels and `q_own` is better estimated still, so the deployed
    effect would be smaller than measured here, not larger. The negative is therefore conservative.

    **The contrast with screen-kNN is the transferable lesson.** Borrowing from similar training
    screens at the **gene** level works (+0.033 to +0.145 held out, §4.14); borrowing from them at
    the **direction** level does not. Gene-level borrowing needs only "similar screens share hits",
    which is true. Direction-level borrowing needs "similar screens share a profile-space geometry",
    which the cosine measurements say is false. Future "borrow across screens" ideas should be asked
    which of those two they assume.

    Code: `build_dir_bank.py` (bank builder, prints the pairwise-cosine spread), `screen_shrinkage.py`
    (the alpha sweep and all three controls), and `codep_seed_acq.py` gained `dir_bank`,
    `shrink_alpha`, `shrink_k`, `shrink_drop_top`, `shrink_min_hits` and `shrink_shuffle`. **The
    acquisition code is kept and is harmless: `shrink_alpha=0` is the default and is bit-identical to
    the shipped behaviour, test-pinned.** 11 new tests (34 in `tests/test_codep_seed.py`). Per this
    file's own rule an offline-negative candidate does not get a campaign, so no sweeps exist for it.


17. **Perturb-seq as a second profile atlas: real signal, but 2-4x weaker than DepMap on the
    calibrated instrument. Does not clear the bar alone. The methodological lesson is bigger than
    the result.** Prompted by the question of whether other gene embeddings (ESM-2, Borzoi,
    Geneformer) could help. The structural argument from §4.9 says a fixed per-gene embedding
    competes with `v_g` at its own job and should be expected to land at ~0 — which is what the four
    previous features did. The exception is the form that works: not a COLUMN but a KERNEL, where
    the atlas supplies a similarity used to propagate this screen's confirmed hits. Replogle 2022
    genome-scale Perturb-seq is the best candidate of that kind — a measured perturbation phenotype
    like DepMap, but on a transcriptional rather than a fitness axis.

    Data: `K562_gwps_normalized_bulk_01.h5ad`, figshare+ article **20029387** file id **35773217**,
    375 MB, free and public, via `api.figshare.com/v2/file/download/<id>` (the `ndownloader` route
    hangs, as with DepMap). `perturbseq_to_npz.py` converts it to the same `.npz` layout every
    consumer here already reads: **8248 dims x 9866 genes** from 11258 perturbation rows (514
    non-targeting controls dropped, 788 genes averaged over multiple protospacers). Only **5191 of
    9866** genes show a detectable transcriptional response (energy test p < 0.05), which bounds how
    much of the atlas is usable.

    | instrument | DepMap | Perturb-seq |
    |---|---|---|
    | **EF oracle, diff-of-means, validation** | **+0.342 … +0.365** | **+0.153 [−0.052, +0.396]** |
    | **EF oracle, diff-of-means, test** | **+0.559 … +0.656** | **+0.131 [−0.164, +0.457]** |
    | EF oracle, centroid, validation / test | +0.126…+0.158 / +0.444…+0.467 | −0.052 / +0.088 |
    | kernel AUROC tail, validation | 0.552 | **0.588** |
    | kernel AUROC tail, test | **0.586** | 0.562 |
    | kernel tail top-1k hit rate, validation (base 4.30%) | 7.07% | **11.96%** |

    **Verdict: not worth a campaign on its own.** Both EF-oracle intervals include zero, against
    DepMap's which exclude it, and `centroid` is negative on validation. It is a real feature — the
    controls below establish that — but a weaker one.

    **The instrument lesson, which is the part to carry forward.** I built
    `screen_kernels.py` to pre-screen candidate spaces cheaply, and it **measures the wrong
    quantity**: standalone ranking ability inside a space, not incremental value over `span(V)`.
    Those two diverge completely — §4.11 is the precedent (LLM score, AUROC 0.58 standalone, 0.00
    incremental) and this is the second instance, with the two instruments ranking DepMap and
    Perturb-seq in **opposite orders**. Including `v_g` as a floor row does not fix it, and the
    DepMap reference row actively refutes the decision rule I proposed: **DepMap scores below `v_g`
    on this measure (0.552 vs 0.634) and is still worth +0.175 EF deployed**, because a kernel needs
    to be complementary to `v_g`, not better than it. Use `augment_u.py` / `screen_seed_aggs.py`,
    which fit the added column jointly with û and whose scale is calibrated against a deployed
    result. `screen_kernels.py` is kept and its docstring now carries this warning; it is a smoke
    test for "does this space carry gene information at all", nothing more.

    **Two controls that earned their keep.**
    - **A positive control on known complexes caught a silent bug.** The first run scored
      Perturb-seq at AUROC 0.492 — indistinguishable from its own shuffle, which reads as a clean
      null. The cause was **6832 infinities** in the atlas (genes with zero control variance, so the
      z-score divides by 0) and `np.nan_to_num`'s default mapping of `inf` to 3.4e38, which makes
      the norm of any row holding one astronomical and its cosine with every query ~0. Checking that
      known complexes cluster (mediator +11.3 sd, ribosome +9.3 sd, RNA pol II +4.9 sd above random
      pairs — all at or above DepMap's figures) exposed it immediately. **`posinf=0.0, neginf=0.0`
      is now set in `codep_seed_acq.py`, `screen_kernels.py` and `screen_seed_aggs.py`**; DepMap has
      no infinities so nothing previously measured is affected, but the failure mode is silent and
      would recur with any z-scored atlas.
    - **The `[shuf]` control** (permuting which gene owns which profile row) sits at 0.49-0.51 for
      every space on both sets, so the signal above it is gene correspondence rather than marginal
      distribution.

    **What is still open, and it is the only configuration that could justify a campaign:**
    complementarity. DepMap measures fitness across 808 cell lines; Perturb-seq measures
    transcriptional response in one. The axes are different, so the pair may add as codep and
    screen-kNN did (§4.14). Testing it needs a 2-column block drawn from two different `.npz` files,
    which `screen_seed_aggs.py` cannot currently do (it takes one `--profile-npz`), and deployment
    would need a second profile source in `codep_seed_acq.py`. Both are small; neither is done.
    Treat the combos in the Perturb-seq screen output (`cen+diff+shrunk` +0.371/+0.378) with the
    §4.15 scepticism: multi-direction blocks buy capacity, not signal.

    Code: `perturbseq_to_npz.py`, `screen_kernels.py` (with the limitation documented in its
    docstring). Raw data at `~/perturbseq/`, converted atlas at `~/perturbseq/k562_gwps.npz`.
    `rpe1_normalized_bulk_01.h5ad` (95 MB, file id 35775581) is a second cell line if wanted.

    **On the other embeddings asked about**, none screened, priors from this result: ESM-2 and
    Borzoi are *sequence-derived* and so cannot be conditioned on the screen at all — they would
    enter as a kernel too, but their geometry reflects protein family and structure rather than
    measured perturbation response, which is one step further from the target than Perturb-seq
    already is. Geneformer / scGPT gene embeddings and a co-expression atlas (ARCHS4, GTEx) are the
    cheaper next candidates. Screen any of them with `screen_seed_aggs.py --profile-npz`, not with
    `screen_kernels.py`.


18. **DepMap + Perturb-seq, and the triple with screen-kNN: consistently positive in sign on both
    screen sets, and too small to deploy or even to confirm. STOP HERE for this line of work.**
    The §4.17 follow-up, run with `screen_pair.py` (multi-atlas blocks, difference-of-means
    directions, every block against a width-matched 5-draw Gaussian control).

    | feature block | cols | `paper_validation` | `paper_test` |
    |---|---|---|---|
    | depmap (deployed single) | 1 | +0.375 [+0.071, +0.804] | +0.606 [+0.278, +0.965] |
    | pseq | 1 | +0.149 [−0.049, +0.383] | +0.165 [−0.145, +0.499] |
    | knn | 1 | +0.136 [−0.001, +0.305] | +0.238 [−0.136, +0.734] |
    | depmap+pseq | 2 | +0.444 [+0.080, +0.900] | +0.709 [+0.340, +1.128] |
    | **depmap+knn (DEPLOYED)** | 2 | +0.340 [+0.059, +0.709] | **+0.874 [+0.441, +1.373]** |
    | depmap+pseq+knn | 3 | +0.384 [+0.092, +0.740] | +0.941 [+0.475, +1.484] |

    **The block means are not the answer — the paired contrasts are**, because a third column adds
    fitting capacity whatever is in it. Both sides of each contrast are already net of their own
    width-matched noise:

    | contrast | `paper_validation` | `paper_test` |
    |---|---|---|
    | depmap+pseq **vs** depmap | +0.069 [−0.052, +0.200] 11/20 | +0.103 [−0.017, +0.255] 13/19 |
    | depmap+pseq+knn **vs** depmap+knn (the deployed config) | +0.044 [−0.089, +0.214] 8/20 | +0.066 [−0.076, +0.253] 9/19 |
    | depmap+knn **vs** depmap (for calibration) | −0.035 [−0.173, +0.093] 8/20 | +0.269 [+0.010, +0.678] 11/19 |

    **Verdict: do not build it.** Adding Perturb-seq to the deployed pair is worth **+0.044 / +0.066
    offline**, neither interval excluding zero, with win rates at or below half on both sets. Effects
    shrink ~3x off these curated screens (§1), so the expected deployed gain is **~+0.02 EF against
    the current +0.175** — and the 198-screen held-out fold resolves to about ±0.1 EF, so **the
    effect could not be demonstrated even if it is real.** That is the decisive argument: not that it
    is zero, but that it is below the resolution of the best confirmation set available. Deploying it
    would mean carrying a second 375 MB atlas and a second profile source in the acquisition for a
    gain nothing in this project could measure.

    **A hypothesis of mine, refuted inside this diagnostic.** I expected Perturb-seq and screen-kNN
    to be substitutes — both are "what co-occurs with these hits elsewhere" signals — and that this
    would explain the flat triple. The per-screen correlations between the three feature columns on
    library genes are **+0.113 (depmap|pseq), +0.114 (depmap|knn), +0.119 (knn|pseq)**: uniformly
    low and essentially equal. They are **near-orthogonal, not redundant.** So the triple is flat
    because Perturb-seq's own contribution is small relative to this instrument's resolution, not
    because the information is already present. Orthogonality is necessary for additivity, not
    sufficient.

    **Read the null with the instrument's known sensitivity in mind.** §4.14 recorded that on
    `paper_test` the deployed pair was *indistinguishable* from codep alone offline, and yet deployed
    on 198 held-out screens it was worth +0.04 (+0.145 against +0.106). So this offline instrument at
    n ≈ 20 has already missed one real effect of exactly the size now being measured for Perturb-seq.
    A null here is therefore weak evidence of absence — which is precisely why the argument above
    rests on unmeasurability rather than on the null.

    Code: `screen_pair.py` (`--atlas LABEL=PATH`, repeatable; prints the inter-column correlations
    and the paired-contrast table). No change was made to `codep_seed_acq.py`: it still takes one
    `profile_npz`, and a second source was deliberately not built.


19. **Joint weight sweep on 198 held-out screens: (0.5, 0.5) IS the optimum, on both axes, and
    every variant is worse. The deployed configuration is confirmed, not improved.** Five new arms
    on the `validation` fold (218 screens, 198 fresh), all `query=diffmeans fuse=rrf`, scored with
    **paired contrasts** against the deployed arm so the screen-to-screen variance cancels. Run
    because the weight surface had only ever been explored as a cross with one point in the
    interior: every arm with `knn_weight > 0` used `weight = 0.5`, and every arm with
    `weight != 0.5` used `knn_weight = 0`.

    | arm | dEF/screen | 95% CI | dEF/family | 95% CI |
    |---|---|---|---|---|
    | **DEPLOYED c0.5 / k0.5** | **+0.175** | [+0.074, +0.291] | **+0.331** | [+0.180, +0.515] |
    | c0.5 / k0.25 | +0.143 | [+0.041, +0.253] | +0.296 | [+0.150, +0.460] |
    | c0.5 / k0.75 | +0.103 | [+0.014, +0.199] | +0.191 | [+0.072, +0.333] |
    | c0.5 / k1.0 | +0.043 | [−0.073, +0.166] | +0.150 | [−0.005, +0.335] |
    | c0.25 / k0.5 | +0.090 | [+0.023, +0.161] | +0.179 | [+0.080, +0.292] |
    | c0.75 / k0.5 | +0.106 | [−0.024, +0.244] | +0.279 | [+0.098, +0.481] |

    Paired contrasts against the deployed arm (negative = worse than 0.5/0.5):

    | arm | per screen | 95% CI | per family | 95% CI | wins |
    |---|---|---|---|---|---|
    | c0.5 / k0.25 | −0.032 | [−0.091, +0.026] | −0.036 | [−0.121, +0.046] | 64/198 |
    | c0.5 / k0.75 | **−0.072** | **[−0.147, −0.003]** | **−0.140** | **[−0.235, −0.055]** | 55/198 |
    | c0.5 / k1.0 | **−0.132** | **[−0.229, −0.040]** | **−0.181** | **[−0.309, −0.065]** | 54/198 |
    | c0.25 / k0.5 | **−0.086** | **[−0.166, −0.019]** | **−0.152** | **[−0.272, −0.056]** | 71/198 |
    | c0.75 / k0.5 | −0.070 | [−0.150, +0.014] | −0.053 | [−0.147, +0.048] | 53/198 |

    **Both axes peak at exactly 0.5 and fall away on either side.** Three of the five variants are
    significantly worse on both weightings; the other two are negative but not significant. Nothing
    is better. This is the first time the fusion weights have been explored in the interior on
    held-out data, and the answer is that the setting chosen from the curated 20s **replicates as
    optimal on 198 independent screens** — which is a genuine strengthening of §4.15 even though no
    improvement was found.

    **My prediction was wrong and the direction matters.** I predicted `knn_weight = 0.75` would be
    the best cell at +0.02 to +0.04, reasoning that the axis had only two measured points and both
    trended upward (0.25 → 0.5 was +0.024 on this fold). It is **−0.072 / −0.140, significantly
    negative.** The two-point trend was an artefact of 0.25 being below the optimum, not evidence of
    a slope continuing past it. Related: the §4.13 "dominance threshold" argument (a scorer gains the
    ability to seat its own top pick at the head of a batch between weight 1.0 and 2.0) correctly
    predicted that higher is worse but put the turn in the wrong place — the decline begins
    immediately above 0.5, well before any threshold is crossed.

    **The headline number is NOT contaminated by this tuning.** Measuring five more arms on the same
    fold would introduce selection only if the best-looking one were promoted. The sweep confirmed the
    incumbent, so nothing was selected and **+0.175 remains a held-out measurement of a
    configuration chosen on different data** (`paper_test`, §4.15).

    **Instrument note, since the first run of the scorer printed a wrong interval.**
    `score_ctx.boot_ci` resamples FAMILIES and returns the family-weighted mean, so its CI belongs
    beside the per-family figure. The first version printed it beside the per-screen mean, producing
    an interval that did not contain the number next to it ([+0.180, +0.515] against +0.175). The
    tell was that it did not match the recorded +0.175 [+0.075, +0.293]. `score_weight_sweep.py` now
    has a separate `boot_screens` that resamples screens, and the deployed row reproduces §4.15's
    recorded figures on both weightings — which is what makes the rest of the table trustworthy.
    **Check a new scorer against a known row before reading its new ones.**

    Code: `run_weight_sweep.sh` (five arms, skips existing), `score_weight_sweep.py`
    (`--arm LABEL=sweep-id`, paired contrasts, `--exclude-set` drops the curated screens).
    Two operational fixes in the runner, both from §8's list: it calls `.venv/bin/assayloop`
    directly because **`uv run` exited 120 with empty logs on the first three arms** — inside an
    unattended sweep that is indistinguishable from an arm that ran and produced nothing — and its
    skip-guard matches any index width, since `${id}-[0-9][0-9]-*` only covers indices 00-99 and
    would silently skip a half-finished arm on a 218-screen fold.


20. **`min_seeds = 2` is already optimal. Raising it is monotonically worse, and the reason is
    instructive.** Audited free from the existing deployed runs with `codep_trace_audit.py` on the
    198 fresh screens (`--off sweep-dm-validation-w0 --on sweep-dm-validation-dm05knn`). The idea was
    that the feature fires as soon as 2 hits exist, when the query direction is at its noisiest, so
    delaying activation should raise quality.

    | rule | fires on | mean dEF on those screens | mean dEF overall |
    |---|---|---|---|
    | **current: min_seeds 2** | 155/198 | +0.225 | **+0.176** |
    | min_seeds 5 | 115/198 | +0.232 | +0.135 |
    | min_seeds 10 | 81/198 | +0.301 | +0.123 |
    | min_seeds 20 | 54/198 | +0.222 | +0.060 |
    | min_seeds 10 + observed rate <= 0.40 | 69/198 | +0.352 | +0.123 |
    | min_seeds 10 + observed rate <= 0.25 | 66/198 | +0.368 | +0.123 |

    **The premise is right and the conclusion is wrong.** Delaying activation genuinely improves the
    CONDITIONAL effect (+0.225 -> +0.301 at min_seeds 10), so the direction really is noisier with
    few seeds. But the **117 screens that never reach 10 confirmed hits still average +0.089**
    (29 up / 41 down), so silencing them discards more than the noise costs. Every threshold is worse
    overall, monotonically. This is the same shape as the saturated-screen gate (§4.13 correction):
    a subgroup where the effect is weaker is not a subgroup where switching off helps.

    Caveat on the instrument: `min_seeds` cannot be scored exactly from finished runs, because
    raising it changes which genes are picked in the rounds that still fire, and a sequential
    campaign cannot be spliced. The overall column assumes a silenced screen scores exactly its
    weight-0 run, which is exact for a per-screen gate and approximate here. The decline is monotone
    over four values, so a rerun is unlikely to reverse it, but a campaign would be needed to state
    it as a measurement rather than a prediction.

    **"Optimal" OVERSTATES IT: only values >= 2 were tested.** The sweep above goes upward from the
    default. `min_seeds = 1` has never been run, and §4.23 found a measured reason to try it — 42 of
    198 fold screens never fire the feature at all because they never reach 2 confirmed hits.

    **Process note, recorded because it has now happened twice in one session.** Both this and the
    `knn_weight = 0.75` prediction (§4.19) were reasoned from a plausible mechanism and both were
    wrong in the same direction: I inferred that a parameter should move without measuring it, when
    the measurement was free and took under a minute. `codep_trace_audit.py` and
    `score_weight_sweep.py` both answer this class of question from runs that already exist.
    **Screen a parameter intuition before proposing it, not after.**


21. **GENE-GENE CO-HIT: the first candidate since the deployed method to beat it on a paired
    contrast. Worth a campaign — the only open item that is.** Prompted by asking what MODEL, rather
    than what data source, could add information.

    **The principle first, because it rules most of the question out.** For a FIXED feature column,
    model class is not the constraint and this is already measured: `augment_u.py` fits the column
    jointly with û against the full answer key, so it reports the upper bound achievable by *any*
    model using that column additively at an optimal weight. Perturb-seq's bound was +0.15; no
    estimator recovers more than that from Perturb-seq. §4.9 said it directly — the binding
    constraint is feature quality, not capacity. So "a better model over the same information" is a
    dead end BY CONSTRUCTION. The gap in that argument is that it bounds a fixed column: a different
    model over the same raw data producing a DIFFERENT column is not bounded by the old one.

    **The idea that exploits the gap.** screen-kNN is screen-screen — find training screens
    resembling this one, vote on genes. Co-hit is **gene-gene** — score a gene by how often it is a
    hit *in the same screens* as this screen's confirmed hits. Same 1349 x 20874 matrix, orthogonal
    direction through it, and the exact analogue of what codep does with DepMap using the benchmark's
    own labels. Motivation: `v_g` IS approximately a rank-10 factorisation of this matrix and only
    3.5% of a test-fold hit vector lies in `span(V, b, 1)`, so rank-10 gene-gene structure is already
    in the model and the other 96.5% is not. Per-gene centring removes marginal hit frequency, which
    §10 records as useless alone, leaving pure co-occurrence.

    | feature block | `paper_validation` | `paper_test` |
    |---|---|---|
    | depmap (deployed single) | +0.373 [+0.080, +0.778] | +0.607 [+0.280, +0.962] |
    | **cohit** | **+0.233 [−0.001, +0.482]** | **+0.522 [+0.134, +1.045]** |
    | knn | +0.132 [−0.007, +0.302] | +0.234 [−0.131, +0.728] |
    | depmap+cohit | +0.452 [+0.168, +0.817] | +0.950 [+0.449, +1.512] |
    | depmap+knn (DEPLOYED) | +0.332 [+0.072, +0.697] | +0.879 [+0.445, +1.378] |
    | **depmap+cohit+knn** | +0.376 [+0.065, +0.787] | **+1.076 [+0.600, +1.628]** |

    Paired contrast that matters, against the deployed pair:

    | contrast | `paper_validation` | `paper_test` |
    |---|---|---|
    | depmap+cohit+knn **vs** depmap+knn | +0.044 [−0.061, +0.154] 9/20 | **+0.198 [+0.064, +0.339] 14/19** |
    | depmap+cohit **vs** depmap | +0.079 [−0.044, +0.220] | +0.343 [+0.032, +0.807] |

    **Why this is a stronger case than Perturb-seq (§4.18) despite the same split-decision shape.**
    Perturb-seq gave +0.044 / +0.066 with BOTH CIs spanning zero, i.e. the whole plausible range sat
    below what any confirmation set resolves. Co-hit gives +0.198 with a CI excluding zero on the
    set with 19 independent families. Applying the observed ~3x shrinkage, the expected deployed gain
    is **+0.015 (from validation) to +0.07 (from test)** on top of +0.175, and the 198-screen fold
    resolves about ±0.05 on a paired contrast — so the upper half of that range is demonstrable.

    **The non-replication has a visible mechanism, which is why it is a split and not a refutation.**
    `cohit|knn` column correlation is **+0.559 on validation against +0.283 on test**. On validation
    the two are substitutes (above the 0.3 threshold), so they cannot add there; on test they are
    complements. `paper_validation` is full of near-duplicate families, so screen-screen and
    gene-gene similarity converge on it — which is a property of that screen set, not of the method.
    `cohit|depmap` is +0.137/+0.140, near-orthogonal on both.

    **CAMPAIGN RUN 2026-10-22 — CLOSED NEGATIVE. Co-hit adds nothing deployed.** Three arms on the
    218-screen `validation` fold, paired against the deployed `sweep-dm-validation-dm05knn` on the
    198 fresh screens:

    | arm (depmap / cohit / knn) | dEF/screen vs deployed | 95% CI | dEF/family | 95% CI | wins |
    |---|---|---|---|---|---|
    | 0.5 / 0.25 / 0.5 (total 1.25) | **+0.006** | [−0.024, +0.039] | **+0.030** | [−0.020, +0.085] | 52/198 |
    | 0.5 / 0.5 / 0.5 (total 1.5) | −0.006 | [−0.064, +0.048] | +0.020 | [−0.054, +0.098] | 53/198 |
    | 1/3 each (total 1.0) | −0.067 | [−0.145, +0.001] | **−0.110** | **[−0.227, −0.019]** | 76/198 |

    **The offline margin overstated the deployed gain by 7-30x** (+0.198 offline on `paper_test`
    against +0.006 per screen / +0.030 per family deployed). That is the FIFTH instance of this
    pattern and the worst of them, and it happened on the candidate with the strongest offline case —
    a CI excluding zero on the 19-family set. Revise the rule accordingly: **an offline paired
    contrast with a CI excluding zero still does not predict a deployed gain.** It licenses the
    campaign and nothing more.

    **The one transferable finding: fusion weights are ABSOLUTE DOSES, not shares of a budget.** The
    `eq` arm carries total promotion pressure 1.0, identical to the deployed 0.5 + 0.5, and is
    significantly worse per family (−0.110). So diluting the two good scorers to 1/3 each to make
    room for a third costs real ground, while bolting a third on top at 0.25-0.5 (total 1.25-1.5) is
    merely useless. This agrees with §4.19 finding both axes peak at 0.5 INDEPENDENTLY, and it is the
    design rule for any future term: **add, never reallocate.**

    Code: `codep_seed_acq.py` gained `profile_npz2` + `weight2`, loaded only when `weight2 != 0`;
    `weight2=0` is bit-identical to single-source and test-pinned. Profile loading is factored into
    `_load_profiles` and `_seed_query` is parameterised by source, so ANY atlas in the shared
    `(dims x genes)` layout can be a second term — the extension is not co-hit-specific and is worth
    keeping. 7 new tests (41 in `tests/test_codep_seed.py`). Runner `run_cohit.sh`, atlas at
    `output/analysis/cohit.npz` (dense 1349 x 20874 from `train_hits.mat.npz`).


22. **LEAKAGE AUDIT of the training-label components, and the headline result PASSES.** Both
    screen-kNN and co-hit read the training hit matrix at inference, justified in §4.14 because "the
    folds are disjoint". **Disjoint folds are not the same as no near-duplicate assays across
    folds** — the split is `yearfold0`, i.e. by year, and §6 documents hit-set Jaccard up to 1.00
    between screens *within* one curated set. If an evaluation screen has a twin in the training
    fold, screen-kNN's top-50 would be weighting that twin heavily and effectively reading its
    answer key. This had never been checked.

    Maximum hit-set Jaccard against any of the 1349 training screens:

    | screen set | median | max | above 0.3 | above 0.5 |
    |---|---|---|---|---|
    | `paper_test` | 0.064 | 0.351 | 1 of 20 | 0 |
    | `paper_validation` | 0.106 | 0.567 | **8 of 20** | **6 of 20** |
    | the 198 fresh fold screens | 0.067 | 0.543 | 18 of 160 | 4 of 160 |

    Contamination is real but concentrated in `paper_validation`'s PC-9 fitness family — the same 8
    screens behind the saturated-screen artefact (§4.13 correction). `U_1795_dec` shares **1291 of
    its 1629 hits** with training screen 1458. Those screens are also budget-saturated, where the
    feature cannot help anyway.

    **The decisive check: the deployed gain is strongest on the CLEANEST screens.**

    | max Jaccard | n | dEF/screen | dEF/family |
    |---|---|---|---|
    | **[0.0, 0.1)** | **110** | **+0.193** | **+0.337 [+0.166, +0.527]** |
    | [0.1, 0.2) | 25 | +0.358 | +0.416 [−0.104, +1.186] |
    | [0.2, 0.3) | 7 | −0.001 | −0.010 [−0.200, +0.226] |
    | [0.3, 1.0] | 18 | +0.054 | +0.298 [+0.006, +0.796] |

    On the 110 screens with Jaccard below 0.1 the effect is **+0.193 / +0.337**, at or above the
    overall +0.175 / +0.331. So **+0.175 is not an artefact of borrowing from near-twins**, and the
    obvious reviewer objection to using training labels at inference is answered with data rather
    than with the fold-disjointness argument alone. Quote this if anyone asks.

    Note co-hit is structurally LESS exposed than screen-kNN to this channel: it aggregates
    gene-gene co-occurrence over all 1349 screens, so one twin carries little leverage, whereas
    screen-kNN's top-50 weighting can concentrate on it.


23. **A rate-adaptive weight is refuted; and the 48 "unchanged" screens turn out to be a measured
    21% COVERAGE GAP, which is a different problem from the diffuse harm.** Both checked free from
    runs that already exist (`scratchpad/new_ideas.py` pattern: the `c075k05` arm plus the deployed
    arm plus the baseline).

    **(a) Rate-adaptive weight — REFUTED, no campaign needed.** §4.19 found weight 0.75 globally
    worse than 0.5, but "globally" pools screens where AssayFormer is strong — where any promotion
    displaces a ~40%-reliable pick (§2a) — with screens where it is near-random and there is nothing
    to protect. §4.13's mechanism predicts a higher weight should pay off exactly on the latter. It
    does not:

    | AssayFormer early hit rate | n | dEF at w=0.5 | dEF at w=0.75 | difference |
    |---|---|---|---|---|
    | < 2% | 133 | +0.143 | +0.049 | **−0.121 [−0.229, −0.002]** |
    | 2-5% | 27 | +0.323 | +0.340 | −0.016 [−0.240, +0.224] |
    | 5-20% | 19 | +0.365 | +0.511 | +0.153 [−0.093, +0.417] |
    | > 20% | 19 | +0.004 | −0.236 | −0.132 [−0.289, +0.025] |

    **Spearman(early hit rate, w=0.75 advantage) = −0.045, p = 0.53** — no relationship, and the sign
    is opposite to the hypothesis: on the 133 weakest screens the higher weight is significantly
    worse. So there is no per-screen weight schedule to be had from the observable AssayFormer is
    weak/strong signal, and this closes the continuous version of the gate idea as well as the
    binary one (§4.13 correction).

    **SCOPE CORRECTION (2026-10-23, from §12 item 5).** That sentence overreaches. What is closed is
    every weight rule keyed on **how well AssayFormer is doing** — a property of the model, and the
    only covariate tested. A weight keyed on **how well the kernel's own direction is estimated** —
    bootstrap the seed set, recompute the difference-of-means direction, take the mean pairwise
    cosine across replicates — is a different quantity and is untested. §4.20's own table is the
    signature of exactly that: the conditional effect rises monotonically with seed count
    (+0.225 at min_seeds 2 -> +0.301 at 10) while the screens below the threshold still average
    +0.089, which is a monotone confidence relationship being forced through a binary cut.

    **(b) The coverage gap: the feature never fires on 21% of screens.** The deployed split is
    75 up / 75 down / 48 same, and the 48 decompose as:

    | | count |
    |---|---|
    | never fired — never reached `min_seeds = 2` | **42** |
    | fired but changed no pick | 0 |
    | fired and changed picks, EF still exactly equal | 6 |

    Those 42 screens have a **mean total hit count of 11** and AssayFormer's early hit rate on them is
    **0.8%**: there are essentially no hits to seed from. Across the whole fold the feature fires in
    only **5.9 of 10 rounds** on average. This had never been measured, and it is a distinct problem
    from the diffuse harm on the 75 losing screens — those screens get a bad direction, these get no
    direction at all.

    **(c) `min_seeds = 1` RUN AND CLOSED NEGATIVE — keep 2.** `sweep-ms1-validation`, the deployed
    configuration with `min_seeds=1`, scored with `score_min_seeds1.py` (which splits by whether the
    feature could fire at all under the default, so the overall number cannot hide a mismatch between
    "it improved" and "it improved where predicted"). Firing rose from **5.93 to 7.06 rounds of 10**
    and silent screens fell from 42 to 21. Per-family means:

    | group | n | deployed | min_seeds=1 | difference |
    |---|---|---|---|---|
    | **ALL screens** | 198 | +0.331 | +0.294 | **−0.038 [−0.147, +0.047]** |
    | never fired at min_seeds=2 | 42 | +0.000 | +0.077 | +0.077 [−0.063, +0.281] |
    | did fire at min_seeds=2 | 156 | +0.440 | +0.369 | −0.071 [−0.201, +0.006] |
    | — of those, <= 50 total hits | 84 | +0.639 | +0.476 | **−0.163 [−0.451, +0.006]** |

    **Directionally right, net negative.** The targeted group does gain, and it is paid for by the
    screens that already worked — concentrated on sparse ones where each pick matters most. Mechanism:
    `min_seeds=1` makes the feature fire one round EARLIER using a single gene's profile as the
    direction, and that noisy extra round changes the whole downstream trajectory. A conditional rule
    (fire at one seed only once two look unreachable, e.g. after round 5) is implementable and
    leak-free, but its best case is +0.077 confined to a fifth of the screens, about +0.015 fold-wide
    — below what the fold resolves. Not worth building.

    **TWO FINDINGS WORTH MORE THAN THE ARM.**
    - **11% of screens are permanently beyond ANY seed-based method.** Only 21 of the 42 silent
      screens could be activated at all; on the other 21 **AssayFormer finds zero hits in 1000
      picks**, so there is no seed to be had at any threshold. That bounds this entire family of
      approaches — co-essentiality, co-hit, screen-kNN, pathway propagation all need a confirmed hit
      to propagate from — and it is a clean thing to state in the write-up.
    - **The newly activated screens are too volatile to ship even when positive.** Their mean change
      is +0.134 but the split is 11 down / 6 up, with a mean of 5 total hits each. On a 5-hit screen
      in a 20k library, adjusted EF swings enormously on one pick, so the positive mean is carried by
      a couple of lucky screens. **Check the up/down split before believing a mean on sparse
      screens** — the same lesson as §4.13's correction, where an 8-screen family drove a subgroup
      claim.

    **Still untested, and the only non-weight variance lever left: batch diversity in the FEATURE's
    space.** Each round takes the top 100 of the fused score, and genes promoted by the co-dependency
    term are by construction mutually correlated in the space that term uses, so one batch can spend
    30 slots inside a single module — costless if the module is right, 30 wasted slots if not.
    Deduplicating within the promoted set trades expected hits for variance, and variance is this
    method's measured weakness (75 of 198 screens worsen). Note §10's "do not rebuild" list excludes
    diversity selection in **BPMF** space, because that is AssayFormer's own geometry; diversity in
    **DepMap profile** space is a different object and is not covered by that entry.


24. **DepMap EXPRESSION as a profile kernel is significantly HARMFUL, and the campaign had been
    sitting unscored on disk since 2026-10-21.** §10 item 2 claimed the expression variant "screened
    positive offline (+0.87 test) but was never run as a campaign" — **wrong on the second half**:
    `sweep-cf-validation-expr-w05` has all 218 runs. Nobody had scored it. Scored now on the 198
    fresh screens against its own `sweep-cf-validation-w0` baseline (both arms `query=centroid`,
    weight 0.5, no kNN, so it is directly comparable to the essentiality arm):

    | arm | dEF/screen | 95% CI | dEF/family | up / down / same |
    |---|---|---|---|---|
    | essentiality kernel | +0.106 | [+0.016, +0.202] | +0.245 | 69 / 80 / 49 |
    | **expression kernel** | **−0.109** | **[−0.173, −0.046]** | −0.047 | **36 / 115 / 47** |
    | paired: expression **vs** essentiality | **−0.215** | **[−0.317, −0.120]** | −0.292 | — |

    **Offline +0.87 on test, deployed −0.109.** This is the sixth instance of an offline margin
    overstating the deployed result and **the first where the sign flipped**, which retires the last
    defence of offline screening as anything but a filter.

    **Why it matters beyond this one feature.** Expression is the closest structural analogue to
    co-essentiality available — same matrix shape (1673 lines x 19193 genes), same kernel machinery,
    same conditioning on the screen's confirmed hits — and it is actively harmful. So **"has the same
    shape as the dependency matrix" is not the property that makes a kernel work.** That substantially
    weakens the case for the obvious next atlases in this form (GTEx tissue expression, HPA
    localisation, STRING neighbourhoods, protein-language-model embeddings): all of them are closer to
    expression than to dependency in what they measure. Plausible mechanism: co-expression tracks
    lineage identity, so genes co-expressed with the hits are largely genes expressed in the same cell
    type — a big, phenotype-nonspecific set — whereas a dependency profile tracks functional
    requirement, which is far nearer to "is a hit in this screen".

    Note this also explains §4.10's scalar expression result flipping sign between screen sets: the
    signal was never there to begin with, in either form.

    **Operational lesson: score every arm you run.** This campaign existed for a day while §10 listed
    "run the expression variant" as an open task worth an hour. Before proposing an experiment, check
    whether its runs are already in `output/runs/` — §8's inventory rule applies to SCORING, not just
    to launching.


25. **WHAT KIND OF GENE IS A HIT ON THE SCREENS NOTHING CAN TOUCH — and a correction to §4.23.**
    `characterize_hard_hits.py`, free, on the 198 fresh screens. Groups them by how AssayFormer does,
    pools each group's hit genes, and measures those genes on two axes defined OUTSIDE the screen:
    how many of the 1349 training screens call the gene a hit, and the spread of its DepMap knockout
    effect across 808 cell lines (a near-flat profile carries no co-dependency signal even when
    present, so amplitude bounds what the working method could ever extract).

    | group | screens | hits | % of hits NEVER a hit in training | median training screens hit | % in atlas | median DepMap amplitude |
    |---|---|---|---|---|---|---|
    | **A. AssayFormer finds ZERO hits** | 20 | 53 | **0.0%** | **27** | 94.3% | 0.169 |
    | **B. very hard (EF < 2)** | 39 | 22837 | 1.5% | **9** | 91.1% | **0.138** |
    | C. middling (2 <= EF < 5) | 51 | 5116 | 0.3% | 39 | 95.6% | 0.172 |
    | D. AssayFormer does well (EF >= 5) | 88 | 28701 | 0.1% | **484** | 95.7% | 0.205 |
    | *(reference) 4000 random genes* | — | 4000 | *4.2%* | *8* | *82.6%* | *0.137* |

    **CORRECTION to §4.23: the zero-hit screens are DEGENERATE, not unreachable, and they bound
    nothing.** §4.23 said 21 such screens "bound this whole approach". Wrong. Their hits are *better*
    represented in training than average (median 27 training screens against 8 for random genes),
    fully present in DepMap, with above-average amplitude — so it is neither a coverage nor a
    representation failure. The actual reason is that **these screens have 1-4 total hits each**: ten
    have exactly ONE hit, in a 681-gene focused library evaluated against a ~20k full-genome candidate
    pool. All 20 together hold **53 of the fold's ~56,700 hits, i.e. 0.09%.** Because adjusted EF on a
    1-hit screen leaps if you happen to find it, they are **noise amplifiers**, which is also why
    `min_seeds=1` produced 11 down / 6 up on them with a positive mean (§4.23).

    **Do not target them.** Any method aimed at the zero-hit screens — seedless direction prediction,
    a reserved seed-acquisition quota, lowering `min_seeds` — is buying variance, not EF. This
    diagnostic retires that whole line of work for an hour's compute.

    **The data-coverage story is REAL but it belongs to group B.** 39 screens holding **22,837 hits**
    whose genes are barely distinguishable from random on either axis: median training-hit frequency
    **9** against 8 for random genes and 484 for the easy screens, and median DepMap amplitude
    **0.138** against 0.137 for random and 0.205 for the easy ones. Two independent measurements
    agree: on the genuinely hard screens the hits are genes the training matrix has almost never seen
    as hits AND whose dependency profiles are nearly flat. That is why the working method extracts
    little there, and it is the one regime where **more training screens** (BioGRID ORCS; the paper's
    Fig. 4C shows no plateau at 1349) would plausibly help more than a better model would.

    **The practical reading for prioritisation.** The easy screens are easy because their hits recur
    constantly across training (median 484 of 1349 screens). The hard screens are hard because theirs
    do not recur at all. That is a statement about the DATA, not about AssayFormer, and it predicts
    that no reweighting of the existing 1349 screens — which is what every scorer here does — can fix
    group B. Sharpen the claim in §4.12 accordingly: the ~10 EF of task headroom sits on screens whose
    hits the training corpus has essentially never observed.


26. **TWO FINDINGS THAT REDIRECT THE FEATURE HUNT: the richest text field is never embedded, and
    drug screens carry 80% of the hard-screen hits.** Both from reading the data rather than running
    anything.

    **(a) `screen_description_text` omits `notes`, which is where the experiment is actually
    described.** `amortized/text_embed.py:~80` renders exactly these: `phenotype`,
    `cleaned_phenotype`, `cell_line`, `cell_type`, `organism`, `library_methodology`,
    `condition_clause`, `contrast_label`, `description`. Its own docstring notes that public screens
    "often leave `description` empty". Meanwhile the ScreenRecord also carries **`notes`** (a full
    paragraph: *"A genome-wide CRISPR screen was combined with high-speed image-enabled cell sorting
    to quantify the nuclear translocation of fluorescently tagged..."*), plus `ranking_rationale`,
    `screen_rationale`, `significance_criteria`, `library_type`, `time_point` and `direction_str` —
    **none of which are embedded.** §4.7 flagged only `direction_str`; the real gap is far larger.

    **CORRECTION 2026-10-23 (§4.45): `direction_str` and `time_point` are EMPTY on every
    screen in both folds, as are the already-rendered `contrast_label` and `description`.**
    This section implied `direction_str` was available and merely unrendered. It is not
    populated at all, so it cannot be added — the direction has to be CLASSIFIED from
    `phenotype` + `notes` (100%/98% coverage). `notes` itself is the real find: 100% and
    96.4% populated at 233-355 chars. See §4.45's audit table before planning anything here.
    The description channel is worth ~0.9 EF on test (3.017 -> 3.930) and is the one input that can
    change without retraining, so this is the largest untried lever in the project. Bounded by the
    oracle-û ceiling (+0.41 / +0.84) since it feeds û, cheap to try (~$0.001 of
    `text-embedding-3-small`), and off-distribution risk is real because the checkpoint was trained
    on the current renderer.

    **(b) The hard screens are mostly DRUG-RESPONSE screens, and the compounds are named.** Of the 36
    hard (EF < 2) non-degenerate held-out screens, **25 are drug/compound screens holding ~18,000 of
    the 22,582 hard-screen hits**, and the compound appears verbatim in `phenotype`: gemcitabine,
    irinotecan, oxaliplatin, 5-FU, in PANC-1 and Bx-PC3 backgrounds. Several are CRISPRa
    ("Activation of the hit gene confers resistance"), which matters because DepMap is a knockout
    atlas. Only 6 hard screens are virus/pathogen (1024 hits).

    **The feature this implies, and why it is different from the six that failed.** For a screen
    testing drug *d*, score gene *g* by `corr(DepMap dependency profile of g, PRISM sensitivity
    profile of d)` across cell lines — genes whose knockout reproduces the drug's cell-line pattern
    are in its mechanism-of-action pathway, which is where resistance genes sit. PRISM Repurposing
    (~4500 compounds x ~900 lines) is free, on the same figshare article as DepMap 24Q4. It is
    **gene x condition**, the only architecturally non-redundant form (§4.9), and unlike the LLM score
    it is **not fame-biased**.

    **A PRE-SCREEN that would have killed four of the six previous nulls, worth applying before
    building anything:**
    1. **`corr(feature, cross-screen hit frequency)` must be near zero.** A feature correlated with
       how often a gene is a hit in general is probably inside `span(V)` — "gene is often a hit" is
       what `v_g` encodes. That is the §4.11 LLM-score failure stated mechanically.
    2. **It must have values for group B's hit genes** — those appearing in a median of 9 of 1349
       training screens (§4.25). A feature covering only well-studied genes cannot reach the
       headroom by construction.
    3. **Then `augment_u.py`**, against its noise floor, as a filter and never an estimate.

    Honest prior: additive features are **0 for 6** here, and PRISM shares the LLM score's structural
    weakness — it is a fixed prior, not conditioned on the screen's revealed labels, so it cannot
    adapt the way the working method does. It is nonetheless the first candidate to pass the
    fame-bias and coverage tests on paper, aimed at the screens holding 80% of the hard-group hits.


27. **PRISM gene x DRUG correlation: PRE-SCREENED AND KILLED in half an hour, and the failure mode
    yields the sharpest principle in the project.** The §4.26 candidate, built and gated by
    `prescreen_prism.py`. Data: Repurposing Public 24Q2, figshare article **25917643**, file id
    **46630984** (`Extended_Primary_Data_Matrix`, 72 MB, compounds x 919 ACH cell lines) plus the
    compound list (46630981) and cell-line metadata (46630978) — free, and the ACH ids align with
    `achilles.npz` rows out of the box. Feature: for a screen testing drug *d*,
    `corr(dependency profile of g, sensitivity to d)` over the 437 shared cell lines.

    | gate | result |
    |---|---|
    | 1. fame bias, `corr(feature, cross-screen hit frequency)` | **PASSED: −0.001 to −0.141** |
    | 2. coverage of each screen's hits | **PASSED: 70-95%** |
    | 3. incremental value at the oracle | **FAILED: AUROC 0.490-0.525, mean margin −0.019** |

    AUROC is **0.5 on all 15 matching screens**, and the margin over width-matched noise is −0.019
    (median +0.007, 8/15 positive) against the +0.5 build threshold. Note those screens are ids
    1639-1654 — one study, two cell lines x four drugs x KO/activation — so effectively **one or two
    families, not 15 observations.**

    **This is a DIFFERENT failure from §4.11 and that is the informative part.** The LLM score had
    AUROC 0.553-0.584 — real standalone signal, zero incremental value, i.e. redundancy with
    `span(V)`. PRISM has **no signal at all**. The mechanism: PRISM correlation measures variation
    ACROSS cell lines (which lineages the drug kills), while these screens run WITHIN one fixed line
    (PANC-1, Bx-PC3) asking which genetic modifiers change the response in that background.
    Cross-lineage sensitivity says almost nothing about within-line modifiers — the same lesson as
    the expression kernel (§4.24), where cross-line covariation turned out to encode lineage identity.

    **THE PRINCIPLE, now supported four independent ways** (expression kernel §4.24, LLM score
    §4.11, Perturb-seq §4.17, PRISM here):

    > **Gene x condition is NECESSARY BUT NOT SUFFICIENT. Every feature that is a fixed prior per
    > screen has failed. The only thing that has ever worked conditions on the screen's OWN REVEALED
    > LABELS** — DepMap co-essentiality uses the cross-line matrix only to define gene-gene
    > similarity and then points it with this screen's confirmed hits.

    This predicts the failure of the whole remaining class of "identify the condition, look it up in
    a database" features: virus-host PPI for the pathogen screens, drug-target annotation, pathway
    membership of the named condition. **Do not build them.** A candidate feature now has to answer:
    what does it do with the hits this screen has already revealed? If the answer is "nothing", the
    four results above predict zero.

    **The gates earned their keep.** Gates 1 and 2 cost minutes and passed, which is what made the
    gate-3 null interpretable as "no signal" rather than "wrong instrument" — and the whole
    pre-screen, download included, took about half an hour against the days a build would have cost.
    Keep this three-gate pattern for any future feature.


28. **BioPlex protein-interaction co-membership as a kernel: the SECOND-BEST kernel found, and still
    not worth building.** Closes the gap §4.1 left open (the PPI half of that diagnostic never ran
    because `reactome_interactions.txt` had no gene symbols). `prescreen_ppi.py`.

    Data: BioPlex 3.0, free and symbol-level —
    `bioplex.hms.harvard.edu/data/BioPlex_293T_Network_10K_Dec_2019.tsv` (118k interactions) and
    `..._HCT116_Network_5.5K_Dec_2019.tsv` (71k). 14551 proteins, **12872 of 21147 universe genes
    (61%)**, against DepMap's 83%; only 347 genes are in BioPlex but not DepMap, so it adds a
    different QUANTITY rather than reach. A gene's profile is its row of the interaction adjacency,
    kept sparse because the dense form is 21213 x 14551.

    **Why it was the best remaining candidate.** §4.25 showed the hard screens' hits have DepMap
    amplitude 0.138 against 0.137 for random genes — flat profiles — so a better estimate of the SAME
    quantity (DEMETER2 RNAi, Sanger Project Score) cannot reach them: averaging two noisy estimates of
    a flat profile gives a flat profile. Physical complex co-membership is a different estimand, and
    two genes can share a complex with uncorrelated dependency profiles. It also passes §4.27, since
    the direction is built from this screen's own confirmed hits.

    | gate | result |
    |---|---|
    | 1. fame bias, `corr(interaction degree, cross-screen hit frequency)` | **+0.203 — partial pass**, hub proteins are better studied |
    | 2. coverage of each screen's hits | pass: median 77.3%, min 55.3% |
    | 3. margin alone vs width-matched noise | **+0.272 [+0.017, +0.588] validation / +0.134 [−0.008, +0.282] test** |
    | DepMap on the same instrument, for scale | +0.334 validation / +0.625 test |
    | **DECISIVE: DepMap + PPI vs DepMap alone** | **+0.186 [+0.004, +0.390] 14/19 validation / +0.033 [−0.062, +0.140] 8/19 test** |

    **It is the only kernel besides DepMap to come out positive on BOTH screen sets** — Perturb-seq's
    intervals both spanned zero, co-hit was split. That is a real finding and makes it the best
    candidate for anyone continuing. But the decisive complementarity contrast is **one set out of
    two**, and the magnitude (+0.03 to +0.19 offline) maps to ~0 deployed given this project's
    measured shrinkage: co-hit went +0.198 offline -> +0.006 deployed, expression +0.87 -> −0.109.

    **HELD OUT AND SPLIT BY DIFFICULTY — the measurement that settles it.** Re-run on the
    `validation` fold excluding the curated 20 and excluding metric-degenerate screens
    (`--exclude-set paper_validation --min-library 1001`): **81 screens, 57 families**, far better
    powered than the curated 19.

    | group | n | PPI alone | DepMap | pair | **pair − DepMap** |
    |---|---|---|---|---|---|
    | **hard (EF < 2)** | **28** | +0.152 | +0.166 | +0.190 | **+0.024 [−0.059, +0.103]** |
    | middling (2-5) | 18 | +0.219 | +0.300 | +0.606 | +0.306 [−0.062, +0.704] |
    | does well (EF >= 5) | 35 | +0.525 | +0.618 | +0.908 | +0.290 [−0.163, +0.740] |
    | overall | 81 | +0.289 | +0.356 | +0.546 | +0.190 [+0.002, +0.387] |

    **PPI adds nothing where the headroom is.** On the 28 hard screens it contributes **+0.024**,
    interval spanning zero; its entire apparent overall value (+0.190) comes from the middling and
    easy screens, where AssayFormer is already strong. This is the one result that could have
    overturned the verdict, and it is measured on a properly powered sample.

    **CORRECTION 2026-10-23 — that paragraph's premise is one §4.12 already retracted, and it is
    the load-bearing half of this section's verdict.** "Adds nothing where the headroom is" assumes
    the headroom is on the hard screens. §4.12's correction measured the opposite and said so in
    terms: the screens where AssayFormer already does well still leave ~7 EF unclaimed and carry
    **28% of the total headroom**, mid carries 27%, and its own words are "**a method does not have
    to specialise in hard screens to find room**". Mid and easy together hold ~55% of the headroom,
    and that is exactly where BioPlex's +0.306 and +0.290 sit. The honest reading of the difficulty
    split is therefore **not** "it helps only where nothing is needed" but "it helps on the 55% of
    the headroom that is easiest to reach, and not on the 37% that is hardest".

    Two readings worth keeping:
    - **It is consistent with the partial fame-bias failure.** Interaction degree correlates +0.203
      with cross-screen hit frequency, and §4.25 found the easy screens' hits recur in a median of
      **484** of 1349 training screens against **9** for the hard ones. So PPI is largely picking up
      the hub/fame signal that helps where hits are already predictable.
    - **BOTH kernels are weak where the headroom is.** DepMap falls from +0.618 on easy screens to
      **+0.166** on hard ones — §4.25's flat-profile finding showing up directly in the margin. The
      hard screens resist both data sources for the reason already measured, which is why no further
      atlas is likely to help.

    **Verdict: do not build.** Rejected on unmeasurability, as Perturb-seq was (§4.18), with one extra
    cost: unlike Perturb-seq this is NOT a config change. The dense profile form is ~2.5 GB, so
    `codep_seed_acq.py` would need a new sparse-profile path — half a day of work for an effect that
    is +0.024 where it would need to matter.

    **CORRECTION 2026-10-23 — THE COST ESTIMATE IS WRONG, AND THE SPARSITY PROBLEM DISSOLVES
    ANALYTICALLY.** Nothing needs densifying and nothing is approximated. `_load_profiles` centres
    each row and unit-normalises it, and centring an adjacency row is what makes it dense — but the
    centring term CANCELS in the only operation that uses it. The query `q` is a difference of two
    means of centred rows, and every centred row sums to zero, so `sum(q) = 0` and therefore

        (A[g] - m_g * 1) . q  =  A[g] . q  -  m_g * (1 . q)  =  A[g] . q      exactly,

    while the norm is closed form for a binary row with `k_g` edges over `d` dimensions:

        || A[g] - m_g * 1 ||^2  =  k_g - k_g^2 / d.

    So the kernel score is `(A_sparse @ q) / sqrt(k_g - k_g**2 / d)` — one sparse matvec over
    ~190k edges, a few MB, exact rather than an approximation of the dense form. **Revised cost:
    1-2 hours including tests, not half a day, with no memory issue at all.** Both stated reasons
    for not building this are now retracted; what remains is the magnitude argument (+0.190 offline,
    ~3x shrinkage measured, so ~+0.06 deployed against a fold that resolves to +-0.1) and the
    partial fame-bias failure.

    **SETTLED 2026-10-23 BY DEPLOYMENT — §4.43.** It was built (2 hours, exact, not half a day) and
    run on 218 screens in three arms. Paired against the deployed configuration: **+0.000 per screen
    [−0.031, +0.038], −0.007 per family, 48/198 wins**, with the higher dose actively harmful and
    the degree-residualised control identical. **The conclusion of this section was right and both
    of its reasons were wrong** — which is worth keeping as a pattern, since a verdict reached for
    bad reasons is not reliable and this one had to be re-earned. The difficulty profile predicted
    here (+0.306 mid, +0.290 easy) did not appear in any cell.

    **A measurement note worth carrying.** The PPI-alone margin read +0.166 and then +0.134 on the
    same screens, changed only by adding noise draws EARLIER in the RNG order. So these margins carry
    **+-0.03 run-to-run even with 5 draws** — §4.15 said the same and this is a second instance. It
    puts the +0.033 test contrast inside its own noise, and it is the right error bar for comparing
    any two rows in a screening table.


29. **SYSTEMATIC RE-CHECK of every discounted candidate on the HARD screens — all stay discounted,
    and three new facts fall out.** Prompted by §4.28: BioPlex looked worth +0.190 overall and was
    +0.024 where the headroom is, so benefit is not evenly distributed across difficulty and
    anything discounted on a pooled average deserved re-examination. `score_by_difficulty.py` does
    this from campaigns ALREADY RUN — no new runs — on the 100 non-degenerate held-out screens
    (36 hard / 24 mid / 40 easy), bucketed by AssayFormer's own EF on the shared weight-0 baseline.

    **PAIRED contrasts against the deployed arm, hard screens only (n=36).** The unpaired
    arm-vs-own-baseline table has intervals too wide to rank arms by; differencing per screen is
    what gives power.

    | arm | delta on hard screens | 95% CI | better on |
    |---|---|---|---|
    | co-hit added (0.25) | +0.026 | [−0.003, +0.064] | 17/36 |
    | min_seeds 1 | +0.022 | [−0.003, +0.060] | **7/36** |
    | co-hit added (0.5) | +0.016 | [−0.037, +0.089] | 13/36 |
    | codep weight 0.75 | +0.012 | [−0.070, +0.093] | 17/36 |
    | codep weight 1.0 | +0.009 | [−0.098, +0.107] | 17/36 |
    | knn weight 0.25 | +0.008 | [−0.082, +0.084] | 16/36 |
    | codep alone (no kNN) | +0.004 | [−0.084, +0.075] | 18/36 |
    | knn weight 0.75 | −0.039 | [−0.111, +0.030] | 12/36 |
    | codep weight 0.25 | **−0.067** | **[−0.136, −0.004]** | 13/36 |
    | all three at 1/3 | **−0.067** | **[−0.134, −0.008]** | 15/36 |
    | centroid+kNN | **−0.064** | **[−0.130, −0.007]** | 11/36 |
    | screen-kNN alone | **−0.096** | **[−0.165, −0.037]** | 13/36 |
    | EXPRESSION kernel | **−0.156** | **[−0.274, −0.066]** | 6/36 |

    **Every positive delta has a CI including zero**; the largest is +0.026. Five arms are
    significantly worse on hard screens. **No discounted candidate is rescued.** Note `min_seeds 1`
    posts +0.022 on only **7 of 36** screens — 29 unchanged or worse with the mean pulled up by a
    couple, the volatility pattern §4.23 found and exactly what §4.25's "check the up/down split"
    rule is for.

    **THREE NEW FACTS.**
    - **On genome-scale screens the deployed method is worth +0.485 per family, not +0.331.** The 98
      metric-degenerate screens (library below the 1000-pick budget, §4.12 correction) dilute the
      headline by about a third. **+0.485 is the honest figure for screens where the metric means
      anything**, and it should be quoted alongside +0.175/+0.331 with that caveat stated.
    - **screen-kNN contributes ~nothing on hard screens**: +0.013 alone, and removing it costs
      +0.004 [−0.084, +0.075]. Its entire value is mid/easy (+0.224/+0.161). This follows directly
      from §4.25 — hard screens' hits appear in a median of 9 of 1349 training screens, so
      screen-similarity voting has nothing to surface them with. The deployed arm's hard-screen
      performance is essentially all the DepMap component.
    - **A higher codep weight preserves hard-screen performance and destroys easy-screen
      performance**: at weight 1.0, hard is +0.117 against the deployed +0.109 while easy collapses
      from +0.658 to +0.075. §4.13's displacement mechanism from a new angle, and it closes the
      difficulty-adaptive-weight idea for good — there is no hard-screen gain available from a
      higher dose, only easy-screen damage to avoid.

    Two bugs worth noting in passing: the arm-spec parser split labels at the FIRST `=`, so an arm
    labelled `min_seeds=1` silently resolved to a nonexistent sweep and reported "runs missing"
    rather than failing — now `rpartition`. And the degenerate-screen filter matters enormously
    here: without `--min-library 1001` every figure in this section is diluted by the 98 screens
    whose metric is meaningless.


30. **WHAT THE HARD SCREENS ARE, AND A MODALITY MISMATCH THAT EXPLAINS PART OF THE FAILURE.** The
    first structural characterisation of the hard screens as EXPERIMENTS, on the 100 non-degenerate
    held-out screens (36 hard / 24 mid / 40 easy).

    **The sharpest discriminator is the perturbation modality, and it is not knockout-vs-other but
    LOSS- vs GAIN-of-function.** DepMap is a CRISPR KNOCKOUT atlas, so the codep feature asks which
    genes, when knocked out, give the same cell-line fitness pattern as the genes that worked here.

    | modality | direction | n | dEF all | dEF hard only | up |
    |---|---|---|---|---|---|
    | Knockout | loss of function | 86 | +0.546 [+0.318, +0.837] | +0.103 [+0.012, +0.209] | 48/86 |
    | Inhibition (CRISPRi) | loss of function | 6 | +0.397 [+0.084, +0.734] | **+0.380 [+0.084, +0.734]** | 5/6 |
    | **Activation (CRISPRa)** | **gain of function** | 8 | **+0.022 [−0.017, +0.080]** | +0.022 | **2/8** |

    CRISPRi is knockdown — the same direction as knockout — and it works, in fact **best of any
    subgroup on the hard screens**. CRISPRa hits are genes whose OVEREXPRESSION confers resistance,
    and gain/loss-of-function phenotypes are asymmetric, so a knockout atlas has little to say. The
    exposure is total: **100% of CRISPRa screens are hard (8/8)**, 83% of CRISPRi, 27% of knockout.
    **So part of the hard-screen weakness is a modality mismatch rather than a representation
    failure.** Caveat: CRISPRi is 6 screens over only 3 families, so +0.380 rests on 3 independent
    observations — directionally clear, statistically thin.

    **Three further structural differences from the easy screens.**
    - **The readout is often not survival.** "increases/decreases cell-surface TFRC/CD71 levels",
      "increased/decreased protein binding" (**5 hard / 1 mid / 0 easy**), nuclear translocation.
      DepMap measures proliferation; a gene regulating surface TFRC need have no fitness phenotype.
    - **Genetic-background modifier screens**: "proliferation in a RIT1 M90I mutant background", "in
      the TAPR1/HAPSTR1 knockout background". Epistasis, context-specific by construction, and
      DepMap's lines are not that background.
    - **The easy screens are precisely what DepMap measures**: knockout, proliferation readout,
      essential genes. "Cell-essential genes and growth restrictive" is 13 easy / 2 mid / 3 hard.
      The kernel working there is close to tautological.
    - **Hit rate is NOT the difference**: hard median 2.06% against easy 1.36%. Hard screens are not
      sparser, they are differently shaped.

    Cell lines are mostly in DepMap (PANC-1 x8, BxPC-3 x7, Calu-3 x5, PC-9 x5), so line coverage is
    not the main issue — the exceptions are HEK-293, STAT1−/− fibroblasts and a cardiomyocyte line.

    **FULL INVENTORY of the 36 hard screens (cells and treatments).**

    | cell line | tissue | screens | hits |
    |---|---|---|---|
    | **PANC-1** | pancreatic cancer | **8** | **13597** |
    | **BxPC-3** | pancreatic adenocarcinoma | **7** | **4207** |
    | Calu-3 | lung adenocarcinoma | 5 | 1082 |
    | PC-9 | lung adenocarcinoma | 5 | 971 |
    | NALM-6 | pre-B ALL | 3 | 568 |
    | HuH-7 | hepatoma | 2 | 520 |
    | K-562 | CML | 2 | 208 |
    | HEK-293 / HeLa / STAT1−/− fibroblasts / cardiomyocyte line | — | 1 each | 965 / 293 / 2 / 169 |

    | treatment | screens | hits |
    |---|---|---|
    | Gemcitabine | 4 | 4273 |
    | Oxaliplatin | 4 | 4703 |
    | 5-fluorouracil | 4 | 4623 |
    | Irinotecan | 3 | 4205 |
    | SARS-CoV-2 Spike-RBD binding | 3 | 502 |
    | RIT1 M90I background (one also + erlotinib) | 4 | 968 |
    | Transferrin receptor (TFRC/CD71) levels | 2 | 580 |
    | HAPSTR1/TAPR1 KO background | 2 | 566 |
    | none (growth suppressors, `_inc` arms) | 3 | 501 |
    | H2O2 / doxorubicin / MERS-CoV / SARS-CoV-2 / Zika / KRAS-G12V / KH-CB19 | 1 each | — |

    **Three consequences.**
    - **It is not 36 independent problems but about EIGHT experimental programmes.** The dominant one
      is a single factorial: **4 chemotherapeutics x 2 pancreatic lines x (CRISPRa + knockout) = 15
      screens holding 17804 of 21500 hard-screen hits (83%)**. Then Calu-3 molecular readouts (5),
      PC-9 oncogene-background modifiers (5), NALM-6 (3), HuH-7 virus resistance (2), K-562 growth
      suppressors (2), four singletons.
    - **BY HITS the gain-of-function problem is DOMINANT, not marginal.** The 8 CRISPRa screens are
      22% of hard screens but hold **15478 of 21500 hits = 72%**. This section previously framed the
      modality mismatch as "a fifth of the hard screens"; weighted by hits it is nearly three
      quarters, which makes the missing gain-of-function atlas a larger gap than stated.
    - **It sharpens §4.34's condition-transfer failure.** The same four drugs run in BOTH PANC-1 and
      BxPC-3 — two pancreatic adenocarcinoma lines, about as similar a pair as exists — and hit
      overlap is still 0-7%. Background specificity is severe even between closely related lines.

    **WHAT TO TRY, and the reasoning.** The missing information is a **gain-of-function perturbation
    atlas**. The concrete candidate is **LINCS L1000 overexpression signatures** (~3.5k
    gene-overexpression perturbagens, transcriptional readout): gene-gene similarity there means
    "overexpressing these two genes produces the same transcriptional change", the GoF analogue of
    co-essentiality, and it fits the label-conditioned form so it passes §4.27. Sizing: it targets 8
    screens over 8 families, so lifting CRISPRa from +0.022 to ~+0.30 would move the hard bucket by
    about +0.07 — at the edge of measurable. Coverage is the risk at ~3.5k of 20k genes.

    **Cheaper and first: re-screen Perturb-seq on the HARD and NON-FITNESS-READOUT subsets.** Its
    readout is transcriptional response, far closer in kind to "surface protein level" or "protein
    binding" than fitness is, and it was rejected (§4.17, §4.18) on a POOLED average — which §4.28
    and §4.29 twice showed can hide difficulty-specific effects. Data already downloaded, machinery
    already written.


31. **DIRECTION-MATCHED co-hit: §4.30's PRESCRIPTION IS REFUTED, though its diagnosis stands. The
    distinction is the lesson.** §4.30 observed that every hard-screen subgroup asks an ENRICHMENT
    question while DepMap is a DEPLETION assay, and prescribed enrichment-direction training data
    (BioGRID ORCS filtered to enrichment screens). `prescreen_cohit_direction.py` tests the
    prescription with data in hand: `screen_direction.py` classifies the 1349 training screens
    (**991 dropout, 240 enrichment, 74 both, 44 unknown — the corpus is 73% depletion**), and a
    co-hit kernel is built from each subset.

    Paired on the 28 hard held-out screens, each variant ADDED to DepMap:

    | co-hit variant | training screens | delta vs DepMap alone | better on |
    |---|---|---|---|
    | **pooled (all 1349)** | 1349 | **+0.222 [+0.069, +0.420]** | 19/28 |
    | dropout only | 991 | +0.075 [−0.029, +0.183] | 18/28 |
    | **ENRICHMENT only** | 240 | **−0.061 [−0.132, +0.011]** | **8/28** |
    | enrichment + both | 314 | −0.062 [−0.136, +0.006] | 13/28 |

    **Exactly backwards.** Direction-matching to enrichment is WORSE than DepMap alone; pooled is
    best, beating even dropout-only, so it is not that depletion data helps — the corpus's DIVERSITY
    does the work. The mechanism was already in this file: §3 found enrichment-screen hits are
    condition-specific and barely recur, so a co-hit matrix built from 240 enrichment screens has
    almost no co-occurrence structure to learn. Restricting to enrichment discards the recurring
    functional-module structure that makes co-occurrence informative and gains nothing, because
    there is no recurring enrichment structure to find.

    **THE LESSON, and it generalises past this experiment.** §4.30's diagnosis is still supported —
    hard screens ask enrichment questions, DepMap is a depletion assay, their hits' DepMap
    amplitudes are flat (0.138 vs 0.137 random). What does NOT follow is "therefore use
    enrichment-direction data". **A correct diagnosis of why a data source fails does not imply that
    the matched data source will work** — the matched source has to carry structure of its own, and
    here it does not. Keep the two claims separate in the write-up.

    **Do not misread the +0.222.** That is an OFFLINE oracle margin; §4.29's DEPLOYED hard-screen
    contrast for co-hit is **+0.026 [−0.003, +0.064]**, the usual ~8x overstatement. The deployed
    figure stands.

    **WHAT IS DEAD AND WHAT IS NOT — read this carefully, the first version of this section
    overstated it.** Dead: the FILTER. "Use only enrichment-direction screens" is refuted, as is any
    other direction-matched subsetting. **NOT dead, and in fact SUPPORTED by this very table: adding
    MORE screens.** Read the variants as corpus sizes —

    | corpus | screens | delta vs DepMap on hard screens |
    |---|---|---|
    | pooled | **1349** | **+0.222** |
    | dropout only | 991 | +0.075 |
    | enrichment only | 240 | −0.061 |

    Removing 358 screens costs +0.147; removing 1109 costs +0.283. **Every subset is worse than the
    whole, monotonically in size**, and the enrichment screens contribute real value as PART of the
    corpus even though they are useless alone. Nothing here suggests saturation.

    **SUPERSEDED BY §4.33: BioGRID ORCS contains no screens the benchmark does not already use.**
    All 1901 AssayBench screens resolve to ORCS screen ids and human ORCS holds 1952, so the corpus
    cannot be grown. The reasoning below was right that corpus SIZE is the limit and wrong that the
    limit could be lifted. Kept as the record of why size is the binding constraint.

    ~~So BioGRID ORCS (unfiltered) is the strongest remaining idea~~, with three independent lines
    converging: §4.25 (hard screens' hits appear in a median of 9 of 1349 training screens, against
    8 for random genes — the corpus barely knows them), this table (co-hit value scales with corpus
    size), and the paper's own Fig. 4C (no plateau at 1349).

    **And it needs NO retraining.** `v_g` was trained on the 1349, but the co-hit and screen-kNN
    kernels merely READ a hit matrix at inference, so extending `train_hits.mat.npz` extends them
    against the frozen checkpoint. `build_train_hits.py` is 47 lines.

    **The real cost is data engineering, not modelling.** Hit definitions are heterogeneous — these
    1349 screens ship harmonised `significance_criteria`, while ORCS aggregates thousands of papers
    each with its own thresholds, so you would be CHOOSING hit calls and that choice propagates into
    every downstream number. Gene-symbol harmonisation across decades of publications. And §4.22's
    cross-fold near-duplicate audit becomes mandatory rather than optional, since ORCS may contain
    the evaluation screens themselves or near-twins of them. Expect days of work and, per §4.29's
    deployed figure for pooled co-hit, a gain near +0.026 on hard screens.

    A gain-of-function atlas for the 8 CRISPRa screens (LINCS L1000 overexpression, §4.30) is also
    NOT refuted here — that is a different modality rather than a different direction — but it
    addresses 8 of 36 hard screens and should be sized accordingly.


32. **A BETTER ESTIMATOR of co-occurrence does NOT substitute for more screens — and the way it
    fails points AT more screens.** `prescreen_cohit_ppmi.py`. Motivation: §4.25 found the hard
    screens' hits appear in a median of 9 of 1349 training screens, and §4.21's kernel used the
    crudest estimator in the family (raw binary co-occurrence, per-gene mean-centring). "Genes
    co-occurring across screens" is formally "words co-occurring across documents", where the
    rare-token problem is solved with PPMI plus context-distribution smoothing and a truncated SVD.
    The hypothesis: a better estimator could do cheaply what BioGRID ORCS would do expensively.

    Paired on the 28 hard held-out screens, each variant added to DepMap:

    | variant | delta vs DepMap alone | better on |
    |---|---|---|
    | **raw co-hit (binary, FULL rank)** | **+0.227 [+0.053, +0.437]** | **22/28** |
    | PPMI (alpha 0.75) + SVD rank 200 | +0.158 [+0.008, +0.333] | 18/28 |
    | PPMI + SVD rank 500 | +0.152 [+0.007, +0.326] | 16/28 |
    | *[ctl] plain SVD rank 200, no PPMI* | *+0.147 [−0.004, +0.330]* | 18/28 |
    | PPMI + SVD rank 50 | +0.096 [−0.023, +0.241] | 19/28 |

    - **The crudest estimator wins.** Raw full-rank binary co-occurrence beats every sophisticated
      variant.
    - **PPMI is indistinguishable from plain SVD** (+0.158 vs +0.147), so the control earns its
      place: whatever the transform does is low-rank denoising, not the PPMI correction, and the
      word-embedding framing is irrelevant here.
    - **Rank truncation LOSES signal rather than removing noise** — ranks 50/200/500, with or
      without PPMI, all sit below full rank.

    **THE INFERENCE FOR BIOGRID, and it is the reverse of what this test was pitched as.** The test
    was justified as "if a better estimator does nothing on 1349 screens, that is informative about
    whether 5000 would". It is informative — in the opposite direction. **If the matrix were
    noise-dominated for rare genes, denoising would have helped; it hurt.** So the information in
    each observation is real and usable and there are simply too few observations. Better statistics
    cannot manufacture observations; more screens genuinely can. **This does not weaken the ORCS
    case (§4.31) — it locates the limit at observation COUNT rather than at processing.** An
    inference, not a measurement of what 5000 screens would deliver, so hold it loosely.

    **Closed by this:** the whole "better estimator on the existing hit matrix" line — GloVe,
    node2vec, graph embeddings, further factorisations. Three variants plus a control all landed
    below the raw version.

    Instrument check: raw co-hit reads +0.227 here against +0.222 in §4.31, within the +-0.03
    run-to-run noise (§4.28), so the screening instrument is stable across scripts.


33. **AssayBench IS BioGRID ORCS — "add more screens" was never available.** Checked by
    downloading ORCS 2.0.18 human (`downloads.thebiogrid.org/Download/BioGRID-ORCS/Release-Archive/
    BIOGRID-ORCS-2.0.18/BIOGRID-ORCS-ALL-homo_sapiens-2.0.18.screens.tar.gz`, 102 MB gz, 668 MB
    extracted; the CDN rate-limits probing with Cloudflare error 1015, so back off and retry).

    | | |
    |---|---|
    | ORCS human screens | **1952** |
    | AssayBench train + validation + test | **1901** |
    | AssayBench screens resolving to an ORCS SCREEN_ID | **1901 (100%)** |
    | ORCS screens the benchmark does NOT use | **~51 (2.6%)** |

    The schema makes it obvious in hindsight: the ORCS index columns are `PHENOTYPE`, `CELL_LINE`,
    `LIBRARY_METHODOLOGY`, `LIBRARY_TYPE`, `SIGNIFICANCE_CRITERIA`, `CONDITION_NAME` — AssayBench's
    record fields verbatim. The numeric `dataset_name`s (1454, 1639, 2462...) ARE ORCS screen ids;
    the `U_`/`TR_` prefixed ones are merged/derived screens.

    **So the §4.31 and §4.32 conclusions were right about the mechanism and useless as a
    prescription.** Both located the limit at corpus SIZE. The corpus cannot be grown: the screens
    one would add are already in the training fold, or they are the evaluation folds themselves.
    This also reframes the paper's Fig. 4C — the curve not plateauing at 1349 means more data would
    help, but 1349 is the training FOLD, the whole corpus is 1901, and essentially nothing exists
    beyond it.

    **What this cost, and the process lesson.** Three sections (§4.30-§4.32) and several hours were
    spent on "get more/better screen data", with three converging lines of evidence for it — all of
    which addressed whether more screens WOULD help and none of which asked whether more screens
    EXISTED. **Check that a data source contains something you do not already have before building
    an argument for acquiring it.** A five-minute comparison of screen ids would have pre-empted all
    of it.

    **The one useful positive: ORCS ingestion is trivial if new screens ever appear.** Each
    `BIOGRID-ORCS-SCREEN_<id>-2.0.18.screen.tab.txt` carries `OFFICIAL_SYMBOL` and an explicit
    **`HIT` column (`YES`/`NO`)** over the screen's FULL library — exactly the `genes` and `hits`
    fields `build_train_hits.py` needs, with no threshold decisions. Earlier claims in this file that
    hit-calling would be "days of judgement" were wrong; the curators have done it. Hit counts per
    screen range from 1 to ~2000, so a per-screen cap or normalisation is still worth considering so
    that a permissive screen does not dominate the co-occurrence structure.

    Remaining variants of "more screens", all weak: the 51 unused human screens (2.6%); ORCS screens
    from other organisms (needs cross-species mapping); screens published after 2.0.18 (Sept 2025).


34. **CONDITION-MATCHED screen transfer: REFUTED, and the reason kills a whole family.** 838 of
    ORCS's 1952 human screens are chemogenomic and the same agents recur across cell lines
    (gemcitabine 11, cisplatin 24), so **15 of the 36 hard held-out screens have at least one
    TRAIN-fold screen testing the identical condition** — which screen-kNN never uses, since it
    picks neighbours by hit-set overlap and ignores what was tested. The feature: per-gene fraction
    of condition-matched training screens calling it a hit. `prescreen_condition_match.py`.

    | on the 11 HARD matched screens | paired vs DepMap alone |
    |---|---|
    | condition-matched prior | **+0.014 [−0.039, +0.062]** up 8/11 |
    | no same-cell-line | +0.014 [−0.040, +0.062] up 8/11 |
    | *[ctl] random training screens, same count* | ***+0.065*** up 5/11 |

    **THE PREMISE IS EMPIRICALLY FALSE, and this is the transferable finding.** AUROC on the hard
    screens is **0.492-0.507** — chance — and the prior has values for only **0.0-6.7%** of each
    screen's hits. A train-fold oxaliplatin screen shares **0.0%** of its hits with a held-out
    oxaliplatin screen; gemcitabine 3.9-6.7%. **The same drug in a different cell line yields almost
    entirely different hits.** Drug-modifier biology is CELL-LINE-specific, not drug-specific —
    which is presumably why these studies ran PANC-1 *and* BxPC-3 separately.

    **So the entire "find data about the same condition" family is dead**, including the Olivieri
    chemogenomics proposal (~30 genotoxins x genome-wide CRISPR in RPE1/HAP1): those are the same
    agents in *different* cell lines, so the same 0-7% transfer applies. Do not pursue
    condition/drug/pathogen lookup in any form.

    **Two controls that changed the conclusion.**
    - **Across all 21 matched screens the paired contrast reads +0.335 [−0.051, +0.852], and
      dropping same-cell-line matches takes it to −0.037.** The apparent gain was reading near-twin
      screens: 2309 and 2317 (AZD6738) show AUROC 0.92 with 100% hit coverage, which is the same
      assay in the same line, not transfer. Exactly the §4.22 channel, firing where predicted.
      **Without this control the feature would have been reported as a large win.**
    - **Fame bias +0.595** — the worst of anything screened (PRISM ~0.00, PPI +0.203). A gene that
      is a hit across many conditioned screens is just a frequently-hit gene, which is what `v_g`
      already encodes.


35. **DEPMAP LEAKAGE AUDIT: no leakage, tested observationally AND interventionally.** The channel
    is real in principle and had never been checked. DepMap Achilles is 808 CELL LINES x 18119
    genes — the fitness effect of each knockout in each line — so for a knockout FITNESS screen in
    line L, DepMap's row for L is that screen's answer key in the same assay type. The exploit path:
    the diffmeans direction is built from the confirmed hits' 808-dim profiles, hits essential in L
    all carry strongly negative values in the L coordinate, so `q` acquires a large negative L
    component and `corr(profile(g), q)` then favours genes essential in L. **43 of 100
    non-degenerate held-out screens resolve to an Achilles line**, so the exposure is substantial.

    **Observational split — having your line in DepMap confers no advantage:**

    | split | deployed dEF |
    |---|---|
    | line IN DepMap (any readout) | +0.459 [+0.137, +0.914] |
    | line NOT in DepMap | +0.499 [+0.289, +0.728] |
    | in DepMap & fitness readout & KO | +0.768 [+0.290, +1.463] |
    | NOT in DepMap & fitness readout | +0.726 [+0.402, +1.087] |
    | **hard screens, line IN DepMap** | **+0.049 [−0.010, +0.121]** |
    | **hard screens, line NOT in DepMap** | **+0.223 [+0.062, +0.401]** |

    **Interventional ablation — the decisive test.** For each of the 35 screens whose line is in
    Achilles, recompute the feature with THAT LINE'S COLUMN DELETED from every gene's profile:

    | | oracle margin |
    |---|---|
    | with own line | +0.310 [+0.089, +0.550] |
    | **own line removed** | **+0.320 [+0.104, +0.558]** |
    | **paired effect of removal** | **+0.010 [+0.001, +0.024]** |

    Removing it is very slightly BETTER. A reading-the-column feature would collapse. Why the
    exploit fails: one coordinate of 808, after row-centring and unit-normalisation, contributes
    almost nothing to the direction, and co-essentiality is dominated by structure shared across
    many lines rather than any single line.

    **What actually drives the benefit is the READOUT TYPE, not the cell line:** fitness-readout
    screens +0.766, non-fitness +0.302, regardless of whether the line is in DepMap. That sharpens
    §4.30 — the easy/hard split is about WHAT THE SCREEN MEASURES, not which line it used, and the
    mechanism generalises rather than being a lookup.

    **All three data sources now audited.** DepMap: this section. Training labels (screen-kNN,
    co-hit): §4.22, gain strongest on the LEAST contaminated screens. Text embedding: OpenAI
    `text-embedding-3-small` over the description, no label contact. **Quote this when asked — it is
    the first question a reviewer will have about using public atlases at inference time.**


36. **LINCS L1000 OVEREXPRESSION: killed on COVERAGE for 11 MB, before any matrix download.** This
    was the named candidate for the gain-of-function gap across §4.30 and §4.36 drafts — the 8
    CRISPRa screens hold **72% of hard-screen hits** (§4.30 inventory) and LINCS OE is the only
    gain-of-function perturbation atlas at scale. It fails §4.27's gate 2.

    Cheap route to the answer: `GSE92742_Broad_LINCS_sig_info.txt.gz` from GEO (11 MB) lists every
    CMap signature with `pert_type` and `pert_iname`, so the overexpressed-gene set is readable
    without the ~2 GB level-5 matrix. (The LINCS2020 beta S3 path works for
    `siginfo_beta.txt` (444 MB) and `geneinfo_beta.txt` (1.1 MB) but the bucket refuses listing and
    the level-5 filename encodes a signature count, so guessing it 403s.)

    | atlas | universe coverage | coverage of each screen's hits |
    |---|---|---|
    | DepMap (the one that works) | 83% | — |
    | PRISM | — | 70-95% |
    | BioPlex PPI | 61% | 77% |
    | Perturb-seq | 46% | — |
    | **LINCS OE** | **13.5% (3114 genes)** | **16-20%** |

    22205 `trt_oe` signatures cover only **3114 distinct genes**. Per CRISPRa screen: 1639 17.7%,
    1640 16.0%, 1641 16.9%, 1642 16.8%, 1643 16.7%, 1645 18.6%, 1646 19.6%. Under the neutrality
    rule 80-84% of the hits that matter take the mean contribution, i.e. no information; even with
    perfect signal on the covered sixth the feature can reorder only a sixth of the relevant genes.

    **A second, independent reason.** LINCS OE cell lines are PC3, A375, HCC515, HT29, MCF7, HA1E,
    HEPG2, VCAP — **no pancreatic line**, while the CRISPRa screens are PANC-1 and BxPC-3. §4.34
    measured 0-7% hit transfer between two *pancreatic adenocarcinoma* lines, so transfer from
    prostate/melanoma/breast would be worse.

    **Process note.** This candidate was recommended three times across §4.30 and the chat before
    anyone checked its coverage, which took 11 MB and five minutes. **Check gate 2 (does the atlas
    have values for the hits you need to rank?) BEFORE recommending a feature, not after** — it is
    the cheapest of the three gates and it is structural rather than probabilistic: you cannot rank
    genes you have no values for.

    Residual: JUMP Cell Painting (~15k ORF + ~12.6k CRISPR morphological profiles, AWS open data) is
    the remaining gain-of-function option with plausibly better gene coverage. Its cell context is
    U2OS/A549, still not pancreatic, so §4.34's background-specificity result applies. Check gate 2
    first.


37. **PANCREATIC-SPECIFIC DATA: lineage restriction is WORSE, and this is the third instance of a
    general pattern.** 15 of the 36 hard screens are PANC-1 / BxPC-3 holding 83% of hard-screen hits
    (§4.30), so lineage-specific data is the obvious ask. **Both lines are already in Achilles**,
    along with 32 other pancreatic lines (34 of 808). Restricting the profile matrix to those 34:

    | screen group | all 808 lines | 34 PANCREATIC | *[ctl] 34 random* |
    |---|---|---|---|
    | all screens | **+0.364 [+0.134, +0.608]** | −0.040 | *−0.067* |
    | pancreatic-line screens | **+0.042** | +0.004 | *−0.058* |
    | hard screens | **+0.179** | +0.061 | *−0.005* |

    **The lineage signal is real but tiny.** Pancreatic-34 consistently beats random-34 (+0.004 vs
    −0.058 on pancreatic screens; +0.061 vs −0.005 on hard), so matching the lineage does contribute
    above chance. But both collapse against the full matrix: **the specificity gained is swamped by
    the statistical power lost.**

    **THE PATTERN, now three for three — filtering data to the matched subset always loses:**

    | restriction | result |
    |---|---|
    | §4.31 enrichment-only training screens | −0.061 vs pooled +0.222 |
    | §4.32 PPMI / rank-truncated co-hit | +0.158 vs raw full-rank +0.227 |
    | §4.37 pancreatic-only cell lines | −0.040 vs all-808 +0.364 |

    Use everything. This is a real finding about the problem class and the reverse of the intuition
    that keeps being brought to it — specificity arguments sound compelling and lose to sample size
    every time here.

    **Other pancreatic datasets, each already answered:** more pancreatic CRISPR screens — ORCS
    exhausted (§4.33) and 0-7% hit transfer between PANC-1 and BxPC-3 on the SAME drug (§4.34);
    pancreatic organoid / PDX pharmacogenomics — that is drug-response profiling, i.e. PRISM, AUROC
    0.500 (§4.27); TCGA-PAAD / CPTAC-PDAC — expression and proteomics, harmful as a kernel at −0.109
    (§4.24). And §4.35 already measured that cell-line matching is not load-bearing for this feature
    (+0.459 for in-DepMap lines against +0.499 for absent ones), so lineage-specific data targets an
    axis shown not to matter.


38. **JUMP CELL PAINTING ORF: the first gain-of-function candidate to PASS gate 2. Gates 1 and 3
    still open — this is the one live lead.** §4.36 killed LINCS L1000 OE on coverage (3114 genes,
    16-20% of the CRISPRa screens' hits). JUMP is the alternative gain-of-function atlas:
    ~15k ORF overexpression and ~8k CRISPR perturbations profiled by imaging, free on AWS open data.
    Metadata is a 550 KB CSV on GitHub (`jump-cellpainting/datasets`, `metadata/orf.csv.gz`), so
    gate 2 costs one download.

    | atlas | universe coverage | coverage of the CRISPRa screens' hits |
    |---|---|---|
    | LINCS OE (killed, §4.36) | 13.5% | **16-20%** |
    | **JUMP ORF** | **56.1% (12603 genes)** | **63-66%** |
    | JUMP CRISPR | 36.8% (7976 genes) | 31-86% |
    | *DepMap, the one that works* | *83%* | — |
    | *BioPlex, tested* | *61%* | *77%* |

    **Three to four times LINCS's coverage, on exactly the 8 screens holding 72% of hard-screen
    hits** (§4.30). JUMP CRISPR is weaker and redundant with DepMap anyway (same modality), though
    note it covers 81-86% of the Calu-3 CRISPRi screens' hits.

    **What the kernel would mean:** two genes are similar if OVEREXPRESSING them produces the same
    change in cell morphology — a functional similarity in the gain-of-function direction, which is
    the §4.30 gap.

    **The two conditions it still violates** (§4's synthesis of why DepMap works):
    - **Readout**: morphology, not fitness. §4.35 measured benefit tracking readout match (fitness
      +0.766 vs non-fitness +0.302) — though morphology is arguably BETTER matched than fitness for
      the Calu-3 surface-marker and protein-binding screens.
    - **Contexts**: one or two cell lines against DepMap's 808. This is the condition whose
      violation destroyed pancreatic-34 (+0.364 -> −0.040, §4.37) and weakened Perturb-seq (§4.17),
      and it is the likeliest killer here.

    **Next step if pursued:** download the ORF profiles from the Cellpainting Gallery (a few hundred
    MB: ~15k perturbations x a few thousand CellProfiler features), build the `(dims x genes)` npz,
    then gate 1 (fame bias) and gate 3 (oracle margin split by difficulty). Everything needed is in
    `prescreen_ppi.py`'s pattern.

    **CLARIFICATION — "many contexts matters" vs "matched context does not" are NOT in tension.**
    These were stated loosely enough to look contradictory. The cross-cell-line dimension does two
    different jobs:
    - **As a SAMPLE axis** — more lines gives a more precise estimate of gene-gene co-dependency.
      **Large effect: 808 -> 34 lines costs −0.40 EF** (+0.364 -> −0.040, §4.37).
    - **As a CONTEXT axis** — matching the lineage adds specificity. **Small effect: pancreatic-34
      beat random-34 by only +0.03 to +0.07** (§4.37), and deleting a screen's OWN line changes the
      margin by +0.010 (§4.35).

    So **more lines of ANY kind beats fewer lines of the RIGHT kind, by roughly 6-13x.** Earlier
    phrasing in §4.35/§4.37 that cell-line matching "does not matter" or "is not load-bearing"
    OVERSTATED it — pancreatic-34 consistently beat random-34, so the lineage effect is real and
    second-order, not absent.

    Consequence for few-line atlases: **they are handicapped for the ESTIMATION reason, not because
    their cell types are wrong.** And the sharpened prediction from §4's synthesis: a
    gain-of-function atlas with ~800 contexts should work on the CRISPRa screens *regardless of
    whether any line is pancreatic*. Nobody has built one, so it stays untested.

    **Sack et al. 2018** (Cell 173:499-514) genome-scale ORF PROLIFERATION screens — (Cell 173:499-514) genome-scale ORF
conceptually the best candidate, since it is
    gain-of-function AND a fitness readout, the one combination that would directly test §4's
    synthesis. **GSE109326 is NOT the screen data** — it is 12 RNA-seq samples ("Transcriptional
    changes after overexpression of proliferation drivers in human mammary epithelial cells"), the
    transcriptomic follow-up. The ORF proliferation SCORES are in the Cell supplementary tables
    (`mmc*.xlsx`), which need the journal page. **Doubly handicapped if pursued:** few cell lines
    (the estimation problem above) and, by the authors' own headline finding, proliferation drivers
    act "in an unexpectedly highly tissue-specific manner" — the same non-transfer §4.34 measured at
    0-7% between two pancreatic lines.


39. **JUMP ORF: FAILS, and the POSITIVE CONTROL explains why — gain-of-function may not support
    co-functional clustering at all.** The §4.38 lead, built and gated. `jump_to_npz.py` +
    `prescreen_atlas.py` (a generic three-gate runner, reusable for any atlas).
    Atlas: 790 morphology features x 12603 genes, from
    `cellpainting-gallery.s3.amazonaws.com/cpg0016-jump-assembled/source_all/workspace/
    profiles_assembled/ORF/v1.0b/profiles_wellpos_var_mad_int_featselect_harmony.parquet` (396 MB,
    81660 wells, median 5 wells per gene, aggregated by median).

    | gate | result |
    |---|---|
    | 1. fame bias | **PASS** ~0.000 |
    | 2. coverage of hits | **PASS** 66.1% median, 63.6% on hard screens |
    | **0. positive control: do co-complex genes cluster?** | **FAIL** |
    | 3. paired vs DepMap, hard screens | **−0.019 [−0.093, +0.060]** (gene-centred), −0.030 (dim-centred) |

    **THE STANDALONE DIFFICULTY SPLIT, recovered from `output/analysis/jump_gates.log`
    on 2026-10-23 — it was measured on all 83 screens at the time and never written
    down here, which made the pair look less settled than it is.** This is the row that
    answers "is there value in using JUMP and DepMap together":

    | kernel alone, oracle margin per family | ALL | hard (28) | mid (19) | easy (36) |
    |---|---|---|---|---|
    | DepMap | +0.358 | +0.147 | +0.413 | +0.542 |
    | **JUMP ORF** [gene-centred] | **−0.095** | +0.047 [−0.007, +0.102] | **−0.269** | **−0.062** |
    | JUMP ORF [dim-centred] | −0.107 | +0.046 | −0.313 | −0.073 |
    | *BioPlex for contrast (§4.28)* | *+0.289* | *+0.152* | *+0.219* | *+0.525* |

    **This is the whole argument, and it is a clean one.** BioPlex's overall pair gain of
    +0.190 came from the MID and EASY screens, where it is strongly positive standalone.
    JUMP is **negative** in both of those cells, so it cannot repeat that. Its only
    positive cell is the hard screens — and that is precisely where the paired contrast
    WAS run, giving −0.019 (13/28 up). So the one region where JUMP carries positive
    standalone signal is the one region where the pair has been tested directly, and it
    came out at zero. **The pair is closed on measurement, not on inference.**

    **GATE 0 IS THE RESULT.** Co-complex genes are no more similar than random pairs in JUMP ORF
    space, in BOTH centring modes:

    | complex | JUMP ORF | DepMap | Perturb-seq |
    |---|---|---|---|
    | ribosome | **−0.1 sd** | +4.0 sd | +9.3 sd |
    | RNA pol II | **−0.1 sd** | +10.6 sd | +4.9 sd |
    | proteasome | **+0.0 sd** | +6.1 sd | +5.2 sd |
    | mediator | +0.4 sd | +7.1 sd | +11.3 sd |
    | SF3B | **−1.5 sd** | +4.9 sd | +0.4 sd |

    Gates 1 and 2 passing is what makes this a real null rather than a broken instrument — the
    reason gate 0 exists at all (§4.17's 6832 infinities produced an identical-looking null).

    **THE DEEPER FINDING, and it is biologically sensible.** Knocking out ANY ribosomal subunit kills
    the cell the same way, which is why complexes cluster in DepMap. **Overexpressing RPL3 and RPL4
    need not do anything similar** — gain-of-function phenotypes depend on stoichiometry,
    localisation and dosage, which differ per gene, while loss-of-function phenotypes converge on
    shared function. So **gain-of-function perturbation may not support co-functional clustering at
    all**, and that is the property that makes a loss-of-function atlas work as a kernel. This is a
    deeper obstacle than "wrong readout" or "too few contexts": it would apply to ANY GoF atlas.

    **It weakens §4's falsifiable prediction without settling it.** The synthesis predicted a GoF
    FITNESS atlas should work on the CRISPRa screens. This is evidence against from a GoF MORPHOLOGY
    atlas, and the two differ — fitness collapses everything onto one axis where convergent effects
    would show, morphology may not. **Sack et al. is therefore now MORE worth testing, not less:**
    it is the only experiment that separates "GoF does not cluster" from "morphology does not
    cluster". Its scores are in the Cell supplement (`mmc*.xlsx`); GSE109326 is only the 12-sample
    RNA-seq arm.

    **Implementation caveat:** `prescreen_atlas.py`'s gate 1 computes
    `corr(profile norm, hit frequency)`, not `corr(feature, hit frequency)` as its docstring says,
    and carries a dead `col = ...` line. A reasonable proxy that passed either way, and gate 0 is
    decisive here, but the number is NOT comparable to the gate-1 figures quoted for PRISM (~0.00),
    BioPlex (+0.203) or condition-matching (+0.595). Fix before relying on it.

    **CONTROL-SELECTION CAVEAT (2026-10-23, from §12 item 1) — this may be the weakest joint in the
    whole negative.** The mechanism this section states to EXPLAIN the gate-0 failure also predicts
    it, and therefore disqualifies the control: if overexpression phenotypes depend on stoichiometry,
    then overexpressing one subunit of a complex creates an IMBALANCE, and co-complex subunits should
    diverge rather than converge. Failing a co-complex control on a gain-of-function atlas is then
    the expected result whether or not the atlas carries functional signal, so it shows the atlas
    lacks **complex-membership** signal — not that it lacks function. The conclusion drawn here
    ("gain-of-function may not support co-functional clustering at all") does not follow from this
    control. Re-run gate 0 on JUMP ORF with a GoF-appropriate positive set — paralogue pairs, shared
    Pfam domain families, known same-pathway oncogenes — before treating the null as structural. The
    npz is already built, so this costs hours and no new data. §12 item 1.
    **RESOLVED 2026-10-23 — see §4.41.** The criticism was right: the control was
    mis-specified AND underpowered (six complexes, ~75 pairs), and this section's
    headline claim that co-complex genes are "no more similar than random pairs" is
    false at 254 complexes (AUROC 0.541 [0.525, 0.559]). The conclusion that
    gain-of-function "may not support co-functional clustering at all" is **withdrawn**.
    The build decision does not change: JUMP's best ground truth is 4.5x weaker than
    DepMap's, matching this section's own gate-3 null of −0.019.


40. **SACK et al. 2018 closed WITHOUT the data: it screened THREE cell types, so its profile space
    is 3-dimensional.** This was the last named candidate and the one experiment that would have
    separated "gain-of-function does not cluster" from "morphology does not cluster" (§4.39).

    **Data access, for the record — all three routes fail.** The paper points to
    `elledgelab.med.harvard.edu`, which redirects to `elledge.hms.harvard.edu`; its "Supplemental
    Data" section hosts Davoli 2013, Zhang 2023 and Bruno 2023 only, with **no Sack entry and no
    mention of it anywhere on the site**. GSE109326 is only the 12-sample RNA-seq arm (§4.38). PMC
    (PMC6643283) *does* hold all seven supplementary tables — `efetch db=pmc retmode=xml` lists them
    as `NIHMS1015078-supplement-Table_S1..S7.xlsx` — but **PMC serves supplementary downloads behind
    a reCAPTCHA challenge**, so they need a browser. Direct link pattern:
    `https://pmc.ncbi.nlm.nih.gov/articles/PMC6643283/bin/NIHMS1015078-supplement-Table_S2.xlsx`
    (S1 library design, **S2 the screen results**, S3 pathways, S4 KRTAP, S5 SCNA, S6 aneuploidy,
    S7 amplified GO genes).

    **The article metadata settles it anyway.** Figure 2's caption reads "Genome Scale Proliferation
    Screens in **Three** Human Cell Types" — HMEC (mammary epithelial), **HPNE (pancreatic)** and
    IMR90 (fibroblast). HPNE is genuinely relevant to the PANC-1/BxPC-3 screens, but three cell
    types means **the profile dimension is 3**:

    | atlas | profile dims |
    |---|---|
    | BioPlex | 14551 |
    | Perturb-seq | 8248 |
    | DepMap | 808 |
    | JUMP ORF | 790 |
    | *AssayFormer's own `v_g`* | *10* |
    | **Sack** | **3** |

    **A 10,000-gene similarity space cannot be built from three measurements per gene.** Not noise —
    structural: every gene collapses onto a few directions, and the space is SMALLER than `v_g`'s
    rank 10, which §4.9 measured as already inside the model. As a feature rather than a kernel it
    would be 1-3 gene-only columns, which is the fixed-prior gene-only category that has failed six
    times (§4.27). **Earlier drafts of §4.38 guessed "4-10 cell lines" and called the dimensionality
    argument probabilistic; at 3 dimensions it is decisive.**

    **Residual value, descriptive not predictive.** Table S2 holds per-gene results — "391 GO-like
    drivers of proliferation and 1,269 STOP-like drivers with FDR < 0.05". Worth one check for the
    write-up: **are the CRISPRa pancreatic resistance hits enriched for Sack's HPNE GO genes?** Two
    independent gain-of-function screens in pancreatic context agreeing would corroborate the hit
    calls; disagreement would speak to §4.34's background-specificity result. It cannot improve EF.

    **DESCRIPTIVE CHECK RUN 2026-10-23 (user supplied the tables in `Sack_dataset/`): the two
    independent gain-of-function screens DO NOT AGREE.** `sack_crispra_overlap.py`, Fisher exact of
    hits vs non-hits within each screen's own library, restricted to the 14684 genes Sack measured
    so neither library composition nor his coverage can drive it. Sheet S2D carries per-gene
    enrichment/dropout stats for HMEC and HPNE (GO = overexpression promotes proliferation, STOP =
    suppresses); at FDR < 0.05 that is **HMEC GO 268, HMEC STOP 1023, HPNE GO 77, HPNE STOP 454**.

    | modality | HMEC GO | HMEC STOP | **HPNE GO** | HPNE STOP |
    |---|---|---|---|---|
    | **Activation (CRISPRa, n=7)** | 1.03 | 1.12 | **0.80** | 1.23 |
    | Knockout (n=8) | 0.81 | 0.87 | 0.82 | 0.87 |

    **No corroboration.** The CRISPRa pancreatic resistance hits are NOT enriched for HPNE
    growth-promoting genes — mean OR **0.80**, slightly depleted, with one screen significantly
    depleted (1640, OR 0.45, p=0.03) and none enriched. **A third independent line of evidence for
    background specificity**, alongside §4.34's 0-7% transfer between PANC-1 and BxPC-3 on the same
    drug: HPNE is immortalised NORMAL pancreatic epithelium, PANC-1/BxPC-3 are adenocarcinoma, and
    that difference suffices.

    Two caveats. **HPNE GO is only 77 genes**, so the odds ratios are unstable — weak evidence of
    absence, not strong. And the one suggestive pattern, CRISPRa hits leaning toward HPNE STOP
    (OR 1.23, three of seven screens significant), is **at the chance level**: 15 screens x 4 gene
    sets = 60 tests, 5 significant at p < 0.05 against ~3 expected. Biologically attractive —
    overexpressing a growth suppressor slows cycling and slow-cycling cells survive gemcitabine,
    5-FU and irinotecan, a known resistance mechanism — but not a finding on this evidence. The
    clean version would pre-specify HPNE STOP on the CRISPRa screens as one hypothesis and test it
    on the `test` fold.

    **This closes the candidate list.** §4's falsifiable prediction — that a gain-of-function atlas
    with many contexts would work on the CRISPRa screens — remains untested, because no such atlas
    exists: JUMP has the contexts but gain-of-function does not cluster in it (§4.39), and Sack has
    the right readout with three contexts.


**Two earlier results worth not re-deriving:**

1. **Coverage overlap predicts hits, it does not predict waste.** Hit rate by the fraction of a
   category's genes already inside AssayFormer's top picks: <=5% -> 4.2%, 5-20% -> 10.3%,
   20-50% -> **15.6%**. The novelty reweighting in `LLMCategoryComplementAcquisition` therefore
   shifts slots *toward the LLM's weakest categories*. Either disable it
   (`--no-novelty-reweight`) or invert it to favour moderate overlap.
2. **The filter trades hit rate for novelty at a fixed number of new hits** (~4-5/screen across
   every horizon from 0 to 1000). It removes redundant hits rather than creating new ones.

**Filter dose-response** (validation, 20 LLM genes/batch, replay mode, free to re-run):

| Excluded AF top-N | EF | Gemini hit rate | New hits/screen | New share |
|---|---|---|---|---|
| none | 6.01 | 17.7% | 3.8 | 39% |
| 100 | 6.03 | 17.0% | 4.1 | 43% |
| 300 | 5.99 | 15.1% | 4.4 | 53% |
| 500 | 5.98 | 13.4% | 4.8 | 65% |
| 800-1000 (default) | 5.96 | 9.7% | 4.6 | 84% |

---

41. **GATE 0 RE-RUN WITH A GAIN-OF-FUNCTION GROUND TRUTH: §4.39's control WAS
    mis-specified, its stated conclusion is too strong, and the atlas still does not
    reopen.** §12 item 1, built and run 2026-10-23. `build_gate0_sets.py` +
    `gate0_gof.py`. Free, no new atlas: the JUMP npz was already on disk.

    **The criticism that prompted it is correct.** §4.39 failed JUMP on co-complex
    clustering and then explained the failure by stoichiometry — overexpressing one
    subunit creates an imbalance. That mechanism *predicts* the failure, so the control
    cannot separate "no functional signal" from "wrong ground truth". Co-complex
    membership is the right control for a KNOCKOUT atlas, where removing any subunit
    disables the whole machine, and the wrong one for an overexpression atlas.

    **The fix: three ground truths whose predictions differ by perturbation direction.**
    Built from HGNC `gene_group` (curated families given as symbols, so no id mapping),
    split by name because that column mixes families with complexes, plus MSigDB C2 CP
    for pathways. 301 complexes / 1402 paralogue families / 2323 pathways at size 2-12,
    against the six hand-picked complexes §4.39 used. Paralogues are the discriminating
    set: they are substitutes rather than partners, so overexpressing either should do a
    similar thing, while under knockout they BUFFER each other and their dependency
    profiles are routinely uncorrelated. Statistic is AUROC of within-set pair cosine
    against 200k random covered pairs — scale-free, so it compares across atlases whose
    cosine widths differ — bootstrapped over FAMILIES, not pairs.

    | ground truth | JUMP ORF (overexpression) | DepMap (knockout) |
    |---|---|---|
    | complexes | 0.541 [0.525, 0.559] | **0.750 [0.728, 0.773]** |
    | **paralogue families** | **0.555 [0.547, 0.566]** | 0.590 [0.581, 0.601] |
    | pathways | 0.519 [0.516, 0.522] | 0.551 [0.547, 0.555] |
    | **paralog − complex** | **+0.014 [−0.003, +0.033]** P(>0)=0.94 | **−0.160 [−0.187, −0.134]** P(>0)=0.00 |

    **FINDING 1: the dissociation is real and it reverses sign between the atlases.**
    Each atlas scores highest on the ground truth that matches its perturbation
    direction. In DepMap complexes beat paralogues by 0.160 AUROC, decisively. In JUMP
    the ordering flips — paralogues beat complexes — though JUMP's own margin is
    marginal (CI touches zero at P(>0)=0.94). The sign reversal itself is not marginal:
    the two intervals do not come close to overlapping. Identical in both centring
    modes, so it is not a centring artefact.

    **FINDING 2: §4.39's "no more similar than random pairs" is WRONG, and partly for a
    power reason.** All three JUMP groups are significantly above chance — complexes
    0.541 [0.525, 0.559], paralogues 0.555, pathways 0.519 on 151k pairs. §4.39 read
    −0.1 to +0.4 sd on **six** complexes (~75 pairs); at 254 complexes and 1967 pairs
    the same quantity is clearly positive. The original gate 0 was underpowered as well
    as mis-specified. **Fix in `prescreen_atlas.py`: six hand-picked complexes cannot
    support a gate verdict. Use `gate0_gof.py`'s sets.**

    **FINDING 3: and yet the atlas does not reopen, for a reason the repaired control
    makes quantitative rather than disputed.** Excess over chance: JUMP's BEST ground
    truth is **+0.055**, DepMap's best is **+0.250** — **4.5x weaker**. JUMP at its best
    (0.555) barely matches DepMap at its WORST (pathways, 0.551). §4.39's gate 3 already
    measured the consequence directly: paired against DepMap on the hard screens,
    **−0.019 [−0.093, +0.060]**. The null stands; what changes is that it is now
    explained by a measured 4.5x signal deficit instead of by a disqualified control.

    **The concrete illustration, from the spotlight families.** Overexpressing Src-family
    tyrosine kinases produces a strongly convergent morphology (**+2.26 sd** in JUMP)
    while knocking them out does not produce convergent dependency (+0.10 sd in DepMap);
    E2F transcription factors likewise (**+2.47 sd** JUMP, −0.24 sd DepMap). The Jak
    family runs the other way (+0.21 JUMP, +2.04 DepMap). These are not noise — they are
    the dissociation at the level of individual families, and they are why the atlas is
    weak on average rather than empty.

    **What this does to the §4.39 claim.** "Gain-of-function may not support
    co-functional clustering at all" must be **withdrawn** — it does, weakly, and most
    on the ground truth suited to it. The defensible version is narrower and still
    useful: **a gain-of-function morphology atlas carries functional similarity at about
    a fifth the strength a loss-of-function fitness atlas carries it, and the ground
    truth it carries is paralogy rather than complex membership.** That is a finding
    about perturbation modality, not a dead end, and it is the honest form of the §4
    synthesis's third condition (perturbation direction must match).

    **Caveat that this instrument cannot settle.** JUMP has a median of 5 wells per gene
    against DepMap's 808 cell lines, so some of the 4.5x deficit is measurement noise
    rather than biology. Separating them needs a GoF atlas with comparable depth, which
    §4.40 established does not exist (Sack screened three cell types). The deficit is an
    upper bound on how much signal is there, which is the direction that matters for a
    build/no-build decision.

    **FOLLOW-UP, asked immediately and worth recording: JUMP and DepMap are as
    ORTHOGONAL as two atlases get, and it does not help.** `atlas_complement.py`, free,
    on the 11392 genes both atlases cover. Pearson correlation between the two atlases'
    cosines on the SAME pairs is **+0.009** (complex), **+0.010** (paralog), +0.056
    (pathway), +0.017 (random) — against the +0.113 to +0.119 §4.18 measured among
    DepMap, Perturb-seq and screen-kNN. And JUMP keeps essentially all of its signal on
    the pairs DepMap cannot see (cosine at or below the median random pair):

    | ground truth | % of pairs DepMap is blind to | JUMP AUROC, all pairs | JUMP AUROC given DepMap blind |
    |---|---|---|---|
    | paralogue families | 38.0% | 0.557 | **0.552** |
    | complexes | 21.3% | 0.544 | 0.525 |
    | pathways | 43.9% | 0.519 | 0.515 |

    Symmetrically DepMap loses nothing where JUMP is blind (complexes 0.758 -> 0.752).
    So the fusion STRUCTURE is as good as it could be: the two spaces are independent,
    not redundant, and 38% of paralogue relationships are invisible to the deployed
    atlas while JUMP reads them at 0.552. **It still does not matter, and this is the
    sharpest statement of §4.18's lesson in the file: near-perfect orthogonality times
    near-zero signal is near-zero.** JUMP's conditional excess over chance is +0.052
    against DepMap's +0.250, and Perturb-seq — stronger than JUMP and already
    near-orthogonal at r = 0.11 — delivered +0.044/+0.066, below the fold's resolution.
    Orthogonality was never the binding constraint; amplitude is.

    Code: `build_gate0_sets.py` (HGNC + MSigDB -> `output/analysis/gate0_sets.json`;
    documents its two hand exclusions), `gate0_gof.py` (`--atlas LABEL=PATH` repeatable,
    family bootstrap, prints the dissociation contrast), `atlas_complement.py`
    (conditional AUROC under blinding, for any atlas pair). HGNC lives at
    `/home/jeff/hgnc/hgnc_complete_set.txt` (17 MB). Perturb-seq could not be included:
    no npz is on disk any more.

---

42. **RESOLUTION WORK (§12 item 4): the gain is real but SMALL — about 7%, and nAUC
    makes things WORSE. Plus a calibration observation that materially raises the
    expected value of DepMap+BioPlex.** `score_resolution.py`, free, re-reading the
    deployed campaign already on disk. Instrument validated first: the per-round curve
    reproduces `_adj_ef_from_run` on all 198 screens, and the unadjusted figure
    reproduces the recorded headline (+0.175 per screen [+0.077, +0.290], +0.331 per
    family).

    **A gotcha worth knowing: `_adj_ef_from_run`'s `budget` argument does NOT truncate
    the pick list.** It counts every pick in the run and only rescales the denominator,
    so "EF at budget 500" from a 1000-pick run is not EF of the first 500 picks. A real
    budget curve has to be built from pick PREFIXES; `ef_curve()` does that, one reading
    per round boundary.

    **SNR is the statistic to compare, not interval width** — a narrower interval around
    a diluted point estimate is not an improvement. SNR here is `mean / (CI width/3.92)`.

    | metric | unit | mean | unadjusted SNR | + ANCOVA (strict) | + ANCOVA (full) |
    |---|---|---|---|---|---|
    | dEF at budget 1000 | screens | +0.175 | **3.22** | 3.41 (−6% width) | 3.36 |
    | dEF at budget 1000 | families | +0.331 | **3.85** | 3.94 | 3.92 |
    | dnAUC over 10 rounds | screens | +0.144 | **3.00** | 3.08 | 3.19 |
    | dnAUC over 10 rounds | families | +0.299 | 3.81 | 4.03 | **4.13** |

    **FINDING 1: nAUC alone is WORSE, and the per-round curve says exactly why.** The
    deployed effect ACCUMULATES: dEF by round boundary runs +0.000, +0.125, +0.184,
    +0.293, +0.157, +0.117, +0.132, +0.113, +0.143, **+0.175**. Round 1 is identically
    zero — `min_seeds = 2` means the kernel cannot fire before any labels exist — so
    averaging over rounds dilutes the signal faster than it averages down the noise.
    SNR falls 3.22 -> 3.00 per screen. **nAUC is the wrong summary for an intervention
    whose effect grows over rounds**, which is every seed-driven kernel here. Retract
    that half of §12 item 4.

    **FINDING 2: ANCOVA works, modestly.** Covariates explain only R² 0.10 (strict,
    screen-level) to 0.23 (full, family-level) of the delta's variance, and SE scales as
    sqrt(1−R²), so the ceiling was always a few percent. Measured: 2-8% narrower
    intervals. The best cell is dnAUC + full covariates at family level, SNR 4.13 against
    the 3.85 baseline — **+7%, equivalent to about 15% more screens.** Worth keeping
    because it is free, not because it changes a verdict.

    **Implementation note, because the first version silently did nothing.** The
    regression estimator must centre covariates at the FIXED full-sample mean:
    `est = y_rep.mean() + beta' (x_full − x_rep.mean())`. Re-centring Z inside each
    bootstrap replicate makes the OLS intercept identically the replicate mean, the
    adjustment cancels, and every covariate set reports "0% width change" — which is
    what the first run printed, and it looks exactly like a true null.

    **FINDING 3, the one that matters most, and it is a CALIBRATION observation.** §4.28
    was rejected partly on my estimate that +0.190 offline maps to ~+0.06 deployed, using
    the ~3x shrinkage of §1. **That shrinkage factor is the wrong one to apply**: it was
    measured going from the CURATED 20-screen sets to the fold, not from an offline
    oracle margin to a deployed campaign on the SAME fold. The clean calibration point
    for the latter is DepMap itself:

    | DepMap, same 81-screen held-out fold | value |
    |---|---|
    | offline oracle margin, per family (`jump_gates.log`) | **+0.358** |
    | DEPLOYED, per family, 198 screens (§4.19) | **+0.331** |

    **Ratio 0.92 — the offline oracle instrument on this fold is nearly unbiased for the
    deployed per-family figure.** The catastrophic shrinkages on record (co-hit +0.198 ->
    +0.006; expression +0.87 -> −0.109) were both measured on the CURATED 19-20 screen
    sets, which §4.18 separately showed can miss a real effect of this size. If the
    81-screen instrument is calibrated like DepMap's, BioPlex's **+0.190 per family maps
    to roughly +0.17 deployed, not +0.06** — and against a family-level SE of ~0.072 on
    198 screens that is an SNR near 2.4, i.e. **measurable**. This does not prove the
    effect is real; it says the magnitude argument that rejected it was built on the
    wrong shrinkage constant, and the one applicable calibration point is favourable.
    **Treat RULE 1 as scoped: it is established for the curated-set instrument, and has
    exactly one measurement on the 81-screen one, which came out 0.92.**

    **REFUTED SAME DAY by §4.43 — do not use this argument.** BioPlex was deployed on
    the strength of it and came in at **−0.007 per family [−0.031, +0.016]** against a
    prediction of ~+0.17. The 0.92 ratio was a property of DEPMAP, not of the 81-screen
    instrument. RULE 1 stands unscoped, now with seven instances including one from the
    well-powered instrument. §4.43 has the mechanism: the oracle margin measures whether
    a column CONTAINS information, while a deployed campaign measures whether that
    information sits where the base ranker is wrong, and only the second is worth
    anything.

    Code: `score_resolution.py` (`--baseline`, `--arm LABEL=sweep-id` repeatable,
    `--min-library`; prints the per-round curve, both metrics, both units, three
    covariate sets, SNR and R²).

---

43. **DEPLOYED DepMap + BioPlex + screen-kNN: a TIGHT, WELL-POWERED NULL. §4.28's
    verdict was right, both of its stated reasons were wrong, and the calibration
    argument I built in §4.42 to reopen it is REFUTED.** §12 revised item 2, built and
    run 2026-10-23. 3 arms x 218 screens, 0 failures, all on CPU. The first result in
    this project that closes a candidate on a MEASUREMENT rather than on
    unmeasurability.

    **What was built.** `bioplex_to_npz.py` + `codep_seed_acq.py`'s `ppi_npz` /
    `ppi_weight` / `ppi_degree_residual`. §4.28 said this needed "a new sparse-profile
    path -- half a day" because "the dense profile form is ~2.5 GB". The dense form is
    never needed and the identity is exact: `q` is built from centred rows, every
    centred row sums to zero, so `sum(q) = 0` and

        (A[g] - m_g * 1) . q  =  A[g] . q  -  m_g * (1 . q)  =  A[g] . q,

    with `|| A[g] - m_g * 1 || = sqrt(ssq_g - sum_g^2 / d)` precomputed. The whole
    kernel is one sparse matvec over 334798 nonzeros.
    `test_sparse_path_equals_the_dense_computation` pins this against an explicit dense
    computation at rtol 1e-10. **So the implementation-cost half of §4.28's verdict is
    retracted: it was ~2 hours, not half a day, and it is exact.**

    **The result. Paired against the deployed 0.5/0.5 on the same 198 screens:**

    | arm | per screen | per family | wins |
    |---|---|---|---|
    | **ppi 0.25** | **+0.000 [−0.031, +0.038]** | **−0.007 [−0.031, +0.016]** | 48/198 |
    | ppi 0.50 | −0.061 [−0.128, +0.011] | −0.037 [−0.108, +0.028] | 49/198 |
    | ppi 0.25, degree-residualised | +0.003 [−0.028, +0.042] | −0.008 [−0.031, +0.015] | 48/198 |

    **THIS IS NOT UNMEASURABILITY — READ THE INTERVAL WIDTH.** The paired contrast
    resolves to **+-0.035 per screen**, because pairing cancels the screen variance that
    makes the absolute figures wide. §4.18, §4.21 and §4.28 were all closed on "the
    effect could not be demonstrated even if real". **This one excludes every effect
    above +0.04.** The point estimate is zero to three decimals and the win rate is
    48/198, i.e. the arm is WORSE on three quarters of screens than the reference.

    **The dose-response runs the wrong way.** Doubling to 0.5 is actively harmful
    (−0.061 per screen), which is §4.13's displacement mechanism: every slot the PPI
    term promotes comes out of a pick AssayFormer wanted, and here it buys nothing.

    **The fame control is a no-op, which rules out the interesting explanation.**
    §4.28's partial gate-1 failure (degree correlates +0.203 with cross-screen hit
    frequency) suggested the offline margin might be §4.11 fame redundancy. Residualising
    on log degree changes the deployed result by **+0.003** — so it is not that the gain
    was fame; there is no gain either way.

    **THE DECISIVE COMPARISON, and it is the most important measurement here.** §4.28's
    offline oracle margin gave a specific difficulty PROFILE. Deployment contradicts it
    in every cell:

    | group | n | offline oracle margin (§4.28) | **DEPLOYED, paired** |
    |---|---|---|---|
    | hard (EF < 2) | 59 | +0.024 | **−0.015 [−0.041, +0.006]** |
    | middling (2-5) | 51 | **+0.306** | **+0.003 [−0.072, +0.065]** |
    | easy (EF >= 5) | 88 | **+0.290** | **−0.020 [−0.068, +0.017]** |
    | overall per family | 198 | **+0.190 [+0.002, +0.387]** | **−0.007 [−0.031, +0.016]** |

    The one positive cell anywhere is middling genome-scale screens at +0.028
    [−0.016, +0.086], still spanning zero and one cell of twelve chosen after looking.

    **§4.42 FINDING 3 IS REFUTED, and this is the real lesson.** I argued there that
    RULE 1's catastrophic shrinkages were all measured on the CURATED 19-20 screen
    instrument, that the 81-screen instrument should be better behaved, and that
    DepMap's own +0.358 offline -> +0.331 deployed (ratio 0.92) was the applicable
    calibration point — predicting ~+0.17 per family for BioPlex. **Observed −0.007,
    with a CI that excludes +0.17 many times over.** So the 0.92 was not a property of
    the instrument's sample size; it was a property of DepMap. **RULE 1 now has seven
    instances, one sign flip, and — crucially — one instance from the WELL-POWERED
    offline instrument. Sample size was never the issue.**

    **Why the two instruments disagree, which is the mechanism worth carrying.**
    `augment_u.py` fits the feature's weight against the full answer key, so the oracle
    margin answers "is there information in this column usable additively at an OPTIMAL
    weight". The deployed kernel instead drives the column through a seed-built direction
    and fuses by reciprocal rank at a FIXED weight. Those coincide only when the column's
    information is spread across genes AssayFormer ranks poorly. BioPlex's is concentrated
    in hubs, which AssayFormer already ranks highly (§4.25: the easy screens' hits recur
    in a median of 484 of 1349 training screens) — so an oracle weight can exploit it
    while RRF fusion only reshuffles a list that was already right. **The oracle margin
    measures information CONTENT; the deployed campaign measures information that is
    both present AND positioned where the base ranker is wrong. Only the second one is
    worth anything, and no offline instrument in this project measures it.**

    **Verdict: do not deploy, and this time the verdict is a measurement.** BioPlex was
    the best remaining candidate by the project's own numbers and it is now closed with a
    +-0.035 interval centred on zero. Code and arms are kept (`run_ppi_sweep.sh`,
    `sweep-ppi-validation-{w025,w05,w025dr}`) so no one rebuilds it.

---

44. **ATLAS GEOMETRY (§12 item 3): the PREMISE is measurably false for this release —
    there is no dominant shared component to remove, so the campaign was not run.**
    Premise measured 2026-10-23. `codep_seed_acq.py` gains `drop_pcs` and `amp_power`,
    both nested at their off value and test-pinned bit-identical (6 new tests, suite
    265). **This section records a structural measurement and a decision, not a
    deployed null** — see the bottom of the entry before citing it.

    **The idea, which is sound in general.** All eleven of §4's feature candidates
    varied the DATA SOURCE. The geometry of the space the kernel measures similarity in
    was never varied: cosine on raw 808-dim profiles. A dependency matrix's columns are
    not 808 independent measurements — proliferation rate, screen quality and lineage
    are shared by every gene, so two genes both loading on "this line grows fast" earn a
    high cosine for a reason that is about neither of them. Difference-of-means cancels
    that in the DIRECTION (part of why it beats the centroid, §4.15) but nothing cancels
    it in the SPACE, and those are two separate corrections. Removing leading components
    is standard practice in the co-essentiality literature.

    **MEASURE THE PREMISE BEFORE BUYING IT — and here it does not hold.** Eigenvalue
    shares of the per-gene-centred Achilles matrix (18119 genes x 808 lines):

    | components removed | share of variance | that component alone |
    |---|---|---|
    | top 1 | **2.2%** | 2.2% |
    | top 3 | 5.9% | 1.5% |
    | top 5 | 8.2% | 1.1% |
    | top 10 | 12.9% | 0.8% |
    | top 20 | 19.9% | 0.6% |

    A matrix genuinely "dominated by a few huge shared components" has a leading
    component at 20-50%. This one is at **2.2%**, the spectrum is nearly flat, and the
    mean random-pair cosine barely moves (+0.0024 -> +0.0017) after dropping three.
    **The explanation is almost certainly that the correction has already been applied
    upstream:** Chronos models screen quality, proliferation and library effects
    explicitly and ships a corrected gene-effect matrix. The literature advice is aimed
    at less-processed data. **This is a second instance of the §4.33 lesson — check that
    a source holds something you do not already have before arguing for it — and the
    cheapest way to check a geometry argument is to look at the spectrum.**

    **Amplitude weighting has little to work with either.** Median-scaled amplitude runs
    p5 **0.79** to p95 **1.74**: norms over 808 dimensions concentrate, so there is no
    long tail of genuinely flat profiles to demote. The reason to test it anyway is
    §4.25, which measured the hard screens' hits at amplitude 0.138 against 0.137 for
    random genes. **Note this knob cuts AGAINST the hard screens by construction**, so
    the difficulty split matters more than the overall figure.

    **CLOSED ON THE PREMISE, DELIBERATELY NOT CAMPAIGNED — user decision, 2026-10-23.**
    Five arms were specified (`pc3`, `pc10`, `pc20`, `amp05`, `amp10`, doses set from
    the variance curve rather than the usual 1/3/5 since `pc20` is the first that
    removes a fifth of the variance) and `run_geometry_sweep.sh` is committed and
    working. The sweep was launched, then stopped at 120 of 218 screens of the first arm
    and the partial run dirs deleted — a half-finished arm is worse than none, because
    the launcher's guard would SKIP it on any retry (§8). **Nothing was measured
    deployed, and this section is NOT a null result.** The call was that a 2.2% leading
    component and a p5-p95 amplitude range of 0.79-1.74 do not support a correction
    worth half an hour of compute and a numbered section, and that is a reasonable
    reading of the premise table above.

    **WHAT IT WOULD COST TO FINISH:** `bash scripts/diagnostics/run_geometry_sweep.sh
    /home/jeff/depmap validation` (~30 min on CPU, free, no API calls), then
    `score_weight_sweep.py --baseline sweep-dm-validation-w0 --reference
    sweep-dm-validation-dm05knn --arm ...`. The code, the tests and the launcher are all
    in place, so the only cost is the compute.

    **WHAT WOULD MAKE IT WORTH REOPENING.** Two things, and neither is "try it anyway":
    - **A less-processed dependency matrix.** The premise fails because Chronos already
      corrected this release. The raw CERES/readcount-level data, or Sanger's Project
      SCORE at an earlier processing stage, would have a real leading component. That
      is a different experiment from this one.
    - **A space that is NOT already corrected.** The spectrum check is cheap and should
      be run on any future atlas before assuming its geometry is clean — JUMP ORF, for
      instance, is Harmony-corrected, while a raw co-hit or adjacency matrix is not.
      **Add the eigenvalue-share table to `prescreen_atlas.py` for new atlases.**

    **THE TRANSFERABLE LESSON, which is the actual output of this section.** The geometry
    argument arrived with a confident mechanism attached — "dominated by growth rate,
    screen quality, lineage" — and the mechanism is real in the literature and false for
    this file. **Measuring the premise cost one eigendecomposition of an 808x808 matrix
    and saved five campaigns.** That is the §4.33 lesson a second time (check that a
    source holds something you do not already have) and the §4.41 lesson a second time
    (check that the control measures what you think). Put the cheap structural check
    before the expensive measurement, every time.

---

45. **THE DESCRIPTION CHANNEL: a FIELD AUDIT that corrects §4.26(a), and three arms.**
    Built 2026-10-23 on the user's instruction to pursue the `notes` + `direction_str`
    idea. **`direction_str` IS EMPTY ON EVERY SCREEN** — 0 of 218 validation, 0 of 334
    test. §4.26(a) and §4.7 both described it as present but unrendered; it is not
    populated at all, so the arm this project has recommended three times **cannot be
    run as specified.** Audit first, every time.

    | field | rendered? | validation 218 | test 334 | mean chars |
    |---|---|---|---|---|
    | phenotype, cleaned_phenotype, cell_line, cell_type, organism, library_methodology | YES | 100% | 100% | 6-57 |
    | condition_clause | YES | 89.9% | 87.1% | 45-53 |
    | **contrast_label** | YES | **0%** | **0%** | — |
    | **description** | YES | **0%** | **0%** | — |
    | **notes** | no | **100%** | **96.4%** | **233 / 355** |
    | ranking_rationale | no | 100% | 100% | 36-47 |
    | screen_rationale | no | 100% | 100% | 36-39 |
    | significance_criteria | no | 100% | 100% | 33-43 |
    | library_type | no | 100% | 100% | 7-10 |
    | **direction_str** | no | **0%** | **0%** | — |
    | time_point | no | 0% | 0% | — |

    **Three corrections fall out of that table.** (a) The "nine-field renderer"
    effectively renders **seven**: `contrast_label` and `description` are empty
    everywhere, which is why the module's own docstring says public screens "often leave
    `description` empty". (b) **`notes` is the real prize** — 100%/96.4% populated at
    233-355 chars, roughly DOUBLING the description, and never read. (c) Four more
    fully-populated fields are also unread (`ranking_rationale`, `screen_rationale`,
    `significance_criteria`, `library_type`).

    **`notes` often carries the direction semantics in prose**, which is what
    `direction_str` was supposed to supply. Example (screen 1454): *"Genes with a
    negative CASTLE effect are depleted upon KO in response to oxidative stress. Genes
    with a positive CASTLE effect are enriched upon KO..."*

    **So the direction arm DERIVES the label instead of reading it**, via
    `classify_selection`. Which text to classify from was measured rather than assumed:

    | classified from | validation informative | test informative |
    |---|---|---|
    | phenotype alone | 209/218 (96%) | 312/334 (93%) |
    | notes alone | 161/218 (74%) | 151/334 (45%) |
    | **phenotype + notes** | **217/218 (100%)** | **328/334 (98%)** |

    It uses both, whether or not `notes` is rendered. **Be clear what this arm is: a
    DERIVED label over text the model already sees, not new information.** It tests only
    whether a frozen encoder benefits from an explicit canonical token — a weaker
    hypothesis than `notes`.

    **Implementation.** `ASSAYLOOP_DESC_EXTRA` on `screen_description_text`, values
    `notes` / `dir` / `rationale` / `all` / `none`, unset being **byte-identical** to the
    shipped renderer (test-pinned, and that test matters: one character of drift
    invalidates the shipped cache and yields DIFFERENT vectors for the same screen).
    The renderer is called by `amortized_ranker.py:223`, so the switch reaches the model
    and the cache keys with no plumbing. 6 new tests, suite 265 -> 271.

    **Embeddings bought (one-time, ~$0.002 all in):** test fold 230 new vectors
    ($0.00035); validation `notes` 212, `dir` 201, `all` 212. Cache 1356 -> 2211. Four
    transient `model_not_found` 403s were absorbed mid-run — a freshly granted model
    permission propagates unevenly and **the OpenAI SDK does not retry 403**, which is
    why `embed_screen_set.py` chunks at 32, saves after every chunk and retries.

    **ARMS: measured on AssayFormer ALONE (weight=0 knn_weight=0), not on the deployed
    pair**, because the description acts only through `u` and isolating it is cleaner —
    and because `sweep-dm-validation-w0` already IS AssayFormer alone under the default
    renderer, so the reference is free. `sweep-desc-validation-{notes,dir,all}`,
    `run_desc_sweep.sh`.

    **THE CEILING, which bounds every number this section can produce.** The channel is
    worth ~0.9 EF (b_g-only 3.017 -> 3.930 on `paper_test`), but fitting `u` to the
    complete answer key is worth only **+0.41 / +0.84** (§4.7), and this acts solely
    through `u`. **+0.84 is the hard cap.** Off-distribution risk is the main reason it
    may lose: the checkpoint was TRAINED on the current renderer.

    **RESULTS: see §4.47 — all three arms are NEGATIVE, and round 1 shows why.**
    `+notes` −0.012 [−0.056, +0.034], `+dir` −0.052, `+all` −0.105 [−0.247, −0.022]; in
    round 1, where `u` comes from the text alone, the three cost −0.860, −0.520 and
    −0.781 EF. The off-distribution risk flagged above is the finding, not a caveat.

---

46. **THE SECOND HELD-OUT CONFIRMATION (§10 item 1): THE RESULT REPLICATES AT THE
    ENDPOINT, ALMOST EXACTLY — AND THE BUDGET PROFILE DOES NOT.** 2026-10-23, 314 fresh
    `test`-fold screens that chose nothing: no weight, no parameter, no gate. Both arms
    334/334, `n_failed: 0`, zero cache misses. `run_testfold_confirm.sh`,
    `sweep-tf-test-{w0,dm05knn}`.

    **THE ENDPOINT. Confirmed, independently, on two folds:**

    | fold | n | dEF/screen | dEF/family | up/dn/same |
    |---|---|---|---|---|
    | validation (recorded) | 198 | +0.175 [+0.077, +0.290] | +0.331 [+0.180, +0.515] | 75/75/48 |
    | **test (NEW)** | **314** | **+0.257 [+0.050, +0.476]** | **+0.398 [+0.114, +0.750]** | 124/85/105 |

    Genome-scale only (library > 1001, excluding the metric-degenerate screens, §4.29),
    which is the honest comparison:

    | fold | n | dEF/screen | **dEF/family** |
    |---|---|---|---|
    | validation | 100 | +0.381 [+0.227, +0.563] | **+0.485 [+0.281, +0.730]** |
    | **test** | **248** | +0.330 [+0.070, +0.591] | **+0.484 [+0.130, +0.917]** |

    **+0.485 and +0.484 per family, on two folds that share no screens.** Every interval
    excludes zero. The headline is no longer a single held-out measurement.

    **THE CAVEAT, AND IT IS SERIOUS. The per-round budget curves disagree completely.**
    Both genome-scale, dEF at each round boundary:

    | | r1 | r2 | r3 | r4 | r5 | r6 | r7 | r8 | r9 | r10 |
    |---|---|---|---|---|---|---|---|---|---|---|
    | validation | +0.00 | +0.35 | +0.47 | **+0.53** | +0.41 | +0.38 | +0.34 | +0.29 | +0.36 | +0.38 |
    | **test** | +0.00 | +0.10 | −0.41 | −0.85 | −1.05 | −1.31 | **−1.44** | −0.63 | +0.32 | **+0.33** |

    On validation the method is ahead at every budget. **On test it is BEHIND from 300 to
    800 picks, by as much as −1.44 EF at round 7, and only overtakes in the final two
    rounds.** Confirmed by dnAUC averaged over the ten budgets: validation **+0.144
    [+0.053, +0.241]**, test **−0.397 [−0.855, −0.008]** — both intervals exclude zero,
    with opposite signs.

    **What this does and does not mean.**
    - The headline stands **at the paper's budget**: 10 rounds of 100, EF at 1000. That
      is the benchmark's metric and the number to quote, now twice-confirmed.
    - It is **NOT robust to stopping early.** A practitioner who halted at round 5-7 on
      test-fold-like screens would be worse off than running AssayFormer alone. Any
      claim of the form "this method helps" must say "over a full ten-round campaign".
    - The mechanism is consistent with §4.13's displacement account: the kernel gives up
      picks AssayFormer wanted in order to mine co-essential neighbourhoods, and that
      trade is paid for up front and collected at the end. On validation the collection
      starts immediately; on test it arrives only in rounds 9-10. **Why the two folds
      differ in WHEN the payoff lands is not explained and is the obvious next
      question.**

    **THIS UN-RETRACTS HALF OF §4.42.** That section concluded "nAUC is the wrong summary
    for an intervention whose effect accumulates — retract that half of §12 item 4",
    because on the validation fold nAUC merely diluted a uniformly positive effect and
    cost SNR. On the test fold the same instrument is the ONLY one that reveals a −1.44
    EF trough, and it is significant with the opposite sign to the endpoint. **nAUC is
    not a better summary than EF-at-1000 and not a worse one; it answers a different
    question — "is the gain robust across budgets?" — and here the answer is no.**
    Report both from now on. A single-endpoint metric hid this for the entire project.

    **Instrument notes.** The per-round curve requires pick PREFIXES;
    `_adj_ef_from_run`'s `budget` argument does not truncate the pick list, it only
    rescales the denominator (§4.42). ANCOVA is worth ~0-2% here against 2-8% on
    validation (R² 0.014-0.056 on the endpoint), so the covariates explain even less of
    the test fold's spread. 286 of 334 test-fold descriptions had to be embedded first
    (230 unique texts, $0.00035, §4.45).

---

47. **THE DESCRIPTION CHANNEL IS CLOSED: giving AssayFormer the metadata it never reads
    makes it WORSE, and round 1 identifies the mechanism exactly.** 2026-10-23, the
    §4.45 arms. Three arms x 218 screens, `n_failed: 0`, zero cache misses, measured on
    AssayFormer alone (weight=0 knn_weight=0) against `sweep-dm-validation-w0`, which is
    that same configuration under the default renderer.

    | arm | text added | dEF/screen at budget 1000 | dEF/family | dnAUC |
    |---|---|---|---|---|
    | **+notes** | +186 chars | **−0.012 [−0.056, +0.034]** | +0.015 | −0.122 [−0.248, +0.004] |
    | +dir (derived label) | +29 chars | −0.052 [−0.183, +0.022] | −0.095 | −0.061 [−0.183, +0.051] |
    | **+all unused fields** | +375 chars | **−0.105 [−0.247, −0.022]** | −0.125 | −0.112 [−0.261, +0.042] |

    **ROUND 1 IS THE RESULT.** In round 1 no labels exist yet, so `u` is derived from the
    text embedding **and nothing else** — it is the cleanest isolation of the description
    channel available anywhere in this project. All three variants are heavily negative
    there:

    | | r1 | r2 | r3 | r5 | r8 | r10 |
    |---|---|---|---|---|---|---|
    | +notes | **−0.860** | −0.166 | +0.008 | +0.018 | −0.069 | −0.012 |
    | +dir | **−0.520** | −0.080 | −0.029 | +0.044 | −0.017 | −0.052 |
    | +all | **−0.781** | −0.034 | +0.030 | +0.048 | −0.106 | −0.105 |

    **A richer description produces a WORSE screen latent, by 0.5-0.9 EF, when the
    description is all the model has.** The label feedback then repairs most of it over
    the following rounds, which is why the endpoints are only mildly negative. **The
    off-distribution risk §4.26(a) flagged is not a caveat here — it is the whole
    finding.** The checkpoint was trained on the seven-field renderer, and text outside
    that distribution degrades `u` rather than informing it.

    **THE CLOSURE IS TIGHT, WHICH IS WHAT MAKES THIS WORTH SOMETHING.** `+notes` is the
    strongest hypothesis in the group — genuinely unread text, 100% populated, roughly
    DOUBLING the description — and its interval is **±0.045**, excluding any gain above
    **+0.034**. Against a channel worth ~0.9 EF and an oracle-`u` ceiling of +0.84
    (§4.7), essentially the entire plausible range is excluded. **§4.26(a) called this
    "the largest untried lever in the project". It is now tried, and it is closed.**

    **Three readings worth keeping.**
    - **More text is not monotonically worse, so this is not a length effect.** `+dir`
      adds 29 characters and lands below `+notes`, which adds 186. A novel categorical
      line (`Selection direction: enrichment`) appears to perturb the encoder more per
      character than prose resembling the fields already there — though the two
      intervals overlap, so treat that as a hypothesis, not a measurement.
    - **The derived-label arm was always the weaker hypothesis and it behaved like it.**
      It added no information, only an explicit token over text the model already saw
      (§4.45). The one thing it rules out is that an explicit canonical direction token
      helps a frozen encoder.
    - **This is the only input-side lever that needed no retraining, and it is spent.**
      What remains for the description channel is retraining the encoder on the richer
      renderer, which is out of this project's compute budget (§4.12's conclusion).

    **Cost: ~$0.002 of `text-embedding-3-small`** for 625 new vectors across the three
    variants (§4.45). `run_desc_sweep.sh`, `sweep-desc-validation-{notes,dir,all}`.
    `ASSAYLOOP_DESC_EXTRA` stays in the code, defaulting off and test-pinned
    byte-identical, so the arms can be re-scored or extended without rebuilding.

---

48. **"COULD THE LLM HELP *IN COMBINATION* WITH `codep_seed`?" — the seed-quality
    mechanism it would need is measured, and it is absent.** Free, from runs already on
    disk, 2026-10-23. Asked because the question is natural and recurring: the two
    interventions act on different channels, so §4's LLM nulls do not automatically
    cover the pair.

    **PRECISION FIRST, because "the LLM as a nominator" sounds like a new idea and is
    not.** The LLM's ROLE here is exactly the paper's warm start: occupy early slots,
    propose genes, have them tested. Nothing about that changes. What changes is only
    **which downstream consumer reads the resulting labels.** A hit the LLM finds
    contributes through three channels:

    | channel | paper's design | with `codep_seed` |
    |---|---|---|
    | (a) the hit itself counts toward EF | yes | yes |
    | (b) its label updates AssayFormer's `u` | yes | yes |
    | **(c) it points the kernel's direction in DepMap space** | **no** | **yes** |

    (a) and (b) are shared and both already closed: (b) is capped at +0.41/+0.84 by §4.7
    and §4.11 found the LLM's gene knowledge fully inside span(V) with zero incremental
    value, so in the paper's design a new hit's INFORMATION has almost nowhere to go.
    **(c) is the only new channel, and it is new for one structural reason** — the
    co-essentiality direction lives OUTSIDE span(V, b, 1), where only 3.5% of a test-fold
    hit vector is expressible, which is the same fact that makes `codep_seed` work where
    twelve `u`-based interventions failed. So the hypothesis is not "a new role for the
    LLM" but "the labels it generates now have a consumer that is not capped at +0.84".

    **The COST side is unchanged.** Every slot the LLM takes still displaces an
    AssayFormer pick; §2 measured ~245 new hits against ~1730 given up on the curated 20.
    Adding the kernel does nothing about that trade.

    **Which is why the oracle test is the right one: it isolates channel (c) and sets the
    cost to ZERO.** The kernel is handed unreachable hits for free, no slots spent, and
    the seeds are blocked from being picked so they earn no credit. It therefore measures
    the entire upside of (c) with none of the downside. If the 200-seed arm is flat, (c)
    is empty, and (a) and (b) being already closed, the family goes with it.

    **That story requires seed quality to drive the kernel's benefit. It does not.**

    | hits AssayFormer alone confirms by round 3 | n | mean deployed dEF |
    |---|---|---|
    | 0-1 | 71 | +0.154 |
    | 1-3 | 70 | +0.312 |
    | 3-10 | 53 | +0.066 |
    | 10-246 | 50 | +0.171 |

    **Spearman(seeds by round 3, deployed dEF) = +0.032, p = 0.66.** No relationship, and
    the quartiles are non-monotone. Feeding the kernel more early hits is not what makes
    it work.

    **A CORRECTION worth carrying: seeds are SCARCE on the fold, not abundant.** §10
    item 6 says "AssayFormer reveals ~40 hits per batch on validation" — that is the
    CURATED 20, where its early hit rate is 39.6% (§2a). On the 198-screen fold the
    median is **3 hits by round 3**; **71 of 198 screens (36%) never reach the 2 hits the
    kernel needs in that window**, and 148 of 198 are below 10. The displacement
    arithmetic that sank the LLM warm start was computed on the curated 20 and **does not
    transfer to the fold**, where AssayFormer's own early marginal pick is far weaker.
    That is a genuine argument FOR re-examining the pair, and the strongest one available.

    **But the decisive measurement already exists and it is not about displacement.**
    §3 restricted LLM mixing to the enrichment screens — exactly the subgroup where the
    LLM is at parity with AssayFormer (4.3% against 5.3%), so displacement is close to
    free. Total hits went from −145/−224 to **+13/+9** and EF from −0.01/+0.05 to
    **+0.02/+0.05**: it removed the damage and **created no gain**. Where the LLM costs
    nothing, it also contributes nothing. Its ~12 new hits per screen do not convert into
    EF on their own, and §4.48 now shows they do not convert through the kernel either.

    **§4.46 adds a new argument against.** The deployed method's payoff already arrives
    LATE on the test fold — it is behind from 300 to 800 picks. An LLM warm start spends
    rounds 1-3, which deepens exactly that early deficit.

    **THE RIGHT EXPERIMENT IS FREE, AND IT SHOULD BE RUN BEFORE ANY BILLED ONE.**
    Substitute an ORACLE for the LLM: seed the kernel with hits taken from the answer key
    and chosen to be ones AssayFormer ranks poorly, then measure the deployed gain. That
    is an upper bound on "better seeds help the kernel", it costs nothing, and it follows
    the project's own pattern of bounding a family with an oracle before building a
    member of it (§4.7, §4.9, §4.11). **If the oracle-seeded kernel is flat, the whole
    LLM-plus-kernel family is closed for free.** If it is large, the question becomes
    whether any real nominator can approach it — and only then is an LLM run worth money.

    **Expected outcome if someone runs the LLM pair anyway:** parity at best. The one
    LLM configuration ever positive was the decline-triggered handoff (+0.08 on
    `paper_test`, §1), and that is one cell of many chosen after looking.

    **ANSWERED THE SAME DAY by §4.49: the oracle bound is NEGATIVE.** Channel (c) does
    not merely fail to help — 12 free unreachable true hits per screen cost **−0.140
    [−0.221, −0.070]** against the deployed configuration, monotonically worse at 50 and
    200, with the mechanism explained by §4.25's flat-profile measurement. **Do not
    spend on an LLM nominator.**

---

49. **THE ORACLE-SEED BOUND IS NEGATIVE, MONOTONICALLY: perfect seeds make the kernel
    WORSE, and §4.25 says why. The LLM-as-nominator family is closed, for free.**
    2026-10-23, the §4.48 test. 4 arms x 218 screens, CPU, no API calls.
    `run_oracle_seeds.py` + `oracle_seeds` on `CodepSeedAcquisition`.

    Paired against the deployed configuration on the same 198 fresh screens:

    | arm | seeds | per screen | per family | wins |
    |---|---|---|---|---|
    | **oracle 12** | 12 unreachable true hits | **−0.140 [−0.221, −0.070]** | −0.176 | 56/198 |
    | oracle 50 | 50 | −0.221 [−0.311, −0.145] | −0.304 | 49/198 |
    | oracle 200 | 200 | **−0.263 [−0.353, −0.186]** | −0.365 | 42/198 |
    | *CONTROL: 50 random non-hits* | 50 | *−0.127 [−0.351, +0.061]* | *−0.438* | *84/198* |

    **Every oracle interval excludes zero, and the dose-response is monotone in the
    WRONG direction.** Handing the kernel more perfectly-correct hits it could not
    otherwise reach makes it progressively worse. There is no upper bound to chase: the
    ceiling on "a better nominator" is **negative**.

    **WHY, and this is the satisfying part — it follows from §4.25.** The query direction
    is `mean(hit profiles) − mean(non-hit profiles)`. Oracle seeds are by construction
    hits **AssayFormer cannot reach**, and §4.25 measured exactly what those genes look
    like: on the hard screens the hits have DepMap amplitude **0.138 against 0.137 for
    random genes** — near-flat profiles. **So the hits a nominator could contribute are
    precisely the genes whose co-essentiality profiles carry no signal, and averaging
    them into the query DILUTES it.** The seeds are correct about biology and useless
    about geometry. §4.25 and §4.48 were measuring two halves of the same fact.

    **The control confirms the reading.** 50 random NON-hits score −0.127, statistically
    indistinguishable from 12 true unreachable hits at −0.140. **In profile space a true
    unreachable hit is worth about as much as a random gene** — which is the same
    statement as the amplitude figures above, arrived at from the other direction.

    **ONE CONFOUND, STATED PLAINLY, AND WHY IT DOES NOT RESCUE THE HYPOTHESIS.** Oracle
    seeds are removed from the candidate pool so they cannot be banked as free hits
    (§4.48). That necessarily denies the arm any credit for finding those genes itself,
    and the handicap grows with the dose — **so the orc50 and orc200 slopes are partly
    the block, not only the seeds, and the dose-response should not be read as a pure
    dilution curve.** But the decisive arm is **orc12**, where only 12 genes per screen
    are blocked (median pool 17), none of which either comparison arm ever acquired, out
    of 1000 picks — a handicap far too small to produce −0.140 [−0.221, −0.070] at 56/198
    wins. **orc12 alone refutes the synergy hypothesis**, and the mechanism above
    explains the direction. Anyone wanting the dose-response clean would re-run with the
    seeds left pickable and the free hits subtracted out afterwards.

    **WHAT THIS CLOSES.** §4.48 laid out the three channels a nominated hit can act
    through: (a) counting as a hit, (b) updating `u`, (c) pointing the kernel. (a) is
    the displacement trade, already measured at ~245 new hits against ~1730 given up.
    (b) is capped at +0.41/+0.84 and §4.11 found the LLM's knowledge fully inside
    span(V). **(c) is now measured and NEGATIVE.** All three channels are closed, so
    **the LLM-plus-kernel pair needs no LLM spend to rule out** — which was the point of
    bounding the family with an oracle rather than buying the experiment.

    **CONFIRMED BY A BETTER-POSED TEST, §4.50**, which removes this section's one
    confound entirely: injected hits enter the campaign as real picks and the free-hit
    credit cancels in a 2x2. Same answer, more sharply — the kernel's benefit goes from
    **+0.331 to −0.074** per family when the campaign is handed perfect early hits.
    **Prefer §4.50 when citing this.**

    **A GENERAL LESSON WORTH MORE THAN THE RESULT.** "Find hits the model misses" sounds
    like an unambiguous good and is the implicit goal of every warm start in the paper
    and in §4. It is not: a hit is only useful to a similarity-based method if the gene
    also carries usable structure in that similarity space. **Coverage of the answer and
    usefulness as a seed are different properties, and on the hard screens they are
    close to anti-correlated.**

---

50. **THE 2x2 INJECTION TEST — the best-posed version of §4.48, and it settles the
    question: a perfect nominator makes the kernel REDUNDANT, not better. The two
    interventions are SUBSTITUTES.** 2026-10-23, user's design, which is strictly better
    than §4.49's and worth understanding as a piece of experimental method.

    **WHY THIS DESIGN BEATS §4.49's.** §4.49 handed the scorers seeds from OUTSIDE the
    campaign, which forced those genes out of the candidate pool so they could not be
    banked as free hits — and that denied the arm any credit for finding them itself, a
    handicap growing with dose (its one confound). Here the injected genes enter the
    campaign as **real picks**, taking the lowest-ranked slots of each of the first three
    batches: budget stays exactly 10 x 100, nothing is confiscated, and it is literally
    what a nominator does. The free-hit credit is handled **structurally** — the same
    genes go to a kernel-ON and a kernel-OFF arm, so the credit is identical and cancels
    in a difference of differences. No accounting correction.

    | kernel effect, per family | value |
    |---|---|
    | with ordinary seeds (the deployed result) | **+0.331 [+0.180, +0.515]** |
    | with 4 injected hits/round x 3 | **−0.074 [−0.290, +0.112]** |
    | with 12 injected hits/round x 3 | **−0.086 [−0.284, +0.087]** |
    | **SYNERGY (difference of differences), 4/round** | **−0.405 [−0.641, −0.209]**, up 57/198 |
    | **SYNERGY, 12/round** | **−0.418 [−0.640, −0.228]**, up 59/198 |

    **Hand the campaign perfect early hits and the kernel's entire benefit evaporates**
    — +0.331 becomes −0.074, both injected intervals now spanning zero, and the change
    is significant at both doses. The two doses agree (−0.405 vs −0.418), which is
    reassuring given the high dose is capped on about half the fold (only 92 of 218
    screens hold 36 unreachable hits; 122 hold 12).

    **For scale, what the injection itself bought, kernel OFF: +5.357 and +6.868 EF**
    (absolute EF 4.688 -> 9.553 -> 10.511). That is not a result and nobody can deploy
    it; it is the arithmetic of being handed 12 to 36 free hits, and it is reported only
    to show the arms were doing what they were told.

    **THE CONCLUSION, and it is stronger than §4.49's.** `codep_seed` and a better
    nominator are **substitutes, not complements**. Both attack the same deficit — too
    few and too uninformative early hits — and once either one fixes it, the other has
    nothing left to add. This is why no LLM-plus-kernel configuration can work: it is
    not that the LLM is too weak (§4 refuted that directly), nor that its hits are
    redundant with span(V) (§4.11, true but about a different channel). **It is that the
    kernel's value is CONTINGENT on the campaign being in a low-hit regime, and a
    nominator good enough to help would remove the regime the kernel needs.**

    **TWO MECHANISMS ARE CONSISTENT WITH THIS AND THE TEST DOES NOT SEPARATE THEM.**
    (a) Dilution: the injected hits are by construction genes AssayFormer cannot reach,
    which §4.25 measured as having near-flat DepMap profiles (amplitude 0.138 vs 0.137
    random), so averaging them into the query direction degrades it — the §4.49 reading,
    supported there by random non-hits scoring as well as true hits. (b) Headroom: with
    EF at 9.5-10.5 against perfect play's 16.77 the campaign is at ~60% of ceiling rather
    than ~28%, AssayFormer's own picks are better because its history is better, and the
    kernel's promotions now displace picks that were already good — the §4.13
    displacement mechanism. **Both point the same way for the practical question, so the
    answer does not depend on which it is.** Separating them would need a budget-neutral
    arm injecting random NON-hits into the same slots, which isolates displacement from
    seed quality.

    Code: `inject_genes` / `inject_per_round` / `inject_rounds` on
    `CodepSeedAcquisition` (+ `_apply_injection`, applied on EVERY return path — it was
    initially skipped on the `weight=0` shortcut, which is exactly the kernel-OFF arm,
    and a test caught it), `run_inject_2x2.py`, `sweep-inj-validation-d{4,12}k{0,1}`.
    5 new tests, suite 276 -> 281.

---

51. **LLM WARM START + `codep_seed`, RUN DIRECTLY: the two effects are ADDITIVE AND
    INDEPENDENT. The kernel rescues the warm start; the warm start still costs more than
    the kernel adds.** 2026-10-27, user's request. Free -- the Flash-Lite warm-start
    rounds were already on disk (`sweep-gem35val`, `sweep-gem35test`, 3 steps each with
    `llm_calls.jsonl`), so the LLM rounds are REPLAYED as history and no API call was
    made. `handoff_fullgenome.py --inner-acq codep_seed` (new flag, defaults to `greedy`
    so the shipped behaviour is unchanged). Arms `sweep-hybrid-{pv,pt}-gem35codep`.

    | | paper_validation | paper_test |
    |---|---|---|
    | AssayFormer alone | 6.100 | 4.828 |
    | Flash-Lite warm start -> AssayFormer *(the paper's design)* | 5.170 | 4.569 |
    | **Flash-Lite warm start -> OUR METHOD** *(new)* | **5.346** | **5.118** |
    | OUR METHOD, no LLM | **6.217** | **5.389** |

    Paired, family-bootstrapped, same 20 screens:

    | contrast | paper_validation | paper_test |
    |---|---|---|
    | combined **vs the paper's warm start** | **+0.257 [+0.036, +0.521]** 15/20 | **+0.549 [+0.303, +0.832]** 18/20 |
    | combined **vs our method alone** | **−0.461 [−0.801, −0.173]** 4/20 | −0.271 [−0.932, +0.335] 10/20 |
    | combined vs AssayFormer alone | **−0.264 [−0.620, −0.014]** 5/20 | +0.290 [−0.400, +0.927] 11/20 |

    **THE DECOMPOSITION IS THE RESULT.** The kernel contributes almost exactly the same
    amount whichever start it is given:

    | kernel's gain, measured from... | paper_validation | paper_test |
    |---|---|---|
    | an LLM-seeded start (rounds 1-3 from Flash-Lite) | +0.177 | **+0.549** |
    | a cold start (AssayFormer's own rounds 1-3) | +0.118 | **+0.562** |

    So `combined ≈ AssayFormer alone + (LLM cost) + (kernel gain)`, with the two terms
    independent. **The kernel does not repair the warm start's displacement cost; it adds
    its own gain on top of whatever it is handed.** On `paper_test` that nets out to
    +0.29 over AssayFormer alone; on `paper_validation`, where AssayFormer's early hit
    rate is 40% and the LLM's is 15%, the displacement cost is −0.93 and the combination
    ends up **significantly WORSE than doing nothing at all**.

    **Answering the question as asked.** Does the combination beat the best LLM method?
    - **Against the paper's own design, yes and significantly, on both sets.** Adding the
      kernel to the Flash-Lite warm start is worth +0.257 and +0.549.
    - **Against the BEST LLM arm** (decline-triggered handoff, 4.904 on `paper_test`):
      +0.213 [−0.136, +0.617], 10/20. **Cannot be claimed on 20 screens.**
    - **Against our method with no LLM at all: no.** Worse on both sets, significantly so
      on `paper_validation`. **There is no configuration in which adding the LLM helps.**

    **Consistent with §4.50 and arrived at independently.** §4.50 showed the kernel and a
    better gene-nominator are SUBSTITUTES, by giving the campaign perfect early hits. This
    gives the campaign WORSE early hits and finds the same structure from the other side:
    the kernel's contribution is insensitive to what it is seeded with (§4.48, Spearman
    +0.032 between seed count and benefit), so the LLM's cost passes straight through to
    the total.

    **Caveats.** 20 curated screens each, which chose this project's settings, so the
    intervals are wide and the absolute figures run high. A fold-scale version is not
    possible without new LLM spend: warm-start traces exist only for the two curated 20s.

---

52. **THE PAPER'S OWN PER-ROUND DATA ARRIVED, AND OUR BASELINE REPRODUCES IT TO THE
    HIT.** 2026-10-27, user supplied `paper_per_round_data.xlsx`: 59 methods x 10 rounds
    x the paper's 20 held-out test screens. Read with `read_paper_xlsx.py` (stdlib only —
    `openpyxl` is not in this environment and an .xlsx is a zip of XML).

    **TWO IDENTITY CHECKS, BOTH PASSED — do not plot against this file without them.**
    (a) The workbook's 20 screens ARE `paper_test`: every name matches and every
    `total_hits` / `library_size` agrees. (b) **Our AssayFormer-alone run reproduces the
    paper's `AssayFormer (BPMF + GRPO)` row to the hit on all 20 screens** — mean
    |difference| **0.000**, correlation 1.000, at every one of the ten rounds. The
    reproduction is exact, not approximate.

    **A TRAP THAT WOULD HAVE COST TWO EF.** The workbook ALSO has a row plainly labelled
    `AssayFormer`, and it is NOT the released model — it scores EF ~2.4 against the
    released checkpoint's 4.83. The released model is the **BPMF + GRPO** variant. Match
    by per-screen hit counts, never by the name.

    **EF from this file is a BAND, not a point.** Our adjusted EF forgives picks that are
    real genes outside a screen's own library; `By_screen_long` records `n_acquired` and
    `out_of_universe` but not that forgiven count, which exists only as a cross-screen
    mean in `Mean_long`. So paper arms get `[forgive nothing, forgive the mean]`. The
    bracket is validated: for the arm we can check exactly, the true 4.828 lies inside
    its band 4.407–4.896.

    **THE COMPARISON, on the paper's own 20 test screens:**

    | arm | EF | source |
    |---|---|---|
    | Paper: Gemini 3.1 Pro + AssayLoop | **5.482** [5.265, 5.699] | paper |
    | Paper: GPT-5.6 Sol + AssayLoop | 5.471 [5.250, 5.692] | paper |
    | **OUR METHOD, no LLM at all** | **5.389** | this project |
    | + BioPlex | 5.348 | this project |
    | LLM warm start -> our method (§4.51) | 5.118 | this project |
    | Paper: GLM-5.1 + AssayLoop | 5.094 [4.871, 5.317] | paper |
    | best LLM arm we ran (decline-triggered) | 4.904 | this project |
    | AssayFormer alone | 4.828 | both, identical |
    | **Paper: Gemini 3.1 Pro ALONE** | **4.774** [4.757, 4.791] | paper |

    **TWO FINDINGS WORTH THE WHOLE EXERCISE.**
    - **Our method lands INSIDE the paper's Gemini-3.1-Pro + AssayLoop band, using no
      language model at all.** 5.389 against [5.265, 5.699]. The paper's headline result
      and ours are not distinguishable on this set — and ours needs no API calls, no
      Gemini access and no retraining.
    - **Gemini 3.1 Pro ALONE scores 4.774, BELOW AssayFormer's 4.828.** The LLM is not a
      better picker than the model; the paper's gain comes entirely from the handoff
      structure, not from the language model being good at this task.

    **WHY OUR LLM REPRODUCTION FELL SHORT, now visible per round.** Gemini 3.1 Pro hits
    **14.1%** in round 1 and its AssayLoop handoff **17.3%**, where our Flash-Lite manages
    **10.7%** and AssayFormer **7.2%**. On `paper_test` AssayFormer's early rate is low
    enough that a strong LLM genuinely helps — §2a's displacement arithmetic is a
    `paper_validation` phenomenon, where AssayFormer opens at 40%. **The paper's claim is
    credible on its own test set; what we could not reproduce was the MODEL, not the
    method.** §10 item 8 (run Gemini 3.1 Pro) is the remaining gap, and this file
    quantifies exactly what it would be worth.

    Code: `read_paper_xlsx.py`; notebook figure 10 and 5 extra CSV rows carrying a
    `source` column. Arms plotted dashed with the EF band shaded.

---

53. **TOP-UP FIXES THE SHORTFALL AND NOT THE SCORE, CONFIRMED ON `paper_test` WITH A PAID
    RUN. Our handoff still beats the paper's by +0.422 when the warm start is COMPLETE, so
    the gain was never GLM-5 leaving slots empty.** 2026-10-27, user's explicit request.
    **Cost $0.36** (135 Bedrock calls, 190k in / 52k out, measured with `bedrock_cost.py`) --
    the only paid run of the session. `sweep-glm5testtopup` (`llm_single_topup`, n=3,
    `collect-glm-5-bedrock.yaml`), then both handoffs replayed off it for free:
    `sweep-glm5tu-pt-af` (greedy) and `sweep-hybrid-pt-glm5tucodep` (`codep_seed`).

    Top-up did its job: `paper_test` shortfall **22.1% -> 2.9%**, and round 3 -- where the
    plain run collapses to 54.5 genes of 100 -- recovers to **93.7**.

    | `paper_test`, 20 screens | EF |
    |---|---|
    | GLM-5 no top-up -> AssayFormer *(the paper's design)* | 4.647 |
    | GLM-5 **+ top-up** -> AssayFormer | 4.805 |
    | GLM-5 no top-up -> OUR METHOD | 5.182 |
    | **GLM-5 + top-up -> OUR METHOD** | **5.227** |
    | AssayFormer alone *(reference)* | 4.828 |
    | **OUR METHOD, no LLM** *(reference)* | **5.389** |

    Paired, family-bootstrapped, same 20 screens:

    | contrast | delta |
    |---|---|
    | our handoff **vs the paper's**, both with top-up | **+0.422 [+0.271, +0.577]** 16/20 |
    | top-up vs none, **under AssayFormer** | +0.158 [−0.026, +0.350] 13/20 |
    | top-up vs none, **under our method** | +0.045 [−0.183, +0.246] 12/20 |
    | LLM+top-up+our method vs **our method with NO LLM** | −0.163 [−0.799, +0.407] 11/20 |

    **FOUR THINGS THIS SETTLES.**

    (a) **§2c's validation finding holds on test.** Filling every slot moves the paper's
    design +0.158, interval crossing zero, and leaves it at 4.805 -- still below
    **AssayFormer alone at 4.828**. Three rounds of a frontier LLM, fully spent, is worth
    slightly less than not running it. The LLM's problem was never unfilled slots.

    (b) **Our advantage is not an artefact of the LLM's shortfall.** The obvious objection
    to §4.51 was that `codep_seed` only looked good because it had more rounds of real
    budget to work with. With the warm start at 97.1% filled the gap is **+0.422
    [+0.271, +0.577]**, barely moved from +0.535 on the incomplete warm start.

    (c) **Top-up and our method are near-substitutes, a third confirmation of §4.48-51.**
    Top-up buys +0.158 under AssayFormer but only **+0.045** under `codep_seed`: once the
    kernel is doing the picking, the extra LLM genes stop paying. Same pattern as oracle
    seeds and as the better nominator.

    (d) **The best LLM arm still loses to no LLM at all.** 5.227 against 5.389, point
    estimate negative. The interval is wide (11/20 screens), so this is "no evidence the
    LLM helps", not "proof it hurts" -- but after top-up there is no remaining story in
    which the LLM is being held back by its plumbing.

    **NAMING.** Bedrock hosts `zai.glm-5`; the paper's GLM-5.1 is the teacher model and is
    not available on this account. Every row above is **GLM-5**.

---

54. **THE DEPMAP ATLAS SURVIVES LOSING EVERY EVALUATION CELL LINE, AND THE
    CHROMOSOMAL-PROXIMITY CORRECTION THE LITERATURE ASKS FOR IS A WASH.** 2026-10-27,
    user's request, free (CPU only). Two independent interventions on the atlas, both
    scored paired and family-bootstrapped against CPU baselines whose round-1 batches
    are byte-identical to the variant's.

    **(a) GLOBAL CELL-LINE ABLATION.** `achilles_noeval.npz` drops every Achilles line
    used by ANY of the four evaluation sets -- **40 of 808 lines (5%), touching 154
    screens** -- so for every screen we report on, its own line is absent BY
    CONSTRUCTION rather than measured unexploited. 18119 genes keep a non-degenerate
    profile. Resolver: HGNC-free, `Model.csv` name -> ModelID, the same one §4.35 used
    (its 43-of-100 figure reproduces exactly).

    | set | n | full atlas | 40 lines removed | paired delta |
    |---|---|---|---|---|
    | paper_validation | 20 | 6.217 | 6.231 | +0.014 [−0.014, +0.041] |
    | paper_test | 20 | 5.388 | 5.362 | −0.026 [−0.056, −0.001] |
    | validation | 100 | 5.054 | 5.017 | **−0.036 [−0.070, −0.008]** |
    | test | 248 | 6.896 | 6.922 | **+0.026 [+0.001, +0.056]** |

    **(b) PROXIMITY CORRECTION.** `gene_locus=` + `prox_level=arm` excludes seeds that
    share a candidate's chromosome arm, which is what the co-essentiality literature
    does before calling a network and what this project never did.

    | set | n | uncorrected | arm-level correction | paired delta |
    |---|---|---|---|---|
    | paper_validation | 20 | 6.217 | 6.210 | −0.007 [−0.038, +0.020] |
    | paper_test | 20 | 5.388 | 5.430 | +0.042 [−0.049, +0.172] |
    | validation | 100 | 5.054 | 4.990 | **−0.064 [−0.122, −0.016]** |
    | test | 248 | 6.896 | 6.961 | **+0.065 [+0.015, +0.141]** |

    **NEITHER IS WORTH ADOPTING FOR THE SCORE, AND THE REASON IS THE SAME PATTERN.**
    Both interventions remove information from the atlas; both help on `test` and hurt
    on `validation`, with intervals excluding zero in OPPOSITE directions on the two
    held-out folds, and both at about 0.03-0.07 against a method effect of +0.485. This
    is the §4.46 fold disagreement again (endpoint replicates, budget curve does not) in
    a third guise, and it is the open question in §10 worth the most.

    **The ablation is still worth shipping, for the CLAIM rather than the number.** It
    converts "we audited for leakage and found none" into "the atlas cannot contain this
    screen's answer". Keep §4.35's leave-one-out as the leakage RESULT (+0.010
    [+0.001, +0.024], no sign disagreement): a global ablation conflates removing the
    leaky column with shrinking the atlas 5% for every screen, and the split below shows
    the confound is real -- the loss lands on the right screens but flips sign by fold,
    which no leakage mechanism explains.

    | set | own line removed | own line never in the atlas |
    |---|---|---|
    | validation | −0.070 (n=43) | −0.011 (n=57) |
    | paper_test | −0.049 (n=7) | −0.014 (n=13) |
    | test | **+0.094** (n=48) | +0.010 (n=200) |

    **THREE PROCESS FAILURES IN ONE SESSION, ALL THE SAME SHAPE: an arm that differs
    from its baseline in more than the thing under test.**
    1. **`assayloop run` needs `--full-genome --min-screen-freq 2`** or it scores against
       each screen's own library, not the shared universe. The first ablation run looked
       like a −2.94 catastrophe; it was a different candidate pool. **Check
       `config.json`'s `universe_genes` length against the baseline's before scoring.**
    2. **Two baselines ran on `cuda`** (`sweep-dm-pt-dm05knn`, `sweep-dm-validation-dm05knn`)
       and GPU0 is now dead, so every new arm is CPU. Measured, the device is worth
       **−0.001 / +0.003** -- small, but it breaks byte-identity. `sweep-cpubase-{pt,validation}`
       are the CPU twins; use them for anything compared against a CPU arm.
    3. **ROUND 1 IS THE CHECK.** It has no seeds, so the kernel is inactive and the batch
       must be byte-identical to the baseline's. It caught (1) and (2) before either
       reached a conclusion, and it costs one comparison. **Run it before scoring any new
       arm.**

    Code: `codep_seed_acq.py` gains `gene_locus` + `prox_level` (band/arm/chrom), default
    OFF and inert; `_seed_sets()` extracted so the plain and corrected scorers share one
    definition of which labels feed the query. The correction is EXACT, not sampled -- the
    mean over seeds is a sum, so excluding a locus is one matrix-vector product per
    distinct locus (42 at arm level) rather than one per gene. 6 new tests, suite
    `test_codep_seed.py` 70 -> 76. Arms `sweep-nolk-*`, `sweep-prox-*`, `sweep-cpubase-*`.

---

55. **THE PAPER SHIPS BOTH kNN BASELINES, AND NEITHER SURVIVES FUSION. "We added a
    screen-kNN" WAS WRONG -- the idea is theirs; the estimator is what is ours.**
    2026-10-27/28, user's request. Free (CPU) apart from a 3.4 GB download.

    **THE CORRECTION FIRST, because it was asserted repeatedly in this session and in
    the slide deck.** `github.com/Genentech/AssayLoop` ships TWO kNN baselines:
    `models/knn_gene_embedding.py` (gene-kNN, distance-weighted vote over the k=10
    nearest already-tested genes in GenePT space) and **`models/screen_knn.py`
    (`Screen-kNN`, the workbook's `classical` row at mean e_t 3.302)**. The claim that
    the paper's kNN was over GENES and ours over SCREENS is false: the paper has both,
    `screen_knn` is registered at `runner.py:105`, and `full_genome_table.py:185`
    already lists it as a table row. Nobody here had ever run either.

    **Standalone, on `paper_test`, in OUR adjusted-EF metric:**

    | baseline | EF | note |
    |---|---|---|
    | paper gene-kNN (GenePT) | 1.929 | no matching workbook row exists |
    | **paper Screen-kNN** | **3.401** | workbook e_t 3.302 -- consistent, ours forgives out-of-library picks |
    | AssayFormer alone | 4.828 | |

    **FUSED, on the 100-screen `validation` fold, DepMap at 0.5, paired against the
    same arm at weight 0, round-1 batches byte-identical:**

    | weight | OURS (§4.19) | paper Screen-kNN | paper gene-kNN |
    |---|---|---|---|
    | 0.25 | +0.143 | +0.011 [−0.027, +0.045] | -- |
    | **0.50** | **+0.175** | −0.064 [−0.159, +0.029] | −0.094 [−0.399, +0.226] |
    | 0.75 | +0.103 | **−0.499 [−0.811, −0.231]** | **−0.552 [−0.879, −0.213]** |
    | 1.00 | +0.043 | **−0.678 [−1.018, −0.366]** | **−0.841 [−1.193, −0.464]** |

    Ours is positive across the whole grid and peaks at 0.5. **Theirs never leaves
    zero** -- +0.011 with a ±0.04 interval and 45/100 wins is a coin flip, a tight null
    rather than an inconclusive one -- **and both of the paper's collapse as soon as
    they carry weight.** The difference is the estimator: theirs shrinks toward a
    global hit prior by `lam = n/(n+100)`, so early in a campaign it mostly reports
    "how often is this gene a hit in general", which is the fame axis `b_g` already
    encodes. Ours is hit-overlap over `sqrt(|H_s|)`, hard top-50, no prior.

    **THE DEFENSIBLE CLAIM, narrower than the old one:** the benchmark's nearest-screen
    baseline fused as a ranker is worth +0.011 [−0.027, +0.045] on 100 held-out
    screens; replacing its Bayesian-shrinkage estimator with hit-overlap plus a size
    correction makes the same idea worth **+0.175**. Measured from both directions.

    **THE CURATED 20s LIED AGAIN, and this time the trap was baited.** On `paper_test`
    the gene-kNN fused at `w=0.5` read **+0.155 [+0.028, +0.332], interval excluding
    zero, 13/20 screens** -- and the 100-screen fold reversed it to −0.094. I had also
    argued for trying HIGHER weights off that n=20 trend; `w=0.75` costs **−0.638**
    against the deployed method. Deploying on the curated result would have been a
    half-point regression.

    **AND THE ORTHOGONALITY ARGUMENT DID NOT PREDICT THE OUTCOME.** The gene-kNN picks
    an almost disjoint gene set -- **Jaccard 0.040** against AssayFormer's 1000 picks,
    where our method sits at 0.743 -- and **66% of the hits it finds are ones
    AssayFormer never reaches**, on 20 of 20 screens. That profile (weak but
    orthogonal) is exactly what justifies a fusion slot, and it still failed held out.
    **Finding genes the model misses is not the same as finding genes worth picking**
    -- the §4.49 oracle-seeds lesson, reproduced by a different mechanism. Do not treat
    low pick-overlap as evidence a ranker will help.

    Code: `codep_seed_acq.py` gains `gpt_npz`/`gpt_weight`/`gpt_k` (their gene-kNN) and
    `pknn_weight`/`pknn_train`/`pknn_tau`/`pknn_kappa` (their Screen-kNN), all default
    OFF and inert. `scripts/diagnostics/genept_to_npz.py` dumps GenePT to a 139 MB npz
    (22,641 x 1536, rows pre-normalised). `_GPT_CACHE` and `_PKNN_CACHE` are
    process-level and single-flight behind a lock: one acquisition is built PER SCREEN,
    so without them a 218-screen sweep at `--parallel 4` reads a 139 MB file 218 times
    and builds the 228 MB Screen-kNN state once per worker. Tests 76 -> 88. Arms
    `sweep-paperknn-pt`, `sweep-gpt*`, `sweep-pk*`, `sweep-gptv-*`, `sweep-pkv-*`.

    **FOUR TRAPS, all of the same kind: an arm that is not the thing it claims to be.**
    1. **`max_labels` must NOT apply to a kNN.** It caps the DepMap scorer's cost of
       averaging profiles; in a kNN it silently drops the earliest rounds from the
       NEIGHBOUR POOL, so the first `_gpt_scores` was a measurement of *a* GenePT kNN,
       not of theirs. (The corrected re-run moved +0.152 -> +0.155, so it cost nothing
       here -- but it was not knowable in advance.)
    2. **Pin a reimplementation against the class it reimplements.** `_pknn_scores` was
       tested against `ScreenKNNModel.predict` from the start; `_gpt_scores` was only
       tested against hand arithmetic until the gap was noticed. Both now have
       exactness tests.
    3. **Their gene-kNN cannot be matched bit-for-bit, and that is inherent.** With
       `epsilon = 1e-6`, a gene that is its own nearest neighbour takes weight
       `1/(1-1+eps) ~ 1e6`, so a float-level difference in a similarity near 1.0 is
       amplified a millionfold. The test asserts the RANKING (all RRF consumes) plus
       scores to 2e-4.
    4. **`ScreenKNNModel.predict` is unusably slow** -- 21k candidates x 1349 screens
       in Python, 0 of 20 screens completed before a session ended. Use the vectorised
       form in `plot_lopo_results.py` (`_build_screen_knn_matrices`,
       `_screen_knn_rollout`), which that file records as reproducing it exactly at
       ~800x. `_screen_knn_rollout` returns adjusted EF directly, so the standalone
       baseline needs no sweep at all.

    **DATA.** PRESAGE cache at `/home/jeff/presage/cache` (3.4 GB, Zenodo 15587986).
    `scripts/fetch_presage_cache.sh` defaults to a path INSIDE the repo that is not
    gitignored -- set `ASSAYLOOP_PRESAGE_CACHE`. One file fails to extract on this
    machine (`WeightedDeepset...Node2Vec...pkl`, filename too long for the encrypted
    home's eCryptfs limit) and `set -e` then aborts the script before it moves the
    subdirs up, which is why the cache root is `.../presage/cache` and not
    `.../presage`. GenePT is unaffected; coverage 97.9% of the gene universe.

    ---

    **ADDENDA, 2026-10-28: five measurements made after the above was written.**

    **(a) k = 50 IS CONFIRMED, and "all neighbours decayed" is just the far end of that
    curve.** Neighbourhood size swept on the 100-screen `validation` fold, paired against
    the deployed k=50, arms `sweep-knnk-k{010,025,100,200}`:

    | k | paired delta vs k=50 |
    |---|---|
    | 10 | −0.009 [−0.061, +0.043] |
    | 25 | −0.004 [−0.040, +0.031] |
    | 100 | **−0.100 [−0.188, −0.021]** |
    | 200 | **−0.119 [−0.223, −0.028]** |

    Flat from 10 to 50, significantly worse at 100 and 200. The deployed value was never
    tuned -- it was the first number tried -- and the sweep says that was luck in the
    sense that anything from 10 to 50 would have served equally.

    **(b) DECOMPOSITION: the neighbourhood rule carries the whole effect, and the two
    mechanisms I argued for are worth nothing.** Their estimator differs from ours in
    three places. Swapping each one ALONE into our arm, 100-screen fold:

    | swapped | theirs | ours | worth |
    |---|---|---|---|
    | global prior | blended by `lam = n/(n+100)` | none | +0.004 |
    | similarity | squared error over shared genes (misses count) | hit overlap / `sqrt(|H_s|)` | −0.057 |
    | **neighbourhood** | all 1349, exponentially decayed | **hard top-50** | **+0.204** |

    **This overturned my own explanation, twice stated above and in the deck.** I argued
    the difference was (i) their distance counting agreement on MISSES, which is cheap
    when hits are ~5% of genes, and (ii) their shrinkage toward a global hit prior
    repeating what `b_g` already encodes. Both sounded mechanistic. Removing the prior
    moves the result by four thousandths; switching the similarity to hits-only makes it
    slightly WORSE. Restricting the neighbourhood is the entire effect.

    **(c) THE DECOMPOSITION IS NOT ADDITIVE -- the transfer test refutes it.** If the
    neighbourhood rule were a transferable component, giving THEIR estimator a hard
    top-50 should recover most of the gap. Predicted +0.140; **measured −0.017**. So
    "restrict the neighbourhood" is not a part that can be bolted onto their machinery:
    it only pays in combination with our normalisation. **A residual +0.169 is still
    unexplained** and lives in the normalisation (`sum(w*Y)/sum(w)` against
    `(1-lam)*prior + lam*evidence`) plus the weight shape.

    **(d) THE TAU PREMISE WAS WRONG, and measuring it took two minutes.** I proposed
    tuning their `tau` on the theory that `exp(-dist/0.1)` is far too peaked and starves
    the vote. Measured across 21 screens: the distance spread over the top 50 neighbours
    is **0.120**, so `tau = 0.1` gives a **3.3x** weight ratio between the nearest and
    the 50th -- against **3.4x** for our linear weighting. **Their tau already matches
    ours.** Any `tau` above 0.1 would make their weighting FLATTER than ours, in a
    direction nothing supports. The weight shape is not the problem.

    **(e) SHRINKING TOWARD THE PRIOR WOULD PROMOTE MEASURABLY WORSE GENES, so the one
    principled fix left points the wrong way.** The remaining idea was per-gene
    shrinkage, `(sum(w*Y) + m*p0) / (sum(w*measured) + m)`, on the premise that a gene
    measured by few of the 50 neighbours is UNCERTAIN and our estimator treats it too
    harshly by shrinking it toward 0. Training density is 71.5%; of the 100 genes their
    estimator ranks highest, 20.2% are measured by fewer than 10 of 50 neighbours and
    4.2% score exactly 1.0. The premise is then one measurement, over 43 screens, on
    genes in their top 200:

    | coverage | hit rate |
    |---|---|
    | thin (< 10 of 50 neighbours measured it) | **7.3%** |
    | good (>= 10) | **29.8%** |
    | the screen's own base rate | 5.8% |

    **Thin coverage does not signal uncertainty, it signals a BAD PICK** -- 7.3% against
    a 5.8% base rate, four times worse than well-covered picks. Shrinking those toward
    the prior is well calibrated, which for a RANKING is worse than burying them: our
    harsher shrink-to-zero is the better policy and is presumably part of the +0.169.
    **Neither remaining variant was run**, and that is the finding: the normalisation
    swap converges to our arm by construction, and the shrinkage variant is now expected
    to underperform it. Two minutes of measurement replaced three arms.

    **NAMING, asked for explicitly.** Ours is **user-based collaborative filtering**
    (screens are users, genes are items, hits are ratings, cosine on binary vectors --
    the `sqrt(|H|)` of the current screen cancels, being equal for every candidate
    neighbour). Theirs is **`Screen-kNN`**, kernel-weighted regression with
    Bayesian shrinkage. Neither half is a new estimator; the contribution is the
    composition and the measurement.

    **A CORRECTION TO THE DECK'S OWN TABLE.** The `ours` column on the screen-kNN slide
    was quoted from §4.19 -- **198 screens including ones whose library is smaller than
    the budget, measured against AssayFormer ALONE**. Theirs was measured on 100
    genome-scale screens against AssayFormer + DepMap. Those are different comparisons
    and the mismatch flattered us. Rescored on the identical 100 screens against the
    identical baseline: **+0.061 / +0.087 / −0.093 / −0.138** at weights
    0.25 / 0.50 / 0.75 / 1.00. Ours is still better at every weight, but at 0.5 the
    interval crosses zero and the honest headline is "ours adds a little where theirs
    adds nothing", not the clean win the old table implied. **See §4.56(d): on the
    ablated atlas, which is now the default, that +0.087 becomes −0.022 (§4.59).**

    **A NAMING HAZARD THAT ALMOST COST A RESULT.** `sweep-nb-agree-k50` is a PREFIX of
    `sweep-nb-agree-k50-prior`, so `glob("sweep-nb-agree-k50-*")` matched 436 directories
    across both arms. The contaminated names parsed to an index of `00-1454` and were
    silently skipped, so the numbers came out right BY LUCK. **No sweep id may be a
    prefix of another.**

---

56. **THE ABLATED ATLAS IS NOW THE PROJECT DEFAULT. THE SWAP COSTS NOTHING OVERALL -- AND
    IT TAKES THE SCREEN-kNN'S INCREMENT TO ZERO ON THE FOLD THAT DECIDES THINGS.**
    2026-10-28, user's request. Free, CPU only, 654 campaigns in 22 minutes.

    **THE DECISION.** `achilles_noeval.npz` -- 768 of DepMap's 808 Achilles cell lines,
    every line used by ANY of the four evaluation sets dropped -- replaces
    `achilles.npz` as the default everywhere. §4.54 built it and measured it on one arm;
    this section adopts it. The reason is the CLAIM, not the score: "no screen's own cell
    line is in the atlas" stops being an audit result (§4.35) and becomes true by
    construction, so the leakage question cannot be reopened against any arm that uses
    it. The cost, measured below, is within a few hundredths on every set.

    **A REPRODUCIBILITY HOLE, CLOSED FIRST.** `achilles_noeval.npz` was about to become
    the default that every campaign reads, and **it had no builder in the repo** -- it
    was written by a script in a session scratchpad that no longer exists, so the default
    atlas could not be rebuilt or checked. `scripts/diagnostics/build_noeval_atlas.py`
    now rebuilds it, importing §4.35's resolver from `screen_depmap.py` rather than
    restating it so the two cannot drift. `--verify` reproduces the shipped file and
    asserts equality: **IDENTICAL**, 40 lines dropped, 154 screen slots touched
    (10 `paper_validation`, 7 `paper_test`, 56 `validation`, 81 `test`), 18119 genes
    non-degenerate. Two traps in writing it, both in the CHECK rather than the build:
    `np.array_equal` without `equal_nan` reports DIFFERS on an identical file, and plain
    `std` instead of `nanstd` flagged 1936 genes as constant -- the matrix carries 28922
    NaNs that the acquisition `nan_to_num`s at load.

    **ARMS RE-RUN.** Only the DepMap users inside the notebook's ten-arm focus, as
    asked. Our method was already done (`sweep-nolk-*`, §4.54).

    | method | new arms | sets |
    |---|---|---|
    | DepMap + screen-kNN | `sweep-nolk-{pv,pt,validation,test}` *(already on disk)* | all four |
    | DepMap co-essentiality only | `sweep-dmne-cd-{pv,pt,validation}` | three |
    | + BioPlex | `sweep-dmne-ppi-{pv,pt,validation}` | three |
    | GLM-5 warm start -> ours | `sweep-dmne-hyb-glm5-{pv,pt}` | two |
    | GLM-5 + top-up -> ours | `sweep-dmne-hyb-glm5tu-{pv,pt}` | two |
    | *full-atlas CPU twins, for the paired baseline* | `sweep-dmcdbase-{pv,pt,validation}` | three |

    `sweep-dmcdbase-*` exist because the registered full-atlas codep-only arms were
    unusable as a paired baseline: `sweep-dm-pt-dm05` and `sweep-dm-validation-dm05` ran
    on **cuda**, and `sweep-cdv-w05` ran an **older recipe** (no rrf, no `train_hits`, no
    `query`) with a `profile_npz` pointing into a dead scratchpad. Without the twins the
    delta would have been confounded by device and by recipe, not just by atlas.

    **ROUND 1 IS STILL THE CHECK (§4.54 failure 3), and it passed 14 pairs of 14** --
    every new arm byte-identical to its twin in round 1, no universe-size or device
    mismatch anywhere. The hybrids pass trivially (round 1 is the replayed LLM batch).

    **(a) OUR METHOD -- reproduces §4.54 exactly, which is the consistency check:**

    | set | n | full | ablated | paired delta *(re-scored, §4.59)* |
    |---|---|---|---|---|
    | paper_validation | 20 | 6.217 | 6.231 | +0.008 [−0.035, +0.052] |
    | paper_test | 20 | 5.388 | 5.362 | −0.026 [−0.056, −0.001] |
    | validation | 100 | 5.054 | 5.017 | **−0.043 [−0.087, −0.006]** |
    | test | 248 | 6.896 | 6.922 | +0.019 [−0.010, +0.050] |

    **(b) THE OTHER THREE METHODS -- nothing of consequence moves:**

    | method | set | n | full | ablated | paired delta |
    |---|---|---|---|---|---|
    | DepMap only | paper_validation | 20 | 6.176 | 6.184 | +0.019 [−0.021, +0.064] |
    | DepMap only | paper_test | 20 | 5.312 | 5.532 | **+0.220 [+0.008, +0.539]** |
    | DepMap only | validation | 100 | 4.967 | 5.014 | **+0.059 [+0.019, +0.110]** |
    | + BioPlex | paper_validation | 20 | 6.292 | 6.295 | +0.005 [−0.061, +0.073] |
    | + BioPlex | paper_test | 20 | 5.348 | 5.395 | +0.048 [−0.026, +0.162] |
    | + BioPlex | validation | 100 | 5.045 | 5.038 | −0.014 [−0.065, +0.024] |
    | GLM-5 -> ours | paper_validation | 20 | 5.359 | 5.358 | +0.010 [−0.017, +0.045] |
    | GLM-5 -> ours | paper_test | 20 | 5.182 | 5.200 | +0.018 [−0.035, +0.063] |
    | GLM-5 +tu -> ours | paper_validation | 20 | 5.431 | 5.442 | +0.011 [−0.029, +0.053] |
    | GLM-5 +tu -> ours | paper_test | 20 | 5.227 | 5.251 | +0.025 [−0.028, +0.095] |

    **The one eye-catching number, `DepMap only` +0.220 on `paper_test`, is a curated-20
    artefact and should not be quoted.** 10/20 wins -- the median screen does not move --
    and the same contrast on the 100-screen fold reads +0.047. This is the §4.55(b)
    trap in its third appearance: an n=20 interval excluding zero on a mean that a few
    screens carry.

    **(c) THE HEADLINE SURVIVES.** On `validation`, our method against AssayFormer alone
    is **+0.443 [+0.247, +0.687]**; on `test`, **+0.506 [+0.159, +0.945]**. The atlas
    swap does not touch the project's result.

    **(d) WHAT IT DOES TOUCH, AND THIS IS THE REAL FINDING OF THE SECTION: the
    screen-kNN's increment over DepMap alone goes to ZERO on the default fold.** Measured
    with both arms on the SAME atlas, so the contrast is the kNN and nothing else:

    | set | n | full atlas | ablated atlas *(re-scored, §4.59)* |
    |---|---|---|---|
    | paper_validation | 20 | −0.033 [−0.217, +0.114] | −0.044 [−0.220, +0.094] |
    | paper_test | 20 | +0.077 [−0.069, +0.210] | **−0.170 [−0.542, +0.110]** |
    | validation | 100 | **+0.080 [−0.048, +0.253]** | **−0.022 [−0.152, +0.128]** |

    **Both intervals always crossed zero, so nothing here was ever significant** -- the
    honest reading is that the kNN half's increment was never distinguishable from zero
    on this fold and its point estimate is now exactly that. But the deck's screen-kNN
    slide quotes `+0.087` as the number that makes our estimator worth having over the
    paper's, and on the default atlas **that number is +0.003**. §4.55's closing claim
    ("ours adds a little where theirs adds nothing") now rests on +0.003 against +0.011,
    which is no claim at all. **The slide needs a decision, not a patch** -- either
    re-measure the whole weight grid on the ablated atlas (3 arms x 218 screens, ~15
    minutes, free) or restate the comparison as a null on both sides. NOT DONE.

    **(e) THE PAPER'S Screen-kNN IS NOT NESTED IN THE BASELINE, AND OURS IS. The §4.55
    comparison was never symmetric.** Found while costing the re-measurement above.
    `_pknn_scores` has **no `min_seeds` gate** -- line 992 fires it whenever
    `pknn_weight != 0.0` -- so with zero observations every training screen sits at
    distance 0, `w_s` is uniform, `lam = n/(n+kappa) = 0`, and the score is the **pure
    global hit prior**. Verified both ways: the trace reads `pknn_active: True` at
    round 1 where every one of our scorers reads `False`, and the round-1 batch differs
    from the DepMap-only baseline on **218 of 218 screens even on matched CPU**
    (`sweep-pkv-w025` vs `sweep-dmcdbase-validation`, 0/218 identical), while our column
    is **218/218 identical** (`sweep-wt-validation-c05k025` vs `sweep-dm-validation-dm05`).

    So §4.55's two columns do not measure the same thing. Ours is "what the kNN term adds
    given identical early picks". Theirs is a different campaign from pick 1, which spends
    round-1 budget on famous genes before any information exists -- and that is a plausible
    mechanism for its collapse at high weight (−0.499 at 0.75, −0.678 at 1.00) that has
    nothing to do with its estimator. **This is faithful to `ScreenKNNModel`** -- a
    standalone model with no data SHOULD predict the prior -- so it is not a
    reimplementation bug. It is a fusion-harness asymmetry, and it means **the round-1
    invariant cannot validate any `pknn` arm**, which is why none of them was ever caught.

    **A SECOND, SMALLER FINDING FROM THE SAME CHECK: device changes round 1.**
    `sweep-dmcdbase-validation` (cpu) vs `sweep-dm-validation-dm05` (cuda) is round-1
    identical on only **162 of 218** screens. §4.54 recorded the device as worth
    −0.001/+0.003 in EF, which is true and remains the reason not to worry about it, but
    it is NOT byte-identical and 56 screens reorder. **Never pair across devices and call
    the result nested.** This is exactly why today's 14 pairs were all cpu-vs-cpu.

    **Why the two halves stop being additive is open and worth a line in §10.** DepMap
    alone GAINS +0.047 from the ablation while the fused arm LOSES −0.036; the kNN's
    ranking is unchanged by construction, so the fusion is absorbing the atlas change
    differently from the atlas scorer alone. That is the same shape as §4.48-51's
    near-substitutes result -- two sources that each add less once the other is present.

    **WHAT WAS NOT RE-RUN, and it is flagged rather than hidden.** 15 arms still read the
    full atlas: the three richer-description arms, BioPlex-degree-residualised, co-hit,
    expression, both Flash-Lite hybrids, and the oracle/injection bounds. The notebook now
    carries an `atlas` column per arm, `list_bars()` marks the full-atlas ones, and the
    master-table cell prints a warning if the current selection mixes matrices. With the
    ten-arm focus it reports **atlas-consistent: every atlas user is 'ablated'**.

    **NOTEBOOK (`output/runs/compare_all_results.ipynb`), executed clean, 22 figures (19 at the time of this section; Figures 11-13 were added in §62).**
    `SELECT` now holds the user's ten arms and feeds `DEFAULT_SELECT`, so every figure
    follows one switch; terms carry no sweep ids so each matches its method on all four
    sets. Three fixes fell out of doing it:
    1. **`pick()` conflated two kinds of miss.** One global `SELECT` printed
       `!! nothing matched` **40+ times per run** -- every figure draws from its own
       subset (Figure 2 drops the paper's arms for having no paired delta, the per-round
       figures for having no curve), so a term legitimately missing from nine panels
       looked identical to a typo. Now a term matching nothing ANYWHERE is still loud; a
       term matching real arms absent from THAT panel is silent. A `_known_term`
       self-check in the cell asserts the typo path still fires, because a check that
       cannot fire reads as a guarantee.
    2. **Two BioPlex arms carried the stem `Our method + BioPlex protein interactions`**,
       which is outside `_STEM_ORDER`, so the same method drew in a different colour on
       different facets -- the one rule a multi-panel chart most often breaks, and the
       one the colour map exists to enforce. All three now share `+ BioPlex protein
       interactions`.
    3. **`Our method + BioPlex, famous-protein effect removed` never matched its
       `DESCRIPTIONS` key**, so the notebook printed `!! no description for:` on every
       run. Pre-existing; renamed to the key.

    **Code:** `build_noeval_atlas.py` (new). `ATLAS=${ATLAS:-achilles_noeval.npz}` in
    `run_{weight,desc,geometry,ppi,cohit}_sweep.sh`,
    `run_codep_{confirm,diffmeans,twofeat}.sh`, `run_curated_backfill.sh`,
    `run_testfold_confirm.sh`. `--depmap` / `--profile-npz` defaults flipped in
    `run_oracle_seeds.py`, `run_inject_2x2.py`, `handoff_fullgenome.py`. **No change to
    `codep_seed_acq.py`**: `profile_npz` has no default there and should not gain one --
    an arm must say which atlas it read. Suite unchanged at 308 passing; the 2 failures
    (`test_cli_full_genome`, `test_gemini_client`) are pre-existing and touch none of
    this.

---

57. **THE SYMMETRIC GRID SETTLES §4.55: ON THE DEFAULT ATLAS NEITHER ESTIMATOR ADDS
    ANYTHING. Ours is consistently the less harmful of the two, and the nesting flaw I
    found and fixed was worth 0.016.** 2026-10-28, user's request ("run the 8 arms").
    Free, CPU only, 8 arms x 218 screens. Resolves §10 item 4b in favour of its CHEAP
    option: restate §4.55 as a null on both sides.

    **THE SETUP.** `knn_weight` (ours) or `pknn_weight` (theirs) over
    {0.25, 0.50, 0.75, 1.00}, DepMap pinned at `weight=0.5`, ABLATED atlas, `validation`
    fold, every arm paired against `sweep-dmne-cd-validation` (AF + DepMap, both kNN
    terms off, EF 5.014) **on the same atlas**, plus the ungated control at 0.50.

    **ROUND 1 PASSED EXACTLY AS SPECIFIED, which is what makes the table readable:** all
    7 gated arms 218/218 identical to the baseline, the ungated control 0/218. The gate
    does what §4.56(e) said it would.

    | weight | OURS | theirs, GATED *(re-scored, §4.59)* |
    |---|---|---|
    | 0.25 | −0.014 [−0.063, +0.039] 47/100 | **−0.059 [−0.108, −0.019]** 38/100 |
    | 0.50 | −0.022 [−0.152, +0.128] 43/100 | **−0.179 [−0.304, −0.070]** 35/100 |
    | 0.75 | **−0.170 [−0.345, −0.003]** 38/100 | **−0.375 [−0.574, −0.197]** 37/100 |
    | 1.00 | **−0.217 [−0.412, −0.037]** 32/100 | **−0.673 [−0.973, −0.400]** 34/100 |

    **(a) OUR HALF IS NOT NEUTRAL ON THE DEFAULT ATLAS, IT IS MILDLY HARMFUL.** Every
    cell is negative: the best is −0.014 at 0.25 with an interval spanning zero, turning
    significantly negative by 0.75. The full-atlas `+0.087` the deck was built on was this
    quantity's high-water mark and it did not survive the atlas change.
    **`+0.087 -> −0.022`: the claim "our formulation adds a little where theirs adds
    nothing" is dead, and so is the weaker "both add nothing" -- fused as a standalone
    ranker on `validation`, ours costs a little and theirs costs more.**

    **(b) THE ORDERING SURVIVES, AND IT IS ALL THAT SURVIVES.** Ours beats theirs at
    every one of the four weights, and the gap widens with dose (0.045 / 0.157 / 0.205 /
    0.456). Theirs is significantly negative at EVERY weight; ours only from 0.75. So the
    defensible sentence is **"neither estimator earns a fusion slot on this fold, and
    theirs degrades faster under weight"** -- an ordering between two negatives, not a
    win.

    **(c) THE GATE WAS WORTH 0.026, AND THAT IS THE MOST USEFUL NUMBER HERE.**
    Ungated −0.205 [−0.337, −0.074] against gated −0.179 [−0.304, −0.070], same atlas,
    same weight, same screens.

    §4.56(e) was a real **VALIDITY** defect -- their arms were not nested, their round-1
    batch differed from the baseline on 218 of 218 screens, and the round-1 invariant
    could not police them. Fixing it moves the result by sixteen thousandths. **So the
    round-1 famous-gene prior was NOT what sank their arm.** Their estimator is simply
    worse under weight.

    **This is the §4.55(b)/(c) pattern for the third time: a mechanistically compelling
    story, carefully argued, worth essentially nothing when measured.** First the global
    prior (+0.004) and the miss-counting similarity (−0.057); then the neighbourhood
    rule, which carried +0.204 alone but failed to transfer (−0.017 against a predicted
    +0.140); now the round-1 gate (+0.026). **The mechanism being real does not make it
    load-bearing, and in this codebase it usually is not.** Fix validity defects because
    they are validity defects, and price them separately from the result.

    **(d) WHAT THIS DOES TO THE DEPLOYED METHOD, and it is the open question now.** On
    the ablated atlas the deployed arm (5.017) and DepMap alone (5.014) are
    indistinguishable on this fold, and §4.56(b) has DepMap alone GAINING +0.047 from
    the ablation. So **on the default atlas the deployed two-term method is, on
    `validation`, DepMap alone with a term that costs a little attached.** Do NOT drop the
    kNN half on this evidence: it reads +0.047 on `paper_validation` and −0.170 on
    `paper_test`. **ANSWERED in (e) below, and the answer is KEEP BOTH:** the `test`
    fold puts the increment at +0.191 and dropping the term would cut that fold's
    headline from +0.506 [+0.159, +0.945] to +0.315 [−0.025, +0.747] -- from significant
    to not.

    **(e) THE TEST FOLD ANSWERS (d) AND REVERSES THE PARSIMONY CASE. KEEP BOTH TERMS.**
    2026-10-28, same session, user's request. `sweep-dmne-cd-test` (334 dirs, free);
    round 1 **334/334 identical** to `sweep-nolk-test`, same 22656-gene universe, both
    cpu, so the only difference between them is `knn_weight`.

    | contrast, `test` fold, n=248, 128 families | paired delta *(re-scored, §4.59)* |
    |---|---|
    | AF + DepMap vs **AF alone** | +0.315 [−0.025, +0.747] 107/248 |
    | AF + DepMap + screen-kNN vs **AF + DepMap** | **+0.191 [−0.019, +0.415]** 99/248 |

    **The screen-kNN increment across all four sets, every interval crossing zero:**

    | set | n | increment over AF + DepMap |
    |---|---|---|
    | paper_validation | 20 | −0.044 [−0.220, +0.094] |
    | paper_test | 20 | −0.170 [−0.542, +0.110] |
    | validation | 100 | −0.022 [−0.152, +0.128] |
    | **test** | **248** | **+0.191 [−0.019, +0.415]** |

    **So the half is UNRESOLVED, and NEGATIVE on three of four sets** -- point estimates
    from −0.170 to +0.191, no interval excluding zero, and only `test` positive. §4.57(d) read the `validation` null as grounds for calling
    the method one component; **the test fold refutes that**, and the decisive number is
    what each method is worth against AssayFormer alone on that fold:

    | arm, `test` fold | vs AF alone |
    |---|---|
    | AF + DepMap | +0.315 [−0.025, +0.747] *(crosses zero)* |
    | **AF + DepMap + screen-kNN** | **+0.506 [+0.159, +0.945]** |

    **Dropping the kNN half would cut the headline on `test` from +0.506 to +0.315 and
    move the interval from clear of zero to CROSSING it.** On `validation` the two are
    equal (+0.443 vs +0.466). So parsimony buys nothing on one fold and costs the
    result's significance on the other: **ship both terms, and describe the
    second as unresolved rather than as a contributor.**

    **AND IT IS THE §4.46 FOLD DISAGREEMENT AGAIN -- the fourth guise.** Cell-line
    ablation helps `test`, hurts `validation` (§4.54a). Proximity correction, same
    (§4.54b). Now the screen-kNN increment: −0.022 on `validation`, +0.191 on `test`.
    Note also that the DepMap term itself is weaker on `test` (+0.315) than on
    `validation` (+0.466) while the FUSED arm is stronger there (+0.506 vs +0.443) --
    the two terms trade places between folds. **Whatever drives the fold split is now
    visible in the headline method and not just in diagnostics, which promotes §10's
    open question from interesting to load-bearing for the write-up.**

    **PROCESS. Three failures, and the first is the one to remember.**
    1. **THE GATE BROKE THE UNIVERSE, and memory was only the symptom.**
       `_pknn_universe` is pinned at the scorer's FIRST call from `candidates`, which is
       the UNTESTED pool and shrinks 100 genes a round. That was safe only because the
       first call was always round 1. Gating moved it to round 2-4, where the acquired
       genes differ per screen -- so every screen got its own cache key, its own ~230 MB
       state (OOM-killed at 7.4 GB after ~32 of 218 screens), **and a global hit prior
       built over "whatever was still untested when the gate opened" rather than over a
       fixed universe.** On a bigger box it would have completed and produced a
       clean-looking table that was not the measurement claimed. The OOM was the only
       thing that caused anyone to look. Fixed by reconstructing the universe as
       `candidates | history`, which rests on "a gene leaves `candidates` only by being
       acquired"; `_PKNN_CACHE` is now capped at 2 with a warning naming this invariant.
    2. **I COSTED MEMORY BY ESTIMATE TWICE AND WAS WRONG BOTH TIMES.** "~228 MB per
       worker" ignored that `-j` is a THREAD pool, then a 4 GB alert threshold picked
       without a baseline fired on a legitimate startup transient and I killed a healthy
       run, calling it a persisting leak. **Measured:** a one-time ~1.6 GB when the first
       pknn call loads 1349 training screens and builds two 114 MB matrices, a transient
       to ~5.0 GB while N threads each load a model and the atlas before that shared
       state is warm, then a FLAT plateau at ~3.4 GB; driven in-process over 4 screens
       the scorer is flat at 2129 MB. There is no leak. `-j 3` for any pknn arm.
    3. **A PARTIAL ARM SCORED AS A FINISHED ONE.** `run()` skipped on "any run dirs
       exist", so `sweep-g8-pk-w025` passed with **18 of 218** directories left by the
       killed run and would have entered the table as an n=18 number with nothing
       complaining. The launcher now requires the exact expected count and clears
       anything short. **This is the §4.54 "an arm that is not what it claims to be"
       family again -- and the only member of it so far that fails SILENTLY.**

    Code: `pknn_gate` on `CodepSeedAcquisition`, **default `True`** -- a deliberate
    default CHANGE, since ungated is a measurement artefact for a fused term rather than
    a feature; `pknn_gate=False` reproduces the pre-2026-10-28 `sweep-pkv-*` arms
    exactly. `_pknn_universe` reconstruction, `_PKNN_CACHE_MAX`. Tests 96 -> 104 in
    `test_codep_seed.py` (the gate: round 1 greedy, `False` restores, opens at
    `min_seeds`, honours a custom `min_seeds`, leaves `_pknn_scores` untouched; the
    universe: identical whichever round it first runs in, cache stays at one entry).
    Suite 313 passing, the 2 failures pre-existing and unrelated. Arms
    `sweep-g8-ours-k{025,075,100}`, `sweep-g8-pk-w{025,050,075,100}`,
    `sweep-g8-pkung-w050`.

---

58. **BIOPLEX ON THE REDUCED ATLAS: SMALL, CONSISTENTLY POSITIVE, AND IT PASSES THE FAME
    CONTROL THAT KILLED EVERYTHING ELSE. The first candidate in a while that deserves a
    test-fold run.** 2026-10-28, user's request. Free, CPU, 476 campaigns. Round 1
    identical to `sweep-nolk-*` on all four new arms (20/20, 20/20, 218/218, 218/218),
    so the only moving part against that baseline is the BioPlex term.

    Plain BioPlex at 0.25 was already on the reduced atlas from the §4.56 SELECT re-runs
    (`sweep-dmne-ppi-*`). What was still full-atlas and is re-run here: the
    degree-residualised control on all three sets, and the 0.50 dose on `validation`.

    **All seven contrasts, over AF + DepMap + screen-kNN, reduced atlas:**

    | arm | paper_validation (20) | paper_test (20) | validation (100) |
    |---|---|---|---|
    | BioPlex 0.25 | **+0.064 [+0.008, +0.143]** | +0.033 [−0.019, +0.088] | +0.021 [−0.015, +0.060] |
    | BioPlex 0.25, degree-residualised | **+0.063 [+0.006, +0.143]** | +0.045 [−0.013, +0.107] | +0.027 [−0.010, +0.067] |
    | BioPlex 0.50 | -- | -- | +0.031 [−0.018, +0.082] |

    **(a) THE FAME CONTROL PASSES, AND THAT IS THE RESULT.** §4.11's fame redundancy --
    gene-level knowledge already living in the BPMF factors -- is the single mechanism
    that explained every other null in this project, and BioPlex's own gate 1 reads
    +0.203 because hub proteins are better studied. So the obvious reading was that any
    BioPlex gain is fame. **It is not.** Residualising on log degree leaves the effect
    untouched: +0.024 -> +0.030 on `validation`, +0.030 -> +0.036 on `test`, +0.090 ->
    +0.090 on `paper_validation`, and priced directly on `test`, dr against plain is
    **+0.006 [−0.003, +0.020]**. Whatever the term contributes, it is not protein fame,
    **and the control passes on BOTH held-out folds. This is the first candidate data
    source to survive that control at all.**

    **(b) EVERY CONTRAST IS POSITIVE, AND IT DOES NOT FLIP BETWEEN FOLDS.** Nine of nine
    positive: +0.024, +0.030, +0.030, +0.031, +0.033, +0.036, +0.045, +0.090, +0.090.
    **The two held-out folds AGREE -- +0.024 on `validation` against +0.030 on `test`** --
    which nothing else in this project manages: the screen-kNN half reads −0.022 and
    +0.191 on those same two folds, and the atlas ablation flips sign across them. On
    `validation`
    the absolute ordering runs BioPlex 0.50 (5.048) > dr (5.044) > plain (5.038) > no
    BioPlex (5.017) > DepMap alone (5.014). **Read that consistency carefully: these
    arms are nested and share screens, so nine positive contrasts are not nine
    independent confirmations** -- they are one small effect measured nine overlapping
    ways. **No interval excludes zero on any set.** The honest statement is a consistent
    sign at 0.02-0.09, agreeing across both held-out folds, every interval crossing zero.

    **(c) `paper_validation` IS THREE TIMES THE FOLD ESTIMATE AND SHOULD BE DISCOUNTED.**
    +0.090 [−0.002, +0.230] on 20 screens / 11 families is the same shape as §4.55(b)'s
    gene-kNN trap (+0.155 excluding zero on the curated 20, reversing to −0.094 on the
    fold). Lead with the folds, which agree at +0.024 / +0.030.

    **(d) THE DECIDING ARM RAN THE SAME DAY, AND IT AGREED.** `sweep-dmne-ppi-test` and
    `sweep-dmne-ppidr-test`, 334 dirs each, free, round 1 identical to `sweep-nolk-test`
    on both. It was the project's cheapest open question because §4.57(e) had the
    screen-kNN half going −0.022 on `validation` to +0.191 on `test`, and §4.54(a) has
    the atlas ablation flipping sign across the same pair. **BioPlex does not flip.**
    Against AssayFormer alone on `test`, the BioPlex arms read **+0.536 [+0.175, +0.985]**
    plain and **+0.542 [+0.182, +0.991]** residualised, against the two-term method's
    +0.506 [+0.159, +0.945].

    **WHAT TO DO WITH IT, stated carefully, because it is easy to oversell.** This is not
    a result: no interval excludes zero anywhere, and +0.03 against a method effect of
    +0.44 to +0.51 is a rounding error on the headline. It is the **best-behaved
    candidate the project has** -- positive on all four sets, agreeing across both
    held-out folds, and immune to the one mechanism (§4.11 fame redundancy) that explains
    every other null here. The next question is not another weight or another set; it is
    whether a fold-consistent +0.03 justifies a third data dependency in the deployed
    method. **That is a judgement for the write-up, not a measurement.**

    Arms `sweep-dmne-ppidr-{pv,pt,validation}`, `sweep-dmne-ppi05-validation`. The
    notebook's BioPlex-dr registry entry is repointed from `sweep-ppi-validation-w025dr`
    (full atlas) to the reduced-atlas arm.

---

59. **A BUG IN `score_pair.py` MADE EVERY "FAMILY-BOOTSTRAPPED" DELTA A PLAIN SCREEN MEAN
    WITH TOO-NARROW INTERVALS. All of §4.54, §4.56, §4.57 and §4.58 are re-scored; two
    conclusions moved, none reversed.** 2026-10-28, found by reconciling a figure against
    the notebook.

    **THE BUG.** `score_ctx.families()` returns a **dict**, `{screen_name: family_index}`.
    `score_pair.boot()` consumed it as if it were a list of lists:

        fams = families(load_screens(target_set=sset)); fam = {}
        for i, f in enumerate(fams):        # yields dict KEYS -- screen name strings
            for s in f:                     # iterates each NAME'S CHARACTERS
                fam[s if isinstance(s, str) else s.dataset_name] = i

    So `fam` came out with **22 single-character keys** (`{'1': 325, '7': 333, 'U': 331,
    ...}`), every real screen name missed the lookup, and `fam.get(s, s)` fell back to the
    screen's own name -- **making every screen its own family**. Two effects on every
    number the script produced:
    1. the delta was a **plain per-screen mean**, never family-averaged;
    2. the bootstrap resampled **248 screens where 128 independent families exist**, so
       **the intervals were too narrow.**

    `paper_test` is unaffected -- 20 screens, 20 families, nothing to cluster -- which is
    part of why it went unnoticed. `paper_validation` (11 families / 20 screens),
    `validation` (76 / 100) and `test` (128 / 248) all move.

    **HOW IT SURFACED, and this is the transferable part.** Not from inspection -- from a
    figure that would not reconcile. The test-fold BioPlex headline came out +0.356 from
    `score_pair` and +0.506 from the notebook, for **the same arms on the same screens
    against the same baseline**. The notebook was right: `boot_families()` uses
    `families_for()` correctly. **Two independent implementations of the same quantity is
    what caught this, and nothing else would have** -- the wrong numbers were individually
    plausible, internally consistent across a whole session, and reproduced §4.54's
    recorded table exactly, because that table came from the same bug.

    `boot()` now consumes the dict directly, `n` is 20000 to match the notebook, the
    family count is printed in every row, and a guard aborts if the map resolves fewer
    than 90% of screens.

    **WHAT MOVED. Corrected figures, all on the reduced atlas unless stated.**

    (a) **The screen-kNN increment over AF + DepMap flipped sign on two sets** and the
    fold disagreement sharpened:

    | set | n | was | **corrected** |
    |---|---|---|---|
    | paper_validation | 20 | +0.047 | **−0.044 [−0.220, +0.094]** |
    | paper_test | 20 | −0.170 | −0.170 [−0.542, +0.110] |
    | validation | 100 | +0.003 | **−0.022 [−0.152, +0.128]** |
    | test | 248 | +0.145 | **+0.191 [−0.019, +0.415]** |

    (b) **§4.57(e)'s "keep both terms" SURVIVES and strengthens.** Against AssayFormer
    alone on `test`: DepMap alone **+0.315 [−0.025, +0.747]**, which now CROSSES ZERO,
    against the two-term method's **+0.506 [+0.159, +0.945]**. Dropping the kNN half
    takes that fold's headline from significant to not.

    (c) **§4.57's grid is uniformly worse, and ours is now negative at every weight on
    `validation`:**

    | weight | OURS | theirs, gated |
    |---|---|---|
    | 0.25 | −0.014 [−0.063, +0.039] | −0.059 [−0.108, −0.019] |
    | 0.50 | −0.022 [−0.152, +0.128] | −0.179 [−0.304, −0.070] |
    | 0.75 | −0.170 [−0.345, −0.003] | −0.375 [−0.574, −0.197] |
    | 1.00 | −0.217 [−0.412, −0.037] | −0.673 [−0.973, −0.400] |

    §4.57's conclusion is unchanged in kind and firmer in degree: **neither estimator adds
    anything, ours is better at every weight, and on `validation` ours is now mildly
    HARMFUL rather than neutral.** The gate's worth moves +0.016 -> **+0.026**
    (ungated −0.205, gated −0.179) -- still negligible, still not what sank their arm.

    (d) **§4.56's atlas-swap table.** `test` no longer excludes zero:

    | set | n | was | **corrected** |
    |---|---|---|---|
    | paper_validation | 20 | +0.014 | +0.008 [−0.035, +0.052] |
    | paper_test | 20 | −0.026 | −0.026 [−0.056, −0.001] |
    | validation | 100 | −0.036 | **−0.043 [−0.087, −0.006]** |
    | test | 248 | +0.026 *(excluded zero)* | **+0.019 [−0.010, +0.050]** *(crosses)* |

    **This weakens §4.54(a)'s headline claim**, which was that the ablation's intervals
    exclude zero in OPPOSITE directions on the two held-out folds. Only the `validation`
    side still excludes zero. The fold DISAGREEMENT in sign survives; the claim that both
    sides are individually significant does not. **§4.54(b)'s proximity-correction table
    is from the same scorer and has NOT been re-scored** -- treat its intervals as too
    narrow until it is.

    (e) **§4.58's BioPlex result survives in full, and the pattern I was about to report
    was an artefact of the bug.** With the buggy scorer `test` read +0.002 and the four
    sets fell in a tidy monotone sequence with sample size (+0.064, +0.033, +0.021,
    +0.002) -- which I was one message away from writing up as the signature of a
    small-sample artefact. **Corrected, there is no such sequence and `test` does not go
    to zero.** See §4.58's own table, now rewritten.

    **THE LESSON, and it is not "be careful".** Every number this scorer produced for two
    sessions was plausible, self-consistent, and agreed with the project's own recorded
    history -- because that history came from the same code. The only thing that exposed
    it was computing one quantity two independent ways and insisting they reconcile.
    **Where a figure carries a conclusion, compute it twice from different code.** The
    notebook and `score_pair` must now be cross-checked on at least one shared contrast
    whenever either changes.

---

60. **THE DEPMAP WEIGHT IS RE-CONFIRMED AT 0.5 ON REDUCED-ATLAS EVIDENCE ALONE, with a
    coherent inverted-U dose response.** 2026-10-28, user's request. Free, CPU, 654
    campaigns. Round 1 **218/218** identical to `sweep-nolk-validation` on all three arms,
    same 21213-gene universe.

    **WHY IT NEEDED RUNNING AT ALL.** Every reduced-atlas arm pinned `weight` at 0.5, and
    that value's only support was §4.19/§4.20 **on the full atlas** -- which §8's
    reduced-atlas-only rule puts out of bounds. It is also the term carrying the effect,
    and §4.56(b) had shown the ablation redistributing value between terms (DepMap-only
    GAINED +0.059 on `validation` from the swap while the fused arm LOST −0.043), so the
    optimum could not be assumed to have stayed put. `knn_weight` held at 0.5, no BioPlex,
    so only the DepMap dose moves. Full-atlas counterparts, recorded and NOT compared
    against: `sweep-wt-validation-c025k05`, `sweep-wt-validation-c075k05`.

    | DepMap weight | vs AssayFormer alone | vs the deployed 0.50 |
    |---|---|---|
    | 0.25 | +0.245 [+0.114, +0.400] | **−0.198 [−0.362, −0.076]** |
    | **0.50** *(deployed)* | **+0.443 [+0.253, +0.689]** | -- |
    | 0.75 | +0.467 [+0.233, +0.732] | +0.023 [−0.112, +0.162] |
    | 1.00 | +0.253 [+0.018, +0.505] | **−0.190 [−0.378, −0.014]** |

    **0.5 STANDS, and the reasoning is deliberately not "it scored highest" -- 0.75 did,
    by 0.024.** Three things decide it: 0.75 is statistically indistinguishable from 0.5
    (+0.023, interval spanning zero); BOTH neighbours, 0.25 and 1.00, are significantly
    WORSE than 0.5 with intervals excluding zero; and the pre-stated rule for this sweep
    was to move off 0.5 only on an interval that excludes it. **Picking 0.75 off a 0.024
    lead across a four-point grid with ±0.14 intervals is fitting noise** -- the §4.55(b)
    trap in its cheapest form.

    **THE CURVE IS COHERENT, WHICH IS ITSELF EVIDENCE.** An inverted U peaking in
    0.5-0.75, falling off by about 0.2 in both directions. §4.19 counted "coherent
    dose-response with the same optimum on both screen sets" among the reasons the
    positive holds; this reproduces that shape on the reduced atlas, which is a
    replication rather than a new result. Reading the two tails: at 0.25 the DepMap term
    is too quiet to move the ranking much (+0.245, roughly half the deployed gain), and at
    1.00 it starts overwhelming AssayFormer's own ordering (+0.253). On THIS fold the
    deployed value sits at the top of a broad, flat plateau -- **but (b) below shows that
    plateau is an artefact of `validation`, and the optimum is sharper than it looks.**

    **WHAT IS STILL UNSWEPT ON THE REDUCED ATLAS**, so the next person does not assume
    otherwise:
    - **`weight` at 0.25 and 1.00 on `test`** -- only 0.75 was run there, see (b). Not
      needed to confirm 0.5, since 0.75 was the only live alternative, but the shape of
      the curve on `test` is known at one point only.
    - **`knn_weight` and `ppi_weight` grids are also `validation` only.** Deliberately not
      extended (§4.58): both terms' increments are within noise on `test`, so a weight
      optimum fitted there would be fitting noise.
    - **`gpt_weight`** (the paper's gene-kNN) has never been run on the reduced atlas at
      any weight. It was refuted on the full atlas (§4.55) and is not worth revisiting.

    **(b) 0.75 ON THE `test` FOLD IS SIGNIFICANTLY WORSE, AND THIS IS THE CLEANEST
    VINDICATION OF NOT CHASING POINT ESTIMATES THE PROJECT HAS.** 2026-10-28, user's
    request, `sweep-dmwtest-c075`, 334 dirs, free. Round 1 **334/334** identical to
    `sweep-nolk-test`, same 22656-gene universe. 0.75 was the only dose on `validation`
    that 0.5 did not beat, so it was the only live alternative worth testing on the
    second fold.

    | weight | vs the deployed 0.50 | vs AssayFormer alone |
    |---|---|---|
    | 0.50 *(deployed)* | -- | **+0.506 [+0.159, +0.945]** |
    | 0.75 | **−0.374 [−0.722, −0.050]** | +0.132 [−0.427, +0.696] *(crosses zero)* |

    **At 0.75 the method stops being distinguishable from AssayFormer alone on `test`.**
    EF falls 6.922 -> 6.486, the paired loss against 0.5 excludes zero, and the absolute
    gain collapses from +0.506 to +0.132 with the interval spanning zero. Note the win
    rate barely moves (112/248 against 0.5's 111/248) while the mean falls by 0.374 --
    **too high a DepMap dose does not lose slightly on many screens, it occasionally
    wrecks one.**

    **THE DISCIPLINE PAID, AND THE AMOUNT IS MEASURABLE.** On `validation`, 0.75 scored
    HIGHER than 0.5 (+0.467 vs +0.443) and the pre-stated rule -- move only on an
    interval excluding the incumbent -- is the only reason it was not adopted. Had that
    0.024 lead been taken at face value, the deployed method would have **lost 0.374 on
    the second fold and ceased to be significantly better than AssayFormer alone there.**
    This is the §4.55(b) curated-20 trap in a new place: not a small screen set this
    time, but a single fold and a flat-looking optimum.

    **AND IT RETIRES THE "BROAD PLATEAU" READING.** `validation` makes 0.5-0.75 look
    like a flat top; `test` says the function is much sharper. **A dose-response measured
    on one fold does not establish the width of an optimum, only its location.**

    **(c) AND THE SAME 0.5 IS OPTIMAL FOR THE ONE-TERM METHOD TOO, WHICH REFUTES THE
    OBVIOUS OBJECTION TO (a) AND (b).** 2026-10-28, user's catch. Both sweeps above held
    `knn_weight=0.5` throughout, so they establish the best DepMap dose INSIDE the
    two-term stack -- not for AssayFormer + DepMap alone. That mattered: under RRF each
    added term contributes `w/(60+rank)` on top of AssayFormer's list, so total added
    weight is 1.0 in the two-term stack and only 0.5 with DepMap alone, and halving the
    displacement from AssayFormer's ordering could plausibly have moved the one-term
    optimum HIGHER. If so, every "DepMap alone vs two-term" comparison in §4.57 and
    §4.58 had been handicapping DepMap alone with a weight tuned for a different stack.
    `sweep-dmonly-c{025,075,100}`, `knn_weight=0.0`, round 1 **218/218** identical to
    `sweep-dmne-cd-validation` and the scorer asserted the kNN term off on every arm.

    | DepMap weight, kNN OFF | vs AssayFormer alone | vs the 0.50 point |
    |---|---|---|
    | 0.25 | +0.077 [+0.023, +0.145] | **−0.389 [−0.596, −0.207]** |
    | **0.50** | **+0.466 [+0.277, +0.683]** | -- |
    | 0.75 | +0.306 [+0.067, +0.571] | **−0.160 [−0.296, −0.032]** |
    | 1.00 | +0.260 [+0.015, +0.530] | **−0.205 [−0.348, −0.070]** |

    **The hypothesis is dead: 0.5 wins, and more decisively than in the two-term stack --
    all three alternatives are significantly worse here, where in the two-term stack 0.75
    was a tie.** So DepMap alone was never handicapped; it was already running at its own
    optimum in every comparison, and §4.57(e)'s conclusion is unaffected. On `validation`
    DepMap alone at 0.5 ties the two-term stack (the kNN increment is −0.022); on `test`
    it is +0.315 crossing zero against +0.506. **Nothing changes, and that is the value of
    the run -- it closes a live loophole rather than opening one.**

    **ONE OBSERVATION WORTH FOLLOWING UP, flagged as speculation.** The one-term dose
    curve is much SHARPER: at 0.25 it loses −0.389 where the two-term stack loses −0.198,
    and at 0.75 it loses −0.160 where the two-term stack was flat at +0.023. A plausible
    reading is that the second term dilutes a misdosed first one, so the two-term stack is
    more forgiving of weight error than it is better at the optimum. That would fit
    §4.57(e)'s finding that the kNN half makes the result hold up rather than raising it,
    and it is testable by sweeping `knn_weight` at a deliberately wrong DepMap dose. Not
    run; not needed for any current claim.

    Arms `sweep-dmw-c{025,075,100}` (validation, kNN on), `sweep-dmwtest-c075` (test, kNN
    on), `sweep-dmonly-c{025,075,100}` (validation, kNN OFF).

---

61. **GLM-5 + TOP-UP HANDED TO AF + DEPMAP ALONE. DepMap is worth +0.362 on top of the
    paper's own design, the LLM warm start still costs more than it buys, and the
    screen-kNN half FLIPS SIGN from a cold start to a warm one.** 2026-10-28, user's
    request. **Free** -- the top-up warm-start rounds were replayed from disk and
    **zero `llm_calls.jsonl` were written**, verified after the fact. Round 1 is the
    replayed LLM batch and matched the two-term arm off the same warm start 20/20 on
    both sets, which is what makes the replay comparable at all.

    `sweep-dmcd-hyb-glm5tu-{pv,pt}`: `weight=0.5`, `knn_weight=0.0`, reduced atlas. Its
    twin `sweep-dmne-hyb-glm5tu-{pv,pt}` hands the same warm start to the two-term
    method, so the LLM question and the one-vs-two-component question stop being
    entangled.

    | arm | paper_validation | paper_test |
    |---|---|---|
    | warm start -> AssayFormer alone *(the paper's design)* | -- | 4.805 |
    | warm start -> AF + DepMap *(new)* | 5.345 | 5.166 |
    | warm start -> AF + DepMap + screen-kNN | 5.442 | 5.251 |
    | **AF + DepMap, NO language model** | **6.184** | **5.532** |

    **(a) DEPMAP IS WORTH +0.362 [+0.188, +0.553], 16/20, ON TOP OF THE PAPER'S OWN
    DESIGN.** 4.805 -> 5.166 on `paper_test`. The paper's handoff is not a dead end to be
    replaced; bolting the atlas onto it recovers a third of a point, with an interval
    clear of zero. That is the most directly useful sentence in this section for anyone
    who wants to keep the LLM stage.

    **(b) THE WARM START STILL COSTS MORE THAN IT BUYS, which reproduces §4.53 for the
    one-term method.** Against AF + DepMap with no language model at all: **−0.516
    [−0.859, −0.199], 3 of 20 screens** on `paper_validation`, and −0.365 [−0.976,
    +0.184] on `paper_test`. Three screens out of twenty is about as one-sided as this
    project gets. The LLM stage is a cost whichever method receives the handoff.

    **(c) THE SCREEN-kNN HALF FLIPS SIGN BETWEEN A COLD START AND A WARM ONE, and this is
    the first mechanism offered for its instability.**

    | set | increment over AF + DepMap, COLD start | same, after the GLM-5 top-up warm start |
    |---|---|---|
    | paper_validation | −0.044 [−0.220, +0.094] | **+0.087 [+0.003, +0.189]** |
    | paper_test | −0.170 [−0.542, +0.110] | +0.085 [−0.076, +0.242] |

    Negative cold on both curated sets, positive warm on both, and on
    `paper_validation` the warm figure's interval excludes zero. A plausible mechanism:
    cold, the term switches on at round 2 with a handful of labels; after a three-round
    warm start it switches on at round 4 with ~300, so it is being asked to vote from a
    far better-populated neighbourhood. **That would also explain why this term keeps
    reading as "unresolved" -- its value may depend on how many labels it has rather than
    on which fold it is measured on.**

    **DO NOT PROMOTE (c) TO A FINDING YET, and the reason is this session's own record.**
    It rests entirely on two 20-screen curated sets, 11 and 20 families, and
    `paper_validation`'s interval-excluding-zero is the exact shape that lied in §4.55(b)
    (gene-kNN, +0.155 excluding zero on the curated 20, reversing to −0.094 on the fold)
    and again in §4.58(c). **It cannot currently be checked on a fold at all**, because
    the LLM warm starts only exist for `paper_validation` and `paper_test` -- producing
    them for `validation` or `test` would mean paid API calls across 100-248 screens.
    The cheap partial test is the label-count hypothesis itself, which needs no LLM: give
    the kernel a non-LLM warm start (random or oracle seeds, both already built in
    §4.49/§4.50) and see whether the kNN increment turns positive as the seed count
    rises. **That is the arm to run before (c) is believed.**

    Arms `sweep-dmcd-hyb-glm5tu-{pv,pt}`.

62. **SUPPLY, AS THREE FIGURES (11-13) IN THE COMPARISON NOTEBOOK. Every LLM on this
    benchmark fails first by not filling the slots, and that failure is already inside
    every EF number in this file.** 2026-10-29, user's request. **Free, no runs** --
    reads the paper's workbook (`Mean_long.n_acquired`, already a mean over its 20 test
    screens) and our own pick lists.

    The quantity is `n_acquired`: genes the policy supplied and the benchmark banked. It
    is **not** the budget. A round answered with 83 symbols, four of them not gene names,
    spends its 100 slots and banks 79; `ef_curve` charges the intended budget and not the
    picks delivered, so the shortfall is an EF penalty already paid (notebook trap (b),
    §4.14).

    **(a) FIGURE 11 -- the 21 policies that ran all ten rounds themselves.** Workbook
    families `llm` (15), `agent` (3) and `assayllm` (3); the handoffs are excluded,
    because after round 2 or 3 the picks are AssayFormer's. **Exactly one of the 21
    delivers all 1000 picks (LLMNN).** The range is 382 (Claude Haiku-4.5) to 1000, with
    the strongest prompted models clustered at 948-969 (GLM-5.2, GPT-5.5, GPT-5.6 Sol,
    Gemini-3.1-pro) and GLM-5.1 at 894. A quarter of the field banks under 710.

    **(b) FIGURE 12 -- the warm start's own supply, cut at round 3, and the k=2/k=3 split
    that makes a naive cut wrong.** The paper's six handoffs do **not** all hand over at
    the same round: GLM-5.1, AssayLLM and AssayLLM (handoff-trained) run the LLM for
    **two** rounds, Gemini-3.1-Pro, GPT-5.5 and GPT-5.6 Sol for **three**. AssayFormer
    fills every slot it is given, so a bar cut at round 3 credits a k=2 arm with 100
    picks its LLM never produced. Each bar is therefore split at k and only the left
    segment is the LLM's. `warm_rounds()` infers k as the last round that fell short,
    **asserted** against the requirement that every later round delivered exactly 100 --
    which is what makes it an inference rather than a guess.

    | warm start | k | supplied by round 3, of 300 |
    |---|---|---|
    | Flash-Lite *(ours)* | 3 | 217.2 |
    | GLM-5 *(ours)* | 3 | 233.8 |
    | GLM-5.1 *(paper)* | 2 | 290.1 *(190.1 of the LLM's own 200)* |
    | **GLM-5 + top-up** *(ours)* | 3 | **291.4** |
    | AssayLLM, both variants *(paper)* | 2 | 291.8 / 291.9 |
    | GPT-5.6 Sol / GPT-5.5 / Gemini-3.1-Pro *(paper)* | 3 | 295.9 / 296.1 / 296.6 |
    | **Flash-Lite + top-up** *(ours)* | 3 | **296.9** |

    **The top-up arms put our warm starts inside the paper's own range and the
    un-topped-up ones are 60-80 picks short of it.** That is the first side-by-side
    confirmation that §4.53's top-up closed a gap against the PUBLISHED runs rather than
    only against our own earlier ones.

    Our rows are deduplicated **on the warm start, not on the downstream method**:
    `sweep-glm5tu-pt-af`, `sweep-dmne-hyb-glm5tu-pt` and `sweep-dmcd-hyb-glm5tu-pt`
    replay ONE GLM-5 warm start into three different scorers, and three identical bars
    would read as three measurements. The dedup asserts that invariant (same picks in
    rounds 1..k) rather than assuming it. This figure reads `REG`, not `DF`, on purpose:
    `REDUCED_ATLAS_ONLY` correctly hides full-atlas arms from every figure that compares
    METHODS, but a warm start happens before any atlas is read.

    **(c) FIGURE 13 -- the paper's GLM-5.1 against our GLM-5, same 20 screens, same 100
    slots a round.** The closest like-for-like pair in the project. **The version differs
    -- ours is `zai.glm-5`, the paper's is GLM-5.1 -- so this is a comparison of supply
    discipline between two builds, not of capability.** Rounds 1-3 only, and not as a
    crop: our GLM runs carry `n_steps: 3` because they were commissioned as warm starts.

    | run | R1 | R2 | R3 | total of 300 |
    |---|---|---|---|---|
    | GLM-5.1, the paper's | 94.85 | 95.25 | 92.65 | 282.8 *(94.2%)* |
    | GLM-5, ours, no top-up | 91.75 | 87.60 | **54.50** | 233.8 *(78.0%)* |
    | **GLM-5 + top-up, ours** | 100.00 | 97.75 | 93.70 | **291.4 *(97.2%)*** |

    **Our un-topped-up GLM degrades with the round and the paper's does not** -- 91.8 ->
    87.6 -> 54.5 against a flat 94.9 -> 95.3 -> 92.7. The top-up run is flat too, and
    slightly ahead of the paper's. So the shortfall was ours, in how the run was driven,
    not a property of the model family; and it was a growing one, which is exactly the
    shape that a single-round pilot would have missed. (The paper's GLM-5.1 went on to
    finish all ten rounds at 893.7 of 1000.)

    **REGISTRY GAP CLOSED.** `sweep-glm5testtopup` -- the top-up warm start on
    `paper_test`, scored alone -- had never been registered although its
    `paper_validation` twin `sweep-glm5valtopup` had been there all along. 76 arms now.

    **ONE DELIBERATE PALETTE DEPARTURE, recorded so it is not read as drift.** Everywhere
    else in that notebook a hue means one method. 21 methods cannot have 21 hues, so in
    Figures 11 and 12 **the hue means the FAMILY, or which policy supplied the pick**, and
    the method is carried by the axis label. Figure 13 is back to one hue per run.
    Figure 1's other convention -- published rows neutral AND hatched -- is for charts
    that mix sources; in Figure 12 the hue is already carrying the supplier, so only the
    hatch marks the paper's rows, and in Figure 13's right panel the hatch had to be
    drawn in the surface colour because grey-on-grey rendered as a plain solid bar and
    the source marking simply was not there.

    **(d) TWO DEFECTS IN THE PER-ROUND SECTION, FOUND BY THE NOTEBOOK'S OWN WARNINGS
    WHILE ADDING THESE FIGURES, BOTH FIXED.**

    1. **`_our_specs` built its series from `REG`, not `DF`, so `REDUCED_ATLAS_ONLY` was
       never applied to it.** The per-round figures were the one place in the notebook
       that could still draw a full-atlas arm beside an ablated one -- the single
       comparison that flag exists to prevent -- and
       `Flash-Lite warm start -> AF + DepMap + screen-kNN` (full atlas) was being drawn
       that way in the EF comparison. Now filtered at the source; that term was removed
       from the cell, since it is one of the 8 orphaned methods with no ablated twin.
    2. **`resolve_series` never got `pick()`'s exact-method-name stage**, so after the
       arms were renamed by component its substring fallback turned a 4-series request
       into 9, and the 8-hue cap then silently dropped one. Same ladder as `pick()` now:
       exact label, exact group, exact method name, then substring. A third term,
       `"Paper's method: GLM-5 rounds"`, had not existed since the renaming and matched
       nothing, so that figure had been drawing one series short with a printed warning
       nobody read.

    **Both were visible in the notebook's own `!!` output on every run and neither was
    noticed until a new figure made me read it.** The warnings are only worth printing if
    the executed output is read; an assertion would have stopped the run.

---

## 5. Code added this session

**`src/assayloop/acquisitions/codep_seed_acq.py`** (NEW 2026-10-21, extended 2026-10-22) — the
project's only working intervention. `--acq codep_seed`, registered in `experiment/runner.py`.
AssayFormer is frozen; only the acquisition changes. Full parameter list:

| `--acq-param` | default | what it does |
|---|---|---|
| `profile_npz` | required | DepMap matrix. **Use `achilles_noeval.npz`, the default atlas since §4.56** (`build_noeval_atlas.py`); `achilles.npz` is the superseded full 808-line matrix and `Expression.npz` was refuted in §4.24. Deliberately has NO default here, so every arm's `config.json` records which atlas it read |
| `weight` | 1.0 | the DepMap scorer's fusion weight; **deployed at 0.5**. **0 reproduces plain greedy exactly.** Confirmed on the full atlas (§4.19/§4.20) and RE-confirmed on the reduced atlas (§4.60) on BOTH folds: on `validation` 0.25 and 1.00 are significantly worse and 0.75 is indistinguishable (+0.023); on `test` 0.75 is **significantly worse (−0.374)** and the method stops beating AssayFormer alone there. **Do not raise it.** Also optimal at 0.5 with the screen-kNN term OFF (§4.60c), so DepMap-alone comparisons were never handicapped |
| `query` | `diffmeans` | `diffmeans` = hit mean minus non-hit mean (§4.15); `centroid` = hit mean only, the pre-2026-10-22 behaviour |
| `min_nonhits` | 10 | below this many non-hits **present in the atlas**, `diffmeans` falls back to the centroid |
| `max_labels` | 600 | cap on labels per mean, for speed |
| `fuse` | `rrf` | `rrf` (rank-based, scale-free) or `zsum` (standardised sum) |
| `min_seeds` | 2 | hits required before anything fires. **Confirmed optimal (§4.20); raising it is monotonically worse overall** |
| `train_hits` | none | prefix from `build_train_hits.py`, enables the screen-kNN scorer |
| `knn_weight` | 0.0 | that scorer's weight. **0 is bit-identical to DepMap-only** |
| `knn_k` | 50 | training screens averaged |
| `gene_locus` | none | HGNC `hgnc_complete_set.txt`; turns the chromosomal-proximity exclusion ON (§4.54) |
| `prox_level` | `arm` | `band`/`arm`/`chrom`; only read when `gene_locus` is set |
| `gpt_npz` / `gpt_weight` / `gpt_k` | none / 0 / 10 | the PAPER's gene-kNN as a ranker (§4.55); `genept_to_npz.py` builds the npz |
| `pknn_weight` / `pknn_train` / `pknn_tau` / `pknn_kappa` | 0 / `train` / 0.1 / 100 | the PAPER's `Screen-kNN` as a ranker (§4.55), its own estimator, not ours |
| `pknn_gate` | `True` | **holds their scorer silent below `min_seeds`, as ours are.** `False` is the ungated behaviour the pre-2026-10-28 `sweep-pkv-*` arms ran, where with no labels the score is the pure global hit prior and round 1 is reordered — so those arms were NOT nested in the baseline (§4.56(e)). The default CHANGED; `False` reproduces them exactly |
| `dir_bank` | none | direction bank from `build_dir_bank.py`; only read when `shrink_alpha > 0` |
| `shrink_alpha` | 0.0 | shrinkage toward the borrowed prior direction. **0 is bit-identical to unshrunk.** REFUTED at every value tested (§4.16) — kept because the code is nested and harmless, not because it is recommended |
| `shrink_k` | 25 | training screens whose directions form the prior |
| `shrink_drop_top` | 0 | nearest neighbours discarded first, the near-duplicate control |
| `shrink_min_hits` | 50 | atlas-covered hits a training screen needs to enter the bank |
| `shrink_shuffle` | −1 | >= 0 permutes the bank: **a control arm**, tagged `_SHUF` in the acquisition name so it cannot be mistaken for the method in a sweep listing |

Scoring is a three-component rank fusion:
`fused(g) = 1/(k + r_AF) + w_codep * phi_codep(g) + w_knn * phi_knn(g)`, k = 60.

**Four invariants, all test-pinned, all load-bearing:**

- `weight=0` and `knn_weight=0` reproduce AssayFormer **pick for pick** — verified end to end, not
  just in a unit test (all 1000 picks and 47 hits matched on screen 1561). This is what makes every
  sweep contain its own baseline and the worst case parity rather than regression.
- **Round 1 is always plain greedy**, as is every round of a screen that never reveals a hit. That
  bounds the downside precisely on the hard screens.
- **A gene missing from a data source is made NEUTRAL, not unranked.** Under RRF, merely *having* a
  score is worth up to `weight/60`, so omitting the 3.4k genes absent from DepMap would have
  systematically outranked them — the coverage confound re-entering through the fusion rule, right
  after §4.10 measured it. They take the MEAN contribution of the covered genes instead. Shared
  `_z` and `_rrf_term` helpers apply this to both scorers.
- Either scorer fires **independently** of the other, so the additivity §4.14 measured offline can
  actually be realised.

**Supporting scripts** (all free, no API calls):

| script | what it does |
|---|---|
| `depmap_to_npz.py` | caches the DepMap CSVs as `.npz`; the 18k-column CSV parse was the entire runtime |
| `build_train_hits.py` | caches the 1349 x 20874 training hit matrix (1.5M entries, 5.32% dense, 1.3 MB) |
| `build_dir_bank.py` | one offline-optimal profile direction per training screen, from its COMPLETE labels (1289 usable, median 1294 atlas-covered hits); prints the pairwise-cosine spread that says whether neighbour choice can matter |
| `screen_shrinkage.py` | the §4.16 alpha sweep, with all three controls (shuffled bank, drop-nearest-N, width-matched noise) |
| `perturbseq_to_npz.py` | converts Replogle 2022 genome-scale Perturb-seq pseudobulk to the shared `.npz` layout; **zeroes the atlas's 6832 infinities**, which otherwise silently null the kernel (§4.17) |
| `screen_kernels.py` | cheap smoke test that a candidate gene space carries gene information, against a shuffled control. **Measures standalone, not incremental, ability — see the warning in its docstring and §4.17. Do not screen features with it.** |
| `screen_pair.py` | screens COMBINATIONS of atlases plus screen-kNN (§4.18). `--atlas LABEL=PATH` repeatable. Prints inter-column correlations (substitutes vs complements) and a paired-contrast table, which is the output to read — a block mean rewards column width, the contrast does not |
| `run_weight_sweep.sh` / `score_weight_sweep.py` | the §4.19 joint weight sweep and its scorer. The scorer reports per-screen and per-family means **each with its own bootstrap CI** (`boot_screens` resamples screens, `score_ctx.boot_ci` resamples families — do not print the latter beside the former) and PAIRED contrasts against a reference arm, which is where the power is |
| `screen_seed_aggs.py` | screens candidate query directions and screen-kNN offline, 1-column, 5-draw noise controls |
| `screen_profiles.py` | screens multi-column profile blocks with width-matched noise controls |
| `run_codep_campaigns.sh` | the original weight sweep |
| `run_codep_twofeat.sh` | adds the screen-kNN scorer |
| `run_codep_diffmeans.sh` | the current best configuration |
| `run_codep_confirm.sh` / `score_codep_confirm.sh` | held-out fold confirmation, excluding the curated screens |
| `codep_gated.py` | scores a deployable gate; `--exclude-set` drops the screens that chose a setting |
| `codep_trace_audit.py` | reads `n_seeds` per round to predict what a different threshold would suppress |
| `task_ceiling.py` | perfect-play ceiling per screen (§4.12) |
| `run_cohit.sh` | the §4.21 three-arm co-hit campaign |
| `score_min_seeds1.py` | scores `min_seeds=1` SPLIT by whether the feature could fire at all under the default, so the overall number cannot hide a mismatch between "it improved" and "it improved where predicted" (§4.23) |
| `score_by_difficulty.py` | re-checks every discounted arm on the HARD screens from campaigns already run, with paired contrasts. `--arm LABEL=sweep:baseline`, `--min-library 1001` to drop metric-degenerate screens (§4.29) |
| `characterize_hard_hits.py` | what KIND of gene is a hit on the unreachable screens: training-label coverage and DepMap amplitude by screen group (§4.25) |
| `prescreen_prism.py` | the three gates on a gene x DRUG feature from DepMap PRISM (§4.27) |
| `prescreen_ppi.py` | BioPlex protein-interaction kernel, sparse, with the difficulty split and the paired-vs-DepMap contrast (§4.28) |
| `prescreen_cohit_direction.py` | co-hit restricted by TRAINING-screen selection direction (§4.31) |
| `prescreen_cohit_ppmi.py` | PPMI + SVD variants of co-hit, with a plain-SVD control (§4.32) |
| `prescreen_condition_match.py` | condition-matched screen transfer, with the random-screen and no-same-cell-line controls that decided it (§4.34) |
| `codep_seed_acq.py` ppi path | NEW (§4.43). `ppi_npz`, `ppi_weight`, `ppi_degree_residual`, all defaulting to off and bit-identical to the shipped behaviour when off. `_ppi_scores` scores the sparse adjacency **without densifying**, exactly: the centring term cancels in the numerator because `q` is built from centred rows and so sums to zero, and the denominator is precomputed. `test_sparse_path_equals_the_dense_computation` pins the identity against an explicit dense computation |
| `prescreen_atlas.py` | **the generic three-gate runner — use this for any new atlas.** Adds GATE 0, a positive control on known complexes. Both centring modes always reported. **TWO KNOWN DEFECTS, both from §4.41: (a) its GATE 0 uses six hand-picked complexes (~75 pairs), too few to support a verdict, and complexes are the WRONG ground truth for a gain-of-function atlas — take the sets from `gate0_gof.py` instead; (b) its gate 1 computes `corr(profile norm, hit frequency)`, not `corr(feature, ...)` as documented, so its gate-1 numbers are not comparable to the PRISM/BioPlex figures** |
| `read_paper_xlsx.py` | NEW (§4.52). Stdlib reader for `paper_per_round_data.xlsx`, the paper's own 59 methods x 10 rounds. `per_round()` joins display names onto the per-screen sheet; `ef_bounds()` returns the EF bracket the file can support. **Match paper methods by per-screen hit counts, not by name — its `AssayFormer` row is a weaker ablation and the released model is `AssayFormer (BPMF + GRPO)`** |
| `handoff_fullgenome.py --inner-acq` | NEW flag (§4.51). `greedy` (default, unchanged) or `codep_seed`: what runs AFTER the warm-start rounds. With `codep_seed` the LLM seeds rounds 1-n and the fused ranker drives the rest off those seeds. **Replays warm-start rounds from disk, so it costs no API calls** -- `--warm-dir output/runs/baseline/gemini --warm-prefix sweep-gem35test-` |
| `output/runs/compare_all_results.ipynb` | **SELECTION CONVENTION: all 14 figure functions take `select=`, routed through one `pick()` helper, and `DEFAULT_SELECT` in the setup cell changes every figure at once.** **`list_bars()` prints every arm as a copy-pasteable Python list** (full label, EF and group per line, highest first; `list_bars(sset)` for one panel, `include_truncated=True` to add the warm-start-only arms the hit-rate figures can draw) -- the intended workflow is copy the block, delete lines, pass the rest as `select=`. **Three selection vocabularies exist and they are NOT interchangeable:** 8 figures select ARMS through `pick()` (listed by `list_bars`); figures 10 and the `fig_rounds` comparisons mix our arms with the paper's workbook and are listed by `list_series(sset)`; **Figures 11-13 (§62, supply) take no `select=` at all** -- they are fixed sets defined by the question, not by this notebook's arm registry; figure 5 selects injection CONDITIONS (`select=["ordinary", "4 free"]`). **`pick()` resolves a term as exact label, then exact GROUP, then the `ours`/`paper` shorthands, then substring -- in that order.** The group check MUST come first: there is a group literally named `ours`, the shorthand meant "everything not the paper's", and these frames hold no paper rows at all, so `select=["ours"]` silently expanded to EVERY arm and a 3-arm request drew 10. **THE PAPER'S PUBLISHED ARMS ARE ROWS OF `DF`**, not just of the exported CSV, so the absolute-EF figures can draw them beside our runs (`"Paper: GLM-5.1 + AssayLoop"` = 5.094, band [4.871, 5.317]). The workbook is loaded with the master table rather than in the per-round section, because a figure cannot plot an arm that does not exist yet and the figures run first. **Where the paper PUBLISHES an EF, that number is plotted, not a reconstruction** (`PAPER_REPORTED`: AssayFormer 4.83, Gemini 3.1 Pro + AssayLoop 5.67, GLM-5.1 + AssayLoop 5.27). This also settles which estimator the reconstruction should use: each published figure lands just UNDER our band's high end (4.83 vs 4.896, 5.27 vs 5.317, 5.67 vs 5.699), so the midpoint was wrong and arms the paper does not report now fall back to `hi`. **Our EF formula is NOT different from the paper's** -- our AssayFormer run gives 4.828 against their reported 4.83 on screens whose hit counts match exactly, so the old band was a MISSING COLUMN, not a metric disagreement: rebuilding EF needs the per-screen count of picks that were in-universe but outside that screen's own library, and the workbook carries it only as a cross-screen mean (`Mean_long.in_universe_not_lib`; `By_screen_long` has `out_of_universe` but not it). **The workbook's own `e_t` column is a THIRD metric and must not be mistaken for the paper's headline EF:** the mean of `e_t` over rounds 1-10 is algebraically identical to `cum_hits / (1000 * h_rand)`, which forgives nothing and equals our band's LOW end exactly (4.4065 vs 4.407). They carry `external=True` and `band_lo`/`band_hi`: their EF is a BAND (the workbook records picks outside the universe but not picks outside each screen's own library, which our metric forgives), so `ef_screen` is only the midpoint. **Every figure built on `curves` or on `delta_vs_af` must filter `~DF.external`** -- they have no per-screen pick lists, so there is no paired delta and no per-round curve. They are drawn hatched, in the NEUTRAL, with the band as a rule through the bar: with 10 groups the palette modulo wrapped the paper group onto slot 1, the same orange as "Our method", and the legend showed two different things in one colour. **Figures 8 and 9 (screen-type cuts) now take `select=` as well**, via `type_rows()`/`type_arms()`: `type_table()` was always arm-generic and only `TT` hardcoded the deployed arm, so per-arm tables are built on demand and cached rather than eagerly for all 50-odd arms. Both cap at 7 arms plus the AssayFormer reference, because the palette holds 8 hues and a 9th is a repeat that makes the chart lie about identity. A term is a full label, a case-insensitive substring, a group name, or `ours`/`paper`/`all`; an unmatched term prints `!! nothing matched` rather than drawing an empty panel. **Any figure added later must take `select=None` and pass it through `pick()`.** `fig_absolute` sorts by value (lowest at top; `sort="value_desc"` or `"group"` to change). NEW. **The comparison notebook: 76 arms x 4 screen sets, 14 figure functions rendering 22 figures, three CSVs** (`compare_all_results.csv`, `compare_hit_rates.csv`, `compare_by_screen_type.csv`). Beyond EF it carries per-round and cumulative HIT RATE (both denominators: in-library picks, and slots offered -- the gap between them is the LLM shortfall, reproducing §2c's 20.0% -> 2.1% exactly), a perfect-play ceiling on every rate plot, drawn as a SEPARATE with-ceiling figure rather than a broken axis (`show_max=True/False`), because the ceiling sits so far above the arms that sharing a panel destroys the granularity the lower lines exist to show, and **per-screen-type breakdowns that independently reproduce §4.30 on BOTH folds: knockout +0.55/+0.42 against activation +0.02/+0.11**. Two flags at the top govern everything: `GENOME_SCALE_ONLY` (drop libraries below the 1000-pick budget) and `EXCLUDE_CURATED_FROM_FOLDS` (keep the 20 curated screens out of the fold panels, which is what makes the fold panels mean the 198/314 FRESH screens this file quotes). **Three traps it documents, all of which bit during construction:** (a) EF averaged over SCREENS (`ef_screen`, the paper's convention, 6.10 on paper_validation) differs sharply from EF averaged over FAMILIES (`ef`, 4.30) because 20 screens collapse to 11 families -- absolute figures and the paper's 5.67 line use the screen mean, deltas use the family mean; (b) the per-round budget must be the INTENDED budget, not the picks delivered, or LLM shortfall is forgiven and `sweep-glm5val` reads 1.74 instead of 0.50 on screen 1561; (c) arms handed free hits are held off the comparison figures entirely. **LEGEND CONVENTION: every legend sits OUTSIDE the axes, below the panels, with a reserved band in figure coordinates (`fig.legend(..., loc="lower center")` plus a `rect` bottom in `tight_layout`).** The bar charts sort descending, so the longest bars land at the bottom right -- exactly where an in-axes `loc="lower right"` legend sits, and it covered the value labels on the best arms in `fig_absolute`, `fig_pct` and `fig_type_absolute`. Value labels also carry a background `bbox` so the reference lines (theoretical maximum, the paper's 5.67) cannot run through the digits. **SERIALIZATION TRAP: write notebook `source` lists with trailing `\n` on every line but the last.** A builder that used `text.split("\n")` stored 44 of 45 cells without newlines; each cell then collapses to ONE line, so everything after the leading `#` comment becomes part of that comment and the cell silently executes as a no-op -- it does not raise, it just does nothing. Verify any notebook edit by executing it: `.venv/bin/jupyter nbconvert --to notebook --execute --output-dir <tmp> --output executed.ipynb compare_all_results.ipynb`, then assert zero `output_type == "error"` outputs and 22 `image/png` outputs |
| `bioplex_to_npz.py` | NEW (§4.43). Saves BioPlex as a SPARSE adjacency plus the two row statistics the kernel needs (`row_sum`, `row_ssq`, `norm`). **The dense form is never built.** 189128 interactions, 14551 proteins, 334798 nonzeros (0.16% dense) |
| `run_ppi_sweep.sh` | NEW (§4.43). Three pre-specified arms on the validation fold: ppi_weight 0.25, 0.5, and 0.25 degree-residualised (the fame control). **`DEV=cuda` to override the device — it defaults to cpu because GPU0 on this machine is in a driver error state (`nvidia-smi` reports no devices), and `device=cuda` fails every screen with a `torch.load` error while the runner still exits 0** |
| `score_resolution.py` | NEW (§4.42). Re-measures any paired arm with a per-round EF curve and ANCOVA. Note `_adj_ef_from_run`'s `budget` argument does NOT truncate the pick list, so a real budget curve needs pick prefixes — `ef_curve()` builds them and asserts the last point against `_adj_ef_from_run` |
| `atlas_complement.py` | NEW (§4.41). Conditional AUROC for any atlas pair: can atlas B see the pairs atlas A is blind to? |
| `build_gate0_sets.py` | NEW (§4.41). Builds GATE 0's positive sets from HGNC `gene_group` + MSigDB C2 CP into `output/analysis/gate0_sets.json`: **301 complexes, 1402 paralogue families, 2323 pathways** at size 2-12, plus a spotlight of 12 families with textbook convergent overexpression phenotypes. Its docstring states the only two hand choices in the file (families not expressed in culture; families above 12 genes) |
| `gate0_gof.py` | NEW (§4.41). Runs GATE 0 across several atlases and all three ground truths at once (`--atlas LABEL=PATH`, repeatable), reports AUROC of within-set against 200k random covered pairs with a **bootstrap over families rather than pairs**, and prints the paralog-minus-complex dissociation contrast. **Use this, not `prescreen_atlas.py`'s gate 0, whenever the candidate atlas is not a loss-of-function fitness atlas.** Pair cosines are chunked: the direct form needs 5 GB on the pathway group and gets OOM-killed |
| `jump_to_npz.py` | JUMP Cell Painting ORF profiles -> the shared `.npz` layout (§4.38, §4.39) |
| `sack_crispra_overlap.py` | Fisher test of CRISPRa resistance hits against Sack's HPNE GO/STOP gene sets (§4.40) |
| `run_profile_screens.sh` | multi-column profile-block screens (pre-dates the §4.15 reframing) |
| `oracle_u.py` / `augment_u.py` | the ceiling and added-feature machinery (§4.7, §4.9) |


**`output/runs/analyze_results_full.ipynb` section 7** (NEW 2026-10-22, cells `codep7md` /
`codep7`) — visualises the `codep_seed` campaigns. Auto-discovers whichever sweeps exist and skips
missing ones, so it renders against partial results. Produces a summary table (dEF, 95% CI over
screen families, gated dEF, screens improved), a **dose-response plot** with error bars for both
screen sets overlaid, and **per-screen bars** coloured red where a screen holds more hits than the
1000-pick budget and grey where the gate disables the feature. Folds with more than 40 screens show
only the two tails. Needs section 6's imports (`HAVE_ASSAYLOOP`). Its `CODEP_EXPERIMENTS` list is
where new sweep families get registered.

**ARM IDENTITY: the LLM sweeps carry no distinguishing tag**, and three were mislabelled on a
first pass through the notebook. `sweep-glm5val` / `sweep-glm5valtopup` / `sweep-gem35val` /
`sweep-gem35test` / `sweep-glm5test` are `llm_single`, the LLM driving **all ten rounds** (EF
0.63-1.38) -- NOT warm starts. The paper's actual 3-round handoff arms record `acq: greedy`,
because the LLM rounds were precomputed and injected, so the acquisition on record is
AssayFormer's. Identity was recovered by matching each arm's EF to §1's table:
**paper_validation** `sweep-8f95fe0c` 5.220 = GLM-5+top-up, `sweep-7d0a513b` 5.170 = Flash-Lite,
`sweep-d7630374` 5.148 = GLM-5, `sweep-afalone-val` / `sweep-fd2017be` 6.100 = AssayFormer alone;
**paper_test** `sweep-gem35decline` 4.904, `sweep-a7dd71ea` 4.696, `sweep-0d295c6a` 4.678,
`sweep-da9a8c7e` 4.647, `sweep-d8d10355` 4.623, `sweep-84fa85f8` 4.569, `sweep-cdt-w0` 4.828.
**Do not infer an arm's identity from its sweep id or its recorded `acq`.**

**`methods_summary.csv`** (repo root) + **`scripts/diagnostics/make_methods_summary.py`** (NEW
2026-10-21) — every method tried in this project as one table, written for a reader who has skimmed
the AssayLoop paper but does not know this codebase. **CURRENT as of 2026-10-27: 116 rows, and all 53
of §4's diagnostics are cited by at least one row** (checked mechanically, not by eye — `handoff_ref`
holds multiple refs on 11 rows, so verify coverage by regex over `S4\.(\d+)`, not by set equality).
Regenerate with `python3 scripts/diagnostics/make_methods_summary.py`; **the generator is the source
of truth and the CSV is a build artefact, so edit the generator.** Its numbers are transcribed from
§4 of this file.

Columns: `group`, `method`, `description`, `screen_set`, `metric`, `value`, `vs_baseline`, `ci_95`,
`verdict`, `handoff_ref`. Groups: **14 "0. Working method"** (numbered 0 so sorting by `group`
puts the positive result first), 15 LLM-as-picker, 11 AssayFormer surgery, 26 additive features,
16 ceiling measurements, 22 refuted hypotheses, and **9 "6. Method and gates"** — the three
pre-screen gates and four rules (offline never predicts deployed; filtering to the matched subset
always lost; a feature must be pointed by the screen's own results; and §4.41's, that a positive
control must suit the data's perturbation direction and be powered). Those seven are not methods
tried: they decide whether any other row's number means anything, so read them first.

Superseded detail: Groups were **9 "0. Working method"** (numbered 0 so sorting by `group` puts
the positive result first — the deployed figure, its components, the two settings confirmed on
held-out data, and the two checks it survived), 15 LLM-as-picker, 11 AssayFormer surgery, 16
additive features, 10 ceiling measurements, 10 refuted hypotheses.

**`metric` is the main way to misread the sheet.** Seven constants in the generator carry most of
the rows: absolute EF; change vs AssayFormer alone; the same change averaged over FAMILIES (about 2x
the per-screen figure here); change vs the **deployed fused method**, paired — where a positive
number means better than the best known configuration, not better than AssayFormer; the ORACLE
change, an upper bound with the weight fitted to all answers and never an achieved result; `GATE`,
a structural pre-screen property rather than any EF; and `STRUCT`, a fact about the data rather than
a method result. The remaining rows carry a one-off metric spelled out in the cell. **No row has a
blank `metric`** — four narrative refutation rows did until 2026-10-23, and they now carry the
measurement that refuted them.

Two things to preserve if it is edited:

- **`metric` is load-bearing, because the groups do not share one.** Group 1 is adjusted EF against
  AssayFormer alone; group 3 is the *oracle* delta (weight fitted to all answers — an upper bound,
  never an achieved result); group 4 is absolute EF. Reading a group-3 number as an achieved gain is
  the obvious misreading.
- **Three rows are controls, not methods**, labelled as such in `verdict`: oracle relevance
  (+9.59/+13.10) and Gaussian noise (-0.003/-0.047) bracket every group-3 result, and the per-gene
  bias row (+0.013) confirms the test returns ~0 for information already in the model. Without them
  a +0.05 is uninterpretable.

Every description is self-contained — someone who sorts the sheet and reads one line should still
understand the method and the number. Rows appear twice where validation and test disagree in sign
(the denoising dose-response, DepMap expression) rather than being averaged into one misleading
figure. Note the repo root is **not** git-ignored, unlike `output/`, so this file shows up in
`git status`.

**`scripts/diagnostics/llm_score_features.py` + `screen_llm_feature.py`** (NEW 2026-10-21) — the LLM
as an *additive feature* rather than a picker, the one role it has never had here. It contributes a
per-gene number to `a(g) = û·v_g + b_g + w_llm·llm_score(g)`, so `w_llm = 0` recovers AssayFormer
exactly and the LLM cannot displace a pick. Unlike every annotation screened in §4 diagnostic 10, the
score is conditioned on the screen's phenotype text, so it is not a fixed property of the gene.

- **Collector** (`llm_score_features.py`, costs money): prompts for an integer 0-10 per gene in
  batches of 100 over AssayFormer's top-N opening ranking. **Select screens with `--only
  1817,1818,...`, not `--screens N`.** `--screens` takes the first N in manifest order, which is
  alphabetical and says nothing about difficulty: the first pilot ran on 1561/1581/1584, which are
  ranks 7, 8 and 9 of the 10 enrichment screens by AssayFormer EF (2.04, 2.07, 4.74) — the *easiest*
  of the hard tier — and I initially described them as the hardest. The six genuinely hard
  enrichment screens are 1817 (0.98), 1818 (1.11), 1815 (1.18), 1644 (1.39), 2104 (1.39) and
  2048 (1.77). `--only` preserves the given order and errors on an unknown name rather than
  quietly scoring fewer screens. Caches every reply to
  `--out-dir/<tag>-<set>/<screen>__<offset>.json`, so a crash or a re-run costs nothing for work
  already done. **Nothing is called without `--go`**; the default prints the plan, token counts and a
  price estimate, and a dry run deliberately works with no credentials resolved.
- **`--pool` must exceed the 1000-gene budget, and the default is 3000.** A feature defined only
  inside AssayFormer's top 1000 can reorder those picks but can never change *which* genes are
  picked, so EF cannot move; `pool - 1000` is how far down the ranking the LLM may promote a gene
  from. The script warns below 1200. This was caught by the self-test, not by reasoning.
- **Scorer** (`screen_llm_feature.py`, free): oracle ceiling (w fitted to all labels), a deployable
  variant (w fitted to rounds 1-k of a real run, then applied to the rest), the AUROC of the raw
  score inside the pool, and a per-screen noise-floor column. Prints a verdict against the +0.5
  threshold from §4 diagnostic 9.
- **Self-test** (`selftest_llm_feature.py`, free, no credentials): synthesises a perfect feature and
  a random one over the real pool. Measured 2026-10-21 — perfect: AUROC 1.000, oracle **+3.174**,
  deployable **+2.145**, verdict WORTH BUILDING; random: AUROC 0.496, oracle **-0.087**, deployable
  -0.887, verdict below threshold. **Run this first if the scorer is ever touched**: it is what makes
  a null result on real scores trustworthy rather than ambiguous.
- Reply parsing and screen selection are tested (`tests/test_llm_score_feature.py`, 16 tests),
  including the two cases that have already bitten: a reasoning model truncated mid-JSON (the parser
  salvages the pairs that arrived instead of discarding the batch) and positional screen selection.
- Cost and wall clock measured by dry run at `--pool 3000` (the dry run also prints a pacing ETA
  from the config's `requests_per_minute`, which the client does enforce — `client.py:572`):

  | config | pilot (3 screens) | full set (20 screens) |
  |---|---|---|
  | `collect-gemini-3.5-flash-lite.yaml` ($0.30/$2.50 per 1M) | 90 calls, ~$0.19, ~6 min | 600 calls, ~$1.27, ~40 min |
  | `collect-glm-5-bedrock.yaml` ($1.00/$3.20 per 1M) | 90 calls, ~$0.27 | 600 calls, ~$1.77 |

  The 40-minute wall clock is the config's 15 rpm cap, not the model. Flash-Lite needs
  `GEMINI_API_KEY`, which is **not** in `.env` (that file has `AWS_BEARER_TOKEN_BEDROCK` and the
  `VLLM_*` block only). Flash-Lite prices are from a third-party aggregator, not Google's own page.
  `reasoning_effort: "minimal"` in that config is right for this job — the task is emitting 100
  integers, not reasoning.


**`src/assayloop/models/amortized_ranker.py`** (MODIFIED 2026-10-21) — two test-time interventions
on AssayFormer itself, the first work in this project that touches the model rather than the warm
start. Both are free: no retraining, no API calls, and a 20-screen 10-round validation sweep runs in
~20 s on the laptop GPU. New constructor kwargs, all reachable from the CLI as `--model-param k=v`:

| kwarg | default | what it does |
|---|---|---|
| `ctx_ensemble` | 1 | K context variants per prediction; û is a deterministic function of the observed-gene *set*, so bagging the set gives K different û |
| `ctx_subsample` | 0.8 | fraction of observations kept per variant; variant 0 is always the full context |
| `ctx_fuse` | `mean` | `mean` (average scores) or `rrf` (k=60 over the K rankings) |
| `ctx_weak_hits` | `keep` | `keep` / `drop` / `negative` — what to do with revealed hits in the weak tail |
| `ctx_weak_q` | 0.25 | quantile boundary defining that tail, over this screen's revealed hit relevances |
| `ctx_seed` | 0 | subsampling RNG seed; variants depend on (seed, n_obs) so a round is reproducible |

Two facts that made this possible and that are easy to miss:

- **`relevance_score` is revealed to the agent at acquisition time.** `tasks/gene_batch.py:116` sets
  `label={"hit": ..., "relevance_score": ...}`, and `predict` previously read only `hit`. Acting on
  revealed strength is therefore legitimate, not answer-leakage. This is *not* the retired
  "graded hit labels" idea — that one tried to add label embeddings the model never learned; this
  keeps the two embeddings it has and changes which observations get fed in.
- **The released checkpoint is `variational: false` with `dropout: 0.0`.** There is no posterior head
  over û and no MC-dropout, so perturbing the context is the only available source of û spread.

`net.encode` already takes a batch dimension, so the K variants go through in one forward pass,
padded with `ctx_pad` (True = pad). `ctx_ensemble=1` with `ctx_weak_hits="keep"` keeps the batch at 1
and is **bit-identical** to the released model — the nested baseline is enforced by
`tests/test_ctx_interventions.py::test_defaults_reproduce_the_released_model_exactly` and confirmed
end to end (`sweep-ctx-base` scores adjusted EF 6.100, the recorded AssayFormer-alone number).
Agreement is only to float32 precision once K > 1, because batching changes the reduction order.

**Tests:** `tests/test_ctx_interventions.py` (16). Suite now 195 pass, same 2 pre-existing failures.

**Scripts:** `scripts/diagnostics/run_ctx_sweeps.sh` (the eight validation sweeps, skips any that
already exist) and `scripts/diagnostics/score_ctx.py` (adjusted EF per screen, paired delta vs the
nested baseline with a 95% bootstrap CI over screen families, per-screen table, and the
direction-gate view). The scorer takes `--screen-set`, `--base`, `--base-dir`, `--runs-dir`, `--glob`,
so it scores the test set against `sweep-00866f93` in `baseline/gemini/` too.


All new files; nothing pre-existing was behaviourally changed at default settings.

**Acquisitions** (`src/assayloop/acquisitions/`)

- `llm_category_acq.py` — `llm_category` / `llm_category_blind`. Per round: plan N categories
  (spread over 8 readout stages) -> optional stage-coverage check -> fill each category -> per-category
  top-up. Registered in `experiment/runner.py`, `cli.py`, `acquisitions/__init__.py`.
  Env: `ASSAYLOOP_CATEGORY_N` (8), `ASSAYLOOP_CATEGORY_COVERAGE` (1), `ASSAYLOOP_CATEGORY_MAX_SHARE` (0.25).
  Hooks for subclasses: `_stages_block()`, `_plan_guidance()`, `_fill_guidance()`, `_after_first_fill()`,
  `covered_genes`, `fill_oversample`.
- `decline_handoff_acq.py` — `DeclineHandoffAcquisition`. LLM keeps picking until a round's hits fall
  below `ratio` x the previous round's (default 0.5), or hits zero, or `max_llm_rounds` (default 3);
  then AssayFormer for good. Uses only the LLM's own hit counts — no `|H|`, so deployable.
- `mixed_batch_acq.py` — `MixedBatchAcquisition` (LLM contributes `n_llm` genes/batch, AssayFormer
  fills the rest), `LLMComplementAcquisition` (paper prompt + "genes another model will test anyway"),
  `LLMCategoryComplementAcquisition` (categories aimed at what AssayFormer misses), `replay_proposer`
  (reuse recorded warm starts, zero LLM calls), `live_proposer`, `ENRICHMENT_RULES`.

**Other source**

- `src/assayloop/screen_direction.py` — `classify_selection(phenotype) -> (direction, why)`,
  `should_use_llm(direction, mode, include_ambiguous)`. Constants `ENRICHMENT`/`DROPOUT`/`BOTH`/`UNKNOWN`.
- `src/assayloop/acquisitions/llm_single_topup_acq.py` — only change: the rejection-reason string is
  now the class attribute `_rejected_reason` so subclasses can widen it. Default text unchanged.

**Scripts (repo root)**

- `handoff_mixed.py` — mixed-batch runner. `--llm replay|live`, `--llm-mode genes|categories`,
  `--n-llm`, `--mix-rounds`, `--af-horizon`, `--selection-mode all|enrichment|dropout`,
  `--include-ambiguous`, `--enrichment-prompt auto|always|never`, `--n-categories`,
  `--fill-oversample`, `--stage-check`, `--no-novelty-reweight`, `--sweep-id`, `--resume`.
  Prints every screen's direction classification before any LLM call.
- `handoff_decline.py` — decline-triggered handoff runner. `--ratio`, `--min-llm-rounds`,
  `--max-llm-rounds`, `--no-switch-on-zero`.
- `handoff_fullgenome.py` — pre-existing; gained `--screen`, `--selection-mode`,
  `--include-ambiguous`, `--sweep-id`, `--resume`. The universe is always built from the full screen
  set so subset runs stay comparable.
- `bedrock_cost.py` — cost of a sweep from its own `llm_calls.jsonl`, plus extrapolation.
  `--in-price`, `--out-price`, `--calls-in-full-run`, `--runs-dir`. Warns on `finish_reason: length`.

**Configs** (`configs/lm/`)

- `collect-glm-5-bedrock.yaml` — **working**. GLM-5 via Bedrock Mantle.
- `collect-gpt-oss-120b-bedrock.yaml` — works (verified by a single call), not used in any run.
- `collect-gpt-5.6-sol-bedrock.yaml` — written but **blocked**, see §7.
- `collect-gemini-3.5-flash-lite-high.yaml` — Flash-Lite at `reasoning_effort: high`.

**Tests** (all passing)

`tests/test_llm_category.py` (12), `tests/test_decline_handoff.py` (6), `tests/test_mixed_batch.py`
(11), `tests/test_screen_direction.py` (30), `tests/test_codep_seed.py` (**54** — 41 plus seven for
§4.43's sparse interaction kernel and six for §4.44's geometry knobs). Full suite as of
2026-10-27: **281 pass, 2 fail** (it
was 179 when this section was written), both pre-existing and
unrelated:
- `test_cli_full_genome.py::test_open_vocab_llm_defaults_to_shared_f2_universe` — failed before this work.
- `test_gemini_client.py::test_yaml_loader_maps_new_fields_and_nulls_inherited_extensions` — the test
  expects `reasoning_effort: "low"` but the user's yaml now says `minimal`. Needs updating.

**Notebook** — `output/runs/analyze_results_full.ipynb` gained two cells at the end of section 4,
both driven by one `CURVE_GROUPS` list (empty = all groups; entries match group names by substring):

- a selector cell that prints every available group ready to paste into `CURVE_GROUPS`;
- a per-screen hit-rate cell (table + heatmap) after the hit-recovery plot. Columns are groups, rows
  are screens sorted by `best / base`, and `base %` is the screen's own hit density (total hits /
  pool) so each row can be read against what random picking would score. `HEAT_SCALE` picks the
  heatmap colouring, because screen difficulty spans ~10x and one absolute scale renders every hard
  screen uniformly dark: `"row"` (default) scales each screen to its own min..max, so colour shows
  ranking within a screen; `"rel_ref"` shows % difference from `HEAT_REF` (e.g. `"af-alone"`)
  diverging around 0, which makes equal *proportional* gains look equal regardless of difficulty;
  `"fold"` is enrichment over the screen's base rate; `"absolute"` is the raw rate. The table always
  shows raw hit rates.

**Analysis outputs** — `output/runs/analysis_tables/`: `val_handoff_n3_hits.csv` (every hit in the
20 validation screens, found/missed, round, found_by, hit_rank, relevance_score, in_shared_universe),
plus `_found`/`_missed`/`_summary` splits and `val_hit_rate_by_screen.csv`.

---

## 6. Methodology notes that matter

**Screen families.** `paper_validation` contains near-duplicate screens: `U_1795_dec`/`U_1795_merged`
(hit-set Jaccard 1.00), `U_1793_dec`/`U_1793_merged` (0.99), `1581`/`1584` (0.98),
`U_1797_dec`/`U_1797_merged` (0.91), `1817`/`1818` (0.42), and the eight `U_179x` screens overlap each
other at 0.53-0.72. **There are ~11 independent families, not 20 screens.** Bootstrap over families,
not screens. `paper_test` has no pairs above 0.3, so all 20 are independent there.

**Metric.** Use `_adj_ef_from_run` from `assayloop.scripts.full_genome_table` with budget 1000 — the
authors' own code, so numbers are on the paper's scale. `n_hits_vs_random` in `result.json` is the
harness EF and is *not* the paper's EF.

**Primary outcome.** Total unique hits by round 10, with EF alongside. Interventions here trade
warm-start hit rate for downstream coverage, so those can move in opposite directions.

**`hit_rank` and `relevance_score`** in the CSVs: `hit_rank` was computed by this work — it ranks each
hit within its screen by the screen's own `relevance_scores` (higher = stronger hit). It is ground
truth, not any policy's score. They correlate at r = −0.9985 by construction. Several screens have
large ties at the top (1990: 121 hits tied, 1561: 29, 1725: 12, 1644: 10), where `hit_rank` ordering
within a tie is arbitrary — compare `relevance_score` instead.

---

## 7. Environment and credentials

**Screen sets: the curated 20s are subsets of much larger folds.** `paper_test` and
`paper_validation` are 20 screens each, curated by YAML manifests shipped inside the installed
package (`.venv/.../assaybench/data/screen_sets/assayloop-{test,validation}.yaml`). The manifests
name the parent fold plus the filter: 18-22k genes, >= 50 hits, hit rate < 15%, and **a Gemini
AnDCG@100 >= 0.05 signal floor**. So the paper's headline screens were pre-selected partly because an
LLM shows measurable signal on them — which biases *in favour* of the warm start, and makes the
negative result stronger. The full folds, resolved from the HuggingFace dataset
`Genentech/assaybench` config `biogrid` on the `yearfold0` column:

| `--screen-set` | screens | fresh (not in the curated 20) |
|---|---|---|
| `train` | 1349 | all (AssayFormer trained here) |
| `validation` | 218 | 198 |
| `test` | 334 | 314 |

Cached at `~/.cache/huggingface/hub/datasets--Genentech--assaybench/` plus
`~/.cache/huggingface/datasets/Genentech___assaybench/` (1.1 GB).

**EXPOSED CREDENTIALS IN THE RUN FILES — measured 2026-10-27, not previously recorded here.**
Every LLM run wrote its full `LLMClientConfig` into `result.json`, **api_key included, in plain
text**. `grep -rl` finds the Bedrock token and a Google key in **194 files** under `output/runs/`,
which is far more than the ~67 figure quoted elsewhere in this file. Reading any such run config —
to identify an arm, say — prints the secret. **Both keys must be rotated**, and the files scrubbed;
`output/runs/` is gitignored, so nothing reached git history, but the files are on disk and travel
with any copy of the output directory. The scrub is a one-line rewrite of the `api_key` field in
those 194 files and has NOT been done.

**`OPENAI_API_KEY` is needed to run most of the `test` fold, and for nothing else here.**
AssayFormer reads a 1536-d `text-embedding-3-small` vector of each screen's description. The shipped
cache covers the screens the paper used: all 40 curated ones and, as it happens, all 218 of the
`validation` fold — but only 48 of the 334 `test` fold screens. The rest raise
`RuntimeError: Text not found in the embedding cache`. Embedding the missing 286 descriptions costs
well under a cent and is written to `output/rankers/text_emb_cache/cache.npz`, so it is one-time.
**This is not the LLM**: the codep work uses no language model, and the key only fills a
preprocessing input AssayFormer has always required.

**CACHE COVERAGE, verified 2026-10-23 with `text_embed._text_key` (do not recompute the key by
hand — the separator is `\x00`, and getting it wrong matches 0 of 1356 entries and looks exactly
like an empty cache):** validation **218/218**; test **48/334**. Of those 48, **28 are FRESH**
(not among the 20 curated screens that chose settings) and **18 of those are genome-scale**
(library > 1001). **So a FREE partial test-fold run is possible on 28 screens with no API key at
all** — see §10 item 1 for what that is and is not worth.

**CREDENTIAL STATE, 2026-10-23 — all three gates checked in order, one token spent:**
1. Key valid. A `sk-proj-` key is in `.env` (loaded by `config.py:17`, gitignored). An earlier
   20-character value beginning `key` was a name rather than a secret and gave 401.
2. Model access granted. Initially 403 `Project proj_... does not have access to model
   text-embedding-3-small` with only `gpt-4o-mini` on the project allowlist; the owner enabled the
   model and `models.list()` now shows it.
3. **BLOCKED HERE: `insufficient_quota` — the account has no prepaid balance.** OpenAI's minimum
   top-up is $5, against actual usage of well under a cent for everything this project needs.
   `insufficient_quota` is reachable only past gates 1 and 2, so it confirms those are fine.
**Any future session hitting this should not re-debug the key: it is a billing state, and the user
has to decide whether to fund it.** There is deliberately no fallback — `LocalTextEmbedder` is
384-d MiniLM in a different vector space and would silently corrupt any number computed against
the 1536-d checkpoint, which is why `backend="auto"` will never select it. The cache is keyed on `sha1(model, text)` and
the checkpoint pins the model, so new vectors are consistent with the shipped ones; a local 384-d
alternative exists and the code deliberately refuses to select it implicitly.

**LOCAL DATA INVENTORY** (all free and public; none of it is in git). Every atlas is stored in the
shared `(dims x genes)` `.npz` layout unless noted, so any consumer reads it unchanged.

| path | what | size | used by |
|---|---|---|---|
| `~/depmap/achilles_noeval.npz` | **THE DEFAULT ATLAS since 2026-10-28.** 768 lines x 18119 genes -- `achilles.npz` with every line used by any evaluation set dropped | 49 MB | **the deployed method** (§4.54, §4.56) |
| `~/depmap/achilles.npz` | DepMap 24Q4 CRISPR gene effect, 808 lines x 18119 genes. Superseded as the default; still read by arms not re-run | 54 MB | pre-2026-10-28 arms |
| `~/depmap/Expression.npz` | DepMap expression, 1673 lines x 19193 genes | 68 MB | §4.24 (harmful) |
| `~/depmap/Model.csv` | cell-line metadata, for line matching | 0.6 MB | §4.10, §4.35, §4.37 |
| `output/analysis/train_hits.{npz,mat.npz}` | 1349 x 20874 training hit matrix, 5.32% dense | 1.3 MB | screen-kNN, co-hit |
| `output/analysis/dir_bank_achilles.npz` | per-training-screen optimal directions | 4 MB | §4.16 (refuted) |
| `output/analysis/cohit.npz` | gene-gene co-hit space, 1349 dims x 20874 genes | — | §4.21, §4.31, §4.32 |
| `~/perturbseq/k562_gwps.npz` | Replogle Perturb-seq, 8248 dims x 9866 genes | — | §4.17 (weak) |
| `~/prism/f466309*.csv` | PRISM Repurposing 24Q2, compounds x 919 lines | 72 MB | §4.27 (AUROC 0.50) |
| `~/bioplex/BioPlex_*.tsv` | BioPlex 3.0 AP-MS, 189k interactions, 14551 proteins | 20 MB | §4.28 (weak) |
| `~/orcs/extracted/` | BioGRID ORCS 2.0.18 human, 1952 screens | 102 MB gz | §4.33 (**= AssayBench**) |
| `~/lincs/sig_info.txt.gz` | CMap signature inventory, for the OE gene list | 11 MB | §4.36 (killed on coverage) |
| `~/jump/{orf,crispr}.csv.gz` + `profiles/orf_profiles.parquet` | JUMP Cell Painting ORF, 790 features x 12603 genes | 396 MB | §4.38, §4.39 |
| `output/analysis/jump_orf.npz` | the above, converted | — | §4.39 |
| `Sack_dataset/*.xlsx` | Sack et al. 2018 supplementary tables S1-S7 | 21 MB | §4.40 |

Download notes worth not rediscovering: `api.figshare.com/v2/file/download/<id>` works while
`figshare.com/ndownloader` returns HTTP 202 forever; BioGRID's CDN rate-limits probing with
Cloudflare error 1015, so back off and retry; **PMC serves supplementary files behind a reCAPTCHA**,
so they need a browser; the `cellpainting-gallery` S3 bucket allows `?list-type=2&prefix=...` listing
over plain HTTP but `macchiato.clue.io` does not.

**Gemini 3.5 Flash-Lite** — `configs/lm/collect-gemini-3.5-flash-lite.yaml`, `GEMINI_API_KEY` /
`VLLM_API_KEY` in `.env`. `reasoning_effort: minimal`. Free tier, 15 rpm.

**Bedrock** — `.env` holds `AWS_BEARER_TOKEN_BEDROCK` (a Mantle API key, `ABSK...`) and
`BEDROCK_BASE_URL=https://bedrock-mantle.us-east-1.api.aws/v1`. Account 964028866551, IAM user
`jhouser220`, region us-east-1.

- **Works** via the Mantle `/v1` route with that key: `zai.glm-5`, `openai.gpt-oss-120b`,
  `openai.gpt-oss-20b-1:0`, `moonshotai.kimi-k2.5`, `moonshot.kimi-k2-thinking`, `deepseek.v3.2`,
  `qwen.qwen3-next-80b-a3b`.
- **Blocked**: every proprietary OpenAI model (`openai.gpt-5.6-sol`, `gpt-6-sol`, `gpt-5.5`) returns
  "not available for this account" through *every* route tried — Converse (4 regions, `us.`/`global.`
  prefixes), bedrock-runtime `/openai/v1`, Mantle `/v1` — with both a Bedrock API key and the Mantle
  key. Control models succeed on the same routes with the same credentials, so it is not auth, region,
  request shape, or the inference-profile prefix. The inference profiles exist and are ACTIVE, which
  reflects regional availability, not entitlement.
- **One untested path remains**: the SigV4 route via `openai.BedrockOpenAI`, which currently fails
  for *all* models with a genuine IAM error — `jhouser220` lacks `bedrock-mantle:CreateInference`.
  Attaching `arn:aws:iam::aws:policy/AmazonBedrockMantleInferenceAccess` (needs admin/root) would make
  that testable. If it then still says "not available for this account", the account lacks entitlement
  and the only route left is the console catalog page or AWS Sales.
- **The repo's client has no SigV4 signing** — it builds a plain `OpenAI(base_url, api_key)`. If the
  SigV4 path ever becomes necessary, add a `bedrock` provider to `llm/client.py` that constructs
  `BedrockOpenAI`; the rest of the call path (retries, rate limiting, usage logging) is shared.
- Mantle's `/openai/v1` route is OpenAI-models-only and rejects `zai.glm-5` with "isn't supported on
  this route". Use plain `/v1`.
- GLM-5 cost measured: **$0.0022/call**, so a 60-call warm start is ~$0.13. Bedrock lists GLM-5 at
  ~$1/M input, $3.20/M output (third-party aggregators — verify on the AWS page).
- Bedrock API keys **expire** (short-term, ~12h). A 401 `invalid_api_key` / "security token ... is
  invalid" means regenerate; that is a different failure from `access_denied`.
- Bedrock has `zai.glm-5`, **not** the paper's GLM-5.1. Report it as GLM-5.

---

## 8. Operational rules (learned the hard way)

**`AF + DepMap` ON THE REDUCED ATLAS IS THE STANDARD METHOD. The user's standing
instruction, 2026-10-28.** AssayFormer plus the DepMap co-essentiality ranking at
weight 0.5, with NO screen-kNN and NO BioPlex, is the project's reference configuration
and the thing every other arm is described relative to. Its weight is confirmed at 0.5
both inside the two-term pair (§4.60a/b) and on its own (§4.60c); every term added on top
of it is within noise on the held-out folds (§4.57e, §4.58). In
`compare_all_results.ipynb` it is the global `STANDARD`, it has its own registry group
`standard`, and `is_standard()` draws a heavy dark outline around it in every bar chart.
**Note what it is NOT: the arm with the highest point estimate.** It loses to
`AF + DepMap + screen-kNN` on `test` (+0.315 crossing zero against +0.506) and to the
BioPlex stacks by ~0.03 on three sets. It is standard because it is the simplest
configuration that holds up, not the top of a leaderboard -- and §4.57(e) is the standing
argument for keeping the screen-kNN term shipped alongside it.

**NO "OUR METHOD" IN LABELS OR PROSE, and every arm named by its components.** Also the
user's instruction, same date. `AF` is AssayFormer unmodified, `+ DepMap` /
`+ screen-kNN` / `+ BioPlex` name the fused terms, and `->` marks a handoff, so
`GLM-5 + top-up warm start -> AF + DepMap` says exactly what ran. **One trap this
created:** `stem(label)` is `label.split(" (")[0]`, because the sweep id is the only
thing allowed in parentheses -- so a name containing " (" is silently truncated and
loses its `DESCRIPTIONS` entry. A first pass used
`"GLM-5 warm start -> AF (the paper's design)"` and eight methods lost their
descriptions. Use commas, never parentheses, inside a label stem. **A second trap:** the
names now NEST (`AF + DepMap` is a substring of `AF + DepMap + screen-kNN + BioPlex`),
which is safe only because `pick()` tries an exact method name before any substring.

**REDUCED ATLAS ONLY. The user's standing instruction, 2026-10-28: every future run uses
`achilles_noeval.npz`, and nothing is ever compared against a full-atlas arm again.**
Enforced in three places so it cannot be forgotten: the launchers default to
`ATLAS=achilles_noeval.npz` (§4.56), `compare_all_results.ipynb` has
`REDUCED_ATLAS_ONLY = True` which drops every `atlas == "full"` arm from all 22 figures
and from the CSV, and each arm's `config.json` records which matrix it read so the claim
is checkable after the fact. **The cost is named rather than hidden: 11 methods have NO
reduced-atlas twin and their results are currently unsupported** -- the three
richer-description arms, gene-gene co-hit, cell-line expression, the Flash-Lite hybrid,
and the five oracle/injection bounds. The notebook prints that list on every run. Re-run
any of them before citing it; all are free and CPU-only.

**Compute every load-bearing figure twice, from different code.** §4.59: `score_pair.py`
spent two sessions reporting plain screen means as family-averaged deltas with intervals
too narrow, and nothing caught it -- the numbers were plausible, self-consistent, and
matched this document's own recorded history, because that history came from the same
script. One figure refusing to reconcile with the notebook is what exposed it. Both
implementations now print a family count; cross-check them whenever either changes.

**Never redirect a sweep's stderr to /dev/null.** `run_codep_confirm.sh` originally did, and 286 of
334 `test`-fold screens failed silently on the embedding-cache error above; the runs looked like a
completed sweep with sparse output. Also check the run count against the fold size, not just that
directories exist. A related trap in the same script: a skip-guard of `output/runs/${id}-0*` only
matches indices 00-09, so an arm that had produced one screen was skipped entirely on the retry;
use `${id}-[0-9][0-9]-*`.

**A single Gaussian noise draw is not a control.** When screening an added feature, the margin over a
width-matched random block is the only interpretable number — but ONE draw is itself a noisy estimate
of the capacity baseline. Changing only the order the controls were drawn moved margins by ~0.1 EF.
Average at least 5 draws, seeded per (screen, feature) so results cannot depend on how many features
ran earlier. **§4.9's 50-column numbers (+0.60 / +0.84) predate this fix and carry more uncertainty
than stated there.**

**Inventory `output/runs/baseline/` before running anything.** Finished sweeps live in topic
subfolders (`baseline/gemini/` ~34 sweeps, `baseline/glm_51/` 8, `baseline/categories_8_first_run/` 3),
and folder names carry only the sweep hash plus screen name — identify a sweep by reading
`config.json -> config.tag`, `screen_set` and `n_steps`. On 2026-10-21 I re-ran an AssayFormer-alone
validation campaign that was already there as `sweep-fd2017be`. One-liner:

```bash
for d in output/runs/baseline/*/; do for s in $(ls "$d" | sed 's/-[0-9][0-9]-.*//' | sort -u); do
  f=$(ls -d "$d"$s-*/ | head -1); python3 -c "
import json; c=json.load(open('$f/config.json'))['config']
print('$s', c.get('n_steps'), c.get('screen_set'), c.get('tag'))"; done; done
```


1. **Never run billed model calls unless the user explicitly asked for that specific run.** Hand over
   the command instead. Free actions (reading runs, scoring, replay-mode handoffs, AssayFormer-only
   runs, writing configs) need no permission.
2. **Finished runs live in `output/runs/baseline/`** — the user moves them there. Search recursively.
3. **`--resume` only sees `output/runs/` directly**, not `baseline/`. Keep a sweep's folders there
   until any dependent handoff has run, then move them.
4. **`--warm-prefix` needs its trailing dash** (`sweep-glm5val-`), or it matches sibling sweeps.
5. **`assayloop run` writes the API key into every run folder's `config.json`** (in
   `acq_overrides.llm`). `handoff_*.py` does not. Scrub with `check_secrets(fix=True)` in section 0 of
   the notebook; rotate the key if the folder was shared.
6. **Check `--n-llm` in the recorded tag before comparing mixed runs.** Two runs were made at
   `--n-llm 50` when 20 was intended (`sweep-mixval-livecat-enrich`,
   `sweep-mixval-livecat-llm50-high`), which confounds them against the 20/batch baselines. The
   acquisition name records it: `mixed_llm50_r3_...`.
7. Set `--tag` on every run. Without it the tag defaults to the model name and the notebook groups
   unrelated runs together.

---

## 9. Sweep reference

Baselines:

| Sweep | Set | What |
|---|---|---|
| `sweep-fd2017be` | validation | **AssayFormer alone**, shared f2 universe (adjusted EF **6.100**) — in `baseline/gemini/`, tag `val-af-alone-fullgenome` |
| `sweep-afalone-val` | validation | redundant re-run of the above (2026-10-21); reproduced it digit-for-digit, which is the determinism check |
| `sweep-00866f93` | test | **AssayFormer alone**, shared f2 universe (EF 4.83) |
| `sweep-4bcd8719` | validation | AssayFormer alone, screen-library pool (EF 6.04) |
| `sweep-e388ac40` | test | AssayFormer alone, library pool (EF 4.67) |

`codep_seed` campaigns (2026-10-21/22, all in `output/runs/`, no API calls). Every family's `-w0`
arm is its own nested baseline and reproduces AssayFormer pick for pick.

| Sweep | Set | What |
|---|---|---|
| `sweep-cdv-w{0,05,1,2}`, `-w1-zsum` | paper_validation | the original weight sweep, `query=centroid` |
| `sweep-cdt-w{0,05,1,2}`, `-w1-zsum` | paper_test | same; `w05` is the +0.540 of §4.13 |
| `sweep-cf-test-{w0,ess-w05,expr-w05}` | `test` FOLD (334) | held-out confirmation — **only 48 of 334 screens ran**, see §7 |
| `sweep-cf-validation-{w0,ess-w05,expr-w05}` | `validation` FOLD (218) | held-out confirmation, complete; the 198 fresh screens are the §4.14 headline |
| `sweep-tf-pt-{w0,knn05,both05,both-c5k25}` | paper_test | adds the screen-kNN scorer |
| `sweep-tf-validation-{...}` | `validation` fold | same, held out |
| `sweep-dm-pt-{w0,dm05,dm05knn,dm10}` | paper_test | `query=diffmeans` (§4.15) |
| `sweep-dm-validation-{...}` | `validation` fold | same, held out — **the number to quote** |
| `sweep-ch-validation-{h05,h025,eq}` | `validation` fold | gene-gene co-hit as a THIRD fusion term (§4.21): depmap/cohit/knn at 0.5/0.5/0.5, 0.5/0.25/0.5 and 1/3 each. Baseline + reference are the existing `sweep-dm-validation-{w0,dm05knn}` |
| `sweep-ms1-validation` | `validation` fold | the deployed arm with `min_seeds=1` (§4.20, §4.23) |
| `sweep-wt-validation-{c05k025,c05k075,c05k10,c025k05,c075k05}` | `validation` fold | the §4.19 joint weight sweep, 218 screens each, `query=diffmeans`. Baseline and the 0.5/0.5 reference are the existing `sweep-dm-validation-{w0,dm05knn}` |
| `sweep-codep-smoke0`, `-smokeexpr`, `sweep-cf-diag*` | single screen | smoke tests and the embedding-cache diagnosis; not results |

Score them with `score_ctx.py` (`--base <family>-w0 --base-dir output/runs --runs-dir output/runs
--glob ...`) or, for a fold, `codep_gated.py --exclude-set paper_test|paper_validation`.

**There are no shrinkage sweeps, and that is deliberate** — §4.16 screened it offline, it came out
negative on both screen sets with controls, and an offline-negative candidate does not get a
campaign. Reproduce the screen with `build_dir_bank.py` then `screen_shrinkage.py`; both are free
and take a few minutes.


Test-time interventions on AssayFormer (2026-10-21, all in `output/runs/`, no API calls):

| Sweep | Set | What |
|---|---|---|
| `sweep-ctx-base` | validation | nested-baseline sanity: defaults, scores 6.100 |
| `sweep-ctx-ens{8,16}-s{50,80}[-rrf]` | validation | context ensembling, K and subsample as named |
| `sweep-ctx-drop{1,25,35,50,75}` | validation | drop weakest q of revealed hits (`drop1` = q 0.10) |
| `sweep-ctx-neg25` | validation | weak hits flipped to negative, q 0.25 |
| `sweep-ctxtest-drop{1,25,35,50,75}` | test | the same dose-response, confirmatory |
| `sweep-ctx-noctx` / `sweep-ctxtest-noctx` | val / test | `ignore_context=true` — the context-value ablation |
| `sweep-ctx{,test}-biasonly` / `-nobias` | val / test | **`score_mode=bias_only\|no_bias` — the ablation ladder** |

The ceiling itself is not a sweep: `scripts/diagnostics/oracle_u.py` computes it directly
(`--screen-set`, `--base-sweep`, `--base-dir`, `--l2`).

Baseline for the test-set rows is `sweep-00866f93` in `baseline/gemini/`; score everything with
`scripts/diagnostics/score_ctx.py` (see §5 for its flags).

Warm starts (rounds 1-3 only; feed the handoffs):

| Sweep | Set | Model / method |
|---|---|---|
| `sweep-gem35val` | validation | Flash-Lite, `llm_single` |
| `sweep-gem35test` | test | Flash-Lite, `llm_single` |
| `sweep-glm5val` / `sweep-glm5test` | val / test | GLM-5, `llm_single` |
| `sweep-glm5valtopup` | validation | GLM-5, `llm_single_topup` |
| `sweep-glm5testtopup` | test | GLM-5, `llm_single_topup` (NEW, $0.36) |

Paper baselines, run here for the first time (§4.55):

| Sweep | Set | What |
|---|---|---|
| `sweep-paperknn-pt` | paper_test | the paper's GENE-kNN over GenePT, standalone (EF 1.93) |
| *(no sweep)* | paper_test | the paper's `Screen-kNN`, standalone (EF 3.40) via `_screen_knn_rollout` |
| `sweep-gem35cat3` / `sweep-gem35cat12n3` | test | Flash-Lite, `llm_category` 8 / 12 cats |
| `sweep-gem35topup3` | test | Flash-Lite top-up, 15 screens only |
| `sweep-gem35test10` | test | Flash-Lite, 10 rounds, LLM-only (EF 2.19) |

Handoffs:

| Sweep | Set | What |
|---|---|---|
| `sweep-7d0a513b` | validation | Flash-Lite -> AF, n=3 (5.17) |
| `sweep-d7630374` | validation | GLM-5 -> AF, n=3 (5.15) |
| `sweep-8f95fe0c` | validation | GLM-5 + top-up -> AF, n=3 (5.22) |
| `sweep-valhandoff-enrichonly` | validation | Flash-Lite -> AF, 10 enrichment screens only |
| `sweep-84fa85f8` | test | Flash-Lite -> AF, n=3 (4.57) |
| `sweep-a7dd71ea` | test | Flash-Lite top-up -> AF (4.70) |
| `sweep-0d295c6a` / `sweep-d8d10355` | test | categories 8 / 12 -> AF (4.68 / 4.62) |
| `sweep-gem35decline` | test | decline-triggered handoff (4.90) |
| `sweep-da9a8c7e` | test | GLM-5 -> AF, n=3 (4.65) |
| `sweep-glm5tu-pt-af` | test | GLM-5 + top-up -> AF, n=3 (4.81) |
| `sweep-hybrid-{pv,pt}-gem35codep` | val / test | Flash-Lite -> OUR METHOD (5.35 / 5.12) |
| `sweep-hybrid-{pv,pt}-glm5codep` | val / test | GLM-5 -> OUR METHOD (5.36 / 5.18) |
| `sweep-hybrid-pv-glm5tucodep` | validation | GLM-5 + top-up -> OUR METHOD (5.43) |
| `sweep-hybrid-pt-glm5tucodep` | test | GLM-5 + top-up -> OUR METHOD (5.23) |

**ABLATED-ATLAS ARMS (§4.56). These are the arms the notebook now registers** -- the
`achilles.npz` originals beside each one are superseded but still on disk:

| Sweep | Set | What | supersedes |
|---|---|---|---|
| `sweep-nolk-{pv,pt,validation,test}` | all four | DepMap + screen-kNN (6.23 / 5.36 / 5.02 / 6.92) | `sweep-pv-dm05knn`, `sweep-dm-pt-dm05knn`, `sweep-dm-validation-dm05knn`, `sweep-tf-test-dm05knn` |
| `sweep-dmne-cd-{pv,pt,validation}` | three | DepMap only (6.18 / 5.53 / 5.01) | `sweep-cdv-w05`, `sweep-dm-pt-dm05`, `sweep-dm-validation-dm05` |
| `sweep-dmne-ppi-{pv,pt,validation}` | three | + BioPlex at 0.25 (6.30 / 5.40 / 5.04) | `sweep-ppi-{pv,pt,validation}-w025` |
| `sweep-dmne-hyb-glm5-{pv,pt}` | two | GLM-5 -> OUR METHOD (5.36 / 5.20) | `sweep-hybrid-{pv,pt}-glm5codep` |
| `sweep-dmne-hyb-glm5tu-{pv,pt}` | two | GLM-5 + top-up -> OUR METHOD (5.44 / 5.25) | `sweep-hybrid-{pv,pt}-glm5tucodep` |
| `sweep-dmcdbase-{pv,pt,validation}` | three | **full-atlas CPU twin of DepMap-only**, the paired baseline for row 2 | — |

Mixed-batch (all validation; `llm20`/`llm40` = LLM genes per batch, `h*` = AF horizon):

`sweep-mixval-llm10-h0`, `-llm10-hrem`, `-llm20-h0`, `-llm20-h100`, `-llm20-h300`, `-llm20-h500`,
`-llm20-hrem`, `-llm40-h0`, `-llm40-hrem`, `-llm100-h0` (sanity: reproduces the handoff),
`-llm20-h0-enrichonly`, `sweep-mixval-livecat-llm20` (live categories),
`sweep-mixval-livecat-llm50-high` (**n_llm=50**, high effort),
`sweep-mixval-livecat-enrich` (**n_llm=50**, enrichment prompt).

---

## 10. What I would do next

**Both results are now well supported. The candidate list is NOT closed** — it was closed against new data sources only, and §12 (2026-10-23) reopens two untouched axes, the atlas GEOMETRY and the statistical INSTRUMENT, plus a disqualified positive control. **Read §12's ranked queue alongside item 1 below.** The negative — the
paper's LLM warm start does not reproduce — was established by refuting three hypotheses
experimentally rather than by argument: model capability (GLM-5 ~ Flash-Lite), unfilled slots
(top-up fixed the shortfall, EF unchanged) and prompt engineering (nine variants, none above
parity), with a ceiling argument explaining why none could have worked (§4.7, §4.12). The positive
— `codep_seed` at **+0.175 EF per screen / +0.331 per family** held out, **+0.485 per family** on
genome-scale screens — survives weight confirmation (§4.19), parameter confirmation (§4.20) and
three leakage audits (§4.22 training labels, §4.35 DepMap own-line ablation, text embedding by
construction).

**ONE LIVE ITEM — ~~DONE 2026-10-23, §4.46.~~ The endpoint replicated on 314 fresh test-fold
screens (+0.257/screen, +0.398/family; genome-scale +0.484/family against validation's +0.485) and
the budget curve did NOT: the method is behind from 300 to 800 picks on that fold. See §4.46; the
open question is now WHY the payoff lands in rounds 2-4 on one fold and rounds 9-10 on the other.**

1. ~~A second held-out confirmation on the `test` fold.~~ 314 more screens, turning one independent
   confirmation into two on the headline result. Needs `OPENAI_API_KEY` for 286 screen-description
   embeddings (§7) — well under a cent, one-time, and **not an LLM call**: it fills a preprocessing
   input AssayFormer has always required. Everything else is free. **This is the only action left
   that strengthens what the project has**, and it is worth more than any further candidate.

   **BLOCKED ON BILLING as of 2026-10-23, not on code or credentials** — §7 records all three gates.
   The key and model access are live; the account has no balance.

   **THE FREE PARTIAL, and read its limits before running it.** 28 fresh test-fold screens are
   already in the shipped cache (18 genome-scale), so they need no key. That is a genuinely
   INDEPENDENT sample, but at n=28 the standard error is ~2.7x the 198-screen fold's, so a true
   effect of +0.175 would come back with a CI roughly [−0.11, +0.46] — **it cannot exclude zero and
   must not be reported as a confirmation.** Its honest value is asymmetric: a positive point
   estimate is mild corroboration, while a clearly NEGATIVE one would be a warning worth having.
   Also note the 48 cached screens are whatever the authors happened to embed, plausibly the ones
   in the paper's figures, so they are not a random subset of the fold. Frame any result as a
   consistency check, and do not repeat §4.18's mistake of reading an underpowered null as absence.

**CLOSED, with the measurement that closed each — do not reopen without a genuinely new argument:**

2. ~~Get more of the profile into the campaign~~ — closed from both directions. §4.15: the 50-column
   offline margin meant a better single DIRECTION existed, not that capacity was short, and
   difference-of-means captured it. §4.16: hierarchical shrinkage toward similar training screens is
   refuted as a blend and as a separate column, because `cos(own direction, borrowed prior)` is only
   ~+0.24 — a neighbour's direction is a biased estimate and the bias swamps the variance removed.
   Explicit-penalty shrinkage has no interior to win in: low lambda approaches an online fit (§4.9
   measured a pure confound control at +0.58 EF), high lambda recovers §4.15.
3. ~~Hunt for more additive features~~ — **eleven candidates screened, one winner.** DepMap
   co-essentiality works; expression is harmful (§4.24, −0.109); Perturb-seq (§4.17), co-hit
   (§4.21), BioPlex PPI (§4.28), PRISM (§4.27), condition-matched transfer (§4.34), PPMI/SVD
   variants (§4.32), direction-filtered co-hit (§4.31), lineage-restricted DepMap (§4.37) and JUMP
   ORF (§4.39) all fail, each on a named condition. **"More training screens" is not available:
   AssayBench already IS BioGRID ORCS**, 1901 of 1952 human screens (§4.33).
4. ~~Tuning~~ — every knob is already optimal: both fusion weights peak at exactly 0.5 (§4.19),
   `min_seeds` at 2 in both directions (§4.20, §4.23), the saturated-screen gate is a no-op (§4.13
   correction), and a rate-adaptive weight has no basis (§4.23, Spearman −0.045).

**The one untested idea of any kind**, for completeness: **batch diversity in DepMap profile
space** — pick each round's 100 by a diversity-penalised objective rather than top-k of the fused
score. ~1-2 hours plus tests, ~10 minutes of compute, nested at `diversity=0`. It is the only
remaining lever that changes the SELECTION RULE rather than a score or a data source, so no existing
refutation covers it. Expect flat-to-slightly-negative: it deliberately declines the
highest-scoring genes, which is the §4.13 displacement mechanism, and §4.15 found the deployed
spread between four direction estimators was only ~+0.03.

4b. ~~DECIDE WHAT THE SCREEN-kNN HALF IS WORTH.~~ **ANSWERED 2026-10-28 by the
   symmetric grid, §4.57: it is worth nothing on the default atlas, and neither is
   theirs.** The eight arms ran, round 1 passed 218/218 on all seven gated arms and
   0/218 on the ungated control, and the resolution is the CHEAP option that was on the
   table all along — restate §4.55 as a null on both sides, with an ordering. Ours
   `−0.014 / −0.022 / −0.170 / −0.217` against theirs `−0.059 / −0.179 / −0.375 /
   −0.673` *(re-scored, §4.59)*; ours is better at every weight and the gap widens with
   dose, but **every cell on both sides is negative** — fused as a standalone ranker on
   `validation`, ours costs a little and theirs costs more. The nesting defect was real
   and cost 0.026. **The deck's screen-kNN slide is updated and no longer quotes
   +0.087.**

   **THE FOLLOW-ON RAN THE SAME DAY AND IS ALSO ANSWERED — §4.57(e). KEEP BOTH TERMS.**
   `validation` made the method look like one component (two-term 5.017 vs DepMap alone
   5.014). The `test` fold refuses that: the kNN increment is +0.191 [−0.019, +0.415]
   there, and against AssayFormer alone the two-term arm reads **+0.506 [+0.159, +0.945]
   against DepMap alone's +0.315 [−0.025, +0.747], which crosses zero**. Dropping the
   term costs that fold's significance and buys nothing on the other. The half is
   **unresolved and negative on three of four sets** — point estimates −0.170 to +0.191,
   not one interval excluding zero — so ship it and describe it as unresolved, never as a
   contributor.

   **Related open question, unchanged:** DepMap alone GAINS +0.047 from the ablation
   while the fused arm LOSES −0.036, so the fusion absorbs the atlas change differently
   from the scorer alone — the §4.48-51 near-substitutes shape again.

4c. **CROSS-CHECK EVERY LOAD-BEARING FIGURE AGAINST A SECOND IMPLEMENTATION, AND FIX
   §4.54(b).** §4.59 found that `score_pair.py` had been reporting plain screen means as
   family-averaged deltas, with intervals too narrow, for two sessions — and nothing
   caught it except one figure refusing to reconcile with the notebook. Every number it
   produced was plausible, internally consistent, and matched the project's own recorded
   history, because that history came from the same code. **§4.54(b)'s proximity-
   correction table has NOT been re-scored and its intervals should be treated as too
   narrow until it is** (`sweep-prox-*` vs `sweep-cpubase-*`, free). More generally: the
   notebook and `score_pair` now both print a family count, and any contrast that carries
   a conclusion should be computed from both before it is written down.

5. **Write it up, around the ceiling decomposition, §4.13 and §4.14.** Nothing further is needed for a defensible
   result. Lead with §4 diagnostic 6: the whole value of AssayFormer's label feedback is ~0.9 EF on
   test, which makes the paper's +0.84 warm-start claim implausible on its face and explains all
   eleven null interventions at once. Support it with the within-screen redundancy measurement
   (62-67%), the selection-direction analysis, and the insensitivity result (dropping three-quarters
   of revealed hits from the context moves EF ~0.1). Four interventions that looked promising were
   each refuted by experiment rather than argument (§4 diagnostics 2, 4, 5).
6. **Expand around confirmed hits** — only in EARLY/MID rounds, and only at low dose.
   The terminal-round version is refuted (§4 diagnostic 4): at round 8 the pool is depleted and
   the LLM's picks collide with what AssayFormer already took. If tried at all, use
   `--mix-start` with 10-20 LLM genes per batch from round 4-6, where AssayFormer's marginal
   rate is lower but the pool is not yet exhausted. Expect parity; the value is the clean test.
   Original framing: In rounds 2+, show
   the LLM the confirmed hits from both components and ask for complex partners, pathway neighbours and
   paralogs that AssayFormer does not rank highly. AssayFormer reveals ~40 hits per batch on validation (**the CURATED 20 only -- on the fold the median is 3 by round 3, §4.48**)
   and the LLM currently ignores them. Pair with `--mix-start` to mix in *later* rounds, where
   AssayFormer's marginal rate has fallen to ~29% and the LLM would have hundreds of labels to work from.
7. **AssayLLM** (`Genentech/assayllm-sft`, `Genentech/assayllm-handoff-grpo`, 55.6 GB each, ungated).
   GRPO with a context-delta reward optimizes marginal contribution given history — the objective all
   the prompting was trying to approximate. Needs a rented GPU (~$5; local machine has 4 GB VRAM,
   14 GB RAM, 51 GB free disk — the weights do not fit). Caveat: trained on the same 1349 screens as
   AssayFormer, so it may be redundant for a different reason.
8. ~~**Gemini 3.1 Pro**~~ **LARGELY ANSWERED 2026-10-27 WITHOUT RUNNING IT (§4.52).** The paper's
   own per-round workbook gives its numbers on these exact screens: its Gemini-3.1-Pro handoff
   scores **5.48 [5.265, 5.699]** and **our method scores 5.389 on the same 20 screens with no LLM
   at all** — inside its band. Gemini 3.1 Pro ALONE scores **4.774, below AssayFormer's 4.828**.
   So the paper's +0.84 is reachable, by its model, through the handoff structure rather than
   through the LLM being a good picker — and it is matched here without one. Running it ourselves
   would now test REPRODUCIBILITY of their number, not whether the gain exists. Worth money only
   if that reproduction is the question being asked.

**Do not rebuild:** the LLM as an additive per-gene feature, at any reasoning effort or model size (§4 diagnostic 11: informative but fully redundant with span(V)); DepMap essentiality or expression as additive features, in pan-cancer or cell-line-matched form (§4 diagnostic 10: inside the noise band and sign-flipping across the two screen sets); any scheme for feeding hit strength to AssayFormer — soft hit embeddings, duplicated strong hits, strength-filtered context (§4 diagnostic 8: worth -0.007/-0.003 EF at the oracle, on both sets, under both estimators); context ensembling over û, or denoising the context with the revealed `relevance_score` (§4 diagnostic 5: ensembling is flat to negative at every setting; denoising has no dose-response and its two significant cells point in opposite directions across the two screen sets — the +0.122 on test is one cell of ten, chosen after looking); the pathway-propagation scorer or any RRF fusion of it, at any round (§4 diagnostic 4: it loses to AssayFormer on 16/20 screens, late seeding is weaker than early, and no deployable gate nets positive); a reversed handoff that puts the LLM in the terminal rounds (same diagnostic, plus the already-tested collision penalty); essentiality gating on enrichment screens (measured backwards, §4);
hierarchical shrinkage of the profile direction toward similar training screens, in blended or
separate-column form, and explicit-penalty variants of it (§4.16: the margin declines monotonically
in the shrinkage weight on both screen sets, and the shuffled-bank control scores as high as the
real bank, so the two-column form's apparent lift is capacity);
a do-not-nominate list from historical hit frequency (210 validation genes are
hits in >=8 screens and found in none — excluding them guarantees 0% recovery); any handoff trigger
referencing `|H|` or the random hit rate (not deployable, leaks the answer); diversity selection in
BPMF space (that is AssayFormer's own geometry); graded hit labels into AssayFormer's history tokens
(it learned exactly two label embeddings).

---

## 11. Design spec: pathway-propagation scorer + RRF — RETIRED, kept as a record

> **Do not build this.** §4 diagnostic 4 (2026-10-21) measured the scorer against AssayFormer
> instead of against random and it loses on 16/20 screens at every seeding horizon. The spec is
> kept because the reasoning below is correct as far as it goes — the flaw was the baseline, not
> the mechanism — and because the next person to have this idea should be able to see exactly
> what was tested and why it failed rather than rebuilding it.

This was §10 item 1 written out mechanically. Nothing here was implemented — the 2026-10-21 diagnostic that produced the 3.53x number
was a throwaway script in the session scratchpad and **is no longer on disk**; it has to be rewritten.
The Reactome inputs it used are permanent, in `src/assayloop/data/gene_sets/`.

### Why it can be additive rather than redundant

AssayFormer scores `a(g) = û·v_g + b_g` with `û ∈ R^10`. Appendix F.3: the BPMF gene factors `v_g`
organize **by screen phenotype, not by pathway or complex membership**. So "these two genes sit in
the same small Reactome module" is information the model does not encode. That is the whole
complementarity argument, and it is why this is worth building where nine prompt variants were not.

### The scorer

A re-ranker over untested genes, recomputed each round. No training, no LLM calls, ms per round.

1. `H` = genes revealed as hits in rounds `1..t-1`; `T` = all genes tested so far.
2. For each `g ∈ H`, look up its Reactome modules, restricted to pathways of **2-150 genes**. The
   size cap is load-bearing: without it, generic sets ("Metabolism") link everything and the median
   lift collapses from 2.31x to 1.08x.
3. For each pathway `p`: `c(p) = |H ∩ genes(p)|`. Keep only `c(p) >= 2`. Requiring two independent
   hits in one module does most of the work; a single-hit seed is too weak.
4. Score each untested gene `g`: `score(g) = Σ_{hot p ∋ g} c(p) / |genes(p)|`. The denominator is the
   specificity term — 3 confirmed hits in a 12-gene complex should outrank 2 in a 150-gene pathway.
5. Optionally add PPI adjacency to confirmed hits via `../candidate_retrieval.py` (symbol-level,
   187 cached entries). Complex co-membership is more specific than pathway co-annotation.
   `reactome_interactions.txt` cannot be used directly: UniProt/Ensembl/Entrez ids, no symbols.

Properties: blind in round 1 (no hits yet) and strengthening as labels accumulate — the mirror image
of the LLM, which is strongest at round 1. It therefore acts in **rounds 4-10, where no intervention
tried so far has operated at all.** Failure mode is tunnelling (reinforcing one module, ignoring the
genome), which is the reason to fuse rather than replace, and a reason to consider a per-pathway cap.

### Why RRF and not a weighted sum

The two scores are on incompatible scales — an unbounded bilinear value vs a sparse sum of ratios
that is **exactly zero for most of the genome** (only ~465 genes/screen scored nonzero in the
diagnostic). Calibrating them (z-score, min-max) is fragile because the pathway score's distribution
changes every round as hot modules appear and disappear; a weight tuned at round 4 is wrong at round 8.

### The fusion

```
RRF(g) = Σ_s  w_s / (k + rank_s(g))        k ≈ 60
```

- Ordinal only, so the scorers never need to agree on units.
- Steep at the top, flat in the tail: rank 1 vs 10 matters, rank 5000 vs 8000 does not. Since only
  the top 100 is taken, that is where the resolution belongs.
- `k` sets the temperament: small `k` is winner-take-all, large `k` demands agreement.
- **`w_pathway = 0` reproduces AssayFormer exactly.** The baseline is nested, so the worst case is
  parity rather than regression. Keep this property — it is what makes the experiment safe.

**Sparsity subtlety (the main implementation trap):** most genes score 0, so ranking all ~21k genes
by the pathway score makes one enormous tie block. Genes with `score > 0` get a pathway rank;
everything else contributes **nothing** from that scorer, not `1/(k + 21000)`. The scorer then acts
as a focused promoter of a few hundred genes while AssayFormer keeps ranking the full pool.

### Evaluation

- Enabled on **enrichment screens only** (`screen_direction.py`), per §3.
- **EF on all hits is the headline metric.** The strong-hit view (top quartile by relevance) is
  diagnostic only and must not be substituted for the paper's metric.
- Expect a gain smaller than 3.53x standalone: AssayFormer still fills most of each batch, so this
  reorders 100 slots rather than replacing the policy. 5 of 16 enrichment screens showed no lift
  (1815, 1818, 2054, 2104), so report per-screen, bootstrapped over screen *families* (§6).
- Scoring is free: it can be replayed offline against the existing AssayFormer-alone runs in
  `output/runs/baseline/` with no API spend.

---

## 12. Triage of the 2026-10-23 external proposals

Two independent reviews arrived as lists of candidate improvements — call them **list A** (10 items,
headed "Enhanced Screen Description Engineering") and **list B** (10 items in three tiers, headed
"follows directly from findings you already have"). Both were scored against §4. **List B's figures
were checked and are all transcribed correctly** (r ≈ 0.11 inter-column, §4.18; +0.044/+0.066
Perturb-seq contrast, §4.18; +0.024 hard / +0.190 all for BioPlex, §4.28; the min_seeds conditional
table, §4.20; the 808 -> 34 cell-line cost, §4.37). **List A cites no numbers from this project and
eight of its ten items are already closed here**, several by the measurement that most directly
contradicts them. Keep this section so neither list gets re-proposed.

### The two axes that were genuinely never touched

This is the part worth internalising, because it explains how ten candidates could be spent without
exhausting the space. Every one of §4's eleven feature candidates varied the **data source**. Two
other axes were never varied at all:

1. **The geometry.** `codep_seed` computes cosine on raw 808-dimensional dependency profiles. Those
   are dominated by a few large shared components — growth rate, screen quality, lineage. Difference
   of means cancels some of that in the *direction*; nothing cancels it in the *space*. PC removal
   and variance weighting are standard in the co-essentiality literature and appear nowhere in 97
   rows of `methods_summary.csv`.
2. **The instrument.** Every result is adjusted EF at one budget (1000), unadjusted for per-screen
   covariates, on a 198-screen fold that has served as screening set and confirmation set at once.
   **Three candidates were closed on unmeasurability** (§4.18, §4.21, §4.28) — when three of your
   own verdicts are "below resolution", resolution is the binding constraint, not the candidates.

### Ranked queue

| # | item | from | cost | why it earns the slot |
|---|---|---|---|---|
| 1 | ~~Re-run JUMP gate 0 with a GoF-appropriate positive set~~ **DONE 2026-10-23, §4.41** | B7 | hours, free | Criticism upheld: the control was mis-specified and underpowered, and §4.39's gain-of-function conclusion is withdrawn. Atlas still does not reopen — its best ground truth is 4.5x weaker than DepMap's, matching §4.39's own −0.019 gate-3 null. |
| 2 | ~~Re-render the screen description with `notes` + `direction_str`~~ **DONE 2026-10-23, §4.45 + §4.47: NEGATIVE.** | A1, §4.26 | ~$0.002 | `direction_str` turned out to be EMPTY on every screen. The real candidate, `notes`, reads −0.012 [−0.056, +0.034] — excluding any gain above +0.034 — and in round 1, where the text is all `u` has, it costs **−0.860 EF**. Richer text degrades the latent rather than informing it. **The channel is closed; only retraining the encoder remains.** |
| 3 | ~~PC removal / variance weighting on the dependency profiles~~ **CLOSED 2026-10-23 on the premise, §4.44** | B2 | few lines | The axis was genuinely untouched, but the leading component carries only **2.2%** of the variance — Chronos corrects this upstream. Built and nested; campaign not run. |
| 4 | ~~ANCOVA-adjusted deltas + nAUC over the budget curve~~ **DONE 2026-10-23, §4.42** | B5 | free | ANCOVA gives +7% SNR (R² only 0.10-0.23); **nAUC is WORSE** because the effect accumulates over rounds. The valuable by-product was finding that §4.28's rejection used the wrong shrinkage constant — see §4.42 finding 3. |
| 5 | **Weight continuous in DIRECTION STABILITY** | B4 | campaign | §4.20's table is the signature. Distinct from the refuted rule — see the §4.23 scope correction. |
| 6 | **Refit BPMF at K = 32/64/128, recompute R² of the hit vector** | B9 | ~1 hour | DIAGNOSTIC ONLY. Closes or opens the architectural question behind the 3.5% in §4.9. |

**REVISED ORDER, 2026-10-23, after re-examining §4.28.** Asked whether DepMap and BioPlex had been
tried together, I found the pair HAD been measured offline — **+0.190 [+0.002, +0.387]** over DepMap
alone on 81 held-out screens, the only pair contrast besides the deployed one whose CI excludes zero
on a properly powered sample — but **never deployed**, and rejected for two reasons that both turned
out to be wrong (§4.28 corrections: the headroom premise contradicts §4.12, and the sparse-profile
cost dissolves analytically to 1-2 hours). So the practical order is now:

1. **Item 4, the resolution work** — free, from runs on disk, and it is what decides whether a
   ~+0.06 deployed effect can be seen at all.
2. ~~**DEPLOY DepMap + BioPlex + screen-kNN**~~ **DONE 2026-10-23, §4.43: a tight null.**
   +0.000 per screen [−0.031, +0.038] paired against the deployed configuration, 48/198 wins, the
   degree-residualised control identical, the higher dose harmful. Closed on a measurement rather
   than on unmeasurability — the interval excludes every effect above +0.04. **The best remaining
   candidate in the project is now spent**, and §4.42's calibration argument with it.
3. ~~**Item 3, the atlas geometry**~~ **CLOSED 2026-10-23 on the premise, §4.44.** The leading
   principal component of the Achilles matrix carries **2.2%** of the variance, not the dominant
   share the argument assumed — Chronos already corrects screen quality, proliferation and library
   effects upstream. Code, tests and `run_geometry_sweep.sh` are in place and the campaign was
   deliberately not run. **Not a null result; reopen only with a less-processed matrix.**

**Then, in order, and only after item 4 raises the resolution:** fuse all surviving kernels at once
(B1 — right instinct, but it must be parameterised as *added dose* at the deployed 0.5/0.5, not as
renormalised shares: §4.21 measured three scorers at 1/3 as significantly worse, so list A's
0.5:0.3:0.2 predicts a loss and also drops half the deployed method); batch rather than gene
selection (B10 = §10's one untested idea, better after item 3 makes the diversity metric meaningful);
Sanger Project SCORE columns (B3); CCLE protein co-abundance (B6); a learned stacker (B8, deferred —
highest effort, and §4.19's flat weight surface argues there is little for a learned rule to find).

### Where list B overstates

- **B3 is not "squarely predicted to work by your own measurement."** §4.37's 808 -> 34 result
  establishes that cell lines act as sample size, but §4.25 measured the hard screens' hits at
  DepMap amplitude **0.138 against 0.137 for random genes** — genuinely flat profiles, not noisily
  measured ones, since measurement noise inflates amplitude rather than suppressing it. More columns
  give a better estimate of flat. §4.28 already states this for Project SCORE by name. Expect a
  small overall gain from 808 -> ~1100 and ~0 where the headroom is.
- **B2's variance-weighting half cuts against the hard screens specifically** — the genes whose
  profiles are flat *are* the hard screens' hits (§4.25), so down-weighting them suppresses exactly
  that group. Net it should still be neutral-to-positive there (those genes were receiving noise
  scores either way) and positive elsewhere, but PC removal is the half with the clean prior.
- **B5's third part is already §10 item 1.** Splitting the 198 into screening and confirmation halves
  is strictly worse than using the untouched 314-screen `test` fold for confirmation. The
  multiplicity concern itself is legitimate and conceded: `0.5/0.5 optimal` (§4.19) and
  `min_seeds = 2 optimal` (§4.20) were both selected on the same 198 screens that produced +0.175.

### List A: closed, with the measurement that closes each

| list A item | status |
|---|---|
| 1. Richer screen description | **LIVE, and it is queue item 2** — but the proposed mechanism is wrong. You cannot embed `notes` with MiniLM and concatenate to the latent: AssayFormer consumes a fixed-width embedding from one encoder (§7), and a foreign vector is off-distribution for a frozen checkpoint. Re-render the template and re-embed with the SAME encoder (§4.26a). The linear probe adds nothing — predicting screen type from embeddings does not bear on EF. |
| 2. Relevance-weighted / strength-filtered seed profiles | Dominated by queue item 5, which keys on the better quantity. The weighting half is a fifth member of §4.15's direction-estimator family, whose whole deployed spread was ~+0.03; the discarding half is rule 2 (filtering to a subset: three for three against). §4.5 also found no dose-response for revealed `relevance_score`. |
| 3. Pathway-guided exploration, reserved slots | **Closed twice.** §4.1 + §4.4: pathway propagation's 3.53x lift was measured against random picking, and head to head it loses on 16/20 screens. Reserved slots are the §4.13 displacement mechanism. §4's category-coverage reweighting measured **backwards**. |
| 4. Multi-atlas RRF at 0.5 : 0.3 : 0.2 | Right direction, wrong parameterisation — see B1 above. Also "focus on screens where DepMap underperforms" is the measurement §4.28 already made: **+0.024** on the 28 hard screens. |
| 5. LLM summarises hit themes into the description | Not strictly refuted, and the most novel item in list A — but it works *entirely through* û, so §4.7 caps it at **+0.41/+0.84**, the description channel caps it at ~0.9 EF, it costs billed calls, and it is off-distribution in the same way as item 2 without item 2's advantage of containing real experimental metadata. Do item 2 first; this is strictly the worse way to use the same channel. |
| 6. Latent jittering | §4.5: context ensembling over û is flat to negative at every setting, and random jitter inside a 10-dim space whose *oracle* is +0.41/+0.84 is strictly worse than the oracle. Expect negative. |
| 7. Cell-line-aware eligibility filter (CERES + TPM) | **Closed three times over.** §4.2: essentiality gating is BACKWARDS — essential-proxy genes hit at **1.56x** the base rate. §4.24: expression as a kernel is **−0.109**. Rule 2: filtering to a subset, three for three. It also targets dropout screens, where AssayFormer already reaches 73-87% of ceiling and the budget binds (§1). |
| 8. Hit-streak detection, adaptive batch expansion | Batch size is fixed by the benchmark at 100 x 10, so only reallocation *within* a round is available — which is what the kernel already does continuously. A median+2σ trigger is a binary gate on the quantity §4.13 and §4.20 both showed binary gates fail on, and "50% of the next round" is ~100x the tuned dose. Dominated by queue item 5. |
| 9. JUMP-CP as a tie-breaker | Same atlas as queue item 1, used as a feature before the control is repaired. §4.39 measured it: co-complex genes **0σ** against DepMap's +4 to +11σ, and paired against DepMap on hard screens **−0.019**. Fix the control first; the feature only becomes worth building if item 1 passes. (JUMP ORF is U2OS, not pancreatic.) |
| 10. Benchmark-centric negative mining | This is deployed screen-kNN with the sign flipped, and §10's do-not-rebuild list already covers it: 210 validation genes are hits in >=8 screens and found in none, so excluding them guarantees 0% recovery. §4.21's per-gene centring removes marginal hit frequency precisely because it is useless alone. |

---

## 13. External validation: AF alone vs AF + DepMap on seven published papers

Everything in §1–§12 is measured on AssayBench. This section is the out-of-benchmark test:
seven independent papers, 206 directional screens, run through the deployed standard method
and its own weight-0 control.

### 13.1 The harness

`scripts/external/` — one builder per paper plus a shared driver:

| file | what it does |
|---|---|
| `build_schukken.py` | Schukken 2026 WGL fitness screens; also holds the shared `hgnc_synonyms()` / `remap_library()` that every other builder imports |
| `build_abt.py` | Abt 2026 PARCB small-cell prostate |
| `build_wilson.py` | Wilson 2024 MERCS split-mVenus reporter |
| `build_nguyen.py` | Nguyen 2025 AraC / daunorubicin / etoposide drug response |
| `build_shan.py` | Shan 2025 cholesterol synthetic lethality |
| `build_organoid.py` | Herranz-Ors 2026 Sanger organoid biobank, 178 screens |
| `build_tau.py` | Samelson 2026 tau-proteostasis CRISPRi (§13.5) |
| `dump_benchmark_screens.py` | pickles the 348 benchmark screens so they can go through the same driver |
| `run_external.py` | the driver: injects `ScreenRecord`s into `run_one_screen` and runs three arms |
| `stratify_benchmark.py` | benchmark `dEF` stratified by hit count, reusing `adjusted_ef_value` and the notebook's filters |
| `final_table.py` | the unified table below |
| `analyse_external.py`, `synthesise_external.py` | cluster bootstraps |
| `SOURCES.md` | provenance: the exact file, sheet and hit-call column for all seven papers, plus the routes that worked when publishers blocked |

Every raw source file is collected under **`output/runs/external/sources/`** (gitignored).
`scripts/external/SOURCES.md` is the index.

**Why screens are injected rather than loaded.** `run_sweep` resolves screens through
`load_screens`, which can only select rows of the published HF dataset — its YAML manifest must
name that dataset in its `source` block. `run_one_screen(screen, cfg, run_id=, sweep_id=)` takes a
`ScreenRecord` directly, so that is the injection point. `res.history` holds `StepRecord`
dataclasses: use `st.metrics`, not `st.get(...)`.

**Arms.** Both are `model=amortized_ranker`, `acq=codep_seed`,
`acq_overrides={profile_npz: achilles_noeval.npz, fuse: rrf, query: diffmeans, knn_weight: 0.0}`
and differ only in `weight`: 0.0 is AF alone, 0.5 is AF + DepMap. `profile_npz` is required by the
constructor even at weight 0. A third arm, `acq=greedy`, is the independent check that weight-0
really is AssayFormer untouched; it agrees except where OOV ties break differently (§13.4).

**One metric regime for every row.** `universe_genes=list(s.genes)`, so no pick can land outside
the library, `N_Ū = 0`, and `b_eff = 1000` exactly. Adjusted EF collapses to `hits / (1000·H/D)`,
which is the `n_hits_vs_random` the runner already records — no reconstruction, no disagreement
with the notebook. The 348 benchmark screens were **re-run this way** (`benchlib_results.json`) so
they sit on the same footing.

### 13.2 The table

```
dataset                         n  clus  avail H |  AF/scr  +DM/scr   d/scr |   AF EF  +DM EF     dEF              95% CI |        up
AssayBench folds (held out)   348   204      413 |   100.8    102.4    +1.6 |   3.869   4.099  +0.307 [+0.222,+0.398]* |  144/348
Schukken 2026                  14     7     1495 |   726.1    718.8    -7.4 |  10.275  10.277  +0.010 [-0.123,+0.238]  |    3/14
Abt 2026                        1     1     1398 |   431.0    408.0   -23.0 |   4.458   4.220  -0.238  (n too small)   |    0/1
Wilson 2024                     2     1      312 |    76.0     81.0    +5.0 |   4.423   4.775  +0.352  (n too small)   |    2/2
Nguyen 2025                     6     3      872 |   267.0    261.7    -5.3 |   5.071   5.187  +0.116 [-0.008,+0.288]  |    4/6
Shan 2025                       3     2      318 |    75.0     74.3    -0.7 |   4.093   4.171  +0.046  (n too small)   |    1/3
Organoid 2026                 178   162     1384 |   761.9    753.4    -8.6 |   9.407   9.317  -0.090 [-0.106,-0.075]* |   29/178
Tau CRISPRi 2026                2     1      548 |   101.5    116.5   +15.0 |   3.633   4.379  +0.746  (n too small)   |    2/2
EXTERNAL POOLED (7 papers)    206   177     1343 |   720.4    712.4    -8.0 |                  -0.075 [-0.094,-0.054]* |   41/206
```

`_merged` arms are excluded everywhere — they are unions of their own directional screens, not
independent units. Cluster: hit-set family for the benchmark, cell line for Schukken and Shan, drug
for Nguyen, organoid for the biobank, whole experiment for Abt / Wilson / Tau.

### 13.3 What it means, and what does NOT explain it

The benchmark gain is **real and not a metric artifact**: +0.307 [+0.222, +0.398] with
universe = library, i.e. without any forgiveness. The external pool is
−0.075 [−0.094, −0.054], carried almost entirely by the 178 organoid screens.

Three explanations were tested and **all three failed**:

- **Hit sparsity does not predict the gain.** Spearman −0.705 (p=0.0004) over 23 external screens
  looked decisive, but on 348 benchmark screens it is **+0.035 (p=0.52)**. The external correlation
  is four small datasets, not a property of the method.
- **Forgiveness does not inflate DepMap — it deflates it.** On AssayBench folds the pool is
  `gene_universe(screens, min_screen_freq=2)` (~21,213 genes) against a mean library of 14,166, so
  **562 of every 1,000 picks are forgiven** and `b_eff ≈ 437`. The diagnostic shows AF alone gets
  562.5 picks forgiven and AF + DepMap 555.9: forgiveness runs **against** the method by −1.47%.
- **Headroom does not explain it.** `%` of ceiling is not comparable across setups — it is unbounded
  on the benchmark (max 1212%) because perfect play assumes `N_Ū = 0`.

The one real regularity is the **phenotype class**, which the new CRISPRi row sharpens:
the two sorted-reporter datasets are the two largest positive deltas (Tau +0.746, Wilson +0.352,
4/4 screens up), and the two genome-scale fitness/dependency datasets are the negative ones
(organoid −0.090 over 162 clusters, Abt −0.238). A co-essentiality atlas is measuring the same
quantity as a fitness screen, so its ranking is partly redundant there; on a reporter phenotype it
is orthogonal information. **Untested hypothesis** that remains: hit-label provenance (ORCS
curation of reported-significant genes vs authors' own thresholds, which sweep in the
core-essential background). Testable with the common-essential ablation —
`CRISPRInferredCommonEssentials.csv` from the organoid Figshare item plus the existing `ef_ne`.

### 13.4 Traps recorded so none of these is paid for twice

- **`ACH-*` leakage.** `07_EssMatrix_bDepletionsB2.Rdata` on the organoid Figshare item is NOT
  organoids: 17,647 × 930 with column names `ACH-001421` etc. — the DepMap cell-line reference the
  paper compares against. Building screens from it scores the co-essentiality atlas against the very
  cell lines it is computed from, and manufactures a large fake gain. Organoid calls live only in
  Supplementary Table 5, whose `sample_ID` values are `HCM-SANG-*`.
- **Organoid double counting.** The unit is `(sample_ID, library)`, 178 screens. 16 of the 162
  organoids were screened with both MinLibCas9 and Yusa v1.1; grouping by sample alone doubles their
  genes and pushes hit counts past the paper's stated 297–2160 range. The bootstrap clusters on the
  organoid, so a dual-library model cannot count twice.
- **Duplicate gene symbols corrupt hit labels.** `reveal()` silently collapses duplicate entries, so
  every builder raises `RuntimeError` if its gene list is not unique. The HGNC remap also needs a
  collision guard: skip any remap whose target is already in the library (it created 22 duplicates
  before the guard).
- **Resolve by Ensembl ID where available.** It also repairs Excel date corruption
  (`Mar-01`/`Mar-02` for MARCHF genes), which is otherwise unmappable.
- **R mangles hyphens to `__`.** Nguyen's sheet S2 has `NKX3__1`, `HLA__C`, `ATP5J2__PTCD1`; 34 of
  3,748 calls silently fail to find a library row without the unmangle step.
- **OOV genes all tie.** Genes outside AssayFormer's 24,217-token vocabulary map to `<unk>` (id 0)
  and score identically. `codep_seed` tie-breaks on `(-fused[i], str(candidates[i]))`, alphabetically;
  `greedy` does not. On Wilson merged this explained a 158-vs-156 discrepancy, and on Tau merged a
  194-vs-198 one. Both directional arms agree exactly in both datasets.
- **Not-reproducible published counts: use the authors' own labels.** Wilson's 410/230 is unreachable
  at any threshold (MAGeCK RRA FDRs are tied — 357 genes share 1.4e-05, achievable counts step
  394 → 412); Nguyen's "2 SD" rule reproduces nothing close to its six counts. In both cases the
  authors' explicit call columns reproduce the paper exactly.
- **Publisher blocking.** PNAS, ashpublications, cell.com and sciencedirect return 403; PMC serves a
  browser challenge. Routes that work: Crossref, Europe PMC `fullTextXML` / `supplementaryFiles`,
  GitHub raw, Springer `static-content`, Figshare API, PubMed esummary.
- **Check the process, not the file.** Reading a results file while it is still being written twice
  produced a false "the run failed" conclusion, and a `nohup ... &` wrapper returns exit 0
  immediately, so the completion notification refers to the wrapper and not the python process.

### 13.5 The tau CRISPRi dataset (added 2026-10-30)

Samelson, Ariqat, McKetney et al., "CRISPR screens in iPSC-derived neurons reveal principles of tau
proteostasis", Cell 189:1517–1534.e19 (2026), doi 10.1016/j.cell.2025.12.038, PMID 41610849,
Kampmann lab. The article is not open access and its PMC record (PMC12978015) is inside the NIH
embargo — no free route reaches the data, so the user supplied `mmc2.xlsx` (all genes) and
`mmc3.xlsx` (called hits) in the repo root, as with `Sack_dataset/`.

```
PYTHONPATH=src:scripts/external .venv/bin/python scripts/external/build_tau.py \
  --all mmc2.xlsx --hits mmc3.xlsx \
  --hgnc /home/jeff/hgnc/hgnc_complete_set.txt \
  --vocab <assayformer snapshot>/vocab.json --out $S/tau_screens.pkl
PYTHONPATH=src:scripts/external .venv/bin/python scripts/external/run_external.py \
  --screens $S/tau_screens.pkl --n-steps 10 \
  --arms "AF alone,AF alone (greedy),AF + DepMap" --tag tau \
  --out output/runs/external/tau_results.json
```

- **Only the genome-wide sheet is usable.** The other six sheets are follow-up screens on a focused
  984–986-gene sublibrary. The budget is 1,000 picks, so a run would acquire the whole library and
  both arms would score identically by construction. They are read and printed, not scored.
- **TSS-level library.** 20,436 entries → 18,825 genes; 1,567 genes carry two or more start-site
  entries (`HMBS_P1`, `NDUFA9_ENST00000266544.5`). Collapsed to the base symbol, strongest entry
  keeps the score, hit = any. HGNC symbols never contain an underscore, so the suffix is unambiguous.
- **A further 20,436 rows are `NTC_*` quasi-genes**, 5–7 copies each for the library's sublibraries.
  Dropped from library and hit list alike — 45 of mmc3's 1,142 genome-wide rows are NTCs, which is
  the authors' own false-discovery background.
- **The hit calls are NOT reconstructible from mmc2's columns.** `gene score` is
  `epsilon · −log10(pvalue)`, but no threshold on it reproduces the lists: 1,017 genes clear the
  weakest called negative while only 729 are called, and genes at |score| 1.2 are excluded while
  genes at 0.51 are kept. That is the Kampmann rule, whose discriminant threshold is calibrated
  against the NTC distribution **separately per sublibrary**, and the table does not say which
  sublibrary a gene is in. The authors' lists are used as given, as for every other external set.
- **Direction is fixed by an internal control, not guessed.** MAPT itself is the most negative gene
  in the screen (epsilon −2.454, score −9.51) — knocking down tau must lower tau — so negative
  epsilon is a decrease and positive an increase. CUL5 being the second-strongest positive agrees
  with the abstract's CRL5/SOCS4 finding.
- **mmc2's `V337M_TOC1` sheet is a copy of `V337M_M204`.** mmc3's TOC1 hits are genuinely distinct,
  so the workbook's all-gene TOC1 table is missing. Irrelevant here (that screen is not scored) but
  worth knowing before reusing mmc2.
- **Leakage: none.** January 2026, after every year fold, and no screen in any fold mentions tau,
  MAPT, tauopathy, oligomers, iPSC-derived neurons or this group. The benchmark's only iPSC screens
  are a neural-progenitor viral-infection screen (train, `1172`) and a tunicamycin-resistance screen
  (test, `2454`).
- **Result.** `inc` (366 hits) 79 → 107 found, EF 4.063 → 5.503; `dec` (729 hits) 124 → 126, EF
  3.202 → 3.254. 2/2 up, mean ΔEF +0.746. The largest positive delta of any dataset here, on the
  only CRISPRi one — the train fold has 1,303 Knockout screens against 28 Inhibition.

### 13.6 Correction to the headline figure

The **+0.485 / +0.484** two-fold figures quoted earlier in this document belong to
**AF + DepMap + screen-kNN** (`knn_weight=0.5`), not to AF + DepMap. The standard method alone
gives validation **+0.466 [+0.277, +0.678]** and test **+0.315 [−0.030, +0.746]** — the test
interval **crosses zero**. That is the notebook's own stored number, not a re-derivation.

### 13.7 The team's `external_sources.json` subset (2026-10-30)

A teammate's frozen screen list, 11 entries named by internal experiment id with
only `source_genes`, `source_hits` and a PMID/DOI. All 11 resolve to this
project's seven papers. `scripts/external/team_subset.py` holds the mapping and
prints the comparison; `scripts/external/SOURCES.md` records it too.

| team id | our screen | counts |
|---|---|---|
| `E225_ABT_PARCB_DEPENDENCY` | `ABT_PARCB_dependency` | exact (14459 / 1398) |
| `E236_SHAN_HEPG2_LOVASTATIN` | `SHAN_HepG2_lovastatin` | exact (19063 / 368) |
| `E236_SHAN_HELA_LOVASTATIN` | `SHAN_HeLa_lovastatin` | exact (19063 / 359) |
| `E236_SHAN_HELA_LDLR_KO` | `SHAN_HeLa_LDLRko` | exact (19062 / 227) |
| `E239_ORGANOID_COLO_031` | `ORGANOID_WTSI-COLIVM_031_minLib` | exact (16562 / 984) |
| `E239_ORGANOID_OESO_005` | `ORGANOID_WTSI-OESO_005_minLib` | exact (16562 / 1683) |
| `E239_ORGANOID_OVAR_018` | `ORGANOID_WTSI-OVAR_018_minLib` | exact (16562 / 1473) |
| `E226_WILSON_MERCS_HIGH` | `WILSON_MERCS_inc` | hits exact; genes 18458 vs 18464 |
| `E232_NGUYEN_K562_ARAC_KO_SENSITIZATION` | `NGUYEN_AraC_resistant` | hits exact; genes 19111 vs 19114 |
| `E245_SAMELSON_T22_LOWERING` | `TAU_CRISPRi_dec` | 728 vs 729 hits; 18818 vs 18825 genes |
| `E224_SCHUKKEN_PAIRED_WGL` | **had to be built** | see below |

The ±3–7 gene differences are duplicate-symbol collapse differing by 0.04% of the
library; they cannot move EF. `E245`'s name also settles something left generic in
§13.5: the genome-wide readout is the **T22** oligomer antibody.

**`E232` has a naming contradiction worth resolving with the teammate.** Its label
says `KO_SENSITIZATION`, but 240 hits is the paper's *Resistant* arm (knockout
confers resistance, sgRNA depletion); the sensitising arm is 1800 hits. The count
is unambiguous, the label is not.

**`E224` is the paper's PAIRED contrast, which this project had not built** —
`build_schukken_paired.py`, `WGL` sheet (quantile-normalised, 19,112 genes), not
`Raw WGL` (20,111). It also **corrects §13.1's builder docstring**: that docstring
justifies the invented `beta < −0.75` fitness rule by saying the paired cut
"leaves 15 hits genome-wide", but 15 was computed on the 3,377-gene `DG + WGL`
sheet. Genome-wide the paper's own rule leaves **188**. The paired contrast was
usable all along.

Its hit rule is **not recoverable**: 194 comes from no round threshold — not
p < 0.05 with any 2-dp beta cut (188, or 197 if beta is rounded first), not
BH-FDR at any level (**zero** genes survive multiplicity at q < 0.3 — the paired
t-test over seven pairs has no significant gene), not a z-cut on beta (151 at
−3 SD, 296 at −2.5 SD), not p < 0.01 (206). Both readings were therefore built
and run, and **they disagree in sign**:

```
SCHUKKEN_PAIRED_WGL_paper    188 hits   AF 45 -> +DM 42    EF 4.575 -> 4.270  dEF -0.305
SCHUKKEN_PAIRED_WGL_rank194  194 hits   AF 89 -> +DM 100   EF 8.768 -> 9.852  dEF +1.084
```

Only 84 genes are in both hit sets, so this is a genuinely different screen under
the two rules, not a rounding artifact. **Get the teammate's gene-level hit list
before quoting E224.**

Subset result over all 11 (paper-clustered bootstrap, 7 clusters — weak by
construction at this n):

```
E224 = paper's rule   n=11  hits/scr 295.1 -> 291.5 (-3.6)  mean dEF -0.020  [-0.161,+0.075]  up 5/11
E224 = rank194        n=11  hits/scr 299.1 -> 296.7 (-2.4)  mean dEF +0.107  [-0.062,+0.498]  up 6/11
```

Neither interval excludes zero. The subset is 3 organoid screens and 1 Abt screen
out of 11, so it is weighted toward the fitness/dependency phenotypes where §13.3
says the atlas is redundant, and it drops the `_inc`/`_dec` pairing that gives the
reporter screens their signal.

### 13.8 The teammate's repo: what it settles and what it does not

`github.com/jhouser220/capstone`, branch `Odysseus`,
`experiments/odysseus/assayloop_comparison`. Their records are vendored under
`output/runs/external/teammate/`; `scripts/external/teammate_compare.py` puts
the two experiments side by side.

**Their experiment is a different comparison from ours.** Arms are `Native` (the
original AssayLoop loop: an LLM proposes, pre-batch AssayFormer fills short
batches) against `HSA` = **Historical Subspace Adaptation**, their method —
32 historical response directions appended to the source's 10 features, task
coefficients fitted on measured labels, ranked by availability times integrated
hit probability. Actors are `google.gemma-4-31b` and `qwen.qwen3-32b`, 2 draws
each, so 4 campaigns per screen per arm. **Both arms are LLM-in-the-loop**; ours
has no LLM in either arm.

**Their budget semantics differ.** The pool is a fixed **20,893 genes**, not the
screen's library, and a pick with no measured outcome is charged to the budget:
of 1,000 accepted picks only 773–949 (mean 893) land on a measured gene. Ours
sets `universe_genes` to the library, so all 1,000 are measured. Their own docs
warn the two AUC labels are different metrics: "Requested-slot AUC and source
measured-axis nAUC are different metrics. Never use one label for both."

```
screen                                       H |  Native     HSA      d   meas |    AF   +DM     d
E224 SCHUKKEN PAIRED WGL                   188 |   47.75   50.00  +2.25    909 |    45    42    -3
E225 ABT PARCB DEPENDENCY                 1398 |  322.00  365.50 +43.50    773 |   431   408   -23
E226 WILSON MERCS HIGH                     394 |   97.25  107.75 +10.50    932 |   101   104    +3
E232 NGUYEN K562 ARAC KO SENSITIZATION     240 |   13.50   15.00  +1.50    911 |    20    21    +1
E236 SHAN HEPG2 LOVASTATIN                 368 |   82.00   97.25 +15.25    912 |    78    77    -1
E236 SHAN HELA LOVASTATIN                  359 |  132.00  136.75  +4.75    915 |   133   121   -12
E236 SHAN HELA LDLR KO                     227 |   15.25   26.50 +11.25    903 |    14    25   +11
E239 ORGANOID COLO 031                     984 |  498.75  521.75 +23.00    889 |   564   563    -1
E239 ORGANOID OESO 005                    1683 |  792.00  866.50 +74.50    905 |   942   922   -20
E239 ORGANOID OVAR 018                    1473 |  680.50  747.75 +67.25    892 |   794   797    +3
E245 SAMELSON T22 LOWERING                 729 |  113.00  155.00 +42.00    878 |   124   126    +2
MEAN                                           |  254.00  280.89 +26.89    893 | 295.1 291.5  -3.6
```

HSA beats Native **11/11**; AF + DepMap beats AF alone 5/11 on the same screens.
Both are true, and the reason they can both be true is that the baselines differ:
**AF alone (no LLM) already finds more hits than their HSA arm on 5 of the 11
screens** and more than their Native arm on 8 — most starkly Abt (431 vs 322
Native, 365.5 HSA) and the organoids (942 vs 792 / 866.5 on OESO-005). Their
intervention is measured against LLM-driven AssayLoop; ours against AssayFormer's
own ranking, which is the stronger control. Their per-1000-measured column
narrows but does not close this.

**What the repo does NOT settle: E224's hit rule.** Their own
`EXTERNAL_REPRODUCTION_ASSETS.json` lists
`Experiments/external_aneuploid_e224_2026-10-26/{build_target.py,TARGET.json,FREEZE.json}`
as assets recorded by hash but deliberately excluded from public Git ("No raw
author tables, weights or provider receipts are part of this public payload").
Searches that came back negative: the Schukken repo has no RPE1-WGD screen at all
(`SupData/` is SD1–SD9; `Data/` and `Dependency/` hold no RPE1 clone), so their
`DATASETS.md` phrase "source-derived RPE1-WGD/ESL endpoint" does not point at a
file in that repo; SD7's `WGL.Top10` sheet is the same 19,112-gene paired ranking
as the `WGL` sheet with no 194 flag; SD7's `Top Hit List` sheet unions to 348
genes, not 194; and no round threshold on the paired columns gives 194 (§13.7).
**Ask for `build_target.py` or `TARGET.json` for E224.**

Two of their stated limitations are worth carrying into any write-up: Wilson is
flagged "Source table versus paper hit-count mismatch", which is independently
the MAGeCK RRA FDR tie problem recorded in §13.4, and Samelson is flagged
"CRISPRbrain T22 release used; final-paper annotation/version parity unresolved",
which is why their 728 and our 729 differ by one gene.

### 13.9 Candidate-pool conventions, and the matched comparison

Three conventions for the same 1000-pick budget. They are not interchangeable and
EF is not comparable across them.

| convention | candidates | off-library pick | measured picks of 1000 |
|---|---|---|---|
| **library pool** (ours, §13.1) | the screen's own library | impossible | 1000 |
| **fixed pool, charged** (the teammate's) | one cross-screen pool for every screen | revealed as a miss, costs a slot | 770–970 |
| **fixed pool, forgiven** (AssayBench native) | same pool | advances neither axis | ~437 on the folds |

`run_external.py --universe <genes.json>` reproduces the charged convention
exactly, with no new code in the loop: `GeneBatchTask.reveal` already returns
`hit=False, in_library=False` for a gene the screen never measured
(`src/assayloop/tasks/gene_batch.py:103`), and that observation still occupies a
budget slot. The driver now also records `n_unknown` and `n_measured` per run.

**EF rises when the pool widens, even as raw hits fall.** `hits_auc` sets
`library_size = max(len(screen.genes), pool_size)`
(`src/assayloop/metrics/hits_auc.py:245`), so a wider pool means a worse random
baseline. Abt goes 431 → 349 hits while its EF goes 4.458 → 5.296. **Compare raw
hits in 1000 picks, or hits per 1000 measured picks — never EF across
conventions.** This is the same trap as §13.3's "% of ceiling is unbounded on the
benchmark".

Their 20,893-gene pool is not in the public package (`UNIVERSE_validation.json`
is in their asset list, hash only). Stand-in used here:
`gene_universe(load_screens("validation"), min_screen_freq=2)` = **21,213 genes** —
the same kind of object, independent of all seven external papers, covering
98.7–100% of each library and leaving 98.3–100% of hits reachable. Its
screen-specific unmeasured fraction tracks theirs: Abt has the largest off-library
share in the stand-in (6,759 of 21,213) and the largest unknown count in their
records (225–229 of 1000).

```
                                             library pool  |  fixed pool (charged) |  theirs: fixed pool, LLM
screen                                     H |   AF  +DM    d |   AF  +DM    d  meas |    Nat    HSA      d  meas
E224 SCHUKKEN PAIRED WGL                 188 |   45   42   -3 |   38   39   +1   958 |  47.75  50.00  +2.25   909
E225 ABT PARCB DEPENDENCY               1398 |  431  408  -23 |  349  325  -24   762 | 322.00 365.50 +43.50   773
E226 WILSON MERCS HIGH                   394 |  101  104   +3 |  104  114  +10   885 |  97.25 107.75 +10.50   932
E232 NGUYEN K562 ARAC KO SENSITIZATION   240 |   20   21   +1 |   13   17   +4   922 |  13.50  15.00  +1.50   911
E236 SHAN HEPG2 LOVASTATIN               368 |   78   77   -1 |   85   91   +6   908 |  82.00  97.25 +15.25   912
E236 SHAN HELA LOVASTATIN                359 |  133  121  -12 |  160  148  -12   938 | 132.00 136.75  +4.75   915
E236 SHAN HELA LDLR KO                   227 |   14   25  +11 |   23   27   +4   913 |  15.25  26.50 +11.25   903
E239 ORGANOID COLO 031                   984 |  564  563   -1 |  526  525   -1   882 | 498.75 521.75 +23.00   889
E239 ORGANOID OESO 005                  1683 |  942  922  -20 |  852  836  -16   891 | 792.00 866.50 +74.50   905
E239 ORGANOID OVAR 018                  1473 |  794  797   +3 |  733  734   +1   887 | 680.50 747.75 +67.25   892
E245 SAMELSON T22 LOWERING               729 |  124  126   +2 |  118  137  +19   853 | 113.00 155.00 +42.00   878
MEAN                                         |  295  291 -3.6 |  273  272 -0.7   891 | 254.00 280.89 +26.89   893
```

**What the convention alone is worth:** AF alone, same model, same seed, goes
295.1 → 272.8 hits per screen (−7.5%) purely because 114 of its 1000 picks now
land on genes the screen never measured.

**The matched comparison, and a correction to §13.8.** §13.8 said AF alone beats
their HSA arm on 5 of 11 screens and their Native arm on 8. That was across
mismatched conventions. On the matched footing it is **HSA on 2/11 and Native on
9/11**: our AF alone (272.8) sits *between* their Native (254.0) and their HSA
(280.9). So HSA does real work against a no-LLM AssayFormer control, by about
+8 hits per screen on these 11. The §13.8 point that survives is narrower and
still worth making: their *stated* +26.9 gain is measured against a baseline that
a plain AssayFormer run beats by ~19 hits per screen, so roughly two thirds of
that headline is the LLM-proposal baseline being weak rather than HSA being
strong.

AF + DepMap is ~0 on this subset under either convention: −3.6 hits (5/11 up) on
the library pool, −0.7 hits (7/11 up) on the fixed pool. The convention does not
rescue or explain the flat result.

**To close the remaining gaps**, in order of value: (1) their
`UNIVERSE_validation.json` so the pool is identical rather than a stand-in;
(2) E224's `build_target.py`/`TARGET.json` (§13.7); (3) an HSA arm on top of
AssayFormer-only proposals, or a DepMap arm inside their LLM loop — until one of
those exists, the two interventions have never been measured against the same
control.

Scripts: `scripts/external/pool_convention.py`, and
`run_external.py --universe/--screen-names`. Results:
`output/runs/external/team11_fixedpool_results.json`.

---

## 14. Repository reorganization, 2026-10-30

Jeff asked for three things: outputs sorted into folders that say what method they tested,
anything not relevant to **AF + DepMap** or **top-up** moved out to
`../../AssayLoop_old/`, and the two methods in separate `.py` files. He chose, when asked,
to **keep every output and only sort it** — nothing under `output/` was moved out of the
project or deleted.

### 14.1 `output/runs/` — 20,908 run directories into 28 method folders

Was flat: 179 sweep ids, every `<sweep>-<NN>-<screen>` directory sitting directly in
`output/runs/`. Now each sweep lives in a folder named for the components it fused —
`af_depmap/`, `af_depmap_screenknn/`, `af_alone_baseline/`,
`af_depmap_screenknn_bioplex_ppi/`, `af_description_channel/`, and so on. Full index with
per-sweep run counts: **`output/runs/README.md`**.

Classification read `model`, `acq`, `screen_set` and the two override blocks from one
representative `config.json` per sweep — never the `lm` block, which is where the API keys
still sit (§7). Arms that ran on `achilles.npz` get a `__full_atlas_superseded` folder of
their own, so pooling the two atlases now takes a deliberate act.

**Four families are not classifiable from config and were named from §9 instead.** Their
intervention arrived by environment variable or by a seed list injected straight into the
campaign, so they record *empty* `acq_overrides` and would otherwise have been filed as
plain AssayFormer — which is what they are mechanically and not what they test:
`sweep-desc-*` (`ASSAYLOOP_DESC_EXTRA`), `sweep-orc-*`, `sweep-inj-*`, and the LLM
handoffs `sweep-hybrid-*` / `sweep-dmne-hyb-*` / `sweep-dmcd-hyb-*` / `sweep-glm5tu-pt-af`.
The classifier and its `NAME_OVERRIDE` table are kept in **`scripts/organize_runs.py`**;
re-run it (`--apply`) after a new sweep and it files that sweep too. It is idempotent.
**A new sweep in a family whose intervention is not in `config.json` needs a
`NAME_OVERRIDE` entry first, or it lands in `af_alone_baseline/`.**

### 14.2 The move breaks every flat path, so paths now resolve by search

~25 scorers build `<dir>/<sweep>-<NN>-<screen>/result.json` from a single `--dir`, and a
paired contrast now routinely spans two folders — a family's `weight=0` baseline is in
`af_alone_baseline/` while its treatment arm is in `af_depmap/`. New module
**`scripts/diagnostics/run_paths.py`** indexes every sweep id below a root and resolves it
to whichever folder holds it; `run_dir()` tries the literal path first, so an explicit
`--dir output/runs/af_depmap` still costs no search. 26 call sites in 24 files were
patched to go through it, plus six flat `glob()` sites. **`--dir output/runs` keeps
working, and is now the right thing to pass.**

`compare_all_results.ipynb` resolves the same way: a new `arm_base(root, sid)` helper in
the setup cell, used at the eight sites that had `RUNS if root == "." else RUNS / root`.
**No registry entry changed** — an explicit `root` like `baseline/gemini` still wins and
the default `"."` is searched. The notebook was edited but NOT executed; its 56 cells and
45 stored outputs are intact.

**`--resume` still does not search.** A sweep you intend to resume has to be moved back up
into `output/runs/` first. That was already true of `baseline/`.

Verified after the move: `codep_gated.py` on `validation` reproduces **+0.285
[+0.135, +0.448]** across the two folders, and `scripts/external/final_table.py` reproduces
the external table unchanged (folds +0.307, pooled external −0.075).

### 14.3 `output/analysis/` — logs by method, inputs left addressable

Loose driver logs grouped into `logs_<method>/` (`logs_bioplex_ppi`, `logs_cohit`,
`logs_description_channel`, …), absorbing the pre-existing `*_logs/` directories.
**The `.npz` files stay at the top level on purpose**: each is a feature matrix an
acquisition reads at run time by a hardcoded path (`train_hits=output/analysis/train_hits`,
`ppi_npz=output/analysis/bioplex.npz`), so filing them would break the launchers that name
them. `output/analysis/README.md` records which method each feeds and whether it is live.
`train_hits.npz` is live — **the standard method needs it**.

### 14.4 What moved to `../../AssayLoop_old/`

Code only, with an index in its own `README.md`. Three acquisitions
(`llm_category_acq.py`, `mixed_batch_acq.py`, `decline_handoff_acq.py`), their two root
drivers (`handoff_decline.py`, `handoff_mixed.py`), the 19 diagnostics for the refuted
feature families (BioPlex, co-hit, PRISM, condition-match, Perturb-seq, JUMP, GenePT,
shrinkage, LLM-as-feature, context and description sweeps, Sack overlap), `Sack_dataset/`,
and the four tests covering them.

Dropping `llm_category_acq.py` also required removing its entries from
`acquisitions/__init__.py`, `cli.py` (three acquisition sets plus the `--acq` help text) and
the `llm_category` branch of `experiment/runner.py`. **The `llm_single_topup` entries in
those same places were kept.**

**Kept deliberately, against a first reading of "not relevant":**
`llm/prompt_override.py` (top-up imports it) and `screen_direction.py` (AssayFormer's text
encoder imports it) — removing either breaks a live method. `handoff_fullgenome.py`, which
is how a top-up warm start is handed to AssayFormer. `run_oracle_seeds.py` and
`run_inject_2x2.py`, which are the nominator bound that closed the LLM question for
AF + DepMap. `run_geometry_sweep.sh` (DepMap atlas geometry). The whole upstream package —
Jeff chose "our dead-end additions only", so BPMF, MAML, `animation/` and `deployment/` all
stay and the package still imports.

### 14.5 Separate `.py` files — already true

**AF + DepMap** is `src/assayloop/acquisitions/codep_seed_acq.py` and **top-up** is
`src/assayloop/acquisitions/llm_single_topup_acq.py`. No code coupling: the top-up file has
zero DepMap references, and the only LLM mention in `codep_seed_acq.py` is a comment on
line 275. Nothing needed splitting.

### 14.6 State of the test suite

**270 passed, 2 failed.** Both failures are pre-existing and unrelated to the move —
neither test mentions anything that moved:
`test_cli_full_genome.py::test_open_vocab_llm_defaults_to_shared_f2_universe` asserts the
`screen_set` default is `public` where both HEAD and the working tree say `paper_test`, and
`test_gemini_client.py::test_yaml_loader_maps_new_fields_and_nulls_inherited_extensions`
expects `reasoning_effort == "low"` where the loader now gives `"minimal"`. They come from
the uncommitted work on `cli.py` and `llm/client.py`, not from this reorganization.

### 14.7 Still outstanding

**The Bedrock and Google API keys are still in plaintext in the run `config.json` files
under `output/runs/`, and still need rotating and scrubbing.** The reorganization moved
those files into method folders; it did not touch their contents, so the exposure is
unchanged. `check_secrets(fix=True)` in the notebook now has to walk the method folders.

---

## 15. The HSA reproduction bundle, pulled and verified 2026-10-31

A public GitHub Release on `jhouser220/capstone` carries the frozen inputs for the HSA
comparison. **It was not known to this project until now**, and it supersedes three of the
four "blockers" recorded in §13.9 and §14.

```
tag       odysseus-assayloop-reproduction-2026-10-22
file      assayloop-presentation-reproduction-2026-10-22.tar.gz
bytes     186,856,583          (matches REPRODUCTION_ASSETS.json)
sha256    19d409b2ab048ff1b8f402183d33bd0073efb0995a5ffa8049b6e08dc0e93903   VERIFIED
```

Integrity: 622 archive entries, no absolute paths and no `..`. The bundle's own
`SHA256.json` verifies **403 of 403** files. `bundle/checkpoint/{model_last.pt,vocab.json,
config.json}` are **byte-identical** to our cached `Genentech/assayformer` revision
`4927703e5903984ae0254a9fe8e93ec0b29c27e4`.

### 15.1 The candidate universes are gene-for-gene identical to ours

Settles the open question of §13.9 — their pool needed no reconstruction and never
diverged:

| bundle file | genes | equals | identical? |
|---|---|---|---|
| `UNIVERSE_validation.json` | 20,893 | `gene_universe(paper_validation, min_screen_freq=2)` | **yes** |
| `UNIVERSE_test.json` | 21,147 | `gene_universe(paper_test, min_screen_freq=2)` | **yes** |
| `UNIVERSE_train.json` | 23,717 | `gene_universe(public_train, min_screen_freq=2)` | **yes** |

Not merely the same size — the sorted lists compare equal. So HSA and every arm of ours
that ran `--full-genome --min-screen-freq 2` share one pool construction exactly.

### 15.2 THE 1,568-DIMENSION DESCRIPTION VECTOR DOES NOT EXIST — a correction

§14 and earlier notes record the per-screen `DESCRIPTION*.npy` as 6,272 bytes = **1,568
float32**, and treated the 32 unexplained dimensions as a blocker. That was an arithmetic
error: the `.npy` v1 header is 128 bytes, so

    1536 float32 x 4 bytes = 6144  +  128 header  =  6272

Every description vector in the bundle is `shape (1536,) float32`, L2-normalised to 1.0 —
**AssayFormer's native `text_dim` under `text-embedding-3-small`**, which is exactly what
`output/rankers/text_emb_cache` already holds (1,536-D, same model per the checkpoint's
`config.json`). Nothing about the description channel is unknown, and it is reproducible
on our side for a fraction of a cent.

The 32 SVD directions are real and central to HSA, but they live in `bank.residual`
(n_genes x 32), not in the description.

### 15.3 The banks are complete and pass their own validator

18 banks, one per query screen, loaded through their `assets.HistoricalBank.load`:
**18/18 valid.** Each carries `genes` (23,742), `residual` (23742, 32) float32 with row
norms exactly 1.0 and 25 all-zero rows for uncovered genes, `mean` (32,), `availability`
in [0.0074, 0.8167], and `scalar_c` ~69.3. Leave-one-publication-out is enforced in
`__post_init__`: each bank names 133-134 donor publications with its own query excluded.

`build_bank()` in their `assets.py` is public and shows the construction. The step worth
knowing: `Q,_ = qr(W)` then `centered -= (centered @ Q) @ Q.T`, where `W` is
AssayFormer's `gene_emb.weight` (n, 10). **The 32 historical directions are
orthogonalised against AssayFormer's gene space by construction**, so HSA adds only what
AssayFormer cannot express. That is the cleanest statement of the method, and it explains
why their own "MAP in the same 42D" control (282.0) captures nearly all of full HSA's
gain (284.5).

Two banks key on a DOI / bioRxiv URL rather than a PubMed id (`U_1718_dec`,
`U_1737_merged`) — harmless unless something joins on publication id.

### 15.4 What the bundle does NOT contain

**It is the 18-screen benchmark panel, not the external one.** 342 recorded campaigns
under `bundle/recorded/` (`<screen>__<arm>__<actor>__<draw>.json.gz`, arms `incumbent` and
`native_filled`, actors Gemma 4 31B and Qwen3 32B, draws 0/1). There is no `TARGET.json`,
no `DESCRIPTION*.npy` and no bank for any external experiment — nothing matching E224,
E225, E226, E232, E236, E239 or E245.

So: the ORIGINAL comparison is now fully reproducible here; the **external** HSA figures
(the 280.9 Language+HSA / 284.5 HSA-without-language means over our 11 shared screens)
are **still not reproducible from public artifacts** and still need the 44 assets listed
in `EXTERNAL_REPRODUCTION_ASSETS.json`.

### 15.5 Three of the 18 screens are in AssayFormer's training fold

`1125`, `1628` and `174` are in `public_train`. The bank excludes the query publication,
so the HSA side is leakage-controlled, but the frozen AssayFormer was trained on those
screens. The rest: 6 in the `validation` fold (2 of them curated `paper_validation`),
9 in `test` (2 curated `paper_test`). Their own docs say so plainly — *"All these cohorts
have informed development"*, *"independent confirmation needs a reserved cohort and a
frozen method."* Do not quote any of the 18-screen numbers as held out.

### 15.6 What is now unblocked, all free

With the bundle, these need no provider calls and no money — `run_campaign` exempts
`assayformer_only` and `historical_only` from needing an actor:

1. **HSA with no LLM** (`historical_only`) on the 18 screens.
2. **AssayFormer alone** (`assayformer_only`) on the same, as the matched baseline.
3. **Exact replay** of all 342 campaigns against their receipts (`reproduce.py`), which
   checks 864 language requests and makes zero provider calls.

The 42-D design is fully specified: 10 AssayFormer gene features, then
`0.25 x bank.residual` (32 historical directions); prior `np.r_[z, bank.mean]`; precision
`np.diag(np.r_[ones(10), full(32, 1/bank.scalar_c)])`; damped-Newton MAP; 128 antithetic
draws; score = `availability x P(hit)`.

**Where it is.** Downloaded and extracted under this session's scratchpad
(`.../scratchpad/bundle_dl/bundle.tar.gz` and `.../scratchpad/bundle_extract/`), which is
EPHEMERAL. Move it somewhere durable before relying on it; it is 187 MB compressed and
gitignored territory either way.

---

## 16. HSA and AssayFormer-alone RE-RUN from the frozen assets, 2026-10-31

Not a recount — the model and bank were loaded and the campaigns re-executed through
their own `reproduce.py` in replay mode, which re-runs `run_campaign` and compares every
round's picks AND observations by sha256 digest against the frozen trajectory, then
asserts `unique_hits` equal and `|ef difference| < 1e-10`. Zero provider calls; the two
arms are the ones `run_campaign` exempts from needing an actor.

**Arm names differ between their layers**, which cost some time: the recorded files and
`manifest.json` use `source` and `scalar`; `adapter.py` calls them `assayformer_only` and
`historical_only`; `cohort_means.csv` calls them `AssayFormer_only` and `cold_incumbent`.
`reproduce.py` maps them at line 90. All four names are the same two arms.

### 16.1 Result: 54/54 on every metric, 53/54 byte-exact

| check | result |
|---|---|
| campaigns re-run (18 `source` x1 draw + 18 `scalar` x2 draws) | **54** |
| max abs difference in `unique_hits` vs recorded | **0** |
| max abs difference in `measured` vs recorded | **0** |
| max abs difference in `ef` vs recorded | **0.000e+00** (bit-identical) |
| strict digest replay passed | **53 / 54** |
| published cohort means reproduced (4 cohorts x 2 arms) | **8 / 8 exactly** |

The single strict failure is the tie-ordering artefact their own `REPRODUCING.md` predicts
("numerical libraries on other platforms can affect exact tie ordering"): screen `2456`,
arm `source`, draw 0, round 5 is a **pure transposition of two adjacent tied genes**,
`PABIR1` <-> `KMT2B` at positions 19-20. Identical gene set, identical every metric. Our
numpy 2.4.6 / scipy 1.17.1 match their pinned versions exactly; torch (2.11.0+cu128 vs
2.12.0) and platform (linux x86_64 vs macOS arm64) do not.

Runtime: 54 campaigns in **26 s** on 4 workers, CPU only.

### 16.2 The paired contrast, re-run

Per screen, averaging draws within a screen, hits in 1000 selections:

| fold | n | AF alone | HSA (no LLM) | delta |
|---|---|---|---|---|
| `test` | 9 | 135.7 | 141.0 | **+5.28** |
| `train` | 3 | 43.3 | 56.7 | **+13.33** |
| `validation` | 6 | 155.8 | 161.1 | **+5.25** |
| **all 18** | 18 | **127.6** | **134.2** | **+6.61** [+1.56, +12.39], up on 11/18 |

So HSA's numerical core is worth about **+6.6 hits per screen** over AssayFormer alone on
this panel, interval excluding zero on a screen bootstrap. That is a smaller gain than the
+11.7 inferred in §13.9 from their external table, and it is measured rather than inferred.

### 16.3 TWO FINDINGS THAT QUALIFY THE GAIN

**(a) The largest gain is on the screens AssayFormer was TRAINED on.** `train` is exactly
`1125`, `1628`, `174` (§15.5) and shows +13.33, two and a half times either held-out fold.
The bank excludes the query publication, so HSA is leakage-controlled; AssayFormer is not.
A frozen model that has memorised a screen leaves less headroom in principle, yet HSA
gains MORE there -- the opposite of what a pure ranking improvement predicts. Do not quote
the 18-screen mean without this split.

**(b) Much of the gain is AVAILABILITY, not ranking.** HSA's score is
`availability q x P(hit)`, and `q` steers it toward genes the screen actually measured.
Mean measured slots per campaign: AF alone **672.9**, HSA **736.9** -- HSA converts ~64
more of its 1000 selections into measured outcomes. On the training screens the effect is
extreme (`1125` 575 -> 909, `174` 556 -> 835), which is also where the hit gain is largest.
**Under a screen-library pool, where every pick is measured by construction, this half of
the advantage cannot exist.** Any fair AF + DepMap vs HSA comparison has to either keep the
fixed pool for both arms or report measured-pick-normalised hits.

Related: six `validation` screens barely overlap the 20,893-gene `UNIVERSE_validation`
pool -- `1588` measures 18 of 1000 selections, `U_1508_merged` 56, `1566` 77,
`U_1718_dec` 137 -- so their hit counts are near zero for both arms and contribute almost
no information. That is the pool-mismatch cost of §15.1 made concrete.

### 16.4 How to re-run it

```sh
W=<workdir>                       # holds experiments/odysseus/assayloop_comparison/
export PYTHONPATH="$W:$W/_reproduction/assayloop-src:$W/_reproduction/vendor"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
python -m experiments.odysseus.assayloop_comparison.reproduce \
  --bundle _reproduction/bundle --output out_numeric --workers 4 --arms source scalar
```
`_reproduction/{bundle,assayloop-src,vendor}` come from the §15 tarball; the module files
are fetched from `jhouser220/capstone@Odysseus`. Do NOT run these with `python -I` -- it
discards `PYTHONPATH` and the package will not import. `--output` must be a fresh
directory. Adding `--arms incumbent native native_filled ordinary_full` replays the
language arms from their receipts, still with zero provider calls.

**This covers the 18-screen benchmark panel only.** The external 11-screen HSA figures
(280.9 / 284.5) remain unreproducible from public assets -- see §15.4.

---

## 17. AF + DepMap beside HSA on the bundle's 18 screens, 2026-10-31

The first comparison of the standard method against HSA on ONE pool, ONE panel and ONE
metric. Script: `scripts/external/run_bundle18.py`. Results:
`output/runs/external/bundle18_results.{json,csv}`.

### 17.1 What was held identical

Only the acquisition differs. Same frozen `ScreenRecord`s (`bundle/screens/<id>.json`),
same candidate pool (`bundle/UNIVERSE_{train,validation,test}.json`, chosen per screen by
its manifest `fold`), same 1,536-D description vector injected through
`desc_emb_by_name` (so no text embedding is recomputed and **no API call is made**), same
checkpoint (`bundle/checkpoint`, byte-identical to our cached 4927703e), same budget
(10 x 100, seed 0), same metric (unique MEASURED hits in 1000 selections, unknowns
charged). DepMap atlas `achilles_noeval.npz`, per the standing rule.

**These are NOT the pools our own sweeps used.** `sweep-dmne-cd-{validation,test}` ran on
the 218-/334-screen FOLD universes (21,213 / 22,656); the bundle uses the 20-screen
curated universes (20,893 / 21,147) plus `UNIVERSE_train` (23,717). Hence the re-run.

**BRIDGE CONTROL, and it is exact.** Our `AF alone` (`codep_seed` at weight 0) reproduces
the bundle's own `source` arm on **18/18 screens**, identical on both hits and measured
slots. `codep_seed` at weight 0 also equals plain `greedy` on 18/18. So the DepMap column
is attributable to the DepMap term and nothing else.

### 17.2 The table

Hits in 1000 selections, per screen, averaging HSA's two draws:

| fold | n | AF alone | AF + DepMap | HSA (no LLM) | +DM − AF | HSA − AF |
|---|---|---|---|---|---|---|
| `test` | 9 | 136.8 | 142.0 | 142.1 | **+5.22** | **+5.28** |
| `train` | 3 | 43.3 | 56.7 | 56.7 | **+13.33** | **+13.33** |
| `validation` | 6 | 155.8 | 154.2 | 161.1 | **−1.67** | **+5.25** |
| **all 18** | 18 | **127.6** | **131.8** | **134.2** | **+4.28** | **+6.61** |

Screen bootstrap, all 18: AF + DepMap **+4.28 [−0.89, +10.06]**, up on 9/18 — **crosses
zero**. HSA **+6.61 [+1.56, +12.39]**, up on 11/18 — excludes zero.

### 17.3 The direct contrast is UNRESOLVED, and the mean is misleading

**HSA − AF + DepMap = +2.33 [−2.78, +8.61]** over 18 screens. But the sign test runs the
other way: **AF + DepMap is better on 10 screens, HSA on 7, 1 tied.** The positive mean
comes almost entirely from two high-hit screens — `U_1793_dec` (+41.5) and
`U_1737_merged` (+24.0). Dropping the four screens where both arms find <=10 hits leaves
+2.80 [−3.27, +10.47], still crossing zero. On this panel the two methods are not
distinguishable, and whichever summary you quote should state both the mean and the 10-7
split.

This also deflates §13.9's inference: from their external table HSA looked worth ~+11.7
over AssayFormer alone. Measured on one pool it is +6.61, and AF + DepMap gets +4.28 of it.

### 17.4 The availability gap persists and now has a control

Mean measured slots per campaign: **AF alone 672.9, AF + DepMap 684.7, HSA 736.9.** HSA
converts ~64 more of its 1000 selections into measured outcomes; AF + DepMap only ~12.
So the §16.3(b) mechanism is specific to HSA's `availability q x P(hit)` term and is NOT
something DepMap fusion reproduces. Under a screen-library pool that advantage cannot
exist at all, so the fixed pool is load-bearing for part of HSA's margin.

### 17.5 Caveats

* `train` (3 screens) shows the largest gain for BOTH arms (+13.33 each, by coincidence of
  the means -- the per-screen values differ). Those are the screens AssayFormer was
  trained on (§15.5), and the reduced atlas drops cell lines chosen from the four
  EVALUATION sets, not from `public_train`, so those three screens' own lines may still
  sit in the atlas. The DepMap arm is therefore not leakage-controlled there in the way
  it is on the held-out folds.
* Four `validation` screens barely intersect their pool (`1588` measures 18 of 1000
  selections, `U_1508_merged` 56, `1566` 77, `U_1718_dec` 137) and both arms score near
  zero. They contribute noise, not signal.
* n=18 screens against the 198/314 of our fold confirmations. This does not revise the
  headline `+0.307 EF [+0.222, +0.398]` on the benchmark folds -- different pool,
  different screens, different metric. Read it as a head-to-head against HSA, not as a
  re-measurement of AF + DepMap.

---

## 18. OFFLINE GATE: DepMap and HSA are near-ORTHOGONAL, so build the combined arm

The §17 scope put one gate before any combined arm, following the §4.16 rule that an
offline-negative candidate does not get a campaign: are the DepMap co-essentiality kernel
and HSA's 32 historical directions measuring the same thing? Script:
`scripts/diagnostics/screen_hsa_depmap.py` (see 18.3 for where it lives).

### 18.1 Method

Five bundle screens, spanning cohorts and hit densities. For each, the RECORDED HSA
trajectory (`<sid>__scalar__none__0`) is replayed so both scorers see identical state, and
at every one of the 10 rounds three score vectors are computed over the SAME candidate set:

* `af`  AssayFormer's own score, `model.predict`
* `cd`  our DepMap kernel, read out of a live `CodepSeedAcquisition`'s `_P` / `_row` /
  `_seed_query` -- the same kernel the deployed arm uses, not a reimplementation
* `hsa` HSA's `q*prob`, from THEIR `Recipient(arm='historical_only').predict` -- their code,
  so the comparison cannot be wrong about what HSA scores

### 18.2 Result: the two signals are independent

Means over 10 rounds; `rho` is Spearman over candidates, overlaps are of the top 100:

| screen | hits | rho(cd,hsa) | rho(af,hsa) | rho(af,cd) | cd^hsa | **fused^hsa** | af^hsa |
|---|---|---|---|---|---|---|---|
| `1125` | 127 | **−0.006** | 0.526 | 0.002 | 5.0 | 8.7 | 6.7 |
| `1985` | ~190 | **+0.020** | 0.810 | −0.016 | 7.0 | 18.0 | 16.8 |
| `2346` | 245 | **−0.022** | 0.747 | −0.016 | 6.1 | 23.0 | 23.1 |
| `2462` | 253 | **−0.021** | 0.751 | −0.015 | 9.2 | 24.7 | 24.4 |
| `U_1793_dec` | 845 | **−0.008** | 0.810 | −0.015 | 12.0 | 32.2 | 32.0 |

Three readings, all consistent across every screen:

1. **`rho(cd, hsa)` is zero** (−0.022 to +0.020). The DepMap kernel and HSA's posterior
   score are statistically independent across candidates. **They are not redundant.**
2. **`rho(af, hsa)` is high** (0.53-0.81), and on `1125` it decays 0.861 -> 0.185 across
   the rounds. Expected: AssayFormer enters HSA through the prior mean `z`, the first ten
   design columns `W[ids]` and the bias `b[ids]`, and HSA only drifts from it as measured
   labels accumulate.
3. **`fused^hsa` ~ `af^hsa` on every screen** (23.0 vs 23.1, 24.7 vs 24.4, 32.2 vs 32.0).
   Adding the DepMap term to AssayFormer barely changes how much the selection overlaps
   HSA's -- i.e. DepMap moves the batch in a direction HSA does not travel at all.

The top-100 overlap between what AF + DepMap selects and what HSA selects is only
**18-32 of 100**, so the two arms already spend 70-80% of each batch on different genes.
That is the regime where an ensemble can add, and it is consistent with §17's
10-7 screen split: the two methods win on different screens because they are reading
different structure.

**GATE PASSES -- the combined arm is worth building.** Note the limit of the evidence:
orthogonal SCORES do not guarantee additive HITS. §4.48-4.51 found near-substitute
patterns where two signals each helped alone and not together. This licenses the
experiment; it does not predict the outcome.

### 18.3 One design consequence

Because `hsa` correlates 0.5-0.8 with `af` while `cd` correlates with neither, an RRF
ensemble `1/(K+r_af) + 0.5*rrf(cd) + w*rrf(hsa)` **double-counts AssayFormer** through
both the explicit term and HSA's interior. So the weight grid of §17 should be run
alongside a variant that DROPS the explicit `af` term and fuses `rrf(hsa) + w_dm*rrf(cd)`
only. The nested-at-zero property is lost in that variant, so it needs its own baseline.

Also: `cd` is undefined in round 1 on every screen (`min_seeds=2`, no labels yet), so
DepMap contributes nothing to the first batch while HSA's prior already departs from
AssayFormer. Any round-1 difference between the arms is HSA's alone.

---

## 19. AF + DepMap + HSA: the gains are ADDITIVE, 2026-10-31

The §18 gate said build it, so it was built. Code: a seventh fusion term in
`codep_seed_acq.py` plus an `af_weight` knob; driver `scripts/external/run_bundle18.py`;
tests `tests/test_codep_hsa_term.py`. Results:
`output/runs/external/bundle18_sixarm_results.{json,csv}`.

### 19.1 Implementation

```python
fused  = self.af_weight * (1.0 / (RRF_K + r_af))
fused += self.weight     * _rrf_term(cd)     # DepMap co-essentiality
...
fused += self.hsa_weight * _rrf_term(hsa)    # NEW
```

The HSA math is **pluggable, not vendored**: `hsa_scorer(observations, candidates,
task_context)` is a caller-supplied callable, and the driver wraps their own
`Recipient(arm='historical_only')`. Their bank, source-state encoder and numerics never
become a dependency of this package, and the correctness gate is exact by construction
because their code is what runs. The scorer is called exactly once per round and only
when `hsa_weight != 0`, because `posterior_scores` counts rounds and raises outside 1..10.

`af_weight` defaults to 1.0 = the historical behaviour, bit for bit.

### 19.2 Gates, all 18 screens

| gate | result |
|---|---|
| 1. `AF alone` and `AF + DepMap` unchanged vs §17 | **PASS**, 18 screens x 2 arms |
| 2. weight-0 `codep_seed` == plain `greedy` | **PASS** (§17, 18/18) |
| 3. `HSA only (ours)` == recorded `historical_only` draw 0 | **PASS**, 18/18 on hits AND measured |
| existing `test_codep_seed.py` | 104 passed |
| new `test_codep_hsa_term.py` | 6 passed |

Gate 3 is the load-bearing one: driving their scorer through OUR ensemble at
`af_weight=0, weight=0, hsa_weight=1` reproduces their arm exactly, so the combined
columns differ from HSA only by the DepMap term.

### 19.3 The table, hits in 1000 selections

| arm | af | dm | hsa | test (9) | train (3) | val (6) | **all 18** | measured |
|---|---|---|---|---|---|---|---|---|
| AF alone | 1 | 0 | 0 | 136.8 | 43.3 | 155.8 | **127.6** | 672.9 |
| AF + DepMap | 1 | .5 | 0 | 142.0 | 56.7 | 154.2 | **131.8** | 684.7 |
| HSA only (ours) | 0 | 0 | 1 | 140.9 | 57.3 | 161.5 | **133.8** | 737.6 |
| **DepMap + HSA (no AF term)** | **0** | **.5** | **1** | **146.3** | **70.7** | **158.7** | **137.8** | 762.6 |
| HSA + DepMap eq | 0 | 1 | 1 | 137.0 | 78.7 | 144.7 | **129.8** | 775.9 |
| AF + DepMap + HSA | 1 | .5 | .5 | 147.9 | 46.7 | 159.5 | **134.9** | 682.1 |

Paired over 18 screens, `*` = interval excludes zero:

| contrast | delta | 95% CI | up |
|---|---|---|---|
| **best arm − AF alone** | **+10.28** | **[+2.78, +19.72]*** | 11/18 |
| **best arm − AF + DepMap** | **+6.00** | **[+0.39, +11.89]*** | 11/18 |
| best arm − HSA only | +4.00 | [−1.06, +9.83] | 8/18 |
| HSA only − AF alone | +6.28 | [+0.94, +12.39]* | 10/18 |
| AF + DepMap − AF alone | +4.28 | [−0.94, +9.95] | 9/18 |
| HSA + DepMap eq − best arm | −8.00 | [−19.00, +0.61] | 5/18 |
| AF + DepMap + HSA − AF + DepMap | +3.06 | [−2.22, +8.83] | 10/18 |

### 19.4 THE RESULT: the two gains ADD

    AF + DepMap − AF alone        +4.28
    HSA only    − AF alone        +6.28
                                 ------
    sum of the parts             +10.56
    combined    − AF alone       +10.28   measured

Within noise of perfect additivity. §18 predicted this from near-zero
`rho(cd, hsa)` and it held — which is the first time in this project that an
orthogonality screen has predicted an additive campaign outcome rather than a null
(contrast §4.21 co-hit, §4.28 BioPlex, both orthogonal-ish and both flat).

The combined arm beats **both** single methods with intervals excluding zero against
AF alone and against AF + DepMap. Its increment over HSA alone (+4.00) does not
exclude zero at n=18.

### 19.5 Two design facts the experiment settled

**(a) Weight 0.5 for DepMap, not 1.0.** `HSA + DepMap eq` (dm=1.0) is **−8.00** against
dm=0.5. §4.19 found both fusion weights peak at exactly 0.5 on the benchmark; that
result transfers to this ensemble unchanged.

**(b) Drop AssayFormer's explicit rank.** `AF + DepMap + HSA` keeps it and gains only
+3.06 over AF + DepMap; dropping it gains +6.00. The double-count is real and costs
~3 hits/screen, because at `af_weight=1` the AF term's `1/(60+0)=0.0167` dominates HSA's
`1/(60+0)` share once HSA is down-weighted — and HSA already carries AssayFormer in its
prior mean `z` and its first ten design columns.

### 19.6 A CORRECTION, and the lesson

An interim read of **four** screens (`1125`, `2346`, `2462`, `U_1793_dec`) said DepMap
was *monotonically harmful* on top of HSA. **That was wrong.** Those four happen to
contain the single worst case (`U_1793_dec`, −17) and three flat ones, and they missed
every winner: `174` 80 -> 118 (**+38**), `1985` 85 -> 106 (+21), `U_2082_dec` 88 -> 105
(+17), `2462` 258 -> 273 (+15). On the full panel the sign reverses to **+4.00**. The
monotonic-harm pattern is real only in the WEIGHT (0.5 helps, 1.0 and 2.0 hurt), not in
the screens. **Do not read a 4-of-18 subsample as the panel**, and note the subsample was
not random — it was chosen for high hit counts, which is exactly where DepMap's kernel
has least headroom.

### 19.7 Caveats, unchanged from §17

* `train` (3 screens) carries the largest gains (+70.7 vs 43.3 for AF alone) and those
  are the screens AssayFormer was trained on (§15.5); the reduced atlas drops cell lines
  from the four EVALUATION sets only, so the DepMap arm is not leakage-controlled there.
  `174` alone supplies +38 of the combined arm's margin.
* Four `validation` screens barely intersect their pool and score ~0 for every arm.
* n=18, fixed-pool, unknowns-charged, hits-not-EF. This does not revise the benchmark
  headline (+0.307 EF), and the availability component of HSA's margin cannot exist under
  a screen-library pool at all (§17.4). `measured` rises 672.9 -> 762.6 across the table.
* HSA ran at `draw=0` only. Their recorded arm averages draws 0 and 1; §17's HSA column
  is that mean, this section's `HSA only` is draw 0, which is why the two differ slightly.

---

## 20. HSA adds NO data beyond the training fold — and this reframes §4.7

Asked directly whether HSA informs its gene scores with anything beyond the training
corpus. It does not, and the bank files prove which screens they draw on.

### 20.1 The donor set is provably the training fold

`build_bank` drops the query's own publication (`keep = pubs != query_publication`), so a
bank's recorded `donor_publications` count reveals the pool: it loses one publication iff
the query screen is itself a donor. Checked against `public_train` membership on all 18
bundle banks:

| query screen's fold | in `public_train` | donor publications |
|---|---|---|
| `train` (`1125`, `1628`, `174`) | yes | **133** — own publication excluded |
| `validation` / `test` (15 screens) | no | **134** — nothing to exclude |

**18/18 match.** The donor pool is 134 publications over the 1,349 `public_train`
screens. No validation or test screen contributes to any bank. (The group count is also
independently confirmed by `min(availability) = 1/(2+G)`, §19 method notes.)

### 20.2 Every input to an HSA score

| component | derived from |
|---|---|
| `bank.residual` (n x 32) | `M`, `Y` of the 1,349 TRAINING screens, with AssayFormer's `W` projected out |
| `bank.availability` | `M` of the same screens, grouped by their 134 publications |
| `bank.scalar_c` | singular values of that residual matrix |
| `bank.mean` (32) | ridge map from training-screen description embeddings to their 32-D coords, applied to the query's description |
| `z` (10) | AssayFormer's screen latent: query description + labels so far |
| `W[ids]`, `b[ids]` | AssayFormer's trained gene embedding and bias |
| labels | the campaign's own reveals |

Complete list. **No DepMap, no PPI, no expression, no gene embeddings, no external atlas.**

One qualification: the description embeddings come from `text-embedding-3-small`, a
pretrained external model, so `bank.mean` channels literature knowledge indirectly. But
AssayFormer already depends on that same model for its own screen latent
(`text_model` in the checkpoint config), so it is not a dependency HSA introduces, and it
acts only on the screens' own metadata text -- a similarity metric over descriptions, not
gene-level biology.

### 20.3 The asymmetry with AF + DepMap, and why §18/§19 came out as they did

| | HSA | AF + DepMap |
|---|---|---|
| new biological data | **none** | DepMap, 808 cell lines x 18,119 genes |
| source | AssayBench train fold + AssayFormer's own weights | external CRISPR gene-effect atlas |
| the gain is | better EXTRACTION from data already in the system | genuinely NEW measurements |
| deployment cost | nothing to maintain | atlas held and re-ablated per evaluation set |
| bounded by | the training corpus | atlas coverage (18,119 genes) |

So §18's near-zero `rho(cd, hsa)` and §19's near-perfect additivity (+4.28 and +6.28
giving +10.28) are structural, not luck: the two terms read entirely different sources and
therefore cannot be redundant.

### 20.4 THIS CHANGES HOW §4.7 SHOULD BE READ

§4.7 established a hard ceiling: fitting the screen latent `u` to the COMPLETE answer key
-- a perfect oracle -- is worth only **+0.41 / +0.84 EF**. That bound is why thirteen LLM
warm-start variants were all null, and the project has quoted it as "the ceiling on
anything measured here".

**The correct reading was never "nothing can help" but "nothing acting THROUGH `u` can
help."** HSA does not optimise within the 10-D space; it ADDS 32 dimensions orthogonal to
it, which is the one move §4.7's bound does not cover. Empirically it is worth +6.28
hits/screen on the 18-screen panel.

So HSA is evidence that **AssayFormer is under-parameterised, not under-informed**: its own
training data contains hit structure its 10-D screen-conditioned term cannot represent,
and `build_bank`'s projection step isolates exactly that discarded part.

**Two consequences worth acting on.**

1. **A wider screen-conditioned term is now a candidate intervention.** HSA's gain is in
   principle reachable by retraining AssayFormer with more than 10 gene-embedding
   dimensions, or with a residual head. DepMap's gain is NOT reachable that way -- no
   retraining conjures 808 cell lines of co-essentiality. The two are not substitutes at
   the level of the fix, either.
2. **Any future candidate should be classified by whether it acts through `u`.** If it
   does, §4.7 caps it at +0.84 and it does not deserve a campaign. If it widens the
   representation or fuses outside the ranker, §4.7 is silent and the candidate is live.
   That is a cheaper triage rule than the eleven-feature screen of §4.3.

---

## 21. ATLAS LEAKAGE AUDIT: the by-construction claim FAILS for 61 screens (hTERT-RPE1)

Asked to confirm that no validation/test cell line survives in `achilles_noeval.npz`.
**It does not fully hold.** New reusable check: `scripts/diagnostics/audit_atlas_lines.py`
(exits 1 when it finds a gap).

### 21.1 The build is reproducible; the MATCHING is what leaks

`build_noeval_atlas.py --verify` reports **IDENTICAL** against the shipped file: 808 -> 768
lines, 40 dropped, 18119 genes all non-degenerate. So the artifact is exactly what the
builder produces. That was never the weak link.

The weak link is that the drop set is built by NAME: `resolve(s.cell_line, by)` normalises
to alphanumerics-uppercase, tries the string whole / cut at `(` / cut at `/`, and consults a
4-entry ALIAS table. A screen whose line DepMap carries under an unreachable name resolves
to `None`, nothing is dropped, and the failure is silent.

Over the 552 unique `validation` + `test` screens (93 distinct `cell_line` strings):

| | screens |
|---|---|
| resolve, line in Achilles -> **dropped** | 137 |
| resolve, line not in Achilles -> nothing to drop | 213 |
| **do NOT resolve** -> audited | **202** |

### 21.2 The audit, three passes

1. counts above
2. strip engineering decorations (`ACE2`, `TMPRSS2`, `SunTag`, `Cas9`, `GFP`,
   `deltaSTAT1`, `hTERT` prefix, trailing clone digits) and retry against lines that
   SURVIVE the ablation -> **0 hits, clean**
3. prefix search: any surviving Achilles line whose normalised name EXTENDS the screen's,
   which is how clonal subclones are named -> **one hit, and it is real**

### 21.3 The one real failure

**`hTERT-RPE1`, 61 screens (6 `validation`, 55 `test`).** DepMap holds four single-cell
subclones of that exact line, all still in the ablated atlas:

| ModelID | name | lineage |
|---|---|---|
| ACH-002462 | RPE1-ss48 | Eye |
| ACH-002463 | RPE1-ss77 | Eye |
| ACH-002464 | RPE1-ss6 | Eye |
| ACH-002467 | RPE1-ss51 | Eye |

`norm("hTERT-RPE1") = "HTERTRPE1"` and `norm("RPE1-ss48") = "RPE1SS48"` share no form, and
no ALIAS entry bridges them. It is also **one-to-many**, which the resolver cannot express
at all -- it returns a single ModelID.

So for **61 of 552 screens (11%)** the statement "no screen's own cell line is in the atlas
by construction" is FALSE. For the other 491 it holds.

### 21.4 Two near-misses that are protected only by luck

* **A549** (ACH-000681) is in Achilles, and 16 screens use `A549-ACE2` / `-SunTag` /
  `-deltaSTAT1-SunTag` / `-ACE2-Cas9-GFP1-10`, none of which resolve. It was dropped anyway
  because **2 other screen slots spell it `A-549`**, which does resolve.
* **Huh-7** (ACH-000480) likewise: 13 screens use `Huh-7.5` / `Huh-7.5.1-ACE2`, which do not
  resolve, but **19 slots spell it `HuH-7` or `Huh-7.5.1`** (the latter via
  `ALIAS["HUH751"]`), so the line went.

Both are correct outcomes reached by accident. Had those sibling spellings not existed, 29
more screens would be unprotected. **The resolver's coverage is not a property anyone
checked; it is a coincidence of how curators spelled things.**

### 21.5 What to do, and what NOT to do

The fix is small -- an explicit extra-drop set in `build_noeval_atlas.py`:

```python
# hTERT-RPE1 screens: DepMap carries four clonal subclones of that line under names
# resolve() cannot reach, and the mapping is one-to-many. 61 evaluation screens use it.
EXTRA_DROP = {"ACH-002462", "ACH-002463", "ACH-002464", "ACH-002467"}
```
giving a 764-line atlas.

**DO NOT overwrite `achilles_noeval.npz`.** Every recorded campaign read that file; silently
replacing it would make every stored `result.json` unreproducible and its provenance a lie.
Write `achilles_noeval2.npz`, then re-run the affected arms to MEASURE the impact rather
than assuming it.

Expected size of the effect: small. §4.54 measured the whole 40-line ablation at a few
hundredths of EF, and this adds 4 lines to that. But 55 of the 61 screens are in the `test`
fold, which is where §4.46's replication and §4.57(e)'s kNN increment live, so the arms to
re-score are named and few. **Until that is done, no result should be described as
"leakage-free by construction" without the hTERT-RPE1 exception stated.**

The §4.35 / §4.54 leakage claims in earlier sections are affected and should be read with
this caveat.

### 21.6 MEASURED: the leakage was real and the impact is NEGLIGIBLE (2026-10-31)

`achilles_noeval2.npz` built (**764** lines; `audit_atlas_lines.py` now reports CLEAN), and
AF + DepMap re-run on all **61** affected screens against it, paired with the stored v1
runs. `--legacy` still reproduces the historical 768-line file **IDENTICAL**, so nothing
on disk was invalidated. Script: `scripts/diagnostics/score_atlas_revision.py`. Results:
`output/runs/external/atlas_revision.json`.

Only the 61 can move: for every other screen the two atlases are byte-identical. AF alone
needed no re-run (`weight=0.0` loads the atlas and contributes nothing).

| | value |
|---|---|
| dEF mean | **−0.0428** [−0.1739, +0.0360] — **crosses zero** |
| dEF median | **+0.0000** |
| screens unchanged in hits | **43 / 61** |
| direction | 18 worse, **21 better** — a coin flip |
| total hits over the 61 screens | 3989 -> 3975 (**−14**) |

**The negative mean is one screen.** `2468` (H=22) goes 13 -> 12 hits for **dEF −3.498**,
because EF leverage is enormous when the hit set is tiny. Drop it and the mean **flips
sign** to **+0.0148 [−0.0106, +0.0447]**. Consistent with the split by sparsity: screens
with <50 total hits average −0.092, those with >=50 average +0.032.

**Effect on the headline: none.** The fold mean EF shifts by **−0.0012** on `validation`
(6 of 218 affected) and **−0.0070** on `test` (55 of 334). Against
**+0.307 [+0.222, +0.398]** per screen that is ~2% of the effect and nowhere near the
interval. §4.46's test-fold replication and §4.57(e)'s kNN increment are unaffected at
this resolution.

**Conclusion: a provenance bug, not a results bug.** The by-construction claim WAS false
for 11% of the folds (§21.3) and is now true, but correcting it changes no finding. The
honest statement going forward:

> Every reported AF + DepMap result used `achilles_noeval.npz`, in which four hTERT-RPE1
> subclones survived although 61 evaluation screens use that line. Re-running those 61
> against the corrected 764-line atlas moves the paired delta by −0.0012 (validation) and
> −0.0070 (test), well inside the reported interval.

**Which file to use.** New campaigns: `achilles_noeval2.npz`. Existing results: leave them
on `achilles_noeval.npz` and cite the sentence above -- re-running 552 screens to shift the
third decimal is not worth the compute, and mixing the two atlases inside one table would
be the §4.56 error again. The README's prerequisite block should be updated to build v2
once a re-run of the full folds is actually done; until then it correctly describes what
the recorded results used.

---

## 22. Where everything from §15-§21 now lives (2026-10-31)

The HSA work ran out of this session's scratchpad, which is ephemeral. All of it is now
durable. Nothing is in git: `output/` was already ignored and `_reproduction/` was added
to `.gitignore` (383 MB of third-party data).

### 22.1 `_reproduction/` — the frozen HSA inputs, 383 MB

Kept locally so none of this depends on a 187 MB re-fetch. Has its own
`_reproduction/README.md` covering the layout, the integrity numbers, the **two different
PYTHONPATH recipes** (their pinned upstream for their arms; our `src/` for `codep_seed`),
the four names the same two arms carry across their layers, and the `python -I` trap.

Verified after the move: tarball sha256 `19d409b2…93903` matches `REPRODUCTION_ASSETS.json`,
and `bundle/SHA256.json` re-checks **403 of 403**. A replay of screen `1125` from the new
path passes 3/3.

| path | what |
|---|---|
| `bundle/banks/*.npz` | the 18 frozen banks — the irreplaceable part |
| `bundle/screens/`, `bundle/UNIVERSE_*.json`, `bundle/recorded/` | screens, pools, 342 trajectories |
| `bundle/checkpoint/` | AssayFormer, byte-identical to our 4927703e snapshot |
| `assayloop-src/`, `vendor/` | upstream at `e7251afa…`, AssayBench 0.2.0 |
| `experiments/odysseus/assayloop_comparison/` | their code |
| `git-package/` | their vendored metadata |
| `assayloop-presentation-reproduction-*.tar.gz` | the original archive |

### 22.2 Results, all under `output/runs/external/`

| file | section | what |
|---|---|---|
| `bundle18_results.{json,csv}` | §17 | AF alone / AF + DepMap on the 18 bundle screens |
| `bundle18_sixarm_results.{json,csv}` | §19 | six arms incl. AF + DepMap + HSA |
| `hsa_depmap_orthogonality.json` | §18 | 50 screen-rounds of `rho(cd,hsa)` and top-100 overlaps |
| `hsa_replay/` (56 files) | §16 | the 54 bit-exact replayed campaigns |
| `atlas_revision.json` | §21.6 | the 61-screen atlas-leakage impact |

**One discrepancy to know about.** `hsa_depmap_orthogonality.json` was regenerated against
`achilles_noeval2.npz`; §18's table in this document was computed on
`achilles_noeval.npz`. Values differ in the third decimal (`1125` rho −0.006 -> −0.004,
`fused^hsa` 8.7 -> 8.9) and the conclusion is identical: `rho(cd, hsa)` is zero on all five
screens. The file is the better record; §18's table is the historical one.

`screen_hsa_depmap.py` gained `--out` so this is reproducible rather than scrollback.

### 22.3 Campaign runs and the rest

* `output/runs/<28 method folders>/` — 20,908 run dirs, 8.9 GB, indexed by
  `output/runs/README.md`; `baseline/` holds 794 LLM run dirs
* `output/runs/external/sources/` — 11 vendored raw source files, 191 MB, so the seven
  external datasets rebuild without re-downloading (`scripts/external/SOURCES.md`)
* `output/analysis/` — 115 MB, `logs_<method>/` plus the `.npz` method inputs
* `/home/jeff/depmap/` — `achilles.npz` (808), `achilles_noeval.npz` (768, what every
  recorded result used), `achilles_noeval2.npz` (764, corrected, §21.6)

### 22.4 Still ephemeral, deliberately

The 11 external `ScreenRecord` pickles were written to `/tmp` by the original builds and
are gone. Not worth storing: they rebuild from `output/runs/external/sources/` with the
commands in `scripts/external/SOURCES.md`, and the built pickles are derived data.

### 22.5 Nothing is committed

`git status` shows ~10 modified files under `src/`, plus `HANDOFF.md`,
`scripts/diagnostics/`, `scripts/external/`, `tests/` and the README change all untracked.
The 8.9 GB of runs is correctly ignored and must stay that way.
