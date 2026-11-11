# GPT-5.6 Sol → HSA Matched Validation

## Purpose

This experiment evaluates whether Historical Subspace Adaptation (HSA) improves
the numerical continuation of AssayLoop when both methods receive the same
GPT-5.6 Sol warm start.

This is a matched validation / extension experiment. It does not introduce a
new HSA method or regenerate the GPT-5.6 language-model calls.

## Experimental Design

The experiment replays Jeff's previously recorded GPT-5.6 Sol selections for
rounds 1–3 and then compares two numerical continuation strategies for rounds
4–10:

1. **AssayFormer** — standard AssayFormer continuation.
2. **HSA** — Historical Subspace Adaptation continuation using Odysseus's
   existing HSA implementation.

Both arms therefore receive the same GPT-5.6 warm-start selections before their
numerical policies diverge.

### Configuration

- Dataset: AssayLoop `paper_test`
- GPT-5.6 warm start: rounds 1–3
- Numerical continuation: rounds 4–10
- Batch size: 100
- Seed: 0
- Candidate universe: paper-test f2 universe
- Usable screens: 15 of 20
- AssayFormer checkpoint: frozen `model_last.pt`
- HSA banks: built from the frozen HSA reproduction assets

Five paper-test screens are excluded because the saved GPT-5.6 warm-start
artifact does not contain a complete three-round warm start for those screens.

## Results

Across the 15 usable screens:

| Method | Total Hits | Mean Hits | Mean Recall | Mean EF |
|---|---:|---:|---:|---:|
| AssayFormer | 1,580 | 105.33 | 0.2935 | 6.2058 |
| HSA | 1,652 | 110.13 | 0.2912 | 6.1580 |

HSA produced **72 additional hits overall**, a **4.6% increase in aggregate hit
count** relative to AssayFormer.

However, HSA did not improve mean recall or mean enrichment factor across all
screens. The aggregate hit gain is therefore not evidence that HSA uniformly
outperforms AssayFormer.

## Screen-Density Analysis

The HSA effect varied substantially by screen.

Using the median screen hit density to divide the 15 screens:

| Density Group | Mean Δ Hits | Mean Δ Recall | Mean Δ EF |
|---|---:|---:|---:|
| Denser half | +12.29 | +0.0143 | +0.3030 |
| Sparser half | -1.75 | -0.0168 | -0.3546 |

Screen hit density was positively associated with HSA's improvement in
screen-normalized performance:

- Δ recall vs. hit density: Spearman rho = **+0.575**, p = **0.0249**
- Δ EF vs. hit density: Spearman rho = **+0.575**, p = **0.0249**
- Δ hits vs. hit density: Spearman rho = **+0.456**, p = **0.0878**

These results provide exploratory evidence that HSA's benefit may depend on
screen characteristics, with stronger performance on denser screens and
possible degradation on some sparse screens.

Because this analysis contains only 15 screens, the density relationship should
be treated as a hypothesis-generating result rather than a definitive
generalization.

## Interpretation

The main contribution of this experiment is not the creation of a new method.
It is an independent cross-LLM validation of HSA under a GPT-5.6 Sol warm start,
combined with screen-level analysis of where HSA helps and where it can fail.

The results suggest:

- HSA can increase aggregate hit discovery after a GPT-5.6 warm start.
- The effect is heterogeneous across screens.
- Screen hit density appears to be associated with the direction and magnitude
  of the HSA effect.
- Additional validation across LLM actors and datasets is needed before making
  a general claim about when HSA should be preferred.

## Files

- `run_gpt56_matched.py` — matched experiment runner
- `results/gpt56_matched_methods.csv` — full method-level results for all
  15 screens
- `results/gpt56_screen_level_comparison.csv` — paired HSA vs. AssayFormer
  differences by screen
- `results/gpt56_density_summary.csv` — sparse/dense screen summary

## Provenance

This experiment builds on two existing team artifacts:

- **GPT-5.6 Sol warm start:** Jeff branch
  `3b83fb26e43fa2e42f6d45b383b97ebb82758bf1`
- **HSA implementation and reproduction assets:** Odysseus branch
  `2259d8c846d57a533a0f92ed25cf41917eec71e9`

The saved GPT-5.6 rounds are replayed from disk; this experiment does not make
new GPT-5.6 API calls.

## Scope

These results represent a matched comparison of the complete numerical
continuation policies used here. They should not be interpreted as isolating
the causal contribution of every individual HSA component.
