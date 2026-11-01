"""CodepSeedAcquisition: the nesting guarantee above all.

The whole design rests on weight=0 and the no-seeds case reproducing plain greedy
EXACTLY, so that any measured difference is attributable to the feature and the
worst case is parity rather than regression. These tests pin that down, plus the
coverage rule that an absent gene must not be penalised.
"""
from __future__ import annotations

import numpy as np
import pytest

from assaybench.core.types import ModelPrediction, Observation, StepRecord

from assayloop.acquisitions.codep_seed_acq import CodepSeedAcquisition
from assayloop.acquisitions.greedy_from_model import GreedyFromModel

GENES = [f"G{i:02d}" for i in range(12)]
LINES = 8


@pytest.fixture
def npz(tmp_path):
    """Profiles where G00-G03 covary strongly and G04-G07 covary separately.
    G08-G11 are absent from the atlas entirely."""
    rng = np.random.default_rng(0)
    a, b = rng.standard_normal(LINES), rng.standard_normal(LINES)
    cols, vals = [], []
    for i in range(4):
        cols.append(GENES[i]); vals.append(a + 0.05 * rng.standard_normal(LINES))
    for i in range(4, 8):
        cols.append(GENES[i]); vals.append(b + 0.05 * rng.standard_normal(LINES))
    p = tmp_path / "prof.npz"
    np.savez(p, values=np.array(vals, dtype=np.float32).T,
             rows=np.array([f"L{i}" for i in range(LINES)], dtype=object),
             cols=np.array(cols, dtype=object))
    return p


def _pred(order):
    """Scores descending in the given gene order."""
    return ModelPrediction(scores={g: float(len(order) - i) for i, g in enumerate(order)},
                           metadata={})


def _hist(hits):
    return [StepRecord(step=1, acquired_batch=list(hits), new_observations=[
        Observation(candidate=g, label={"hit": True, "relevance_score": 1.0}) for g in hits])]


# -- the nesting guarantee ---------------------------------------------------

def test_weight_zero_is_exactly_greedy(npz):
    acq = CodepSeedAcquisition(profile_npz=npz, weight=0.0)
    pred = _pred(GENES)
    got = acq.suggest(_hist(["G00", "G01"]), GENES, 4, pred)
    want = GreedyFromModel().suggest([], GENES, 4, pred)
    assert got == want
    assert acq.last_trace()["codep_active"] is False


def test_no_hits_yet_is_exactly_greedy(npz):
    """Round 1 has no seeds, so the acquisition must not deviate at all."""
    acq = CodepSeedAcquisition(profile_npz=npz, weight=5.0)
    pred = _pred(GENES)
    got = acq.suggest([], GENES, 4, pred)
    assert got == GreedyFromModel().suggest([], GENES, 4, pred)
    assert acq.last_trace()["reason"] == "no_seeds"


def test_one_hit_is_below_min_seeds_and_stays_greedy(npz):
    acq = CodepSeedAcquisition(profile_npz=npz, weight=5.0, min_seeds=2)
    pred = _pred(GENES)
    assert acq.suggest(_hist(["G00"]), GENES, 4, pred) == \
        GreedyFromModel().suggest([], GENES, 4, pred)


def test_screens_that_never_hit_never_deviate(npz):
    """A screen where nothing is a hit must behave identically to greedy for all
    ten rounds, which is what bounds the downside on the hardest screens."""
    acq = CodepSeedAcquisition(profile_npz=npz, weight=5.0)
    pred = _pred(GENES)
    hist = [StepRecord(step=i, acquired_batch=[], new_observations=[
        Observation(candidate=g, label={"hit": False, "relevance_score": 0.0})
        for g in GENES[:3]]) for i in range(1, 4)]
    assert acq.suggest(hist, GENES, 4, pred) == GreedyFromModel().suggest([], GENES, 4, pred)


# -- the feature actually fires ---------------------------------------------

def test_covarying_genes_get_promoted(npz):
    """With G00/G01 confirmed, their co-varying partners G02/G03 should be pulled
    up past genes AssayFormer ranked higher."""
    # AssayFormer prefers the G04+ block; codep should rescue G02/G03
    pred = _pred(["G04", "G05", "G06", "G07", "G02", "G03"] + GENES[8:])
    plain = GreedyFromModel().suggest([], list(pred.scores), 4, pred)
    acq = CodepSeedAcquisition(profile_npz=npz, weight=3.0)
    got = acq.suggest(_hist(["G00", "G01"]), list(pred.scores), 4, pred)
    assert "G02" in got or "G03" in got
    assert got != plain
    assert acq.last_trace()["codep_active"] is True


def test_zsum_and_rrf_both_run_and_report(npz):
    pred = _pred(GENES)
    for fuse in ("rrf", "zsum"):
        acq = CodepSeedAcquisition(profile_npz=npz, weight=1.0, fuse=fuse)
        got = acq.suggest(_hist(["G00", "G01"]), GENES, 4, pred)
        assert len(got) == 4 and len(set(got)) == 4
        assert acq.last_trace()["fuse"] == fuse


# -- the coverage rule -------------------------------------------------------

def test_genes_absent_from_the_atlas_are_not_penalised(npz):
    """G08-G11 have no profile. Possessing a co-dependency rank must not itself be
    worth anything, or the fusion reintroduces the DepMap-coverage confound: an
    uncovered gene AssayFormer ranks first has to survive in the batch.

    The already-tested seeds are excluded from the candidate list, as the runner
    does in a real campaign.
    """
    cands = ["G08", "G09"] + GENES[2:8] + ["G10", "G11"]
    pred = _pred(cands)
    acq = CodepSeedAcquisition(profile_npz=npz, weight=1.0, fuse="rrf")
    got = acq.suggest(_hist(["G00", "G01"]), cands, 3, pred)
    assert "G08" in got, "an uncovered top-ranked gene was dropped by the fusion"
    assert acq.last_trace()["n_covered"] == 6


def test_coverage_alone_confers_no_advantage(npz):
    """With exactly ONE covered candidate, the co-dependency term is a constant:
    that gene contributes 1/(k+0) and every uncovered gene contributes the mean of
    the covered contributions, which is the same number. A constant added to every
    candidate cannot reorder them, so the batch must equal plain greedy. If
    uncovered genes were instead left unranked, this gene would jump the queue.
    """
    cands = ["G02", "G08", "G09", "G10", "G11"]      # one covered, four not
    pred = _pred(["G09", "G10", "G02", "G11", "G08"])
    acq = CodepSeedAcquisition(profile_npz=npz, weight=5.0, fuse="rrf")
    got = acq.suggest(_hist(["G00", "G01"]), cands, 3, pred)
    assert got == GreedyFromModel().suggest([], cands, 3, pred)
    assert acq.last_trace()["n_covered"] == 1


def test_trace_reports_overlap_with_greedy(npz):
    pred = _pred(GENES)
    acq = CodepSeedAcquisition(profile_npz=npz, weight=1.0)
    acq.suggest(_hist(["G00", "G01"]), GENES, 5, pred)
    t = acq.last_trace()
    assert 0 <= t["overlap_with_greedy"] <= 5


# -- argument validation -----------------------------------------------------

def test_bad_fuse_rejected(npz):
    with pytest.raises(ValueError):
        CodepSeedAcquisition(profile_npz=npz, fuse="softmax")


def test_batch_larger_than_candidates(npz):
    acq = CodepSeedAcquisition(profile_npz=npz, weight=1.0)
    got = acq.suggest(_hist(["G00", "G01"]), GENES[:3], 10, _pred(GENES[:3]))
    assert len(got) == 3


# -- the screen-kNN second scorer -------------------------------------------
#
# v_g is a 10-number compression of 1349 training screens; this scorer reads the
# raw hit matrix instead. Offline the two scorers are near-perfectly additive, so
# both must be able to fire independently and both must respect the nesting rule.

@pytest.fixture
def trainhits(tmp_path):
    """4 training screens. Screens 0 and 1 both hit {G00,G01,G02,G03}; screen 2
    hits a disjoint set; screen 3 hits everything (so size correction matters)."""
    from scipy.sparse import csr_matrix, save_npz
    genes = GENES
    gi = {g: i for i, g in enumerate(genes)}
    hits = [
        ["G00", "G01", "G02", "G03"],
        ["G00", "G01", "G02", "G03"],
        ["G06", "G07"],
        genes,
    ]
    rows = [r for r, hs in enumerate(hits) for _ in hs]
    cols = [gi[g] for hs in hits for g in hs]
    M = csr_matrix((np.ones(len(rows), dtype=np.float32), (rows, cols)),
                   shape=(len(hits), len(genes)))
    pre = tmp_path / "th"
    save_npz(pre.with_suffix(".mat.npz"), M)
    np.savez(pre.with_suffix(".npz"),
             gene_names=np.array(genes, dtype=object),
             screen_names=np.array([f"S{i}" for i in range(len(hits))], dtype=object),
             n_hits=np.asarray(M.sum(axis=1)).ravel())
    return pre


def test_knn_weight_zero_leaves_behaviour_unchanged(npz, trainhits):
    """Adding the scorer with weight 0 must not perturb the codep-only result."""
    pred = _pred(GENES)
    a = CodepSeedAcquisition(profile_npz=npz, weight=1.0)
    b = CodepSeedAcquisition(profile_npz=npz, weight=1.0,
                             train_hits=trainhits, knn_weight=0.0)
    h = _hist(["G00", "G01"])
    assert a.suggest(h, GENES, 4, pred) == b.suggest(h, GENES, 4, pred)


def test_both_weights_zero_is_exactly_greedy(npz, trainhits):
    pred = _pred(GENES)
    acq = CodepSeedAcquisition(profile_npz=npz, weight=0.0,
                               train_hits=trainhits, knn_weight=0.0)
    assert acq.suggest(_hist(["G00", "G01"]), GENES, 4, pred) == \
        GreedyFromModel().suggest([], GENES, 4, pred)


def test_knn_alone_fires_without_the_profile_scorer(npz, trainhits):
    """knn_weight > 0 with weight = 0 must still deviate from greedy: the two
    scorers are independent, so either should work on its own."""
    cands = ["G02", "G03", "G08", "G09"]
    pred = _pred(["G08", "G09", "G02", "G03"])
    acq = CodepSeedAcquisition(profile_npz=npz, weight=0.0,
                               train_hits=trainhits, knn_weight=5.0, knn_k=2)
    got = acq.suggest(_hist(["G00", "G01"]), cands, 2, pred)
    # seeds G00,G01 match training screens 0 and 1, whose other hits are G02,G03
    assert "G02" in got or "G03" in got
    t = acq.last_trace()
    assert t["knn_active"] is True and t["codep_active"] is False


def test_knn_promotes_genes_from_similar_training_screens(npz, trainhits):
    pred = _pred(["G06", "G07", "G08", "G02", "G03"])
    cands = list(pred.scores)
    plain = GreedyFromModel().suggest([], cands, 2, pred)
    acq = CodepSeedAcquisition(profile_npz=npz, weight=0.0,
                               train_hits=trainhits, knn_weight=8.0, knn_k=2)
    got = acq.suggest(_hist(["G00", "G01"]), cands, 2, pred)
    assert got != plain


def test_knn_needs_min_seeds(npz, trainhits):
    acq = CodepSeedAcquisition(profile_npz=npz, weight=0.0, min_seeds=3,
                               train_hits=trainhits, knn_weight=5.0)
    pred = _pred(GENES)
    got = acq.suggest(_hist(["G00", "G01"]), GENES, 4, pred)   # only 2 seeds
    assert got == GreedyFromModel().suggest([], GENES, 4, pred)


def test_name_reports_both_scorers(npz, trainhits):
    a = CodepSeedAcquisition(profile_npz=npz, weight=0.5)
    b = CodepSeedAcquisition(profile_npz=npz, weight=0.5,
                             train_hits=trainhits, knn_weight=0.5)
    assert "knn" not in a.name() and "knn0.5" in b.name()


# -- the query direction ----------------------------------------------------
#
# `diffmeans` (the default from 2026-10-22) subtracts the non-hit mean from the hit
# mean. Offline it roughly doubles the centroid's margin. These tests pin the
# fallback that keeps it safe when a screen has barely any labelled non-hits, and
# that `centroid` still reproduces every pre-2026-10-22 result.

def _mixed_hist(hits, misses):
    obs = [Observation(candidate=g, label={"hit": True, "relevance_score": 1.0}) for g in hits]
    obs += [Observation(candidate=g, label={"hit": False, "relevance_score": 0.0}) for g in misses]
    return [StepRecord(step=1, acquired_batch=list(hits) + list(misses), new_observations=obs)]


def test_diffmeans_falls_back_to_centroid_with_too_few_nonhits(npz):
    """Under ~10 labelled non-hits the non-hit mean is unusable, so the direction
    must revert to the hit centroid rather than being driven by two genes."""
    pred = _pred(GENES)
    h = _mixed_hist(["G00", "G01"], ["G06"])          # only one non-hit
    a = CodepSeedAcquisition(profile_npz=npz, weight=2.0, query="diffmeans")
    b = CodepSeedAcquisition(profile_npz=npz, weight=2.0, query="centroid")
    assert a.suggest(h, GENES, 4, pred) == b.suggest(h, GENES, 4, pred)


def test_diffmeans_differs_from_centroid_with_enough_nonhits(npz):
    """G00-G03 covary; G04-G07 covary separately. With the G04 block labelled as
    non-hits, subtracting their mean should change the ranking."""
    hits, misses = ["G00", "G01"], ["G04", "G05", "G06", "G07"]
    cands = ["G02", "G03", "G06", "G07"]
    pred = ModelPrediction(scores={g: 1.0 for g in cands}, metadata={})
    h = _mixed_hist(hits, misses)
    # min_nonhits counts non-hits PRESENT IN THE ATLAS, and this fixture only has
    # 8 such genes, so the threshold has to come down for the branch to fire at all
    a = CodepSeedAcquisition(profile_npz=npz, weight=5.0, query="diffmeans", min_nonhits=3)
    got = a.suggest(h, cands, 2, pred)
    assert a.last_trace()["query"] == "diffmeans"
    # subtracting the G04-block mean must leave the hits' own covarying partners on top
    assert set(got) == {"G02", "G03"}


def test_centroid_mode_is_still_available(npz):
    """Every result before 2026-10-22 used the centroid; it must stay reproducible."""
    acq = CodepSeedAcquisition(profile_npz=npz, weight=1.0, query="centroid")
    assert acq.query == "centroid"
    assert "centroid" in acq.name()


def test_query_appears_in_the_name_and_trace(npz):
    acq = CodepSeedAcquisition(profile_npz=npz, weight=1.0, query="diffmeans")
    assert "diffmeans" in acq.name()
    acq.suggest(_mixed_hist(["G00", "G01"], GENES[4:]), GENES[:4], 2, _pred(GENES[:4]))
    assert acq.last_trace()["query"] == "diffmeans"


def test_bad_query_rejected(npz):
    with pytest.raises(ValueError):
        CodepSeedAcquisition(profile_npz=npz, query="fisher")


def test_weight_zero_still_greedy_under_diffmeans(npz):
    """The nesting guarantee must hold for the new default too."""
    pred = _pred(GENES)
    acq = CodepSeedAcquisition(profile_npz=npz, weight=0.0, query="diffmeans")
    assert acq.suggest(_mixed_hist(["G00", "G01"], GENES[4:]), GENES, 4, pred) == \
        GreedyFromModel().suggest([], GENES, 4, pred)


# -- hierarchical shrinkage of the query direction ---------------------------
#
# The live query direction is estimated from ~20-30 hits in 808-dimensional space.
# Shrinkage pulls it toward the average optimal direction of the most similar
# TRAINING screens, which have complete labels. The nesting rule applies here too:
# alpha=0 must leave the unshrunk behaviour untouched bit for bit.


@pytest.fixture
def dirbank(tmp_path, npz):
    """One direction per training screen, in the same space the acquisition scores
    in. Deliberately contrarian: the two NEAREST training screens (S0, S1 -- the
    ones that hit G00-G03, as the live screen does) point at the OTHER block, so a
    test can tell whether the prior is really being followed. The more distant S3
    points at the live screen's own block, which is what `shrink_drop_top` exposes."""
    z = np.load(npz, allow_pickle=True)
    P = np.asarray(z["values"], dtype=np.float64).T
    P -= P.mean(axis=1, keepdims=True)
    P /= np.linalg.norm(P, axis=1, keepdims=True)
    cols = {str(c): i for i, c in enumerate(z["cols"])}
    a, b = P[cols["G00"]], P[cols["G04"]]
    p = tmp_path / "bank.npz"
    np.savez(p, dirs=np.array([b, b, a, a], dtype=np.float32), n_lines=LINES,
             screen_names=np.array([f"S{i}" for i in range(4)], dtype=object),
             n_hits=np.array([100, 100, 100, 100]), n_nonhits=np.array([9] * 4))
    return p


def _shr(npz, trainhits, dirbank, **kw):
    return CodepSeedAcquisition(profile_npz=npz, weight=5.0, query="centroid",
                                train_hits=trainhits, dir_bank=dirbank, **kw)


def test_shrink_alpha_zero_is_bit_identical_to_no_bank(npz, trainhits, dirbank):
    """The nesting guarantee, extended to shrinkage: supplying a bank at alpha=0
    must not perturb anything, so every result measured before this existed stays
    reproducible."""
    pred = _pred(GENES[8:] + ["G02", "G03", "G04", "G05", "G06", "G07", "G00", "G01"])
    hist = _hist(["G00", "G01"])
    with_bank = _shr(npz, trainhits, dirbank, shrink_alpha=0.0)
    without = CodepSeedAcquisition(profile_npz=npz, weight=5.0, query="centroid")
    assert with_bank.suggest(hist, GENES, 4, pred) == without.suggest(hist, GENES, 4, pred)
    assert with_bank._bank is None          # not even loaded


def test_alpha_one_follows_the_borrowed_direction(npz, trainhits, dirbank):
    """alpha=1 discards the live estimate. The two nearest training screens point
    at the G04-G07 block, so that block must be promoted even though the live
    screen's own hits are G00/G01."""
    pred = _pred(GENES[8:] + ["G02", "G03", "G04", "G05", "G06", "G07", "G00", "G01"])
    hist = _hist(["G00", "G01"])
    own = _shr(npz, trainhits, dirbank, shrink_alpha=0.0).suggest(hist, GENES, 4, pred)
    borrowed = _shr(npz, trainhits, dirbank, shrink_alpha=1.0, shrink_k=2
                    ).suggest(hist, GENES, 4, pred)
    b_block = {"G04", "G05", "G06", "G07"}
    assert len(b_block & set(borrowed)) > len(b_block & set(own))


def test_drop_top_changes_which_screens_are_borrowed_from(npz, trainhits, dirbank):
    """Near-duplicate screens recur across this benchmark, so a neighbour can be
    effectively the same assay. Dropping the two nearest leaves only S3, which
    points at the live screen's own block."""
    pred = _pred(GENES[8:] + ["G02", "G03", "G04", "G05", "G06", "G07", "G00", "G01"])
    hist = _hist(["G00", "G01"])
    near = _shr(npz, trainhits, dirbank, shrink_alpha=1.0, shrink_k=2)
    far = _shr(npz, trainhits, dirbank, shrink_alpha=1.0, shrink_k=2, shrink_drop_top=2)
    assert near.suggest(hist, GENES, 4, pred) != far.suggest(hist, GENES, 4, pred)
    assert far.last_trace()["n_prior_screens"] == 1


def test_partial_alpha_lands_between_the_two(npz, trainhits, dirbank):
    """A blend must actually blend: the cosine to the prior has to rise with alpha."""
    pred = _pred(GENES[8:] + ["G02", "G03", "G04", "G05", "G06", "G07", "G00", "G01"])
    hist = _hist(["G00", "G01"])
    acq = _shr(npz, trainhits, dirbank, shrink_alpha=0.5, shrink_k=2)
    acq.suggest(hist, GENES, 4, pred)
    # cos is recorded BEFORE blending, so it measures own-vs-prior disagreement
    assert acq.last_trace()["cos_q_prior"] < 0.5
    assert acq.last_trace()["shrink_alpha"] == 0.5


def test_shrinkage_needs_the_training_matrix_for_neighbours(npz, dirbank):
    """A bank without a similarity source is a silent no-op, so refuse it loudly."""
    with pytest.raises(ValueError, match="train_hits"):
        CodepSeedAcquisition(profile_npz=npz, dir_bank=dirbank, shrink_alpha=0.5)


@pytest.mark.parametrize("bad", [-0.1, 1.5])
def test_alpha_outside_the_unit_interval_raises(npz, trainhits, dirbank, bad):
    with pytest.raises(ValueError, match="shrink_alpha"):
        _shr(npz, trainhits, dirbank, shrink_alpha=bad)


def test_underdetermined_training_screens_are_dropped_not_downweighted(npz, trainhits, dirbank):
    """A training screen with a handful of hits has a direction as noisy as the one
    being shrunk; requiring min_hits above every screen must empty the bank and
    leave the live estimate alone rather than quietly blending in noise."""
    pred = _pred(GENES[8:] + ["G02", "G03", "G04", "G05", "G06", "G07", "G00", "G01"])
    hist = _hist(["G00", "G01"])
    acq = _shr(npz, trainhits, dirbank, shrink_alpha=1.0, shrink_min_hits=10_000)
    plain = _shr(npz, trainhits, dirbank, shrink_alpha=0.0)
    assert acq.suggest(hist, GENES, 4, pred) == plain.suggest(hist, GENES, 4, pred)
    assert acq.last_trace()["n_prior_screens"] == 0


def test_round_one_is_still_exactly_greedy_under_shrinkage(npz, trainhits, dirbank):
    """The invariant that bounds the downside on hard screens: with no hits there is
    no query to shrink, and a borrowed direction must not sneak in."""
    pred = _pred(GENES)
    acq = _shr(npz, trainhits, dirbank, shrink_alpha=1.0)
    assert acq.suggest([], GENES, 4, pred) == GreedyFromModel().suggest([], GENES, 4, pred)
    assert acq.last_trace()["reason"] == "no_seeds"


def test_shuffle_is_a_control_and_is_visible_in_the_name(npz, trainhits, dirbank):
    """Directions share a large common component (pairwise cosine ~ +0.43 on the
    real bank), so shrinking toward an arbitrary training screen may help by
    itself. The shuffled arm isolates that, and must be impossible to mistake for
    the method in a sweep listing."""
    acq = _shr(npz, trainhits, dirbank, shrink_alpha=0.5, shrink_shuffle=0)
    assert "_SHUF" in acq.name()
    assert np.allclose(np.sort(np.linalg.norm(acq._bank, axis=1)),
                       np.sort(np.linalg.norm(
                           _shr(npz, trainhits, dirbank, shrink_alpha=0.5)._bank, axis=1)))


def test_name_records_the_shrinkage_settings(npz, trainhits, dirbank):
    acq = _shr(npz, trainhits, dirbank, shrink_alpha=0.25, shrink_k=10, shrink_drop_top=3)
    assert acq.name() == "codep_seed_centroid_w5_rrf_shr0.25k10d3"


# -- a second profile source (gene-gene co-hit, §4.21) ----------------------
#
# The deployed fusion carries one atlas plus screen-kNN. A second atlas enters as
# its own term, so the same two invariants have to hold for it: weight2=0 must be
# bit-identical to the single-source result, and the scorer must fire independently
# of the first.


@pytest.fixture
def npz2(tmp_path):
    """A SECOND, different space over the same genes: here G00/G04 covary and the
    rest are noise, so it disagrees with the `npz` fixture's blocks. That makes the
    two sources distinguishable in a test rather than interchangeable."""
    rng = np.random.default_rng(7)
    base = rng.standard_normal(LINES)
    cols, vals = [], []
    for g in ("G00", "G04"):
        cols.append(g); vals.append(base + 0.02 * rng.standard_normal(LINES))
    for g in ("G01", "G02", "G03", "G05", "G06", "G07"):
        cols.append(g); vals.append(rng.standard_normal(LINES))
    p = tmp_path / "prof2.npz"
    np.savez(p, values=np.array(vals, dtype=np.float32).T,
             rows=np.array([f"L{i}" for i in range(LINES)], dtype=object),
             cols=np.array(cols, dtype=object))
    return p


def test_weight2_zero_is_bit_identical_to_single_source(npz, npz2):
    """The nesting guarantee for the second source: supplying it at weight 0 must
    not perturb anything, so every result measured before it existed stands."""
    pred = _pred(GENES)
    hist = _hist(["G00", "G01"])
    one = CodepSeedAcquisition(profile_npz=npz, weight=1.0, query="centroid")
    two = CodepSeedAcquisition(profile_npz=npz, weight=1.0, query="centroid",
                               profile_npz2=npz2, weight2=0.0)
    assert two.suggest(hist, GENES, 5, pred) == one.suggest(hist, GENES, 5, pred)
    assert two._P2 is None            # not even loaded


def test_second_source_fires_and_changes_the_order(npz, npz2):
    """With G00 confirmed, source 2 should pull its partner G04 up -- a gene the
    first source groups with the OTHER block, so the effect is attributable."""
    pred = _pred(["G01", "G02", "G03", "G05", "G06", "G07", "G04"] + GENES[8:])
    hist = _hist(["G00", "G01"])
    off = CodepSeedAcquisition(profile_npz=npz, weight=1.0, query="centroid")
    on = CodepSeedAcquisition(profile_npz=npz, weight=1.0, query="centroid",
                              profile_npz2=npz2, weight2=8.0)
    got = on.suggest(hist, GENES, 3, pred)
    assert "G04" in got and "G04" not in off.suggest(hist, GENES, 3, pred)
    assert on.last_trace()["codep2_active"] is True


def test_second_source_works_with_the_first_at_weight_zero(npz, npz2):
    """Either scorer must fire independently, or the additivity measured offline
    cannot be realised."""
    pred = _pred(["G01", "G02", "G03", "G05", "G06", "G07", "G04"] + GENES[8:])
    acq = CodepSeedAcquisition(profile_npz=npz, weight=0.0, query="centroid",
                               profile_npz2=npz2, weight2=8.0)
    got = acq.suggest(_hist(["G00", "G01"]), GENES, 3, pred)
    assert "G04" in got
    assert acq.last_trace()["codep_active"] is False
    assert acq.last_trace()["codep2_active"] is True


def test_all_three_scorers_off_is_exactly_greedy(npz, npz2, trainhits):
    acq = CodepSeedAcquisition(profile_npz=npz, weight=0.0, profile_npz2=npz2,
                               weight2=0.0, train_hits=trainhits, knn_weight=0.0)
    pred = _pred(GENES)
    assert acq.suggest(_hist(["G00", "G01"]), GENES, 4, pred) == \
        GreedyFromModel().suggest([], GENES, 4, pred)


def test_round_one_stays_greedy_with_two_sources(npz, npz2):
    """The invariant that bounds the downside on hard screens, with both atlases."""
    acq = CodepSeedAcquisition(profile_npz=npz, weight=5.0, profile_npz2=npz2,
                               weight2=5.0)
    pred = _pred(GENES)
    assert acq.suggest([], GENES, 4, pred) == GreedyFromModel().suggest([], GENES, 4, pred)


def test_genes_absent_from_the_second_source_are_neutral(npz, npz2):
    """G08-G11 are in neither atlas. Under RRF, merely HAVING a score is worth up to
    weight/60, so they must take the mean contribution rather than go unranked --
    the coverage confound, which must not re-enter through a new term."""
    pred = _pred(GENES)
    acq = CodepSeedAcquisition(profile_npz=npz, weight=1.0, profile_npz2=npz2,
                               weight2=1.0, query="centroid")
    got = acq.suggest(_hist(["G00", "G01"]), GENES, 12, pred)
    assert set(got) == set(GENES)        # all ranked, none dropped


def test_name_records_the_second_source(npz, npz2):
    acq = CodepSeedAcquisition(profile_npz=npz, weight=0.5, profile_npz2=npz2,
                               weight2=0.25, query="diffmeans")
    assert acq.name() == "codep_seed_diffmeans_w0.5_rrf_p20.25"


# -- the sparse interaction kernel (BioPlex, §4.28 / §4.42) ------------------
#
# The point of this block is the EXACTNESS test. §4.28 rejected DepMap+BioPlex partly
# because the dense profile form is ~2.5 GB, and the sparse path only earns that
# rejection's reversal if it computes the identical quantity rather than an
# approximation of it. `test_sparse_path_equals_the_dense_computation` is the one that
# matters; the rest pin nesting and neutrality as elsewhere in this file.

@pytest.fixture
def ppinpz(tmp_path):
    """A weighted symmetric adjacency where G00-G03 interact with one set of partners
    and G04-G07 with a disjoint set, saved in `bioplex_to_npz.py`'s layout."""
    from scipy.sparse import csr_matrix, save_npz
    n = 20                                     # 12 genes + 8 partner-only proteins
    names = GENES + [f"P{i}" for i in range(8)]
    rows, cols, w = [], [], []
    for i in range(4):                         # G00-G03 -> P0..P3
        for j in range(12, 16):
            rows += [i, j]; cols += [j, i]; w += [0.9, 0.9]
    for i in range(4, 8):                      # G04-G07 -> P4..P7
        for j in range(16, 20):
            rows += [i, j]; cols += [j, i]; w += [0.9, 0.9]
    A = csr_matrix((w, (rows, cols)), shape=(n, n), dtype=np.float64)
    A.sum_duplicates()
    rs = np.asarray(A.sum(axis=1)).ravel()
    ssq = np.asarray(A.multiply(A).sum(axis=1)).ravel()
    p = tmp_path / "ppi.npz"
    save_npz(p.with_suffix(".mat.npz"), A)
    np.savez(p, cols=np.array(names, dtype=object), row_sum=rs, row_ssq=ssq,
             norm=np.sqrt(np.maximum(ssq - rs ** 2 / n, 0.0)), n_dims=np.int64(n))
    return p


def test_ppi_weight_zero_leaves_behaviour_unchanged(npz, ppinpz):
    a = CodepSeedAcquisition(profile_npz=npz, weight=0.5)
    b = CodepSeedAcquisition(profile_npz=npz, weight=0.5, ppi_npz=ppinpz, ppi_weight=0.0)
    pred, hist = _pred(GENES), _hist(["G00", "G01"])
    assert a.suggest(hist, GENES, 4, pred) == b.suggest(hist, GENES, 4, pred)


def test_ppi_alone_is_greedy_without_seeds(npz, ppinpz):
    acq = CodepSeedAcquisition(profile_npz=npz, weight=0.0,
                               ppi_npz=ppinpz, ppi_weight=0.5)
    greedy = GreedyFromModel()
    pred = _pred(GENES)
    assert acq.suggest([], GENES, 5, pred) == greedy.suggest([], GENES, 5, pred)
    assert acq.last_trace()["ppi_active"] is False


def test_ppi_promotes_interaction_partners_of_the_hits(npz, ppinpz):
    """G04-G07 share partners with each other; seeding on G04/G05 should pull G06/G07
    up from the bottom of the greedy order."""
    acq = CodepSeedAcquisition(profile_npz=npz, weight=0.0,
                               ppi_npz=ppinpz, ppi_weight=4.0)
    order = ["G08", "G09", "G10", "G11", "G00", "G01", "G02", "G03", "G04", "G05",
             "G06", "G07"]
    picked = acq.suggest(_hist(["G04", "G05"]), GENES, 4, _pred(order))
    assert acq.last_trace()["ppi_active"] is True
    assert {"G06", "G07"} & set(picked)


def test_sparse_path_equals_the_dense_computation(npz, ppinpz):
    """THE test for this block: score every candidate by densifying the adjacency,
    centring each row, unit-normalising and taking the cosine against the same query,
    and require the sparse path's ranking to match. If the centring-cancels identity
    in `_ppi_scores` were wrong, this is what would catch it."""
    from scipy.sparse import load_npz
    acq = CodepSeedAcquisition(profile_npz=npz, weight=0.0,
                               ppi_npz=ppinpz, ppi_weight=1.0)
    hits, misses = ["G00", "G01"], [f"G{i:02d}" for i in range(4, 8)]
    hist = [StepRecord(step=1, acquired_batch=hits + misses, new_observations=[
        Observation(candidate=g, label={"hit": g in hits, "relevance_score": 1.0})
        for g in hits + misses])]

    A = load_npz(str(ppinpz).replace(".npz", ".mat.npz")).toarray()
    meta = np.load(ppinpz, allow_pickle=True)
    names = [str(x) for x in meta["cols"]]
    row = {g: i for i, g in enumerate(names)}
    D = A - A.mean(axis=1, keepdims=True)                 # centre every row
    nrm = np.linalg.norm(D, axis=1, keepdims=True)
    U = D / np.where(nrm > 0, nrm, 1.0)                   # dense unit centred rows
    q = U[[row[g] for g in hits]].mean(axis=0)
    q = q - U[[row[g] for g in misses]].mean(axis=0)      # min_nonhits default is 10
    dense = U[[row[g] for g in GENES]] @ (q / np.linalg.norm(q))

    sparse = acq._ppi_scores(hist, GENES)
    assert sparse is not None
    ok = np.isfinite(sparse)
    # diffmeans needs >= min_nonhits non-hits; with 4 it falls back to the centroid,
    # so compare against the centroid form the acquisition actually used
    qc = U[[row[g] for g in hits]].mean(axis=0)
    dense_c = U[[row[g] for g in GENES]] @ (qc / np.linalg.norm(qc))
    np.testing.assert_allclose(sparse[ok], dense_c[ok], rtol=1e-10, atol=1e-12)
    assert not np.allclose(dense[ok], dense_c[ok])         # the two forms do differ


def test_ppi_uncovered_genes_are_not_penalised(npz, ppinpz):
    """A gene with no interactions has an all-zero row and zero norm; it must land at
    the mean contribution, never below every scored gene."""
    acq = CodepSeedAcquisition(profile_npz=npz, weight=0.0,
                               ppi_npz=ppinpz, ppi_weight=1.0)
    extra = GENES + ["ZZZ1", "ZZZ2"]
    picked = acq.suggest(_hist(["G00", "G01"]), extra, 6,
                         _pred(["ZZZ1", "ZZZ2"] + GENES))
    assert {"ZZZ1", "ZZZ2"} <= set(picked)


def test_ppi_degree_residual_changes_the_ranking_and_the_name(npz, ppinpz):
    plain = CodepSeedAcquisition(profile_npz=npz, weight=0.0,
                                 ppi_npz=ppinpz, ppi_weight=1.0)
    resid = CodepSeedAcquisition(profile_npz=npz, weight=0.0, ppi_npz=ppinpz,
                                 ppi_weight=1.0, ppi_degree_residual=True)
    assert "ppi1" in plain.name() and "ppi1dr" in resid.name()
    hist = _hist(["G00", "G01"])
    a = plain._ppi_scores(hist, GENES)
    b = resid._ppi_scores(hist, GENES)
    ok = np.isfinite(a) & np.isfinite(b)
    assert not np.allclose(a[ok], b[ok])


def test_ppi_name_reports_all_three_scorers(npz, trainhits, ppinpz):
    acq = CodepSeedAcquisition(profile_npz=npz, weight=0.5, train_hits=trainhits,
                               knn_weight=0.5, ppi_npz=ppinpz, ppi_weight=0.25)
    n = acq.name()
    assert "w0.5" in n and "knn0.5" in n and "ppi0.25" in n


# -- atlas geometry: PC removal and amplitude weighting (§4.44) ---------------
#
# Both knobs must be nested at their off value, as every other knob in this file is,
# because that is what makes a measured difference attributable to the knob rather than
# to a rewrite.

def test_drop_pcs_zero_is_bit_identical(npz):
    a = CodepSeedAcquisition(profile_npz=npz, weight=0.5)
    b = CodepSeedAcquisition(profile_npz=npz, weight=0.5, drop_pcs=0)
    np.testing.assert_array_equal(a._P, b._P)
    assert "pc" not in b.name()


def test_amp_power_zero_is_bit_identical(npz):
    a = CodepSeedAcquisition(profile_npz=npz, weight=0.5)
    b = CodepSeedAcquisition(profile_npz=npz, weight=0.5, amp_power=0.0)
    pred, hist = _pred(GENES), _hist(["G00", "G01"])
    assert a.suggest(hist, GENES, 4, pred) == b.suggest(hist, GENES, 4, pred)
    assert "amp" not in b.name()


def test_drop_pcs_removes_the_shared_direction(npz):
    """The fixture's two gene blocks covary within themselves, so the leading component
    is shared structure. Removing it must leave unit rows that are no longer dominated
    by it: the mean absolute cosine between the two BLOCKS should fall."""
    plain = CodepSeedAcquisition(profile_npz=npz, weight=0.5)
    cut = CodepSeedAcquisition(profile_npz=npz, weight=0.5, drop_pcs=1)
    i, j = [plain._row[f"G{k:02d}"] for k in range(4)], \
           [plain._row[f"G{k:02d}"] for k in range(4, 8)]
    before = np.abs(plain._P[i] @ plain._P[j].T).mean()
    after = np.abs(cut._P[i] @ cut._P[j].T).mean()
    assert np.allclose(np.linalg.norm(cut._P, axis=1), 1.0)   # still unit rows
    assert after != before
    assert "pc1" in cut.name()


def test_drop_pcs_cannot_exceed_the_rank(npz):
    """Removing as many components as there are dimensions leaves nothing; the loader
    must still return finite unit-normalised rows rather than nan."""
    acq = CodepSeedAcquisition(profile_npz=npz, weight=0.5, drop_pcs=LINES)
    assert np.isfinite(acq._P).all()


def test_amp_power_demotes_flat_profile_genes(npz, tmp_path):
    """A gene whose profile is nearly flat has a cosine that is noise. With amp_power
    on it must rank below a gene with the same cosine and a real amplitude."""
    rng = np.random.default_rng(1)
    a = rng.standard_normal(LINES)
    cols, vals = [], []
    for i in range(4):                                  # strong, covarying
        cols.append(GENES[i]); vals.append(a + 0.02 * rng.standard_normal(LINES))
    for i in range(4, 8):                               # same shape, 1/50th amplitude
        cols.append(GENES[i]); vals.append(0.02 * a)
    pth = tmp_path / "amp.npz"
    np.savez(pth, values=np.array(vals, dtype=np.float32).T,
             rows=np.array([f"L{i}" for i in range(LINES)], dtype=object),
             cols=np.array(cols, dtype=object))
    off = CodepSeedAcquisition(profile_npz=pth, weight=0.0, amp_power=0.0)
    on = CodepSeedAcquisition(profile_npz=pth, weight=0.0, amp_power=1.0)
    assert on._amp[on._row["G04"]] < 0.1 * on._amp[on._row["G00"]]
    # cosines are near-identical by construction; only the weight separates them
    hist = _hist(["G00", "G01"])
    q, _ = off._seed_query(hist)
    c_strong = float(off._P[off._row["G02"]] @ q)
    c_flat = float(off._P[off._row["G06"]] @ q)
    assert abs(c_strong - c_flat) < 0.2
    w_strong = c_strong * on._amp[on._row["G02"]]
    w_flat = c_flat * on._amp[on._row["G06"]]
    assert w_strong > w_flat
    assert "amp1" in on.name()


def test_geometry_knobs_appear_in_the_trace(npz):
    acq = CodepSeedAcquisition(profile_npz=npz, weight=0.5, drop_pcs=2, amp_power=0.5)
    acq.suggest(_hist(["G00", "G01"]), GENES, 4, _pred(GENES))
    t = acq.last_trace()
    assert t["drop_pcs"] == 2 and t["amp_power"] == 0.5


# -- the description renderer variants (§4.45) --------------------------------
#
# The default MUST be byte-identical to the shipped renderer, because the shipped
# embedding cache is keyed on the exact text and every existing number in the project
# was computed against it. A one-character drift silently invalidates the cache and
# forces fresh API calls that produce DIFFERENT vectors for the same screen.

def _screen_stub():
    class S:
        dataset_name = "X1"
        def context(self):
            return {"phenotype": "increases resistance to MERS-CoV",
                    "cell_line": "HEK-293", "organism": "Homo sapiens",
                    "library_methodology": "Knockout",
                    "notes": "Genes enriched upon KO confer resistance.",
                    "ranking_rationale": "low p-Value",
                    "screen_rationale": "Increased resistance",
                    "significance_criteria": "p < 0.05",
                    "library_type": "CRISPRn",
                    "direction_str": "", "description": "", "contrast_label": ""}
    return S()


def test_desc_extra_unset_is_byte_identical(monkeypatch):
    from assayloop.amortized import text_embed as te
    s = _screen_stub()
    monkeypatch.delenv("ASSAYLOOP_DESC_EXTRA", raising=False)
    base = te.screen_description_text(s)
    for off in ("", "none", "NONE", "  "):
        monkeypatch.setenv("ASSAYLOOP_DESC_EXTRA", off)
        assert te.screen_description_text(s) == base, f"{off!r} changed the text"
    assert "Notes:" not in base and "Selection direction:" not in base


def test_desc_extra_notes_appends_only_notes(monkeypatch):
    from assayloop.amortized import text_embed as te
    s = _screen_stub()
    monkeypatch.delenv("ASSAYLOOP_DESC_EXTRA", raising=False)
    base = te.screen_description_text(s)
    monkeypatch.setenv("ASSAYLOOP_DESC_EXTRA", "notes")
    got = te.screen_description_text(s)
    assert got.startswith(base)                      # purely additive
    assert "Notes: Genes enriched upon KO" in got
    assert "Selection direction:" not in got and "Ranking rationale:" not in got


def test_desc_extra_dir_is_derived_from_phenotype_and_notes(monkeypatch):
    from assayloop.amortized import text_embed as te
    s = _screen_stub()
    monkeypatch.setenv("ASSAYLOOP_DESC_EXTRA", "dir")
    got = te.screen_description_text(s)
    # `direction_str` is empty on every real screen, so the label must be classified
    assert "Selection direction: enrichment" in got
    assert "Notes:" not in got


def test_desc_extra_all_includes_every_group(monkeypatch):
    from assayloop.amortized import text_embed as te
    s = _screen_stub()
    monkeypatch.setenv("ASSAYLOOP_DESC_EXTRA", "all")
    got = te.screen_description_text(s)
    for frag in ("Notes:", "Ranking rationale:", "Screen rationale:",
                 "Significance criteria:", "Library type:", "Selection direction:"):
        assert frag in got, frag


def test_desc_extra_changes_the_cache_key(monkeypatch):
    """The whole reason new embeddings are needed: a different text is a different key,
    so no variant can silently reuse the shipped vectors."""
    from assayloop.amortized import text_embed as te
    s = _screen_stub()
    keys = {}
    for v in ("none", "notes", "dir", "all"):
        monkeypatch.setenv("ASSAYLOOP_DESC_EXTRA", v)
        keys[v] = te._text_key(te.screen_description_text(s),
                               model="text-embedding-3-small")
    assert len(set(keys.values())) == 4, keys


def test_desc_extra_rejects_unknown_values(monkeypatch):
    from assayloop.amortized import text_embed as te
    monkeypatch.setenv("ASSAYLOOP_DESC_EXTRA", "notes,bogus")
    with pytest.raises(ValueError, match="unknown part"):
        te.screen_description_text(_screen_stub())


# -- the oracle-seed control (§4.49) -----------------------------------------
#
# This is a BOUND, never a deployable setting: it hands the scorers genes the answer key
# says are hits. The two properties that make it a bound rather than a cheat are pinned
# here -- the seeds shape the direction, and they cannot be picked for credit.

def test_oracle_seeds_empty_is_bit_identical(npz, trainhits):
    a = CodepSeedAcquisition(profile_npz=npz, weight=0.5, train_hits=trainhits,
                             knn_weight=0.5)
    b = CodepSeedAcquisition(profile_npz=npz, weight=0.5, train_hits=trainhits,
                             knn_weight=0.5, oracle_seeds=[])
    pred, hist = _pred(GENES), _hist(["G00", "G01"])
    assert a.suggest(hist, GENES, 4, pred) == b.suggest(hist, GENES, 4, pred)
    assert "ORC" not in b.name()


def test_oracle_seeds_are_never_picked(npz):
    """The load-bearing property. A gene given away as a known hit must be removed from
    the candidate pool, or the arm banks it as a hit and the bound measures free hits."""
    acq = CodepSeedAcquisition(profile_npz=npz, weight=0.5,
                               oracle_seeds=["G00", "G01", "G02"])
    picked = acq.suggest(_hist(["G04", "G05"]), GENES, 6, _pred(GENES))
    assert not ({"G00", "G01", "G02"} & set(picked)), picked
    assert len(picked) == 6          # still fills the batch from what remains


def test_oracle_seeds_fire_the_kernel_with_no_revealed_hits(npz):
    """36% of fold screens never reach 2 revealed hits in rounds 1-3 (§4.48), so part of
    the bound is the kernel firing at all on those. With oracle seeds it must fire from
    an empty history."""
    acq = CodepSeedAcquisition(profile_npz=npz, weight=4.0,
                               oracle_seeds=["G04", "G05", "G06"])
    order = ["G00", "G01", "G02", "G03", "G08", "G09", "G10", "G11", "G07"]
    picked = acq.suggest([], order, 3, _pred(order))
    assert acq.last_trace()["codep_active"] is True
    # G07 covaries with the oracle seeds G04-G06 and is last in the greedy order
    assert "G07" in picked, picked


def test_oracle_seeds_steer_the_direction(npz):
    """The seeds must change the query, not merely the candidate list."""
    plain = CodepSeedAcquisition(profile_npz=npz, weight=0.5)
    orc = CodepSeedAcquisition(profile_npz=npz, weight=0.5,
                               oracle_seeds=["G04", "G05", "G06"])
    hist = _hist(["G00", "G01"])
    q_plain, _ = plain._seed_query(hist)
    q_orc, _ = orc._seed_query(hist)
    assert float(q_plain @ q_orc) < 0.99      # genuinely a different direction


def test_oracle_seeds_appear_in_the_name(npz):
    acq = CodepSeedAcquisition(profile_npz=npz, weight=0.5,
                               oracle_seeds=["G00", "G01"])
    assert "ORC2" in acq.name()


# -- in-campaign injection (§4.50) -------------------------------------------
#
# The better-posed version of §4.49. Seeds enter the campaign as REAL PICKS in the early
# rounds, so the budget is untouched and nothing is confiscated. The injected hits do
# credit the score, which is why the experiment is only readable as a 2x2 -- these tests
# pin the mechanics the 2x2 relies on.

def test_inject_off_is_bit_identical(npz, trainhits):
    a = CodepSeedAcquisition(profile_npz=npz, weight=0.5, train_hits=trainhits,
                             knn_weight=0.5)
    b = CodepSeedAcquisition(profile_npz=npz, weight=0.5, train_hits=trainhits,
                             knn_weight=0.5, inject_genes=["G08", "G09"],
                             inject_per_round=0)
    pred, hist = _pred(GENES), _hist(["G00", "G01"])
    assert a.suggest(hist, GENES, 4, pred) == b.suggest(hist, GENES, 4, pred)
    assert "INJ" not in b.name()


def test_inject_fills_the_lowest_slots_and_keeps_the_budget(npz):
    """The load-bearing property: the batch is still exactly k picks. If injection grew
    the batch, EF's denominator would change and the 2x2 would not cancel."""
    acq = CodepSeedAcquisition(profile_npz=npz, weight=0.0,
                               inject_genes=["G10", "G11"], inject_per_round=2)
    order = ["G00", "G01", "G02", "G03", "G04", "G05", "G10", "G11"]
    picked = acq.suggest([], order, 4, _pred(order))
    assert len(picked) == 4                       # budget untouched
    assert {"G10", "G11"} <= set(picked)          # injected genes are in
    assert picked[:2] == ["G00", "G01"]           # strongest picks kept
    assert "G02" not in picked and "G03" not in picked   # weakest slots displaced
    assert acq.last_trace()["n_injected"] == 2


def test_inject_stops_after_inject_rounds(npz):
    acq = CodepSeedAcquisition(profile_npz=npz, weight=0.0,
                               inject_genes=["G10", "G11"], inject_per_round=1,
                               inject_rounds=2)
    pred = _pred(GENES)
    acq.suggest([], GENES, 4, pred)
    assert acq.last_trace()["n_injected"] == 1
    acq.suggest(_hist(["G00"]), GENES, 4, pred)                    # round 2
    assert acq.last_trace()["n_injected"] == 1
    acq.suggest(_hist(["G00", "G01"]) * 2, GENES, 4, pred)         # round 3+
    assert acq.last_trace()["n_injected"] == 0


def test_inject_never_repeats_a_gene(npz):
    acq = CodepSeedAcquisition(profile_npz=npz, weight=0.0,
                               inject_genes=["G10", "G11"], inject_per_round=1,
                               inject_rounds=3)
    pred = _pred(GENES)
    first = acq.suggest([], GENES, 4, pred)
    second = acq.suggest(_hist(["G00"]), GENES, 4, pred)
    inj1 = {g for g in first if g in ("G10", "G11")}
    inj2 = {g for g in second if g in ("G10", "G11")}
    assert inj1 and inj2 and not (inj1 & inj2)
    acq.reset()
    assert acq.last_trace() == {}


def test_inject_is_independent_of_the_kernel_weight(npz, trainhits):
    """The 2x2 needs the SAME genes injected whether the kernel is on or off, or the
    free-hit credit does not cancel."""
    kw = dict(profile_npz=npz, train_hits=trainhits,
              inject_genes=["G10", "G11"], inject_per_round=2)
    off = CodepSeedAcquisition(weight=0.0, knn_weight=0.0, **kw)
    on = CodepSeedAcquisition(weight=0.5, knn_weight=0.5, **kw)
    a = off.suggest([], GENES, 4, _pred(GENES))
    b = on.suggest([], GENES, 4, _pred(GENES))
    assert {"G10", "G11"} <= set(a) and {"G10", "G11"} <= set(b)


# --- chromosomal-proximity exclusion ----------------------------------------
# Co-essentiality between neighbouring genes is inflated by the local copy-number
# state, so the literature drops proximal pairs before calling a network. These pin
# that the option is exact and, crucially, OFF by default.

def _hgnc(tmp_path, mapping, name="hgnc.txt"):
    """A minimal HGNC-shaped TSV: symbol and location are the only columns read."""
    p = tmp_path / name
    lines = ["symbol\tlocation"] + [f"{g}\t{loc}" for g, loc in mapping.items()]
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return p


def test_prox_is_off_without_a_locus_file(npz):
    """The default must not move a single pick."""
    base = CodepSeedAcquisition(profile_npz=npz, weight=1.0)
    assert base._locus is None
    assert base.suggest(_hist(["G00", "G01"]), GENES, 4, _pred(GENES)) == \
           CodepSeedAcquisition(profile_npz=npz, weight=1.0).suggest(
               _hist(["G00", "G01"]), GENES, 4, _pred(GENES))


def test_prox_no_shared_locus_leaves_the_ranking_alone(npz, tmp_path):
    """No CANDIDATE sharing a locus with a seed means nothing is excluded, so the
    order is unchanged. The seeds themselves are left out of the pool: a seed does
    share a locus with itself, and dropping its own profile from its own score is
    the correction working, not a regression."""
    far = _hgnc(tmp_path, {g: f"{i + 1}q11.1" for i, g in enumerate(GENES)})
    off = CodepSeedAcquisition(profile_npz=npz, weight=1.0)
    on = CodepSeedAcquisition(profile_npz=npz, weight=1.0,
                              gene_locus=far, prox_level="arm")
    assert (on._locus >= 0).all()
    cands = [g for g in GENES if g not in ("G00", "G01")]
    hist, pred = _hist(["G00", "G01"]), _pred(GENES)
    assert on.suggest(hist, cands, 5, pred) == off.suggest(hist, cands, 5, pred)


def test_prox_excludes_a_seed_sharing_the_candidate_locus(npz, tmp_path):
    """G02 sits with its seeds, so its score drops the seeds it is proximal to."""
    same = _hgnc(tmp_path, {"G00": "1q11.1", "G01": "1q11.1", "G02": "1q22.2",
                            **{g: f"{i + 5}q11.1" for i, g in enumerate(GENES[3:])}})
    acq = CodepSeedAcquisition(profile_npz=npz, weight=1.0,
                               gene_locus=same, prox_level="arm")
    hist = _hist(["G00", "G01"])
    scores, _n, n_corr = acq._codep_prox_scores(hist, GENES)
    # every gene on chromosome 1 is corrected: G00, G01 and G02
    assert n_corr == 3
    # with both seeds on its own arm, G02 has no seed left and scores nan
    assert np.isnan(scores[GENES.index("G02")])
    # a gene on another arm keeps the full, uncorrected mean
    P, row = acq._P, acq._row
    H = [row["G00"], row["G01"]]
    expect = float(P[row["G03"]] @ P[H].sum(axis=0)) / 2
    assert scores[GENES.index("G03")] == pytest.approx(expect)


def test_prox_level_band_is_finer_than_arm(npz, tmp_path):
    """Same arm, different band: excluded at arm level, kept at band level."""
    m = _hgnc(tmp_path, {"G00": "1q11.1", "G01": "1q11.1", "G02": "1q22.2",
                         **{g: f"{i + 5}q11.1" for i, g in enumerate(GENES[3:])}})
    hist = _hist(["G00", "G01"])
    arm = CodepSeedAcquisition(profile_npz=npz, weight=1.0, gene_locus=m, prox_level="arm")
    band = CodepSeedAcquisition(profile_npz=npz, weight=1.0, gene_locus=m, prox_level="band")
    assert np.isnan(arm._codep_prox_scores(hist, GENES)[0][GENES.index("G02")])
    assert np.isfinite(band._codep_prox_scores(hist, GENES)[0][GENES.index("G02")])


def test_unplaced_genes_are_never_excluded(npz, tmp_path):
    """'mitochondria' and 'unplaced' carry no arm, so they get locus -1 and no
    correction -- not a shared bucket that would exclude them from each other."""
    m = _hgnc(tmp_path, {"G00": "mitochondria", "G01": "unplaced",
                         "G02": "reserved", "G03": "1q11.1"})
    acq = CodepSeedAcquisition(profile_npz=npz, weight=1.0, gene_locus=m)
    assert acq._locus[acq._row["G00"]] == -1
    assert acq._locus[acq._row["G01"]] == -1
    assert acq._locus[acq._row["G03"]] >= 0


def test_prox_level_is_validated(npz, tmp_path):
    with pytest.raises(ValueError, match="prox_level"):
        CodepSeedAcquisition(profile_npz=npz, weight=1.0,
                             gene_locus=_hgnc(tmp_path, {"G00": "1q11.1"}),
                             prox_level="chromosome")


# --- the paper's kNN baseline, fused as a fourth ranker ---------------------
# Standalone it is weak (EF 1.93 vs AssayFormer's 4.83) but nearly orthogonal
# (0.040 pick Jaccard), so it is fused rather than substituted. These pin that it
# matches `KNNGeneEmbedding`'s estimator and stays inert by default.

@pytest.fixture
def gptnpz(tmp_path):
    """Unit-norm embeddings: G00-G03 cluster, G04-G07 cluster elsewhere."""
    rng = np.random.default_rng(1)
    a, b = rng.standard_normal(6), rng.standard_normal(6)
    rows, cols = [], []
    for i in range(4):
        rows.append(a + 0.05 * rng.standard_normal(6)); cols.append(GENES[i])
    for i in range(4, 8):
        rows.append(b + 0.05 * rng.standard_normal(6)); cols.append(GENES[i])
    X = np.array(rows, dtype=np.float32)
    X /= np.linalg.norm(X, axis=1, keepdims=True)
    p = tmp_path / "genept.npz"
    np.savez(p, values=X, cols=np.array(cols, dtype=object))
    return p


def _labelled_hist(labels):
    """labels: {gene: bool} revealed in one step."""
    return [StepRecord(step=1, acquired_batch=list(labels), new_observations=[
        Observation(candidate=g, label={"hit": h, "relevance_score": 1.0 if h else 0.0})
        for g, h in labels.items()])]


def test_gpt_is_inert_by_default(npz):
    acq = CodepSeedAcquisition(profile_npz=npz, weight=1.0)
    assert acq._G is None
    assert acq._gpt_scores(_hist(["G00", "G01"]), GENES) is None


def test_gpt_weight_zero_does_not_even_load(npz, gptnpz):
    """weight 0 must skip the file entirely -- the real one is 139 MB per screen."""
    acq = CodepSeedAcquisition(profile_npz=npz, weight=1.0,
                               gpt_npz=gptnpz, gpt_weight=0.0)
    assert acq._G is None


def test_gpt_matches_a_hand_computed_weighted_vote(npz, gptnpz):
    """The estimator is KNNGeneEmbedding's: 1/(1 - cos + eps) weights over the k
    nearest ALREADY-TESTED genes, labels in {0, 1}."""
    acq = CodepSeedAcquisition(profile_npz=npz, weight=1.0, gpt_npz=gptnpz,
                               gpt_weight=0.5, gpt_k=2)
    hist = _labelled_hist({"G00": True, "G04": False, "G05": False})
    got = acq._gpt_scores(hist, GENES)

    G, row = acq._G, acq._G_row
    tested = ["G00", "G04", "G05"]; y = np.array([1.0, 0.0, 0.0])
    g = "G01"
    sims = np.array([float(G[row[g]] @ G[row[t]]) for t in tested])
    top = np.argsort(-sims)[:2]
    w = 1.0 / (1.0 - sims[top] + 1e-6)
    want = float((w * y[top]).sum() / w.sum())
    assert got[GENES.index(g)] == pytest.approx(want, rel=1e-6)


def test_gpt_needs_no_hits_unlike_the_other_scorers(npz, gptnpz):
    """A round of pure misses is still information to a kNN, so its gate is one
    LABEL, not two confirmed hits -- the codep term is off here and it is not."""
    acq = CodepSeedAcquisition(profile_npz=npz, weight=1.0, gpt_npz=gptnpz,
                               gpt_weight=0.5)
    hist = _labelled_hist({"G00": False, "G04": False})
    assert acq._seed_query(hist, acq._P, acq._row)[0] is None      # codep: no seeds
    s = acq._gpt_scores(hist, GENES)
    assert s is not None and np.isfinite(s[:8]).all()              # gpt: scored anyway
    assert np.allclose(s[:8], 0.0)                                 # all labels are 0


def test_gpt_leaves_uncovered_genes_nan(npz, gptnpz):
    acq = CodepSeedAcquisition(profile_npz=npz, weight=1.0, gpt_npz=gptnpz,
                               gpt_weight=0.5)
    s = acq._gpt_scores(_labelled_hist({"G00": True, "G04": False}), GENES)
    assert np.isnan(s[GENES.index("G08")])     # G08-G11 are in no atlas
    assert np.isfinite(s[GENES.index("G01")])


def test_gpt_changes_the_picks_and_is_traced(npz, gptnpz):
    """The model must PREFER the other cluster, or both arms pick the same genes for
    reasons that have nothing to do with the ranker."""
    hist = _labelled_hist({"G00": True, "G04": False, "G05": False})
    pred = _pred(list(reversed(GENES)))          # the model wants G11 .. G08 first
    cands = [g for g in GENES if g != "G00"]     # G00 is already tested
    off = CodepSeedAcquisition(profile_npz=npz, weight=1.0)
    on = CodepSeedAcquisition(profile_npz=npz, weight=1.0, gpt_npz=gptnpz,
                              gpt_weight=5.0)
    got_off = off.suggest(hist, cands, 4, pred)
    got_on = on.suggest(hist, cands, 4, pred)
    assert on.last_trace()["gpt_active"] is True
    assert off.last_trace()["gpt_active"] is False
    assert got_off != got_on
    # G00 was the only hit and G01-G03 are its neighbours, so they get promoted
    assert set(got_on) & {"G01", "G02", "G03"}


# --- the paper's own Screen-kNN, fused as a ranker --------------------------
# `models/screen_knn.py` is the workbook's `Screen-kNN` baseline. Its own loop is
# 21k candidates x 1349 screens per round, so the ranker uses the vectorised form.
# These pin that the fast path IS the model, not an approximation of it.

class _FakeScreen:
    def __init__(self, name, genes, hits):
        self.dataset_name = name
        self.genes = list(genes)
        self.hits = list(hits)
        self.cleaned_phenotype = "x"


def _fake_train():
    return [
        _FakeScreen("t1", ["G00", "G01", "G02", "G03"], [True,  True,  False, False]),
        _FakeScreen("t2", ["G00", "G01", "G04", "G05"], [True,  False, True,  False]),
        _FakeScreen("t3", ["G02", "G03", "G04", "G06"], [False, True,  True,  True]),
        _FakeScreen("t4", ["G00", "G02", "G05", "G06"], [False, False, True,  True]),
    ]


@pytest.fixture
def patched_train(monkeypatch):
    import assayloop.tasks as T
    from assayloop.acquisitions import codep_seed_acq as C
    monkeypatch.setattr(T, "load_screens", lambda *a, **k: _fake_train())
    C._PKNN_CACHE.clear()
    yield
    C._PKNN_CACHE.clear()


def test_paper_screenknn_ranker_is_inert_by_default(npz, patched_train):
    acq = CodepSeedAcquisition(profile_npz=npz, weight=1.0)
    assert acq.pknn_weight == 0.0
    acq.suggest(_hist(["G00", "G01"]), GENES, 3, _pred(GENES))
    assert acq.last_trace()["pknn_active"] is False


def test_paper_screenknn_matches_their_model_exactly(npz, patched_train):
    """The whole point of the vectorised copy: same numbers, 800x the speed."""
    from assayloop.models.screen_knn import ScreenKNNModel
    cands = GENES[:7]
    hist = _labelled_hist({"G00": True, "G02": False, "G05": True})
    obs = [o for s in hist for o in s.new_observations]

    acq = CodepSeedAcquisition(profile_npz=npz, weight=0.0, pknn_weight=1.0)
    fast = acq._pknn_scores(hist, cands, {"dataset_name": "current"})

    slow = ScreenKNNModel().predict(obs, cands, {"dataset_name": "current"}).scores
    for i, g in enumerate(cands):
        assert fast[i] == pytest.approx(slow[g], rel=1e-6, abs=1e-9), g


def test_paper_screenknn_excludes_the_screen_under_test(npz, patched_train):
    """A training screen with the same name must not vote on itself."""
    cands = GENES[:7]
    hist = _labelled_hist({"G00": True, "G02": False})
    acq = CodepSeedAcquisition(profile_npz=npz, weight=0.0, pknn_weight=1.0)
    free = acq._pknn_scores(hist, cands, {"dataset_name": "nothing-matches"})
    held = acq._pknn_scores(hist, cands, {"dataset_name": "t1"})
    assert not np.allclose(free, held)


def test_paper_screenknn_runs_with_no_hits_at_all(npz, patched_train):
    """Its distance counts agreement on MISSES too, so a hitless round still scores
    -- unlike the DepMap term, which needs two confirmed hits to point a direction."""
    acq = CodepSeedAcquisition(profile_npz=npz, weight=1.0, pknn_weight=1.0)
    hist = _labelled_hist({"G00": False, "G02": False})
    assert acq._seed_query(hist, acq._P, acq._row)[0] is None
    s = acq._pknn_scores(hist, GENES[:7], None)
    assert s is not None and np.isfinite(s).all()


# --- the GenePT gene-kNN must BE the paper's, not merely resemble it --------

class _FakeEmb:
    """A GeneEmbeddingProvider over a fixed table, so both implementations see
    byte-identical embeddings and any difference is the estimator."""
    def __init__(self, table):
        self._t = table
    def name(self) -> str:
        return "fake"
    def embed_batch(self, genes):
        d = len(next(iter(self._t.values())))
        X = np.zeros((len(genes), d), dtype=np.float32)
        m = np.zeros(len(genes), dtype=bool)
        for i, g in enumerate(genes):
            v = self._t.get(str(g))
            if v is not None:
                X[i] = v; m[i] = True
        return X, m


def test_genept_knn_reproduces_the_paper_model(npz, gptnpz):
    """Pins `_gpt_scores` against `KNNGeneEmbedding.predict`, the class it
    reimplements. Without this the fused result is a measurement of *a* GenePT kNN,
    not of the paper's.

    Agreement is on the RANKING, which is all the fusion consumes, plus scores to
    about 1e-4. It is not bit-exact, and that is inherent to their estimator rather
    than sloppiness here: with epsilon = 1e-6, a gene that is its own nearest
    neighbour takes weight 1/(1 - 1 + eps) ~ 1e6, so a float-level difference in a
    similarity near 1.0 is amplified a millionfold in the vote. Their normalisation
    happens per call on float32; these embeddings are pre-normalised once. Same
    estimator, different rounding, same order.
    """
    from assayloop.models.knn_gene_embedding import KNNGeneEmbedding

    acq = CodepSeedAcquisition(profile_npz=npz, weight=0.0,
                               gpt_npz=gptnpz, gpt_weight=1.0, gpt_k=3)
    table = {g: acq._G[i] for g, i in acq._G_row.items()}
    cands = GENES[:8]
    # enough labels that a max_labels-style cap would show up as a difference
    hist = _labelled_hist({"G00": True, "G01": False, "G02": True,
                           "G04": False, "G05": True, "G06": False})
    obs = [o for s in hist for o in s.new_observations]

    mine = acq._gpt_scores(hist, cands)
    theirs = KNNGeneEmbedding(provider=_FakeEmb(table), k=3).predict(obs, cands).scores
    tv = np.array([theirs[g] for g in cands])
    # the ranking is the contract -- RRF never sees the scores themselves
    assert list(np.argsort(-mine)) == list(np.argsort(-tv))
    assert np.abs(mine - tv).max() < 2e-4


def test_genept_knn_uses_every_label_not_just_the_last_600(npz, gptnpz):
    """max_labels caps the DepMap scorer; applying it here would drop early rounds
    from the neighbour pool and silently change the estimator."""
    acq = CodepSeedAcquisition(profile_npz=npz, weight=0.0, gpt_npz=gptnpz,
                               gpt_weight=1.0, gpt_k=1, max_labels=2)
    # G00 is the only hit and it is the OLDEST label; with a cap of 2 it would be
    # dropped and every score would fall to 0
    hist = _labelled_hist({"G00": True, "G04": False, "G05": False, "G06": False})
    s = acq._gpt_scores(hist, ["G01"])
    assert s[0] > 0.0


# --- decomposing their estimator: similarity vs aggregation -----------------

def test_pknn_sim_is_validated(npz, patched_train):
    with pytest.raises(ValueError, match="pknn_sim"):
        CodepSeedAcquisition(profile_npz=npz, weight=0.0, pknn_weight=1.0,
                             pknn_sim="pearson")


def test_pknn_kappa_zero_removes_the_prior(npz, patched_train):
    """lambda = n/(n+kappa); kappa=0 gives lambda=1, so the score is pure neighbour
    evidence. This is the switch that separates 'the prior is the problem' from
    'the similarity is'."""
    cands = GENES[:7]
    hist = _labelled_hist({"G00": True, "G02": False, "G05": True})
    shrunk = CodepSeedAcquisition(profile_npz=npz, weight=0.0, pknn_weight=1.0)
    pure = CodepSeedAcquisition(profile_npz=npz, weight=0.0, pknn_weight=1.0,
                                pknn_kappa=0.0)
    a = shrunk._pknn_scores(hist, cands, None)
    b = pure._pknn_scores(hist, cands, None)
    assert not np.allclose(a, b)


def test_pknn_jaccard_is_one_minus_the_jaccard_distance(npz, patched_train):
    """Their squared-error distance, restricted to the union of hits, IS 1 - Jaccard
    on binary vectors. Pinned against the index computed directly from the sets."""
    acq = CodepSeedAcquisition(profile_npz=npz, weight=0.0, pknn_weight=1.0,
                               pknn_sim="jaccard", pknn_kappa=0.0)
    hist = _labelled_hist({"G00": True, "G01": True, "G03": False})
    got = acq._pknn_scores(hist, ["G02"], None)

    H = {"G00", "G01"}
    train = _fake_train()
    sims, hits_of = [], []
    for s in train:
        hs = {g for g, h in zip(s.genes, s.hits) if h}
        j = len(H & hs) / len(H | hs) if (H | hs) else 0.0
        sims.append(np.exp(-(1.0 - j) / 0.1))
        hits_of.append(1.0 if "G02" in hs else 0.0)
    sims = np.array(sims); hits_of = np.array(hits_of)
    # measured-weighted evidence, as the scorer computes it
    meas = np.array([1.0 if "G02" in s.genes else 0.0 for s in train])
    want = (sims * hits_of).sum() / (sims * meas).sum()
    assert got[0] == pytest.approx(want, rel=1e-5)


def test_pknn_cosine_matches_our_similarity_up_to_the_shared_factor(npz, patched_train):
    """`cosine` here divides by sqrt(|H|.|H_s|); our own scorer divides by sqrt(|H_s|)
    alone. The missing sqrt(|H|) is identical for every past screen, so the two rank
    screens the same way -- which is what makes this the like-for-like cell."""
    acq = CodepSeedAcquisition(profile_npz=npz, weight=0.0, pknn_weight=1.0,
                               pknn_sim="cosine", pknn_kappa=0.0)
    H = {"G00", "G01"}
    train = _fake_train()
    full, ours = [], []
    for s in train:
        hs = {g for g, h in zip(s.genes, s.hits) if h}
        inter = len(H & hs)
        full.append(inter / np.sqrt(max(len(hs), 1) * len(H)))
        ours.append(inter / np.sqrt(max(len(hs), 1)))
    assert list(np.argsort(-np.array(full))) == list(np.argsort(-np.array(ours)))


def test_pknn_k_truncates_to_the_top_neighbours(npz, patched_train):
    """k zeroes everything outside the top k by weight, leaving the survivors'
    weights untouched -- so it isolates NEIGHBOURHOOD SIZE from reweighting."""
    cands = GENES[:7]
    hist = _labelled_hist({"G00": True, "G02": False, "G05": True})
    full = CodepSeedAcquisition(profile_npz=npz, weight=0.0, pknn_weight=1.0,
                                pknn_kappa=0.0)
    cut = CodepSeedAcquisition(profile_npz=npz, weight=0.0, pknn_weight=1.0,
                               pknn_kappa=0.0, pknn_k=1)
    a = full._pknn_scores(hist, cands, None)
    b = cut._pknn_scores(hist, cands, None)
    assert not np.allclose(a, b, equal_nan=True)
    # With k=1 one screen votes, so a gene it MEASURED scores exactly its label.
    # A gene it never measured has no evidence at all and correctly falls back to
    # the global prior -- which is why 0/1 alone is not the right assertion.
    train = _fake_train()
    H = {"G00", "G05"}
    best = max(train, key=lambda s: len(H & {g for g, h in zip(s.genes, s.hits) if h}))
    lab = dict(zip(best.genes, best.hits))
    prior = {}
    for g in GENES:
        tot = sum(1 for s in train if g in s.genes)
        hit = sum(1 for s in train if dict(zip(s.genes, s.hits)).get(g))
        prior[g] = hit / tot if tot else 0.0
    for i, g in enumerate(cands):
        if not np.isfinite(b[i]):
            continue
        want = float(lab[g]) if g in lab else prior[g]
        assert b[i] == pytest.approx(want, abs=1e-9), g


def test_pknn_k_zero_and_k_above_the_pool_are_both_no_ops(npz, patched_train):
    cands = GENES[:7]
    hist = _labelled_hist({"G00": True, "G02": False})
    base = CodepSeedAcquisition(profile_npz=npz, weight=0.0, pknn_weight=1.0)
    for k in (0, 999):
        v = CodepSeedAcquisition(profile_npz=npz, weight=0.0, pknn_weight=1.0,
                                 pknn_k=k)
        assert np.allclose(base._pknn_scores(hist, cands, None),
                           v._pknn_scores(hist, cands, None), equal_nan=True)


def test_pknn_shrink_m_zero_is_a_no_op(npz, patched_train):
    cands = GENES[:7]
    hist = _labelled_hist({"G00": True, "G02": False})
    a = CodepSeedAcquisition(profile_npz=npz, weight=0.0, pknn_weight=1.0)
    b = CodepSeedAcquisition(profile_npz=npz, weight=0.0, pknn_weight=1.0,
                             pknn_shrink_m=0.0)
    assert np.allclose(a._pknn_scores(hist, cands, None),
                       b._pknn_scores(hist, cands, None), equal_nan=True)


def test_pknn_shrink_pulls_thin_coverage_toward_the_prior(npz, patched_train):
    """A gene few neighbours measured should move further than a well-covered one."""
    cands = GENES[:7]
    hist = _labelled_hist({"G00": True, "G02": False, "G05": True})
    plain = CodepSeedAcquisition(profile_npz=npz, weight=0.0, pknn_weight=1.0,
                                 pknn_kappa=0.0)
    shrunk = CodepSeedAcquisition(profile_npz=npz, weight=0.0, pknn_weight=1.0,
                                  pknn_kappa=0.0, pknn_shrink_m=5.0)
    a = plain._pknn_scores(hist, cands, None)
    b = shrunk._pknn_scores(hist, cands, None)
    train = _fake_train()
    cover = {g: sum(1 for s in train if g in s.genes) for g in cands}
    moved = {g: abs(b[i] - a[i]) for i, g in enumerate(cands) if np.isfinite(a[i])}
    thin = min(moved, key=lambda g: cover[g])
    thick = max(moved, key=lambda g: cover[g])
    assert moved[thin] >= moved[thick]


def test_pknn_shrink_sends_unmeasured_genes_to_the_prior(npz, patched_train):
    """With a pseudo-count, wM = 0 needs no special case -- it lands on the prior."""
    acq = CodepSeedAcquisition(profile_npz=npz, weight=0.0, pknn_weight=1.0,
                               pknn_kappa=0.0, pknn_shrink_m=5.0, pknn_k=1)
    hist = _labelled_hist({"G00": True})
    s = acq._pknn_scores(hist, GENES[:7], None)
    assert np.isfinite(s).any()


# --- the gate that makes the paper's Screen-kNN nested, like ours (§4.56(e)) ------

def _empty_hist():
    """Round 1: nothing acquired, nothing observed."""
    return []


def test_pknn_is_gated_by_default_so_round_one_is_plain_greedy(npz, patched_train):
    """THE NESTING GUARANTEE, extended to the paper's scorer.

    Every scorer of ours returns None below `min_seeds`, which is what makes round 1
    pure greedy and `weight=0` reproduce AssayFormer pick for pick. `_pknn_scores` had
    no such check, so with no observations it returned the global hit prior and
    reordered round 1 -- its arms were not nested and the round-1 invariant could not
    validate them.
    """
    pred = _pred(GENES)
    gated = CodepSeedAcquisition(profile_npz=npz, weight=0.5, pknn_weight=1.0)
    got = gated.suggest(_empty_hist(), GENES, 4, pred)
    assert gated.last_trace()["pknn_active"] is False
    assert got == GENES[:4], "round 1 must be the greedy order"


def test_pknn_gate_false_restores_the_ungated_behaviour(npz, patched_train):
    """The control arm has to be reproducible, or the gate's worth cannot be priced."""
    pred = _pred(GENES)
    ungated = CodepSeedAcquisition(profile_npz=npz, weight=0.5, pknn_weight=1.0,
                                   pknn_gate=False)
    got = ungated.suggest(_empty_hist(), GENES, 4, pred)
    assert ungated.last_trace()["pknn_active"] is True
    assert got != GENES[:4], ("ungated, the prior must reorder round 1 -- if this ever "
                              "passes trivially the fake training set has a flat prior")


def test_pknn_gate_opens_once_min_seeds_hits_exist(npz, patched_train):
    """Gated is not disabled: past the threshold the two agree exactly."""
    hist = _labelled_hist({"G00": True, "G02": False, "G05": True})
    pred = _pred(GENES)
    a = CodepSeedAcquisition(profile_npz=npz, weight=0.5, pknn_weight=1.0)
    b = CodepSeedAcquisition(profile_npz=npz, weight=0.5, pknn_weight=1.0,
                             pknn_gate=False)
    ga = a.suggest(hist, GENES, 4, pred)
    gb = b.suggest(hist, GENES, 4, pred)
    assert a.last_trace()["pknn_active"] is True
    assert ga == gb


def test_pknn_gate_respects_min_seeds(npz, patched_train):
    """One hit is below the default threshold of two, so the term stays shut."""
    pred = _pred(GENES)
    acq = CodepSeedAcquisition(profile_npz=npz, weight=0.5, pknn_weight=1.0)
    acq.suggest(_labelled_hist({"G00": True, "G02": False}), GENES, 4, pred)
    assert acq.last_trace()["pknn_active"] is False
    acq2 = CodepSeedAcquisition(profile_npz=npz, weight=0.5, pknn_weight=1.0,
                                min_seeds=1)
    acq2.suggest(_labelled_hist({"G00": True, "G02": False}), GENES, 4, pred)
    assert acq2.last_trace()["pknn_active"] is True


def test_pknn_gate_does_not_touch_the_scorer_itself(npz, patched_train):
    """The gate belongs to the FUSION HARNESS, not to the estimator. `_pknn_scores`
    must stay a faithful copy of `ScreenKNNModel`, which with no data correctly
    predicts the prior -- that behaviour is right for a standalone model and only
    wrong for a term fused against a nested baseline."""
    acq = CodepSeedAcquisition(profile_npz=npz, weight=0.5, pknn_weight=1.0)
    assert acq.pknn_gate is True
    s = acq._pknn_scores(_empty_hist(), GENES[:7], None)
    assert s is not None and np.isfinite(s).all()


def test_pknn_universe_is_the_same_whenever_the_gate_opens(npz, patched_train):
    """THE BUG THE GATE INTRODUCED, pinned.

    `candidates` is the UNTESTED pool and shrinks each round. `_pknn_universe` is
    pinned on first call, so it used to be safe only because the first call was always
    round 1. Gating delays that call, and the acquired genes differ per screen -- which
    gave every screen its own ~230 MB cache entry (the real sweep died at 7.4 GB) and
    built the global prior over a screen-specific subset. The universe must be
    reconstructed from candidates PLUS history, so it is identical either way.
    """
    hist1 = _labelled_hist({"G00": True, "G02": True})
    # round 1 equivalent: nothing acquired, the whole pool is on offer
    a = CodepSeedAcquisition(profile_npz=npz, weight=0.0, pknn_weight=1.0)
    a._pknn_scores(hist1, GENES, None)
    early = a._pknn_universe

    # a later round: the same two genes are now OFF the candidate list
    later_cands = [g for g in GENES if g not in ("G00", "G02")]
    b = CodepSeedAcquisition(profile_npz=npz, weight=0.0, pknn_weight=1.0)
    b._pknn_scores(hist1, later_cands, None)
    late = b._pknn_universe

    assert early == late == tuple(sorted(GENES)), (
        "the Screen-kNN universe must not depend on which round it first ran in")


def test_pknn_state_cache_stays_at_one_entry_across_rounds(npz, patched_train):
    """The consequence that actually killed the sweep: one entry, not one per screen.

    The reconstruction rests on one property of the loop -- a gene leaves `candidates`
    only by being ACQUIRED, so it is in `history` -- so each case here keeps the
    acquired set and the dropped set equal, which is what the real loop guarantees.
    """
    from assayloop.acquisitions import codep_seed_acq as C
    rounds = [{},
              {"G00": True, "G02": True},
              {"G00": True, "G02": True, "G05": False, "G06": False},
              {"G00": True, "G02": True, "G05": False, "G06": False, "G07": True}]
    seen_universes = set()
    for labels in rounds:
        acq = CodepSeedAcquisition(profile_npz=npz, weight=0.0, pknn_weight=1.0)
        cands = [g for g in GENES if g not in labels]
        acq._pknn_scores(_labelled_hist(labels), cands, None)
        seen_universes.add(acq._pknn_universe)
    assert seen_universes == {tuple(sorted(GENES))}, (
        f"the universe varied across rounds: {sorted(len(u) for u in seen_universes)}")
    assert len(C._PKNN_CACHE) == 1, (
        f"expected one cached state, got {len(C._PKNN_CACHE)} -- the universe is "
        "varying between screens again")
