"""Test-time context interventions on the amortized ranker.

Two levers, both free of retraining and of API calls:
  * ``ctx_ensemble`` -- bag the observed-gene set to get several screen latents
    and fuse their rankings.
  * ``ctx_weak_hits`` -- act on revealed hits whose ``relevance_score`` puts
    them in the weak tail.

The contract that matters most is the first test: the defaults must reproduce
the released model exactly, so every comparison has a nested baseline.
"""
from __future__ import annotations

import numpy as np
import pytest

torch = pytest.importorskip("torch")

from assaybench.core.types import Observation  # noqa: E402

from assayloop.amortized.data import GeneVocab  # noqa: E402
from assayloop.amortized.model import RankerConfig, RankerNet  # noqa: E402
from assayloop.models.amortized_ranker import AmortizedRankerModel  # noqa: E402

GENES = [f"G{i:03d}" for i in range(60)]
TEXT_DIM = 16


def _model(**kw) -> AmortizedRankerModel:
    vocab = GeneVocab(GENES)
    torch.manual_seed(0)
    net = RankerNet(RankerConfig(
        vocab_size=len(vocab), text_dim=TEXT_DIM, d_model=16, d_gene=10,
        d_hit=8, nhead=2, num_layers=1, dim_feedforward=32, dropout=0.0,
    ))
    m = AmortizedRankerModel(net=net, vocab=vocab, device="cpu", **kw)
    # a fixed description embedding, so no text backend is touched
    m._desc_emb = lambda ctx: np.full(TEXT_DIM, 0.1, dtype=np.float32)
    return m


def _obs(n: int = 20, n_hit: int = 8) -> list[Observation]:
    """n observations; the hits carry increasing relevance 1..n_hit."""
    out = []
    for i in range(n):
        hit = i < n_hit
        out.append(Observation(
            candidate=GENES[i],
            label={"hit": hit, "relevance_score": float(i + 1) if hit else 0.0},
        ))
    return out


def _scores(m, obs, cands=None):
    return m.predict(obs, cands or GENES[30:], {"dataset_name": "S"}).scores


# -- the nested baseline ------------------------------------------------------

def test_defaults_reproduce_the_released_model_exactly():
    """K=1 and ctx_weak_hits='keep' must be bit-identical to the single-variant
    path, or no comparison against the published AssayFormer number is valid."""
    base, same = _model(), _model(ctx_ensemble=1, ctx_weak_hits="keep")
    obs = _obs()
    a, b = _scores(base, obs), _scores(same, obs)
    assert a.keys() == b.keys()
    for g in a:
        assert a[g] == b[g]


def test_subsample_of_one_is_a_no_op_even_with_a_big_ensemble():
    """Every variant is then the full context, so fusion cannot move anything.

    Agreement is to float32 precision rather than bitwise: encoding a batch of 8
    identical contexts reduces in a different order than a batch of 1. The
    defaults path above keeps the batch at 1 and so stays bit-identical.
    """
    a = _scores(_model(), _obs())
    b = _scores(_model(ctx_ensemble=8, ctx_subsample=1.0), _obs())
    for g in a:
        assert a[g] == pytest.approx(b[g], rel=1e-4, abs=1e-8)


# -- ensembling ---------------------------------------------------------------

def test_ensemble_keeps_the_full_context_as_variant_zero():
    m = _model(ctx_ensemble=4, ctx_subsample=0.5)
    ids, hits = list(range(10)), [1] * 10
    variants = m._context_variants(ids, hits)
    assert len(variants) == 4
    assert variants[0] == (ids, hits)
    assert all(len(v[0]) == 5 for v in variants[1:])


def test_ensemble_is_deterministic_across_instances():
    obs = _obs()
    a = _scores(_model(ctx_ensemble=6, ctx_subsample=0.7), obs)
    b = _scores(_model(ctx_ensemble=6, ctx_subsample=0.7), obs)
    assert a == b


def test_ensemble_changes_the_ranking_but_stays_finite():
    obs = _obs()
    a = _scores(_model(), obs)
    b = _scores(_model(ctx_ensemble=8, ctx_subsample=0.5), obs)
    assert all(np.isfinite(v) for v in b.values())
    assert a != b


def test_rrf_fusion_ranks_rather_than_averages():
    obs = _obs()
    mean = _scores(_model(ctx_ensemble=8, ctx_subsample=0.5, ctx_fuse="mean"), obs)
    rrf = _scores(_model(ctx_ensemble=8, ctx_subsample=0.5, ctx_fuse="rrf"), obs)
    # RRF scores live on the 1/(60+rank) scale, so they are small and positive
    assert all(0 < v <= 8 / 60 for v in rrf.values())
    assert set(mean) == set(rrf)


def test_empty_context_survives_an_ensemble():
    s = _scores(_model(ctx_ensemble=4, ctx_subsample=0.5), [])
    assert len(s) == 30 and all(np.isfinite(v) for v in s.values())


# -- denoising ----------------------------------------------------------------

def test_drop_removes_only_the_weak_quartile_of_hits():
    m = _model(ctx_weak_hits="drop", ctx_weak_q=0.25)
    ids = list(range(20))
    hits = [1] * 8 + [0] * 12
    rels = [float(i + 1) for i in range(8)] + [0.0] * 12
    out_ids, out_hits, touched = m._denoise(ids, hits, rels)
    assert touched == 2                      # relevance 1 and 2 are <= quantile
    assert len(out_ids) == 18
    assert sum(out_hits) == 6                # 6 hits survive
    assert out_hits[6:] == [0] * 12          # the 12 non-hits are untouched


def test_negative_keeps_the_gene_and_flips_the_label():
    m = _model(ctx_weak_hits="negative", ctx_weak_q=0.25)
    ids = list(range(20))
    hits = [1] * 8 + [0] * 12
    rels = [float(i + 1) for i in range(8)] + [0.0] * 12
    out_ids, out_hits, touched = m._denoise(ids, hits, rels)
    assert touched == 2
    assert out_ids == ids                    # nothing removed
    assert sum(out_hits) == 6


def test_denoise_is_skipped_when_there_are_too_few_hits():
    """A quantile over 3 points is noise; leave the context alone."""
    m = _model(ctx_weak_hits="drop")
    ids, hits, rels = [0, 1, 2, 3], [1, 1, 1, 0], [1.0, 2.0, 3.0, 0.0]
    assert m._denoise(ids, hits, rels) == (ids, hits, 0)


def test_denoise_ignores_nan_relevance():
    """Out-of-library picks are revealed with relevance NaN; they must not
    silently become the weakest hits."""
    m = _model(ctx_weak_hits="drop", ctx_weak_q=0.25)
    ids = list(range(8))
    hits = [1] * 8
    rels = [float("nan"), 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0]
    out_ids, out_hits, touched = m._denoise(ids, hits, rels)
    assert 0 in out_ids                      # the NaN hit is kept, not dropped
    assert touched >= 1

def test_denoise_reports_through_metadata():
    m = _model(ctx_weak_hits="drop")
    md = m.predict(_obs(), GENES[30:], {"dataset_name": "S"}).metadata
    assert md["ctx_weak_hits"] == "drop" and md["ctx_weak_dropped"] == 2


# -- argument validation ------------------------------------------------------

@pytest.mark.parametrize("kw", [
    {"ctx_fuse": "softmax"},
    {"ctx_weak_hits": "delete"},
    {"ctx_subsample": 0.0},
    {"ctx_subsample": 1.5},
])
def test_bad_arguments_fail_loudly(kw):
    with pytest.raises(ValueError):
        _model(**kw)
