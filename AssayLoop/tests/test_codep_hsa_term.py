"""The seventh fusion term: a caller-supplied per-candidate scorer.

The term exists so an external 42-D historical-subspace model can enter the ensemble
without its bank, encoder or vendored numerics becoming a dependency of this package.
These tests pin the contract rather than any particular scorer: nesting at weight 0,
the once-per-round call discipline, the shape/finiteness guards, and that at
``af_weight=0`` with only the HSA term the ensemble ranks by that scorer alone.
"""
from types import SimpleNamespace

import numpy as np
import pytest

from assayloop.acquisitions.codep_seed_acq import CodepSeedAcquisition


def _profile(tmp_path, genes):
    """Minimal co-essentiality profile npz: rows = genes, cols = cell lines."""
    rng = np.random.default_rng(0)
    p = tmp_path / "atlas.npz"
    np.savez(p, values=rng.normal(size=(len(genes), 8)).astype(np.float32),
             rows=np.array(genes, dtype=object), cols=np.array([f"L{i}" for i in range(8)],
                                                               dtype=object))
    return str(p)


GENES = [f"G{i:03d}" for i in range(40)]


def _pred(cands):
    return SimpleNamespace(scores={c: float(len(cands) - i) for i, c in enumerate(cands)})


def test_hsa_weight_without_scorer_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="requires an hsa_scorer"):
        CodepSeedAcquisition(profile_npz=_profile(tmp_path, GENES), weight=0.5,
                             hsa_weight=0.5)


def test_defaults_do_not_call_the_scorer(tmp_path):
    """hsa_weight=0 must leave the scorer untouched, so existing arms are unchanged."""
    calls = []
    acq = CodepSeedAcquisition(profile_npz=_profile(tmp_path, GENES), weight=0.5,
                               hsa_scorer=lambda *a: calls.append(1) or {})
    acq.suggest(history=[], candidates=list(GENES), batch_size=5,
                model_prediction=_pred(GENES))
    assert calls == []
    assert acq.af_weight == 1.0 and acq.hsa_weight == 0.0


def test_scorer_called_once_per_round(tmp_path):
    """A stateful scorer may count rounds, so it must be invoked exactly once each."""
    steps = []

    def scorer(observations, candidates, task_context):
        steps.append(len(observations))
        return {c: 0.0 for c in candidates}

    acq = CodepSeedAcquisition(profile_npz=_profile(tmp_path, GENES), weight=0.0,
                               hsa_weight=1.0, af_weight=0.0, hsa_scorer=scorer)
    hist = []
    for _ in range(3):
        acq.suggest(history=list(hist), candidates=list(GENES), batch_size=5,
                    model_prediction=_pred(GENES))
        hist.append(SimpleNamespace(new_observations=[
            SimpleNamespace(candidate=GENES[0], label={"hit": False},
                            metadata={"in_library": True})]))
    assert len(steps) == 3, steps
    assert steps == [0, 1, 2], "observations must arrive flattened and cumulative"


def test_hsa_alone_drives_the_ranking(tmp_path):
    """af_weight=0, weight=0, hsa_weight>0: the batch is the scorer's own top-k."""
    want = list(reversed(GENES))
    acq = CodepSeedAcquisition(
        profile_npz=_profile(tmp_path, GENES), weight=0.0, hsa_weight=1.0,
        af_weight=0.0,
        hsa_scorer=lambda o, c, t: {g: float(len(want) - want.index(g)) for g in c})
    picked = acq.suggest(history=[], candidates=list(GENES), batch_size=6,
                         model_prediction=_pred(GENES))
    assert picked == want[:6]
    assert acq._last_trace["hsa_active"] is True


def test_scorer_shape_and_finiteness_are_enforced(tmp_path):
    atlas = _profile(tmp_path, GENES)
    bad_shape = CodepSeedAcquisition(profile_npz=atlas, weight=0.0, hsa_weight=1.0,
                                     af_weight=0.0,
                                     hsa_scorer=lambda o, c, t: np.zeros(len(c) - 1))
    with pytest.raises(ValueError, match="expected"):
        bad_shape.suggest(history=[], candidates=list(GENES), batch_size=5,
                          model_prediction=_pred(GENES))
    nan = CodepSeedAcquisition(profile_npz=atlas, weight=0.0, hsa_weight=1.0,
                               af_weight=0.0,
                               hsa_scorer=lambda o, c, t: np.full(len(c), np.nan))
    with pytest.raises(ValueError, match="non-finite"):
        nan.suggest(history=[], candidates=list(GENES), batch_size=5,
                    model_prediction=_pred(GENES))


def test_name_records_the_new_weights(tmp_path):
    acq = CodepSeedAcquisition(profile_npz=_profile(tmp_path, GENES), weight=0.5,
                               hsa_weight=0.25, af_weight=0.0,
                               hsa_scorer=lambda *a: {})
    n = acq.name()
    assert "hsa0.25" in n and "af0" in n, n
