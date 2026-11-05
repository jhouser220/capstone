"""`GlmHandoffAcquisition(retain=k)`: at most k proposals, the rest from `base`.

The capped protocol is the teammate's ("retain up to 25 ranked language proposals
and numerically fill the remainder"), so the three things worth pinning are that
the cap is respected, that the batch comes out FULL even when the recorded round
is short, and that the kept proposals are the model's top-k of the pool rather
than the LLM's own first k.
"""
from __future__ import annotations

import unittest
from typing import Any

from assaybench.core.types import ModelPrediction
from assayloop.acquisitions.glm_handoff_acq import GlmHandoffAcquisition


class FakeBase:
    """Deterministic stand-in: hands back the candidates in order."""

    def __init__(self) -> None:
        self.calls: list[dict[str, Any]] = []

    def name(self) -> str:
        return "fake"

    def reset(self) -> None:
        self.calls.clear()

    def last_trace(self) -> dict[str, Any]:
        return {"fake": True}

    def suggest(self, history, candidates, batch_size, model_prediction=None,
                task_context=None):
        self.calls.append({"n_cand": len(candidates), "batch_size": batch_size})
        return list(candidates)[:batch_size]


POOL = [f"G{i:03d}" for i in range(50)]
# the recorded round, in LLM reply order: 8 genes, deliberately NOT score order
ROUND1 = ["G010", "G011", "G012", "G013", "G014", "G015", "G016", "G017"]
# model likes the LLM's LAST proposals best
SCORES = {g: float(i) for i, g in enumerate(POOL)}
PRED = ModelPrediction(scores=SCORES)


class TestRetain(unittest.TestCase):
    def test_uncapped_is_unchanged(self):
        acq = GlmHandoffAcquisition(rounds=[ROUND1], n_warm=1, base=FakeBase())
        got = acq.suggest([], POOL, 10, PRED, None)
        self.assertEqual(got, ROUND1)                     # short round stays short
        self.assertEqual(acq.last_trace()["source"], "glm_warm_start")

    def test_cap_then_fill_to_full_batch(self):
        base = FakeBase()
        acq = GlmHandoffAcquisition(rounds=[ROUND1], n_warm=1, base=base, retain=3)
        got = acq.suggest([], POOL, 10, PRED, None)
        tr = acq.last_trace()
        self.assertEqual(len(got), 10)                    # never short
        self.assertEqual(len(set(got)), 10)               # no duplicates
        self.assertEqual(tr["n_retained"], 3)
        self.assertEqual(tr["n_filled"], 7)
        # model-ranked: the three highest-scoring proposals, not ROUND1[:3]
        self.assertEqual(got[:3], ["G017", "G016", "G015"])
        self.assertEqual(base.calls[0]["n_cand"], len(POOL) - 3)
        self.assertEqual(base.calls[0]["batch_size"], 7)

    def test_llm_order_ranking(self):
        acq = GlmHandoffAcquisition(rounds=[ROUND1], n_warm=1, base=FakeBase(),
                                    retain=3, rank="llm")
        got = acq.suggest([], POOL, 10, PRED, None)
        self.assertEqual(got[:3], ROUND1[:3])
        self.assertEqual(acq.last_trace()["ranked_by"], "llm")

    def test_fewer_proposals_than_cap(self):
        acq = GlmHandoffAcquisition(rounds=[["G010", "G011"]], n_warm=1,
                                    base=FakeBase(), retain=25)
        got = acq.suggest([], POOL, 10, PRED, None)
        self.assertEqual(acq.last_trace()["n_retained"], 2)
        self.assertEqual(len(got), 10)

    def test_proposals_outside_the_pool_are_dropped(self):
        acq = GlmHandoffAcquisition(rounds=[["NOPE1", "G010", "NOPE2"]], n_warm=1,
                                    base=FakeBase(), retain=25)
        acq.suggest([], POOL, 10, PRED, None)
        tr = acq.last_trace()
        self.assertEqual((tr["n_recorded"], tr["n_in_pool"]), (3, 1))

    def test_handoff_still_happens_after_n_warm(self):
        base = FakeBase()
        acq = GlmHandoffAcquisition(rounds=[ROUND1], n_warm=1, base=base, retain=3)
        acq.suggest([], POOL, 10, PRED, None)
        got = acq.suggest([], POOL, 10, PRED, None)       # round 2: pure base
        self.assertEqual(got, POOL[:10])
        self.assertEqual(acq.last_trace()["source"], "transformer")

    def test_no_scores_falls_back_to_reply_order(self):
        acq = GlmHandoffAcquisition(rounds=[ROUND1], n_warm=1, base=FakeBase(), retain=3)
        got = acq.suggest([], POOL, 10, None, None)
        self.assertEqual(got[:3], ROUND1[:3])
        self.assertEqual(acq.last_trace()["ranked_by"], "llm_fallback_no_scores")


if __name__ == "__main__":
    unittest.main()
