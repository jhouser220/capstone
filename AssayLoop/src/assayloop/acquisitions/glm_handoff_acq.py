"""GLM -> transformer handoff acquisition.

Replays GLM-5.1's recorded per-round acquisitions for the first ``n`` AL rounds,
then delegates to a base (greedy-from-model) acquisition for the remaining rounds.

The inner loop reveals the GLM genes' *true* labels (from the screen) and feeds
them to the amortized-ranker model as context, so the transformer continues
exactly from GLM's warm start. Because everything runs through the standard
:class:`~assaybench.core.SequentialLoop`, the resulting per-step metrics
(``hits_auc``, ``n_hits_vs_random`` with the requested-budget denominator,
batch hits, ...) match GLM's own dashboard runs and ``eval-ranker``.

``retain`` switches on the *proposal-capped* protocol used by the teammate's HSA
package (``output/runs/external/teammate/README.md``, "Proposal handling"):

    Original AssayLoop retains valid language proposals and fills short batches
    using pre-batch AssayFormer. HSA retains up to 25 ranked language proposals
    and numerically fills the remainder.

With ``retain=k`` a warm round keeps at most ``k`` of the recorded proposals --
ranked by the model's own score within the proposal pool, not by the LLM's
ordering -- and the remaining ``batch_size - k`` slots go to ``base``. So the
numerical method is live from round 1, and the batch is never short, which the
uncapped replay cannot promise (the recorded rounds supply ~97 of 100). The LLM's
share of a 10x100 campaign falls from ~29% to 7.5%. ``retain=None`` (the default)
keeps the uncapped behaviour: the whole recorded round, shortfall included.
"""

from __future__ import annotations

from typing import Any

from assaybench.core.acquisition import AcquisitionFunction
from assaybench.core.types import ModelPrediction, StepRecord


class GlmHandoffAcquisition(AcquisitionFunction):
    """First ``n_warm`` rounds = GLM's recorded picks; the rest = ``base``.

    Args:
        rounds: GLM's per-round acquired gene symbols (round 1 first). Symbols
            not in the live candidate pool are simply dropped by the inner loop
            (counted as shortfall), matching GLM's own under-supplied rounds.
        n_warm: number of leading rounds to take from ``rounds``.
        base: acquisition used once GLM has handed off (typically greedy).
        retain: cap on proposals kept per warm round; ``base`` fills the rest of
            that same round. ``None`` keeps every proposal and fills nothing.
        rank: how the proposal pool is ordered before the cap is applied --
            ``"model"`` by ``model_prediction.scores`` (the teammate's
            source/MAP ranking), ``"llm"`` by the recorded reply order.
    """

    def __init__(
        self,
        rounds: list[list[str]],
        n_warm: int,
        base: AcquisitionFunction,
        retain: int | None = None,
        rank: str = "model",
    ):
        self.rounds = [list(r) for r in (rounds or [])]
        self.n_warm = max(0, int(n_warm))
        self.base = base
        self.retain = None if retain is None else max(0, int(retain))
        if rank not in ("model", "llm"):
            raise ValueError(f"rank must be 'model' or 'llm', got {rank!r}")
        self.rank = rank
        self._step = 0
        self._last_trace: dict[str, Any] = {}

    def name(self) -> str:
        cap = "" if self.retain is None else f"_retain{self.retain}{self.rank[0]}"
        return f"glm_handoff_n{self.n_warm}{cap}+{self.base.name()}"

    def reset(self) -> None:
        self._step = 0
        self._last_trace = {}
        try:
            self.base.reset()
        except NotImplementedError:
            pass

    def last_trace(self) -> dict[str, Any]:
        return self._last_trace

    def suggest(
        self,
        history: list[StepRecord],
        candidates: list[Any],
        batch_size: int,
        model_prediction: ModelPrediction | None = None,
        task_context: dict[str, Any] | None = None,
    ) -> list[Any]:
        self._step += 1
        # GLM phase: replay the recorded round (intersected with the live pool).
        if self._step <= self.n_warm and self._step <= len(self.rounds):
            recorded = self.rounds[self._step - 1]
            cand = set(candidates)
            picked = [g for g in recorded if g in cand]
            if self.retain is None:                 # uncapped: the whole round, shortfall kept
                self._last_trace = {
                    "source": "glm_warm_start",
                    "round": self._step,
                    "n_picked": len(picked),
                    "n_recorded": len(recorded),
                }
                return picked[:batch_size]

            # Capped: rank the proposal pool, keep <= retain, fill the rest with `base`.
            ranked_by = self.rank
            if self.rank == "model":
                scores = getattr(model_prediction, "scores", None) or {}
                scored = [g for g in picked if g in scores]
                if scored:                          # unscored proposals keep reply order, after
                    picked = (sorted(scored, key=lambda g: -float(scores[g]))
                              + [g for g in picked if g not in set(scored)])
                else:
                    ranked_by = "llm_fallback_no_scores"
            kept = picked[: min(self.retain, batch_size)]
            need = batch_size - len(kept)
            keptset = set(kept)
            filled: list[Any] = []
            if need > 0:
                filled = [
                    g for g in self.base.suggest(
                        history=history,
                        candidates=[c for c in candidates if c not in keptset],
                        batch_size=need,
                        model_prediction=model_prediction,
                        task_context=task_context,
                    ) if g not in keptset
                ][:need]
            self._last_trace = {
                "source": "glm_warm_start_capped",
                "round": self._step,
                "n_recorded": len(recorded),
                "n_in_pool": len(picked),
                "n_retained": len(kept),
                "n_filled": len(filled),
                "retain": self.retain,
                "ranked_by": ranked_by,
                "fill": dict(self.base.last_trace() or {}),
            }
            return kept + filled
        # Transformer phase: delegate to the base (greedy) acquisition.
        picked = self.base.suggest(
            history=history,
            candidates=candidates,
            batch_size=batch_size,
            model_prediction=model_prediction,
            task_context=task_context,
        )
        self._last_trace = dict(self.base.last_trace() or {})
        self._last_trace["source"] = "transformer"
        return picked


__all__ = ["GlmHandoffAcquisition"]
