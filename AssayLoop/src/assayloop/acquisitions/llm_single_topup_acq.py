"""``llm_single_topup``: ``llm_single`` plus a top-up loop so each round really delivers ``batch_size`` genes.

The default ``llm_single`` makes ONE call per round, keeps the symbols that are valid and not yet tested, and
accepts whatever number survives. A round can come back nearly empty (for example a reply that mostly repeats
genes it was already shown), and the shortfall is simply lost.

This variant keeps the default behaviour for the first attempt (same system prompt, same user prompt, same
filtering) and only acts when that attempt falls short. It then re-asks, up to ``topup_max_attempts`` times,
with the same prompt asking for a few more genes than are still needed, plus a short note listing

  * the genes already accepted for THIS round, so the model does not repeat them, and
  * the symbols that were rejected (already tested, or not valid symbols).

No hit/non-hit labels are revealed during a round: the model learns nothing about outcomes between attempts,
only which symbols were taken, so the information available to it matches the default method. If the target
is still not reached after the allowed attempts (or attempts stop yielding new genes) the round is left short;
nothing is ever padded with genes the model did not propose.

Settings (constructor arguments, or environment variables so they can be set from the command line):

    ASSAYLOOP_TOPUP_MAX_ATTEMPTS   extra calls allowed per round after the first        (default 4)
    ASSAYLOOP_TOPUP_OVERSAMPLE     ask for ceil(still_needed * this) genes per top-up   (default 1.5)

Compatibility: the trace deliberately has NO ``response_text`` (a round can have several replies, and the
paper-style handoff replay would otherwise try to parse a truncated copy). The final accepted batch is
``acquired_batch``, which both handoff paths replay exactly. Every reply is still recorded in llm_calls.jsonl,
and ``trace["topup"]`` carries per-attempt counts.
"""
from __future__ import annotations

import logging
import math
import os
from typing import Any

from assaybench.core.types import ModelPrediction, StepRecord

from ..llm.client import complete_ex
from ..llm.parse_genes import extract_gene_list
from .llm_single_acq import LLMSingleAcquisition, _llm_trace_dict

log = logging.getLogger("assayloop.acquisitions.llm_single_topup")

_FINAL_MARKER = "=== Final instruction ==="
_STALL_LIMIT = 2          # stop after this many consecutive top-ups that add nothing
_MAX_LISTED_REJECTED = 100


class LLMSingleTopUpAcquisition(LLMSingleAcquisition):
    # why a proposed symbol was not accepted, as stated in the top-up note (subclasses that filter more widen it)
    _rejected_reason = "already tested, or not valid symbols"

    def __init__(
        self,
        *args: Any,
        topup_max_attempts: int | None = None,
        topup_oversample: float | None = None,
        **kwargs: Any,
    ):
        super().__init__(*args, **kwargs)
        if topup_max_attempts is None:
            topup_max_attempts = int(os.getenv("ASSAYLOOP_TOPUP_MAX_ATTEMPTS", "4"))
        if topup_oversample is None:
            topup_oversample = float(os.getenv("ASSAYLOOP_TOPUP_OVERSAMPLE", "1.5"))
        self.topup_max_attempts = max(0, int(topup_max_attempts))
        self.topup_oversample = max(1.0, float(topup_oversample))

    def name(self) -> str:
        cfg = self.llm.resolved()
        suffix = "" if self.reveal_labels else "_blind"
        return f"llm_single_topup{suffix}/{cfg.provider}:{cfg.model}"

    # ------------------------------------------------------------------ prompt for the top-up attempts

    def _build_topup_prompt(
        self,
        history: list[StepRecord],
        task_context: dict[str, Any] | None,
        ask: int,
        accepted: list[Any],
        rejected: list[str],
    ) -> tuple[str, str, int]:
        """The default prompt, asking for ``ask`` genes, with a short exclusion note before the final instruction."""
        system, user, n_requested = self._build_prompt(history, ask, task_context)
        lines = [
            "## Top-up request",
            f"Your previous reply for this round could not be fully used. {len(accepted)} new gene"
            f"{'s' if len(accepted) != 1 else ''} from this round "
            f"{'have' if len(accepted) != 1 else 'has'} already been accepted (their outcomes are not available yet).",
        ]
        if accepted:
            lines.append("Do NOT suggest any of these again: " + ", ".join(str(g) for g in accepted))
        if rejected:
            lines.append(
                f"These symbols were rejected ({self._rejected_reason}): "
                + ", ".join(rejected[:_MAX_LISTED_REJECTED])
            )
        lines.append(
            f"Provide {n_requested} additional, different genes that are hits in this screen, "
            "ranked from strongest to weakest, using the same output format."
        )
        addendum = "\n".join(lines)
        if _FINAL_MARKER in user:
            head, tail = user.rsplit(_FINAL_MARKER, 1)
            user = head.rstrip() + "\n\n" + addendum + "\n\n" + _FINAL_MARKER + tail
        else:                                     # a prompt override removed the marker: append instead
            user = user + "\n\n" + addendum
        return system, user, n_requested

    # ------------------------------------------------------------------ main entry point

    def suggest(
        self,
        history: list[StepRecord],
        candidates: list[Any],
        batch_size: int,
        model_prediction: ModelPrediction | None = None,
        task_context: dict[str, Any] | None = None,
    ) -> list[Any]:
        if not candidates:
            return []

        label = "llm_single_topup_blind" if not self.reveal_labels else "llm_single_topup"
        ctx = task_context or {}
        screen = ctx.get("screen_name") or ctx.get("dataset_name") or "?"
        step = len(history) + 1
        cand_upper = {str(c).upper(): c for c in candidates}

        picked: list[Any] = []
        seen: set[Any] = set()
        rejected: list[str] = []
        rejected_set: set[str] = set()
        attempts: list[dict[str, Any]] = []
        total_raw = total_matched = total_unknown = 0
        first = {}                                  # first attempt's prompt/reply, kept for the trace
        zero_streak = 0

        for attempt in range(1 + self.topup_max_attempts):
            need = batch_size - len(picked)
            if need <= 0:
                break
            if attempt == 0:
                system, user, n_ask = self._build_prompt(history, batch_size, task_context)
            else:
                ask = max(need, int(math.ceil(need * self.topup_oversample)))
                system, user, n_ask = self._build_topup_prompt(history, task_context, ask, picked, rejected)

            try:
                completion = complete_ex(self.llm, prompt=user, system=system)
                text = completion.text
                reasoning = completion.reasoning or ""
            except Exception as e:                  # the client already retried transient errors
                log.warning("[%s step %d] LLM call failed on attempt %d (%s: %s); stopping with %d of %d genes "
                            "(no padding).", screen, step, attempt + 1, type(e).__name__, e, len(picked), batch_size)
                attempts.append({"attempt": attempt + 1, "requested": n_ask, "error": f"{type(e).__name__}: {e}"})
                if attempt == 0:
                    first = {"system": system, "user": user, "preview": "", "reasoning": ""}
                break

            picks_raw = extract_gene_list(text)
            n_matched = n_unknown = n_dup = n_new = 0
            for g in picks_raw:
                c = cand_upper.get(g)
                if c is None:                       # already tested (not in the unrevealed pool) or not a valid symbol
                    n_unknown += 1
                    if g not in rejected_set:
                        rejected_set.add(g)
                        rejected.append(g)
                    continue
                n_matched += 1
                if c in seen:                       # already accepted earlier this round
                    n_dup += 1
                    continue
                if len(picked) >= batch_size:
                    continue
                picked.append(c)
                seen.add(c)
                n_new += 1

            total_raw += len(picks_raw)
            total_matched += n_matched
            total_unknown += n_unknown
            attempts.append({
                "attempt": attempt + 1, "requested": n_ask, "raw": len(picks_raw), "matched": n_matched,
                "new_accepted": n_new, "in_round_dups": n_dup, "unknown": n_unknown,
                "finish_reason": getattr(completion, "finish_reason", None),
            })
            if attempt == 0:
                first = {"system": system, "user": user, "preview": text[:600], "reasoning": reasoning}
            else:
                zero_streak = zero_streak + 1 if n_new == 0 else 0
                if zero_streak >= _STALL_LIMIT:
                    log.warning("[%s step %d] top-up stalled after %d attempts with no new genes.", screen, step, attempt + 1)
                    break

        shortfall = max(0, batch_size - len(picked))
        first_try = attempts[0].get("new_accepted", 0) if attempts and "error" not in attempts[0] else 0
        if shortfall > 0:
            log.warning("[%s step %d] still under-supplied after %d call(s): %d of %d (first attempt gave %d). NOT padding.",
                        screen, step, len(attempts), len(picked), batch_size, first_try)
        elif len(attempts) > 1:
            log.info("[%s step %d] top-up filled the round: %d from the first reply, %d from %d more call(s).",
                     screen, step, first_try, len(picked) - first_try, len(attempts) - 1)

        self._last_trace = {
            "acquisition": label,
            "reveal_labels": self.reveal_labels,
            "llm": _llm_trace_dict(self.llm),
            "prompt_override": self.prompt_override.describe() if self.prompt_override else None,
            "n_candidates_in_prompt": 0,
            "n_requested_from_llm": attempts[0].get("requested") if attempts else batch_size,
            "raw_response_preview": first.get("preview", ""),
            "n_raw_picks": total_raw,
            "n_matched_in_pool": total_matched,
            "n_out_of_pool": total_unknown,
            "n_llm_picked": len(picked),
            "n_already_revealed": 0,
            "shortfall": shortfall,
            "prompt_system": first.get("system", ""),
            "prompt_user": first.get("user", ""),
            "reasoning": first.get("reasoning", ""),
            # NOTE: no "response_text" on purpose (see module docstring); llm_calls.jsonl has every reply.
            "picked_genes": [str(c) for c in picked],
            "topup": {
                "max_attempts": self.topup_max_attempts,
                "oversample": self.topup_oversample,
                "n_llm_calls": len(attempts),
                "first_try_accepted": first_try,
                "topup_accepted": len(picked) - first_try,
                "reached_target": shortfall == 0,
                "attempts": attempts,
            },
        }
        if attempts and "error" in attempts[0]:     # first call failed: same shape as the default's error trace
            self._last_trace["error"] = attempts[0]["error"]
        return picked


__all__ = ["LLMSingleTopUpAcquisition"]
