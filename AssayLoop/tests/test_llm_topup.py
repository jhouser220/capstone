"""llm_single_topup: first attempt == default behaviour; short rounds are topped up; limits, errors, trace, replay."""
from __future__ import annotations

import json
import math
from types import SimpleNamespace

import pytest

pytest.importorskip("assaybench")
from assaybench.core.types import Observation, StepRecord

from assayloop.acquisitions import llm_single_acq as parent_mod
from assayloop.acquisitions import llm_single_topup_acq as topup_mod
from assayloop.acquisitions.llm_single_acq import LLMSingleAcquisition
from assayloop.acquisitions.llm_single_topup_acq import LLMSingleTopUpAcquisition
from assayloop.llm.client import LLMClientConfig

CTX = {"cell_line": "SK-N-DZ", "phenotype": "increases RSL3 resistance.", "organism": "Homo sapiens"}
TESTED = [f"T{i:04d}" for i in range(1, 101)]                 # genes already tested in earlier rounds
POOL = [f"G{i:04d}" for i in range(1, 401)]                   # untested candidates


def _history():
    obs = [Observation(candidate=g, label={"hit": i % 7 == 0}) for i, g in enumerate(TESTED)]
    return [StepRecord(step=1, acquired_batch=list(TESTED), new_observations=obs, acquisition_trace={})]


class FakeLLM:
    """Scripted replacement for complete_ex: pops one reply per call and records every prompt."""

    def __init__(self, replies):
        self.replies, self.calls = list(replies), []

    def __call__(self, llm, prompt, system=None):
        self.calls.append({"prompt": prompt, "system": system})
        item = self.replies.pop(0)
        if isinstance(item, Exception):
            raise item
        return SimpleNamespace(text=item, reasoning="", finish_reason="stop")


def _install(monkeypatch, replies):
    fake = FakeLLM(replies)
    monkeypatch.setattr(topup_mod, "complete_ex", fake)
    monkeypatch.setattr(parent_mod, "complete_ex", fake)
    return fake


def _acq(**kw):
    return LLMSingleTopUpAcquisition(llm=LLMClientConfig(provider="vllm", model="m", api_key="k", base_url="http://x"), **kw)


def _reply(genes):
    return ", ".join(genes)


def _history_block(prompt):
    """The history section only: from its header up to the top-up note (if any) or the final instruction."""
    start = prompt.index("Active-Learning History")
    ends = [prompt.index(m, start) for m in ("## Top-up request", "=== Final instruction") if m in prompt[start:]]
    return prompt[start: min(ends)].rstrip()


# ---------------------------------------------------------------- parity with the default method


def test_good_first_reply_behaves_exactly_like_the_default(monkeypatch):
    reply = _reply(POOL[:100])
    fake = _install(monkeypatch, [reply])
    got = _acq().suggest(_history(), list(POOL), 100, task_context=CTX)
    top_prompt = fake.calls[0]

    fake_p = _install(monkeypatch, [reply])
    ref = LLMSingleAcquisition(llm=LLMClientConfig(provider="vllm", model="m", api_key="k", base_url="http://x")) \
        .suggest(_history(), list(POOL), 100, task_context=CTX)

    assert len(fake.calls) == 1                                  # no extra call when the round is already full
    assert got == ref == POOL[:100]
    assert top_prompt == fake_p.calls[0]                         # same system + user prompt as the default


# ---------------------------------------------------------------- the failure mode this exists for


def test_collapsed_first_reply_is_topped_up_to_the_full_batch(monkeypatch):
    first = _reply(TESTED[:88] + POOL[:9] + ["NOTAGENE"])        # 88 repeats, 9 new, 1 invalid (like the run-2 collapse)
    second = _reply(POOL[9:9 + 137])                             # plenty of new genes
    fake = _install(monkeypatch, [first, second])
    acq = _acq()
    got = acq.suggest(_history(), list(POOL), 100, task_context=CTX)

    assert len(got) == 100 and len(set(got)) == 100
    assert got[:9] == POOL[:9] and got[9:] == POOL[9:100]        # first reply's valid genes kept, top-up fills the rest
    assert len(fake.calls) == 2

    retry = fake.calls[1]["prompt"]
    need = 91
    ask = math.ceil(need * 1.5)                                  # 137
    assert "## Top-up request" in retry
    assert "Do NOT suggest any of these again: " + ", ".join(POOL[:9]) in retry
    assert "NOTAGENE" in retry and "T0001" in retry              # rejected: invalid symbol and an already-tested gene
    assert f"{ask}" in retry.split("=== Final instruction ===")[1]   # final instruction asks for the new count
    assert str(ask) in fake.calls[1]["system"]                   # ... and so does the system prompt
    # the original prompt is otherwise unchanged: same history block, and no outcome for this round's genes
    assert _history_block(retry) == _history_block(fake.calls[0]["prompt"])
    for g in POOL[:9]:
        assert f"{g}: hit" not in retry and f"{g} (hit" not in retry

    t = acq.last_trace()
    assert "response_text" not in t                              # see module docstring (replay compatibility)
    assert t["topup"]["n_llm_calls"] == 2
    assert t["topup"]["first_try_accepted"] == 9 and t["topup"]["topup_accepted"] == 91
    assert t["topup"]["reached_target"] is True and t["shortfall"] == 0 and t["n_llm_picked"] == 100
    assert [a["new_accepted"] for a in t["topup"]["attempts"]] == [9, 91]
    assert t["topup"]["attempts"][0]["unknown"] == 89           # 88 repeats + 1 invalid


def test_several_topups_when_each_adds_a_little(monkeypatch):
    replies = [_reply(POOL[0:40])] + [_reply(POOL[40 + 20 * i: 40 + 20 * (i + 1)]) for i in range(3)]   # 40 + 20 + 20 + 20 = 100
    fake = _install(monkeypatch, replies)
    got = _acq().suggest(_history(), list(POOL), 100, task_context=CTX)
    assert got == POOL[:100] and len(fake.calls) == 4


# ---------------------------------------------------------------- limits: never pads, never loops forever


def test_stops_when_topups_stall_and_returns_a_short_batch(monkeypatch):
    junk = _reply(TESTED[:50])                                   # only repeats: nothing new, ever
    fake = _install(monkeypatch, [_reply(POOL[:30]), junk, junk, junk, junk])
    acq = _acq(topup_max_attempts=4)
    got = acq.suggest(_history(), list(POOL), 100, task_context=CTX)
    assert got == POOL[:30]                                      # nothing invented to fill the gap
    assert len(fake.calls) == 3                                  # first call + 2 consecutive empty top-ups, then stop
    t = acq.last_trace()["topup"]
    assert t["reached_target"] is False and acq.last_trace()["shortfall"] == 70


def test_max_attempts_is_respected(monkeypatch):
    replies = [_reply(POOL[0:10])] + [_reply([POOL[10 + i]]) for i in range(10)]   # each top-up adds exactly one gene
    fake = _install(monkeypatch, replies)
    acq = _acq(topup_max_attempts=3)
    got = acq.suggest(_history(), list(POOL), 100, task_context=CTX)
    assert len(fake.calls) == 4 and len(got) == 13               # 1 + 3 attempts, then it stops
    assert acq.last_trace()["topup"]["n_llm_calls"] == 4


def test_zero_topups_is_the_default_behaviour(monkeypatch):
    fake = _install(monkeypatch, [_reply(POOL[:30])])
    got = _acq(topup_max_attempts=0).suggest(_history(), list(POOL), 100, task_context=CTX)
    assert len(got) == 30 and len(fake.calls) == 1


# ---------------------------------------------------------------- errors


def test_error_on_a_topup_keeps_what_was_already_accepted(monkeypatch):
    fake = _install(monkeypatch, [_reply(POOL[:30]), RuntimeError("quota")])
    acq = _acq()
    got = acq.suggest(_history(), list(POOL), 100, task_context=CTX)
    assert got == POOL[:30] and len(fake.calls) == 2
    attempts = acq.last_trace()["topup"]["attempts"]
    assert "RuntimeError" in attempts[1]["error"] and "error" not in acq.last_trace()


def test_error_on_the_first_call_matches_the_default_error_shape(monkeypatch):
    _install(monkeypatch, [RuntimeError("down")])
    acq = _acq()
    assert acq.suggest(_history(), list(POOL), 100, task_context=CTX) == []
    t = acq.last_trace()
    assert "RuntimeError" in t["error"] and t["shortfall"] == 100 and t["n_llm_picked"] == 0


# ---------------------------------------------------------------- compatibility with the handoff replay paths


def _result_step(acq, picked, step):
    return {"step": step, "acquired_batch": [str(g) for g in picked], "hits": [False] * len(picked),
            "acquisition_trace": dict(acq.last_trace())}


def test_replay_uses_the_accepted_batch_even_with_extra_llm_calls(monkeypatch, tmp_path):
    from assayloop.amortized.warmstart import WarmStart

    # three rounds, each needing a top-up, each with a different (disjoint) accepted batch
    acq, steps, acquired = _acq(), [], set()
    for i in range(3):
        pool = [g for g in POOL if g not in acquired]
        _install(monkeypatch, [_reply(TESTED[:88] + pool[:9]), _reply(pool[9:140])])
        picked = acq.suggest(_history(), pool, 100, task_context=CTX)
        assert len(picked) == 100
        acquired |= set(picked)
        steps.append(_result_step(acq, picked, i + 1))

    d = tmp_path / "sweep-x-00-2054"
    d.mkdir()
    (d / "result.json").write_text(json.dumps({
        "run_id": "sweep-x-00-2054", "task_id": "option2/2054",
        "config": {"acquisition": "llm_single_topup/vllm:m", "batch_size": 100, "n_steps": 3, "model": "null"},
        "final_metrics": {}, "steps": steps}))
    # 5 logged LLM calls for 3 steps: the count mismatch the paper-style replay must tolerate
    (d / "llm_calls.jsonl").write_text("\n".join(json.dumps({"response_text": "A, B", "error": None}) for _ in range(5)))

    universe = set(POOL) | set(TESTED)
    ws = WarmStart.from_run_dirs(str(tmp_path), "sweep-x-", universe_genes=universe, batch_size=100)
    assert [len(r) for r in ws.traces["2054"][0]] == [100, 100, 100]
    legacy = WarmStart.from_run_dirs(str(tmp_path), "sweep-x-")                # what eval-ranker-handoff does
    assert [len(r) for r in legacy.traces["2054"][0]] == [100, 100, 100]


# ---------------------------------------------------------------- configuration and registration


def test_env_vars_configure_the_loop(monkeypatch):
    monkeypatch.setenv("ASSAYLOOP_TOPUP_MAX_ATTEMPTS", "7")
    monkeypatch.setenv("ASSAYLOOP_TOPUP_OVERSAMPLE", "2.0")
    a = _acq()
    assert (a.topup_max_attempts, a.topup_oversample) == (7, 2.0)
    assert _acq(topup_max_attempts=1, topup_oversample=1.2).topup_max_attempts == 1


def test_registered_with_runner_and_cli():
    from assayloop import cli
    from assayloop.experiment.runner import make_acquisition

    acq = make_acquisition("llm_single_topup", seed=3)
    assert isinstance(acq, LLMSingleTopUpAcquisition) and acq.name().startswith("llm_single_topup/")
    for group in (cli._OPEN_VOCAB_LLM_ACQS, cli._LLM_CLIENT_ACQS, cli._ZERO_SHOT_ACQS):
        assert "llm_single_topup" in group
