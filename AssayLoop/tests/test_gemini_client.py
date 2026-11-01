"""Reasoning effort, retry/pacing behaviour, and prompt overrides.

No network and no API key: the OpenAI SDK is replaced by a scripted fake, and
the clock is faked wherever timing matters.
"""
from __future__ import annotations

import sys
import threading
import time
import types
from pathlib import Path
from types import SimpleNamespace

import pytest

from assayloop import config as bl_config
from assayloop.llm import client as client_mod
from assayloop.llm.client import LLMClientConfig, complete_ex
from assayloop.llm.prompt_override import (
    load_prompt_override,
    render_template,
)

REPO = Path(__file__).resolve().parent.parent


# ---------------------------------------------------------------- fakes


def _api_error(name: str, msg: str, status: int, headers: dict | None = None):
    """An exception whose *class name* matches the SDK's (the client classifies
    errors by type name), carrying a status code and response headers."""
    exc = type(name, (Exception,), {})(msg)
    exc.status_code = status
    exc.response = SimpleNamespace(status_code=status, headers=headers or {})
    return exc


def _ok(text: str = "TP53, MYC"):
    msg = SimpleNamespace(content=text, reasoning=None, reasoning_content=None)
    return SimpleNamespace(
        choices=[SimpleNamespace(message=msg, finish_reason="stop")],
        usage=SimpleNamespace(prompt_tokens=1, completion_tokens=2, total_tokens=3),
        model="gemini-3.5-flash-lite",
    )


class FakeOpenAI:
    script: list = []
    instances: list = []

    def __init__(self, **kwargs):
        self.init_kwargs = kwargs
        self.calls: list[dict] = []
        FakeOpenAI.instances.append(self)
        self.chat = SimpleNamespace(completions=SimpleNamespace(create=self._create))

    def _create(self, **kwargs):
        self.calls.append(kwargs)
        item = FakeOpenAI.script.pop(0)
        if isinstance(item, Exception):
            raise item
        return item


@pytest.fixture
def fake_openai(monkeypatch):
    FakeOpenAI.script, FakeOpenAI.instances = [], []
    mod = types.ModuleType("openai")
    mod.OpenAI = FakeOpenAI
    monkeypatch.setitem(sys.modules, "openai", mod)
    monkeypatch.setattr("assayloop.tracing.log_completion_call", lambda **kw: None)
    sleeps: list[float] = []
    monkeypatch.setattr(time, "sleep", lambda s: sleeps.append(s))
    return SimpleNamespace(sleeps=sleeps)


def _cfg(**over) -> LLMClientConfig:
    base = dict(
        provider="vllm",
        model="gemini-3.5-flash-lite",
        base_url="https://example.invalid/v1",
        api_key="k",
        max_retries=3,
        retry_initial_backoff=0.01,
        retry_max_backoff=0.02,
    )
    base.update(over)
    return LLMClientConfig(**base)


# ------------------------------------------------------- reasoning effort


def test_reasoning_effort_is_sent_and_sdk_retries_are_off(fake_openai):
    FakeOpenAI.script = [_ok()]
    out = complete_ex(_cfg(reasoning_effort="low"), "hi")
    assert out.text == "TP53, MYC"
    inst = FakeOpenAI.instances[0]
    assert inst.init_kwargs["max_retries"] == 0
    body = inst.calls[0]["extra_body"]
    assert body == {"reasoning_effort": "low"}


def test_reasoning_effort_supersedes_vllm_thinking_hint(fake_openai):
    FakeOpenAI.script = [_ok()]
    complete_ex(_cfg(reasoning_effort="low", disable_thinking=True), "hi")
    body = FakeOpenAI.instances[0].calls[0]["extra_body"]
    assert "chat_template_kwargs" not in body


def test_disable_thinking_alone_keeps_existing_vllm_behaviour(fake_openai):
    FakeOpenAI.script = [_ok()]
    complete_ex(_cfg(disable_thinking=True), "hi")
    body = FakeOpenAI.instances[0].calls[0]["extra_body"]
    assert body == {"chat_template_kwargs": {"enable_thinking": False}}


def test_nothing_extra_sent_by_default(fake_openai):
    FakeOpenAI.script = [_ok()]
    complete_ex(_cfg(), "hi")
    assert FakeOpenAI.instances[0].calls[0]["extra_body"] is None


def test_env_defaults_fill_reasoning_and_rpm(monkeypatch):
    monkeypatch.setattr(bl_config, "VLLM_REASONING_EFFORT", "low")
    monkeypatch.setattr(bl_config, "VLLM_REQUESTS_PER_MINUTE", 7.0)
    r = LLMClientConfig(provider="vllm").resolved()
    assert (r.reasoning_effort, r.requests_per_minute) == ("low", 7.0)
    # explicit values win over the environment
    r = LLMClientConfig(provider="vllm", reasoning_effort="high").resolved()
    assert r.reasoning_effort == "high"


def test_new_fields_survive_resolved_for_every_provider():
    for provider in ("vllm", "anthropic", "dspy"):
        r = LLMClientConfig(
            provider=provider, model="m", api_key="k",
            reasoning_effort="low", requests_per_minute=5.0,
        ).resolved()
        assert (r.reasoning_effort, r.requests_per_minute) == ("low", 5.0)


# ------------------------------------------------------------------ retries


def test_429_is_retried_and_server_delay_is_honoured(fake_openai):
    FakeOpenAI.script = [
        _api_error("RateLimitError", "slow down", 429, {"retry-after": "7"}),
        _ok(),
    ]
    out = complete_ex(_cfg(), "hi")
    assert out.text == "TP53, MYC"
    assert len(FakeOpenAI.instances[0].calls) == 2
    assert len(fake_openai.sleeps) == 1 and fake_openai.sleeps[0] >= 8.0


def test_google_style_body_delay_is_parsed():
    e = _api_error(
        "RateLimitError",
        "Error code: 429 - [{'error': {'details': [{'retryDelay': '34s'}]}}]",
        429,
    )
    assert client_mod._server_retry_delay(e) == 34.0
    e = _api_error("RateLimitError", "Quota exceeded. Please retry in 12.5s.", 429)
    assert client_mod._server_retry_delay(e) == 12.5
    assert client_mod._server_retry_delay(_api_error("RateLimitError", "x", 429)) is None


def test_daily_quota_is_not_retried(fake_openai):
    err = _api_error(
        "RateLimitError",
        "quotaId: GenerateRequestsPerDayPerProjectPerModel-FreeTier",
        429,
    )
    FakeOpenAI.script = [err]
    with pytest.raises(Exception) as ei:
        complete_ex(_cfg(), "hi")
    assert ei.value is err
    assert len(FakeOpenAI.instances[0].calls) == 1
    assert fake_openai.sleeps == []


def test_bad_request_is_not_retried(fake_openai):
    FakeOpenAI.script = [_api_error("BadRequestError", "unknown field", 400)]
    with pytest.raises(Exception):
        complete_ex(_cfg(), "hi")
    assert len(FakeOpenAI.instances[0].calls) == 1


def test_retries_are_bounded(fake_openai):
    FakeOpenAI.script = [_api_error("RateLimitError", "x", 429) for _ in range(4)]
    with pytest.raises(Exception):
        complete_ex(_cfg(max_retries=3), "hi")
    assert len(FakeOpenAI.instances[0].calls) == 4  # first try + 3 retries


# ------------------------------------------------------------------- pacing


def test_rate_limiter_spaces_slots(monkeypatch):
    monkeypatch.setattr(time, "monotonic", lambda: 100.0)
    waits: list[float] = []
    monkeypatch.setattr(time, "sleep", lambda s: waits.append(s))
    lim = client_mod._RateLimiter(rpm=60)  # one per second
    got = [lim.acquire() for _ in range(4)]
    assert got == pytest.approx([0.0, 1.0, 2.0, 3.0])
    assert waits == pytest.approx([1.0, 2.0, 3.0])


def test_rate_limiter_is_thread_safe(monkeypatch):
    monkeypatch.setattr(time, "monotonic", lambda: 50.0)
    monkeypatch.setattr(time, "sleep", lambda s: None)
    lim = client_mod._RateLimiter(rpm=120)  # 0.5s apart
    results: list[float] = []
    lock = threading.Lock()

    def worker():
        w = lim.acquire()
        with lock:
            results.append(w)

    threads = [threading.Thread(target=worker) for _ in range(8)]
    [t.start() for t in threads]
    [t.join() for t in threads]
    assert sorted(results) == pytest.approx([i * 0.5 for i in range(8)])


def test_limiter_is_shared_per_endpoint_and_off_by_default():
    a = client_mod._get_rate_limiter(_cfg(requests_per_minute=5))
    b = client_mod._get_rate_limiter(_cfg(requests_per_minute=5))
    assert a is b and a is not None
    assert client_mod._get_rate_limiter(_cfg()) is None
    assert client_mod._get_rate_limiter(_cfg(requests_per_minute=0)) is None


def test_every_attempt_is_paced_including_retries(fake_openai, monkeypatch):
    acquired = []
    monkeypatch.setattr(
        client_mod, "_get_rate_limiter",
        lambda cfg: SimpleNamespace(acquire=lambda: acquired.append(1)),
    )
    FakeOpenAI.script = [_api_error("RateLimitError", "x", 429), _ok()]
    complete_ex(_cfg(requests_per_minute=10), "hi")
    assert len(acquired) == 2


# ------------------------------------------------------- YAML -> config


def test_yaml_loader_maps_new_fields_and_nulls_inherited_extensions(tmp_path, monkeypatch):
    from assayloop.llm.assaybench_config import load_lm_config

    monkeypatch.setenv("GEMINI_API_KEY", "secret-key")
    (tmp_path / "collect-predictions.yaml").write_text(
        "lm:\n  provider: local\n  top_k: 20\n  min_p: 0.1\n  api_base: http://localhost:8000/v1\n"
    )
    src = (REPO / "configs/lm/collect-gemini-3.5-flash-lite.yaml").read_text()
    (tmp_path / "collect-gemini-3.5-flash-lite.yaml").write_text(src)
    lm = load_lm_config(tmp_path / "collect-gemini-3.5-flash-lite.yaml").client
    assert lm.reasoning_effort == "low"
    assert lm.requests_per_minute == 10.0
    assert lm.model == "gemini-3.5-flash-lite"
    assert lm.api_key == "secret-key"
    assert lm.top_k is None and lm.min_p is None  # base values nulled out
    assert lm.base_url.startswith("https://generativelanguage.googleapis.com/")


# ------------------------------------------------------ prompt override


def test_render_template_only_touches_known_placeholders():
    t = 'Reply as {"genes": ["A"]} for {n_requested} genes; {unknown} stays.'
    out = render_template(t, {"n_requested": 100})
    assert out == 'Reply as {"genes": ["A"]} for 100 genes; {unknown} stays.'


def test_override_validation(tmp_path):
    def write(text):
        p = tmp_path / "p.yaml"
        p.write_text(text)
        return p

    with pytest.raises(FileNotFoundError):
        load_prompt_override(tmp_path / "missing.yaml")
    with pytest.raises(ValueError, match="unknown key"):
        load_prompt_override(write("sytem: typo\n"))
    with pytest.raises(ValueError, match="screen_block"):
        load_prompt_override(write("user: 'no screen here'\n"))
    with pytest.raises(ValueError, match="regex"):
        load_prompt_override(write("screen_replace:\n  - {pattern: '(', repl: x}\n"))
    ok = load_prompt_override(write("final_instruction: 'go {n_requested}'\n"))
    assert len(ok.sha256) == 64


def _acq_inputs():
    pytest.importorskip("assaybench")
    from assaybench.core.types import Observation, StepRecord

    obs = [
        Observation(candidate="GPX4", label={"hit": True}),
        Observation(candidate="DHODH", label={"hit": False}),
    ]
    warm = StepRecord(
        step=1, acquired_batch=["GPX4", "DHODH"], new_observations=obs,
        acquisition_trace={"warm_start": True},
    )
    ctx = {"cell_line": "SK-N-DZ", "phenotype": "increases RSL3 resistance."}
    return [warm], ctx


def test_example_prompt_file_reproduces_builtin_prompt():
    from assayloop.acquisitions.llm_single_acq import LLMSingleAcquisition

    history, ctx = _acq_inputs()
    builtin = LLMSingleAcquisition(prompt_file=None)._build_prompt(history, 100, ctx)
    example = LLMSingleAcquisition(
        prompt_file=REPO / "configs/prompts/example.yaml"
    )._build_prompt(history, 100, ctx)
    assert example == builtin


def test_override_edits_reach_the_prompt(tmp_path):
    from assayloop.acquisitions.llm_single_acq import LLMSingleAcquisition

    history, ctx = _acq_inputs()
    p = tmp_path / "p.yaml"
    p.write_text(
        "system: 'SYS {batch_size} {convention}'\n"
        "screen_replace:\n  - {pattern: 'ranking genes', repl: 'triaging genes'}\n"
        "final_instruction: 'FINAL {n_requested} | {suffix}'\n"
        "user: |\n  >>> {screen_block}\n  <<< {history_block}\n  {final_instruction}\n"
    )
    system, user, n = LLMSingleAcquisition(prompt_file=p)._build_prompt(history, 50, ctx)
    assert system == "SYS 50 HGNC" and n == 50
    assert "triaging genes" in user and "ranking genes" not in user
    assert user.startswith(">>> ") and "GPX4" in user      # history block present
    assert "FINAL 50 | All gene symbols must use the HGNC" in user


def test_trace_records_provenance_and_never_credentials(monkeypatch, tmp_path):
    from assayloop.acquisitions import llm_single_acq as mod

    history, ctx = _acq_inputs()
    p = tmp_path / "p.yaml"
    p.write_text("final_instruction: 'FINAL {n_requested}'\n")
    monkeypatch.setattr(
        mod, "complete_ex",
        lambda llm, prompt, system: client_mod.Completion(text="TP53, MYC"),
    )
    acq = mod.LLMSingleAcquisition(
        llm=LLMClientConfig(provider="vllm", model="m", api_key="SECRET-KEY",
                            base_url="http://x", extra_headers={"X-Auth": "SECRET-HDR"}),
        prompt_file=p,
    )
    picked = acq.suggest(history, ["TP53", "MYC", "EGFR"], 2, task_context=ctx)
    assert picked == ["TP53", "MYC"]
    trace = acq.last_trace()
    assert trace["prompt_override"]["sha256"] == acq.prompt_override.sha256
    assert "SECRET-KEY" not in str(trace) and "SECRET-HDR" not in str(trace)
