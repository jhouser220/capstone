"""User-supplied prompt overrides for the open-vocabulary LLM acquisitions.

Point ``ASSAYLOOP_PROMPT_FILE`` (or ``LLMSingleAcquisition(prompt_file=...)``)
at a YAML file to change what the LLM is sent, without editing code. With no
file configured nothing changes: the built-in prompt is used byte-for-byte.

Every key is optional. Unset keys fall back to the built-in text.

``system``
    System message. Placeholders: ``{batch_size}`` (alias ``{n_requested}``),
    ``{convention}``, ``{organism}``.

``final_instruction``
    The closing block of the user message, *including* its header line.
    Placeholders: the ones above plus ``{suffix}`` (the gene-symbol-namespace
    sentence).

``screen_replace``
    A list of ``{pattern, repl}`` regex substitutions applied to the rendered
    screen block (AssayBench's Goal / Experimental Context / Hit Definition /
    ... text). This is how to reword that section without touching the
    companion ``assaybench`` package.

``user``
    Full layout of the user message. Must contain ``{screen_block}``; should
    contain ``{history_block}`` (that is where the hit/non-hit feedback goes)
    and ``{final_instruction}``. Without ``user`` the layout is the built-in
    ``screen_block + history_block + final_instruction``.

Only ``{name}`` tokens for the placeholders above are substituted. Other braces
(JSON examples, etc.) pass through untouched. The reply parser accepts a plain
comma-separated list, a JSON object with a ``genes`` key, or text wrapped in
``<Final Answer>...</Final Answer>`` -- keep the requested output format
within those.

The loaded file's path and SHA-256 are recorded in each step's trace, and the
full rendered prompt is already stored there, so a run documents the prompt it
used.
"""

from __future__ import annotations

import hashlib
import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

import yaml

log = logging.getLogger(__name__)

_PLACEHOLDER_RE = re.compile(r"\{([A-Za-z_][A-Za-z0-9_]*)\}")

_BASE_PLACEHOLDERS = frozenset({"batch_size", "n_requested", "convention", "organism"})
_FINAL_PLACEHOLDERS = _BASE_PLACEHOLDERS | {"suffix"}
_USER_PLACEHOLDERS = _FINAL_PLACEHOLDERS | {
    "screen_block",
    "history_block",
    "final_instruction",
}
_ALLOWED_KEYS = frozenset({"system", "user", "final_instruction", "screen_replace"})


def render_template(template: str, values: Mapping[str, Any]) -> str:
    """Substitute ``{name}`` tokens that appear in ``values``; leave the rest."""

    def _sub(m: re.Match) -> str:
        key = m.group(1)
        return str(values[key]) if key in values else m.group(0)

    return _PLACEHOLDER_RE.sub(_sub, template)


def _warn_unknown(field_name: str, template: str, known: frozenset[str], path: Path) -> None:
    unknown = sorted(set(_PLACEHOLDER_RE.findall(template)) - known)
    if unknown:
        log.warning(
            "%s: '%s' contains {%s}, which is not a known placeholder and will "
            "be sent to the model literally (known: %s).",
            path.name, field_name, "}, {".join(unknown), ", ".join(sorted(known)),
        )


@dataclass
class PromptOverride:
    path: Path
    sha256: str
    system: str | None = None
    user: str | None = None
    final_instruction: str | None = None
    screen_replace: list[tuple[re.Pattern, str]] = field(default_factory=list)

    def apply_screen_replace(self, screen_block: str) -> str:
        for pattern, repl in self.screen_replace:
            screen_block = pattern.sub(repl, screen_block)
        return screen_block

    def describe(self) -> dict[str, str]:
        """Small provenance record for the step trace."""
        return {"path": str(self.path), "sha256": self.sha256}


def load_prompt_override(path: str | Path) -> PromptOverride:
    """Load and validate a prompt-override YAML. Raises on anything suspect
    rather than silently falling back to the built-in prompt."""
    p = Path(path).expanduser()
    if not p.is_file():
        raise FileNotFoundError(f"prompt override file not found: {p}")
    raw = p.read_bytes()
    data = yaml.safe_load(raw) or {}
    if not isinstance(data, dict):
        raise ValueError(f"{p}: expected a mapping at the top level")

    unknown_keys = set(data) - _ALLOWED_KEYS
    if unknown_keys:
        raise ValueError(
            f"{p}: unknown key(s) {sorted(unknown_keys)}; allowed: {sorted(_ALLOWED_KEYS)}"
        )

    def _text(key: str) -> str | None:
        v = data.get(key)
        if v is None:
            return None
        if not isinstance(v, str):
            raise ValueError(f"{p}: '{key}' must be a string")
        v = v.strip()
        if not v:
            raise ValueError(f"{p}: '{key}' is empty")
        return v

    system = _text("system")
    user = _text("user")
    final_instruction = _text("final_instruction")

    if system:
        _warn_unknown("system", system, _BASE_PLACEHOLDERS, p)
    if final_instruction:
        _warn_unknown("final_instruction", final_instruction, _FINAL_PLACEHOLDERS, p)
    if user:
        if "{screen_block}" not in user:
            raise ValueError(
                f"{p}: 'user' must contain {{screen_block}}; without it the model "
                "never sees the screen description."
            )
        if "{history_block}" not in user:
            log.warning(
                "%s: 'user' has no {history_block}; the model will receive no "
                "hit/non-hit feedback (a blind ablation, not AssayLoop's setup).",
                p.name,
            )
        _warn_unknown("user", user, _USER_PLACEHOLDERS, p)

    replacements: list[tuple[re.Pattern, str]] = []
    for i, item in enumerate(data.get("screen_replace") or []):
        if not isinstance(item, dict) or "pattern" not in item or "repl" not in item:
            raise ValueError(f"{p}: screen_replace[{i}] needs 'pattern' and 'repl'")
        try:
            replacements.append((re.compile(str(item["pattern"])), str(item["repl"])))
        except re.error as e:
            raise ValueError(f"{p}: screen_replace[{i}].pattern is not a valid regex: {e}")

    return PromptOverride(
        path=p.resolve(),
        sha256=hashlib.sha256(raw).hexdigest(),
        system=system,
        user=user,
        final_instruction=final_instruction,
        screen_replace=replacements,
    )


__all__ = ["PromptOverride", "load_prompt_override", "render_template"]
