"""Find a sweep's run directories after `output/runs` was sorted into method folders.

Run directories used to sit flat in `output/runs`, so every scorer in this project builds
its paths as `<runs-dir>/<sweep>-<NN>-<screen>/result.json`. They now live one level down,
in a folder named for the METHOD the arm tested (`af_depmap`, `af_depmap_screenknn`,
`af_alone_baseline`, ...), the way finished LLM runs have always lived under
`baseline/<subfolder>/`.

`run_dir()` keeps the old call shape working. It tries the literal path first, so an
explicit `--dir output/runs/af_depmap` still resolves with no search at all, then falls
back to an index of every sweep id found below the root. The index is built once per
root and cached.

A paired contrast now routinely spans two folders -- the weight-0 baseline of a family
sits in `af_alone_baseline` while its treatment arm sits in `af_depmap` -- so scorers
that take a single `--dir` only work because of this fallback. Passing the runs ROOT
(`output/runs`) is therefore the right thing to do, not a workaround.
"""
from __future__ import annotations

import re
from functools import lru_cache
from pathlib import Path

_SWEEP_DIR = re.compile(r"^(?P<sid>.+)-(?P<idx>[0-9]+)-(?P<screen>.+)$")


@lru_cache(maxsize=32)
def _index(root: str) -> dict[str, Path]:
    """sweep id -> the directory that holds that sweep's run folders."""
    base = Path(root)
    out: dict[str, Path] = {}
    if not base.is_dir():
        return out
    # the root itself, then one and two levels down (method folders, and
    # baseline/<subfolder> which has nested this way since before the move)
    for depth in ("", "*/", "*/*/"):
        for p in base.glob(depth + "*"):
            if not p.is_dir():
                continue
            m = _SWEEP_DIR.match(p.name)
            if m and (p / "config.json").is_file():
                out.setdefault(m.group("sid"), p.parent)
    return out


def sweep_root(runs_dir: str | Path, sweep: str) -> Path:
    """The directory holding `sweep`'s run folders, searching below `runs_dir`."""
    base = Path(runs_dir)
    hit = _index(str(base)).get(sweep)
    return hit if hit is not None else base


def run_dir(runs_dir: str | Path, sweep: str, index: int, screen: str,
            width: int = 2) -> Path:
    """`<runs_dir>/[method/]<sweep>-<index>-<screen>`, wherever that sweep now lives."""
    name = f"{sweep}-{index:0{width}d}-{screen}"
    direct = Path(runs_dir) / name
    if direct.is_dir():
        return direct
    return sweep_root(runs_dir, sweep) / name


def sweep_ids(runs_dir: str | Path) -> list[str]:
    """Every sweep id found below `runs_dir`."""
    return sorted(_index(str(Path(runs_dir))))
