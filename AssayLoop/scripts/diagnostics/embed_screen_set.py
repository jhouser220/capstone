"""Embed a screen set's descriptions into the writable cache layer, resumably.

Why this exists rather than calling `embed_screens` directly. `embed_texts` embeds
every missing text and saves ONCE at the end, so a failure in the last chunk discards
all the work before it. And OpenAI's own SDK does not retry 403, which is the error a
freshly-granted model permission returns while it propagates across gateway nodes --
observed here flipping between success and `model_not_found` minutes apart on the same
key. So: small chunks, save after every chunk, and retry 403/429 with backoff.

Safe to re-run. Anything already cached is skipped, so an interrupted run resumes
where it stopped and costs nothing for what it already has.

    uv run python scripts/diagnostics/embed_screen_set.py --screen-set test
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, "src")
import assayloop.config  # noqa: F401,E402  -- loads .env
from assayloop.amortized import text_embed as te   # noqa: E402
from assayloop.tasks import load_screens           # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--screen-set", default="test")
    ap.add_argument("--chunk", type=int, default=32)
    ap.add_argument("--retries", type=int, default=12)
    ap.add_argument("--wait", type=float, default=30.0)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    screens = load_screens(target_set=a.screen_set)
    emb = te.get_text_embedder("auto")
    cache = te.EmbeddingCache()
    texts = [te.screen_description_text(s) for s in screens]
    key = lambda t: te._text_key(t, model=emb.model)   # noqa: E731

    todo, seen = [], set()
    for t in texts:
        k = key(t)
        if cache.get(k) is None and k not in seen:
            seen.add(k)
            todo.append(t)
    chars = sum(len(t) for t in todo)
    print(f"{a.screen_set}: {len(screens)} screens, cache holds {len(cache)}, "
          f"{len(todo)} unique texts to embed (~{chars // 4} tokens, "
          f"~${chars / 4 / 1e6 * 0.02:.5f})")
    if a.dry_run or not todo:
        print("nothing to do" if not todo else "dry run, nothing sent")
        return

    done = 0
    for i in range(0, len(todo), a.chunk):
        chunk = todo[i:i + a.chunk]
        for attempt in range(1, a.retries + 1):
            try:
                vecs = emb._embed_uncached(chunk)
                for t, v in zip(chunk, vecs):
                    cache.put(key(t), v)
                cache.save()          # persist after EVERY chunk, so a later
                done += len(chunk)    # failure never discards earlier work
                print(f"  {done}/{len(todo)} embedded", flush=True)
                break
            except Exception as e:
                msg = str(e)
                if "insufficient_quota" in msg:
                    raise SystemExit("STOP: account has no credit; nothing further "
                                     "will succeed. Progress so far is saved.")
                transient = ("does not have access" in msg or "rate_limit" in msg
                             or "429" in msg or "500" in msg or "503" in msg)
                if not transient or attempt == a.retries:
                    raise
                print(f"  chunk {i // a.chunk}: {type(e).__name__} "
                      f"(attempt {attempt}/{a.retries}), waiting {a.wait:.0f}s",
                      flush=True)
                time.sleep(a.wait)

    fresh = te.EmbeddingCache()
    unresolved = [s.dataset_name for s, t in zip(screens, texts)
                  if fresh.get(key(t)) is None]
    v = np.stack([fresh.get(key(t)) for t in texts])
    print(f"\ncache re-read from disk: {len(fresh)} vectors at {fresh.path}")
    print(f"unresolved screens: {len(unresolved)}")
    print(f"all {len(screens)} vectors: shape {v.shape}, dtype {v.dtype}, "
          f"norms {np.linalg.norm(v, axis=1).min():.3f}-"
          f"{np.linalg.norm(v, axis=1).max():.3f}")
    if unresolved:
        raise SystemExit(f"INCOMPLETE: {len(unresolved)} screens still missing; "
                         f"re-run this command to resume.")
    print("COMPLETE: every screen in this set resolves from cache.")


if __name__ == "__main__":
    main()
