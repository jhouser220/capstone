"""Cost of an LLM sweep from its own llm_calls.jsonl, and the extrapolation to a full run.

    uv run python bedrock_cost.py sweep-solpilot --calls-in-full-run 60

Default prices are Bedrock's Aug-2026 promotional rate for GPT-5.6 Sol ($4/M input, $20/M
output); pass --in-price/--out-price for any other model or tier. Reasoning tokens are billed
as output, so a reasoning model's cost is almost entirely the output term.
"""
from __future__ import annotations

import argparse
import glob
import json


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("sweep_id", help="e.g. sweep-solpilot (run folders are matched as <sweep_id>-*)")
    ap.add_argument("--runs-dir", default="output/runs", help="searched recursively")
    ap.add_argument("--in-price", type=float, default=4.0, help="$ per million input tokens")
    ap.add_argument("--out-price", type=float, default=20.0, help="$ per million output tokens")
    ap.add_argument("--calls-in-full-run", type=int, default=60,
                    help="calls the full sweep would make (paper default method: 20 screens x 3 rounds = 60)")
    args = ap.parse_args(argv)

    pin = pout = calls = errors = truncated = 0
    for f in glob.glob(f"{args.runs_dir}/**/{args.sweep_id}-*/llm_calls.jsonl", recursive=True):
        for line in open(f, encoding="utf-8"):
            if not line.strip():
                continue
            rec = json.loads(line)
            if rec.get("error"):
                errors += 1
                continue
            u = rec.get("usage") or {}
            if not u.get("prompt_tokens"):
                continue
            calls += 1
            pin += u["prompt_tokens"]
            pout += u.get("completion_tokens") or 0
            truncated += rec.get("finish_reason") == "length"
    if not calls:
        raise SystemExit(f"no usage records under {args.runs_dir}/**/{args.sweep_id}-*/llm_calls.jsonl")

    cost = pin / 1e6 * args.in_price + pout / 1e6 * args.out_price
    print(f"{calls} successful call(s)" + (f", {errors} error(s)" if errors else ""))
    print(f"  input   {pin:9,d} tokens  mean {pin / calls:8,.0f}/call  ${pin / 1e6 * args.in_price:7.4f}")
    print(f"  output  {pout:9,d} tokens  mean {pout / calls:8,.0f}/call  ${pout / 1e6 * args.out_price:7.4f}")
    print(f"  measured cost ${cost:.4f}  (${cost / calls:.4f} per call)")
    if truncated:
        print(f"  WARNING: {truncated} call(s) stopped at the token cap (finish_reason='length'). Those answers "
              "were cut off but still billed; raise max_tokens or lower reasoning_effort.")
    full = cost / calls * args.calls_in_full_run
    print(f"\nextrapolated to {args.calls_in_full_run} calls: ${full:.2f}")
    print(f"  if the promotional price has lapsed (~2.2x at $8.80/$33): ${full * 2.2:.2f}")
    return cost


if __name__ == "__main__":
    main()
