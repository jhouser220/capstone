# GPT-5.6 Sol warm start on `paper_test`, plus our test-fold text embeddings

Everything needed to reproduce the AssayLoop handoff experiments in this repo without
re-running (or paying for) the LLM. Two independent pieces:

- `warm_start/` - the LLM's picks for rounds 1-3 on the paper's 20 `paper_test` screens,
  with the raw calls behind them.
- `embeddings/` - the `text-embedding-3-small` screen-description vectors for the `test`
  fold that we computed with our own OpenAI key.

## What the warm start is

The paper's default method, unmodified: `--acq llm_single`, k=3 warm rounds, batch 100,
seed 0, **no top-up**. Candidate pool is the paper's f2 universe (genes measured in >= 2
`paper_test` screens, 21,147 genes). `reasoning_effort` was left unset, so OpenAI's
default (medium) applied. Reached through OpenAI direct (`api.openai.com/v1`), not a
gateway and not Bedrock - `openai.gpt-5.6-sol` is listed on Bedrock's Mantle route but
inference returns `401 access_denied` on our account.

Generating command (one screen per `--screen`, or omit it for the whole set):

```bash
assayloop run --model null --acq llm_single --screen-set paper_test \
  --n-steps 3 --batch-size 100 --warm-start 0 --seed 0 \
  --lm-config configs/lm/collect-gpt-5.6-sol-openai.yaml \
  --sweep-id sweep-sol56-pt-cost1 --resume
```

### 15 of 20 screens are usable

`1827`, `U_1863_inc`, `U_1987_merged`, `2060` and `2076` lost warm rounds to OpenAI's
**`bio_policy`** filter, which rejects the prompt with HTTP 400 ("flagged for possible
biological risk"). `llm_single` returns an empty batch and never pads, so those rounds
produced zero picks.

Those five are **named `calls/<screen>_declined.jsonl`** so a refusal is visible from the
file listing alone, and `per_screen.csv` flags them in both `declined` and `usable`.
Four were refused on all three rounds; `2060` was refused on two of three, so one round
of real picks exists but the campaign is still unusable as a k=3 warm start. All five are
excluded from `rounds.json` and their logs are kept so the failure is inspectable.

Any aggregate computed from this warm start therefore covers 15 screens and **is not
comparable to a published 20-screen figure** such as the paper's EF 5.67.

### Files

| file | what |
|---|---|
| `warm_start/rounds.json` | `{screen: [[round1 genes], [round2], [round3]]}` for the 15 usable screens. This is the artifact you replay. |
| `warm_start/per_screen.csv` | one row per screen: picks and hits per round, blocked-call count, token usage, `picks_sha`. |
| `warm_start/manifest.json` | model, route, method, universe, seeds, totals. |
| `warm_start/calls/<screen>.jsonl` | the raw calls: full prompt `messages`, `response_text`, `finish_reason`, and `usage` including `completion_tokens_details.reasoning_tokens`. |
| `warm_start/calls/<screen>_declined.jsonl` | same, for a screen where gpt-5.6-sol **refused at least one call** on `bio_policy` grounds. A refusal row has the rejection in `error` and `usage: null` - nothing was generated, so nothing was billed. |

`picks_sha` is `sha256` of the newline-joined picks across the three rounds, truncated to
16 hex chars - use it to confirm a replay reproduced the same campaign.

### Replaying it

`rounds.json` plugs straight into `GlmHandoffAcquisition`, which replays recorded rounds
for the first `n_warm` rounds and then hands off:

```python
from assayloop.acquisitions.glm_handoff_acq import GlmHandoffAcquisition
rounds = json.loads(Path("warm_start/rounds.json").read_text())
acq = GlmHandoffAcquisition(rounds=rounds[screen_id], n_warm=3, base=your_arm)
# retain=25, rank="model" reproduces the capped ("retain25") protocol instead
```

`scripts/external/run_glm5_handoff_hsa.py` does this across a screen set; point
`--warm-sweep` at the recorded sweep, or adapt it to read `rounds.json` directly.

### Cost, measured

46 billed calls (14 of 60 were rejected by the filter before generating anything):
**47,457 input tokens, 201,750 output tokens, of which 180,375 (89%) were reasoning.**
At the GPT-5.6 Sol promotional rate of $4/$20 per 1M that is **$4.22**.

## The embeddings

`embeddings/test_fold_text_embeddings.npz` holds **230 unique `text-embedding-3-small`
vectors (1536-d, float32)** covering **286 of the `test` fold's 334 screens**. Fewer
vectors than screens because screens with identical descriptions share one cache key.
The other 48 screens are already in the package's shipped cache, so the two together
cover the fold.

Keys are `assayloop.amortized.text_embed._text_key(description, model)`;
`embeddings/key_map.json` maps `dataset_name -> key` so you can see which screen is which
without the description text being republished here.

To use them, merge into the writable cache layer:

```python
import numpy as np, shutil
from assayloop.amortized.text_embed import EmbeddingCache
c = EmbeddingCache()
merged = {**dict(np.load(c.path)), **dict(np.load("embeddings/test_fold_text_embeddings.npz"))}
np.savez_compressed(c.path, **merged)
```

`EmbeddingCache` reads the shipped cache first and this store second, so anything already
shipped is untouched.

## Provenance and safety

The original `output/runs/.../config.json` files are **deliberately not included**: the
runner serialises the whole LLM client config into `acq_overrides.llm`, API key and all.
Everything in this folder was scanned for `sk-`, `ABSK`, `AIza` and `AKIA` patterns before
commit. `result.json` and `llm_calls.jsonl` never contained credentials.

Screen descriptions and hit labels come from the public `Genentech/assaybench` dataset;
nothing here republishes that text. The embeddings are derived values computed from it.
