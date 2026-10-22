# Reproducing the presentation comparison

There are two useful checks:

1. **Exact historical replay:** rerun all ten decisions using the original model responses, with the original prompts and parser. Recompute numerical scores, selections, feedback and final metrics. Compare every decision with the frozen trajectories and every metric with the presentation's campaign table.
2. **Fresh replication:** rerun both main systems with new AWS responses. These stochastic runs assess whether the comparison holds again; their values can differ from the historical presentation.

## Fixed inputs

Download the [companion reproduction bundle](https://github.com/jhouser220/capstone/releases/download/odysseus-assayloop-reproduction-2026-10-22/assayloop-presentation-reproduction-2026-10-22.tar.gz) from the [versioned GitHub Release](https://github.com/jhouser220/capstone/releases/tag/odysseus-assayloop-reproduction-2026-10-22). Its size and SHA-256 are recorded in `REPRODUCTION_ASSETS.json`. It contains the approved 18 screens, three candidate universes, cached description vectors, frozen historical banks, released AssayFormer checkpoint, compact historical trajectories and genuine provider receipts. It also includes the pinned upstream source and its licenses. The bundle is separate from the small Git package. Download this bundle alongside the code to run the commands below.

Extract the bundle into the repository's ignored `_reproduction/` directory. From the repository root:

```sh
mkdir -p _reproduction
tar -xzf /path/to/assayloop-presentation-reproduction-2026-10-22.tar.gz -C _reproduction
```
 The expected layout is `_reproduction/bundle/`, `_reproduction/assayloop-src/` and `_reproduction/vendor/`. `SHA256.json` inside the input bundle checks every input; `RUNTIME.json` records tested package versions and upstream Python source hashes. The frozen banks avoid changes from recomputing an SVD on a different software stack. Each query publication is excluded from its historical donors.

Use Python 3.11 or newer with the versions in `requirements-reproduction.txt`. The measured execution environment is recorded in `RUNTIME.json`; numerical libraries on other platforms can affect exact tie ordering. The upstream source is from Genentech/AssayLoop commit `e7251afacda88c803864896e0d1444196f482471`, with AssayBench 0.2.0. Installation requires ordinary network access for dependencies.

```sh
python -m pip install -r experiments/odysseus/assayloop_comparison/requirements-reproduction.txt
export PYTHONPATH=.:_reproduction/assayloop-src:_reproduction/vendor
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
python -m experiments.odysseus.assayloop_comparison.reproduce \
  --bundle _reproduction/bundle --output _reproduction/replay --workers 2
python -m experiments.odysseus.assayloop_comparison.summarize_run \
  --summary _reproduction/replay/SUMMARY.json \
  --reference experiments/odysseus/assayloop_comparison/recorded_results/campaign_metrics.csv
```

The replay requires 342 campaigns and checks 864 original language requests. Saved selected candidates serve only as comparison targets. Each request is regenerated from the current campaign state and must match its receipt before the stored response is passed through the original parser. A changed history, prompt, request, parser output, selection or measured outcome fails verification. The replay makes zero provider calls.

## Fresh AWS replication

Provide `AWS_BEARER_TOKEN_BEDROCK` through your own secure environment. Check your account's access, quotas, prices and remaining spending allocation first. The runner uses Gemma4 31B through AWS Mantle and Qwen3 32B through Bedrock InvokeModel in us-east-1. It preserves temperature 1, an 8,192-token output limit, and the original language prompts/parser. Prices used for the verified run appear in the verification report.

```sh
python -m experiments.odysseus.assayloop_comparison.live \
  --bundle _reproduction/bundle --output _reproduction/fresh \
  --workers 2 --allow-paid --cap-usd 1
python -m experiments.odysseus.assayloop_comparison.summarize_run \
  --summary _reproduction/fresh/SUMMARY.json
```

The matrix is 18 screens × 2 actors × 2 repetitions × 2 systems = **144 complete campaigns**, normally 432 new language calls. Each campaign has three language decisions and seven numerical decisions, with 100 selections per decision. Each branch receives its own measured feedback and fresh subsequent language calls. Both systems select 1,000 distinct candidates. Available measured outcomes can differ, so retain the measured count, precision and EF alongside unique hits.

`incumbent` is the complete extension; `native_filled` preserves native proposals and fills their shortfalls using pre-batch AssayFormer predictions. Historical replay also covers original native, the same-proposal ordinary scorer, AssayFormer alone and the historical scorer alone. Attribution to the historical subspace specifically remains a separate ablation question.

The runner records requests, responses, trajectories, failures, usage and cost under the chosen output directory. A shared guard reserves costs for in-flight calls, retains uncertain billing, stops on provider failures and requires 10 GiB free disk. Use a new output directory for each attempt. Final `SUMMARY.json` distinguishes complete success from partial execution; retain failed attempts for audit.

All these cohorts have informed development. This is implementation reproduction and a repeat comparison on exposed data. Independent confirmation needs a reserved cohort and a frozen method.
