"""Sort output/runs' flat sweep directories into one folder per METHOD TESTED.

Run it again after a new sweep to file that sweep; it is idempotent and skips a run
directory whose destination already exists. `--apply` moves, no flag dry-runs.

    python3 -I scripts/organize_runs.py            # show the plan
    python3 -I scripts/organize_runs.py --apply    # move

Folders are named by the components the arm actually fused, in the project's own
vocabulary (HANDOFF section 9 and the notebook registry): `af` is AssayFormer, `depmap`
is DepMap co-essentiality via codep_seed, `screenknn` is the screen-kNN scorer. A reader
should be able to tell what an arm tested from the folder name alone.

Classification reads `model`, `acq`, `screen_set` and the two override blocks of one
representative run per sweep -- never the `lm` block, which is where API keys live.

NAME_OVERRIDE exists because four families are NOT distinguishable from config: their
intervention arrived by environment variable (`sweep-desc-*`, ASSAYLOOP_DESC_EXTRA) or by
an injected seed list written straight into the campaign (`sweep-orc-*`, `sweep-inj-*`,
the LLM `sweep-*hyb*` / `sweep-hybrid-*` handoffs). Those record empty `acq_overrides`
and would otherwise be filed as plain AssayFormer, which is what they are mechanically
and not what they test. Their names come from HANDOFF section 9.

Full-atlas arms get a `_full_atlas_superseded` folder of their own. They are the
`achilles.npz` originals of arms re-run on `achilles_noeval.npz`; the project's rule is
that the reduced atlas is the standard and the two are never pooled, so the layout makes
pooling them take a deliberate act.
"""
import json, os, re, shutil, sys
from collections import defaultdict
from pathlib import Path

RUNS = Path("/home/jeff/Documents/Textbooks/DATASCI_210/capstone/AssayLoop/output/runs")
KEEP_AS_IS = {"baseline", "external", "old", "analysis_tables"}
RUN_DIR = re.compile(r"^.+-[0-9]+-.+$")

# sweep-id prefix -> folder, for the families whose method is not in the config
NAME_OVERRIDE = [
    ("sweep-desc-",        "af_description_channel"),
    ("sweep-orc-",         "af_depmap_oracle_nominator"),
    ("sweep-inj-",         "af_depmap_injected_nominator"),
    ("sweep-hybrid-",      "llm_warmstart_to_af_depmap"),
    ("sweep-dmne-hyb-",    "llm_warmstart_to_af_depmap"),
    ("sweep-dmcd-hyb-",    "llm_warmstart_to_af_depmap"),
    ("sweep-glm5tu-pt-af", "llm_topup_warmstart_to_af_alone"),
    ("sweep-glm5testtopup","topup_llm_only"),
    ("sweep-codep-smoke",  "_smoke_tests"),
]


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return 0.0


def classify(sid, c):
    for pre, name in NAME_OVERRIDE:
        if sid.startswith(pre):
            return name, False

    ao, mo = c.get("acq_overrides") or {}, c.get("model_overrides") or {}
    if c.get("acq") == "llm_single_topup":
        return "topup_llm_only", False
    if c.get("model") == "null":
        return "_smoke_tests", False
    if c.get("model") == "knn":
        return "paper_geneknn_alone", False

    if any(k.startswith("ctx_") for k in mo) or "ignore_context" in mo:
        return "af_context_ablations", False
    if "score_mode" in mo:
        return "af_score_mode_ablation", False
    if "prox_level" in ao or "gene_locus" in ao:
        return "af_depmap_screenknn_proximity_correction", False

    dm, knn = num(ao.get("weight")), num(ao.get("knn_weight"))
    atlas = os.path.basename(str(ao.get("profile_npz") or ""))
    extras = [n for n, k in (("cohit", "weight2"), ("bioplex_ppi", "ppi_weight"),
                             ("genept", "gpt_weight"), ("paper_geneknn", "pknn_weight"))
              if num(ao.get(k))]

    if atlas == "Expression.npz" and dm:
        base, full = "af_expression_profile", False
    else:
        parts = ["af"]
        if dm:
            parts.append("depmap")
        if knn:
            parts.append("screenknn")
        parts += extras
        base = "_".join(parts)
        if base == "af":
            base = "af_alone_baseline"
        full = bool(dm) and atlas == "achilles.npz"
    return base, full


def main(apply):
    members, reps = defaultdict(list), {}
    for d in sorted(os.listdir(RUNS)):
        p = RUNS / d
        # Only a RUN directory is ever filed: `<sweep>-<NN>-<screen>` holding a
        # config.json. Anything else at this level is a method folder this script
        # created, or `baseline/`, `external/` and friends -- re-running must leave
        # all of them alone, which is what makes it idempotent.
        if d in KEEP_AS_IS or not p.is_dir():
            continue
        if not RUN_DIR.match(d) or not (p / "config.json").is_file():
            continue
        sid = re.sub(r"-[0-9]+-.*$", "", d)
        members[sid].append(d)
        if sid not in reps:
            reps[sid] = p

    plan, folders = {}, defaultdict(list)
    for sid, dirs in sorted(members.items()):
        rep = reps.get(sid)
        if rep is None:
            folder = "_unclassified"
        else:
            c = json.load(open(rep / "config.json"))["config"]
            name, full = classify(sid, c)
            folder = f"{name}__full_atlas_superseded" if full else name
        plan[sid] = folder
        folders[folder].append(sid)

    total = sum(len(members[s]) for s in plan)
    print(f"{len(plan)} sweep ids / {total} run dirs -> {len(folders)} method folders\n")
    for f in sorted(folders):
        n = sum(len(members[s]) for s in folders[f])
        print(f"  {f:52s} {len(folders[f]):4d} sweeps {n:6d} dirs")

    if os.environ.get("SC"):
        json.dump({s: plan[s] for s in sorted(plan)},
                  open(os.environ["SC"] + "/final_plan.json", "w"), indent=1)
    if not apply:
        print("\n(dry run)")
        return

    moved = 0
    for sid, folder in plan.items():
        dest = RUNS / folder
        dest.mkdir(exist_ok=True)
        for d in members[sid]:
            src, dst = RUNS / d, dest / d
            if dst.exists():
                print(f"  SKIP exists: {folder}/{d}")
                continue
            shutil.move(str(src), str(dst))
            moved += 1
    print(f"\nmoved {moved} run directories")


if __name__ == "__main__":
    main(apply="--apply" in sys.argv)
