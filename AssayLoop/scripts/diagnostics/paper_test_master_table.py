"""Every arm ever run on paper_test, one metric convention. EF / hit rate / total hits.

Canonical adjusted EF (HANDOFF 24.29) throughout:  EF = hits / (eff * H / D),
D = the screen's own library, eff = 1000 - (off-library but in-universe picks).
Full-atlas arms are EXCLUDED by this project's standing rule.
"""
import json, glob, sys, statistics as st
from pathlib import Path
sys.path.insert(0, "src")
from assayloop.metrics.hits_auc import adjusted_ef_value
from assayloop.tasks import gene_universe, load_screens

screens = load_screens(target_set="paper_test")
UNIV = set(gene_universe(screens, min_screen_freq=2))
BY = {s.dataset_name: s for s in screens}
H_AVAIL = sum(s.total_hits for s in screens)
BUDGET = 1000

def score_run(d, s):
    rd = json.loads((Path(d) / "result.json").read_text())
    lib = set(s.genes); n1 = n2 = hits = 0
    for step in rd.get("steps", []):
        for g, h in zip(step.get("acquired_batch", []), step.get("hits", [])):
            if g in lib:
                n1 += 1
                if h: hits += 1
            elif g in UNIV:
                n2 += 1
    ef = adjusted_ef_value(hits, n1, n2, len(lib), s.total_hits, BUDGET)
    return ef, hits, n1, n2

def arm_from_sweep(folder, sid):
    out = []
    for d in sorted(glob.glob(f"output/runs/{folder}/{sid}-*-*")):
        name = json.loads((Path(d) / "config.json").read_text())["screen_name"]
        if name not in BY: return None
        out.append(score_run(d, BY[name]))
    return out if len(out) == 20 else None

SWEEPS = [
 ("af_alone_baseline","sweep-tf-pt-w0","AF alone"),
 ("af_depmap","sweep-dmne-cd-pt","AF + DepMap (dm=0.5)"),
 ("af_depmap_screenknn","sweep-nolk-pt","AF + DepMap + screen-kNN"),
 ("af_depmap_screenknn_bioplex_ppi","sweep-dmne-ppi-pt","AF + DepMap + screen-kNN + BioPlex PPI"),
 ("af_depmap_screenknn_bioplex_ppi","sweep-dmne-ppidr-pt","AF + DepMap + screen-kNN + BioPlex (dir)"),
 ("af_depmap_screenknn_proximity_correction","sweep-prox-pt","AF + DepMap + screen-kNN + proximity corr"),
 ("af_screenknn","sweep-tf-pt-knn05","AF + screen-kNN only"),
 ("af_genept","sweep-gptonly-pt","AF + GenePT only"),
 ("af_paper_geneknn","sweep-pk-only-pt","AF + paper gene-kNN only"),
 ("paper_geneknn_alone","sweep-paperknn-pt","paper gene-kNN alone (no AF score)"),
 ("af_context_ablations","sweep-ctxtest-noctx","AF, no screen context"),
 ("af_context_ablations","sweep-ctxtest-drop50","AF, 50% context dropped"),
 ("af_score_mode_ablation","sweep-ctxtest-nobias","AF, gene bias removed"),
 ("af_score_mode_ablation","sweep-ctxtest-biasonly","AF, gene bias only"),
 ("af_description_channel","sweep-desc-pt-all","AF + description channel (all)"),
 ("af_description_channel","sweep-desc-pt-dir","AF + description channel (direction)"),
 ("af_description_channel","sweep-desc-pt-notes","AF + description channel (notes)"),
 ("llm_topup_warmstart_to_af_alone","sweep-glm5tu-pt-af","GLM-5 +top-up warm start -> AF"),
 # `dmne-hyb-*` are the TWO-TERM receiver (DepMap + screen-kNN), per HANDOFF 4.61:
 # "its twin sweep-dmne-hyb-glm5tu-* hands the same warm start to the two-term method".
 # `dmcd-hyb-*` is the one-term AF + DepMap receiver. These were mislabelled here.
 ("llm_warmstart_to_af_depmap","sweep-dmne-hyb-glm5-pt","GLM-5 warm start -> AF + DepMap + screen-kNN"),
 ("llm_warmstart_to_af_depmap","sweep-dmne-hyb-glm5tu-pt","GLM-5 +top-up -> AF + DepMap + screen-kNN"),
 ("llm_warmstart_to_af_depmap","sweep-dmcd-hyb-glm5tu-pt","GLM-5 +top-up -> AF + DepMap"),
 ("llm_warmstart_to_af_depmap","sweep-hybrid-pt-gem35codep","Flash-Lite warm start -> AF + DepMap"),
 ("topup_llm_only","sweep-glm5testtopup","GLM-5 + top-up, LLM all 10 rounds"),
]

rows = []
for folder, sid, label in SWEEPS:
    r = arm_from_sweep(folder, sid)
    if r is None:
        print(f"  !! skipped {sid}", file=sys.stderr); continue
    rows.append((label, sid, r))

# the HSA arms, from the external JSONs (eff == n_measured there: all picks in-universe, budget filled)
HSA_FILES = {"paper_test_hsa.json": None, "paper_test_hsa_noaf.json": None,
             "paper_test_hsa_realmean.json": " [real mean]"}
for f, tag in HSA_FILES.items():
    data = json.load(open("output/runs/external/" + f))
    arms = {}
    for r in data: arms.setdefault(r["arm"], []).append(r)
    for a, rs in arms.items():
        if a in ("AF alone", "AF + DepMap"): continue      # already have the sweep versions
        per = []
        for r in rs:
            assert r["picks"] == 1000 and r["n_measured"] + r["n_unknown"] == 1000
            per.append((r["hits"] / (r["n_measured"] * r["total_hits"] / r["n_genes"]),
                        r["hits"], r["n_measured"], r["n_unknown"]))
        rows.append((a + (tag or ""), f.replace(".json", ""), per))

print(f"paper_test: 20 screens, {H_AVAIL:,} hits available, pool {len(UNIV):,} genes, "
      f"budget {BUDGET}\nmetric: canonical adjusted EF. Full-atlas arms excluded.\n")
print(f"| arm | EF | hits | % of available | hit rate /1000 slots | hit rate /in-library |")
print(f"|---|---|---|---|---|---|")
for label, sid, per in sorted(rows, key=lambda x: -st.mean(p[0] for p in x[2])):
    ef = st.mean(p[0] for p in per); hits = sum(p[1] for p in per)
    n1 = sum(p[2] for p in per)
    print(f"| {label} | {ef:.3f} | {hits:,} | {hits/H_AVAIL:.1%} | {hits/(20*BUDGET):.2%} | {hits/n1:.2%} |")
