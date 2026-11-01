"""Is there a DEPLOYABLE gate? Does AssayFormer's own observed rounds 5-7 hit rate
predict the screens where pathway propagation beats it in rounds 8-10?
Uses only information available at the end of round 7."""
import json, pathlib
import numpy as np

RUNS = pathlib.Path("output/runs")
rows = {r["screen"]: r for r in json.loads(
    open("/tmp/claude-1000/-home-jeff-Documents-Textbooks-DATASCI-210-capstone/"
         "2599c046-3406-48a5-93fb-83bd0a3f777c/scratchpad/af_seed7.json").read())}

out = []
for d in sorted(RUNS.glob("*/sweep-afalone-val-*")):
    res = json.loads((d / "result.json").read_text())
    name = json.loads((d / "config.json").read_text())["screen_name"]
    if name not in rows:
        continue
    steps = res["steps"]
    r57 = [h for i in (4, 5, 6) for h in steps[i]["hits"]]          # rounds 5-7, observable
    r13 = [h for i in (0, 1, 2) for h in steps[i]["hits"]]
    r = rows[name]
    out.append(dict(screen=name, af_r13=100 * np.mean(r13), af_r57=100 * np.mean(r57),
                    af_late=r["af_late_pct"], path=r["top_k_pct"],
                    delta=r["top_k_pct"] - r["af_late_pct"],
                    extra=r["top_k_hits"] - r["af_late_hits"]))

out.sort(key=lambda r: -r["delta"])
print(f"{'screen':16s} {'AF r1-3%':>8s} {'AF r5-7%':>8s} {'AF r8-10%':>9s} {'path300%':>8s} {'delta':>7s} {'+hits':>6s}")
print("-" * 70)
for r in out:
    print(f"{r['screen']:16s} {r['af_r13']:8.2f} {r['af_r57']:8.2f} {r['af_late']:9.2f} "
          f"{r['path']:8.2f} {r['delta']:+7.2f} {r['extra']:+6d}")

d = np.array([r["delta"] for r in out]); g = np.array([r["af_r57"] for r in out])
print(f"\nSpearman(AF rounds 5-7 hit rate, pathway advantage) = "
      f"{np.corrcoef(np.argsort(np.argsort(g)), np.argsort(np.argsort(d)))[0,1]:+.2f}")
for thr in (1, 2, 3, 5, 10):
    sel = [r for r in out if r["af_r57"] < thr]
    if sel:
        print(f"  gate AF r5-7 < {thr:2d}%: {len(sel):2d} screens, "
              f"net extra hits {sum(r['extra'] for r in sel):+4d}, "
              f"wins {sum(1 for r in sel if r['extra']>0)}/{len(sel)}")
