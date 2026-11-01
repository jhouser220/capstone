import json, numpy as np
B="/tmp/claude-1000/-home-jeff-Documents-Textbooks-DATASCI-210-capstone/2599c046-3406-48a5-93fb-83bd0a3f777c/scratchpad/"
for tag in ("af_seed3","af_seed7"):
    rows=json.loads(open(B+tag+".json").read())
    for lab,sub in (("all 20",rows),("12 non-family",[r for r in rows if not r["screen"].startswith("U_179")])):
        base=np.mean([r["base_pct"] for r in sub])
        path=np.mean([r["top_k_pct"] for r in sub]); af=np.mean([r["af_late_pct"] for r in sub])
        cand=np.mean([r["cand_pct"] for r in sub])
        print(f"{tag}  {lab:14s} base {base:5.2f}%  cand-pool {cand:5.2f}% ({cand/base:4.2f}x)  "
              f"path-top300 {path:5.2f}% ({path/base:4.2f}x)  AF8-10 {af:5.2f}% ({af/base:4.2f}x)")
