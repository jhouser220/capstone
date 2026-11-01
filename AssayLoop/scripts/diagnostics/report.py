import json, sys
import numpy as np

label = sys.argv[2] if len(sys.argv) > 2 else ""
rows = json.loads(open(sys.argv[1]).read())
EASY = lambda n: n.startswith("U_179")          # the near-duplicate fitness family

def lift(r, key="cand_pct", base="base_pct"):
    return r[key] / r[base] if r[base] else float("nan")

print(f"\n=== {label}  ({len(rows)} screens) ===")
hdr = f"{'screen':16s} {'seeds':>5s} {'hot':>5s} {'cands':>6s} {'base%':>6s} {'cand%':>6s} {'lift':>5s} {'top300%':>7s} {'AF8-10%':>7s} {'win':>4s} {'ovlp':>4s}"
print(hdr); print("-" * len(hdr))
for r in sorted(rows, key=lambda r: -lift(r)):
    win = "yes" if r["top_k_pct"] > r["af_late_pct"] else "no"
    print(f"{r['screen']:16s} {r['n_seeds']:5d} {r['n_hot']:5d} {r['n_cands']:6d} "
          f"{r['base_pct']:6.2f} {r['cand_pct']:6.2f} {lift(r):5.2f} "
          f"{r['top_k_pct']:7.2f} {r['af_late_pct']:7.2f} {win:>4s} {r['overlap_af']:4d}")

def agg(sub, name):
    if not sub: return
    L = np.array([lift(r) for r in sub])
    LS = np.array([lift(r, "cand_strong_pct", "base_strong_pct") for r in sub])
    tk = np.array([r["top_k_pct"] for r in sub]); af = np.array([r["af_late_pct"] for r in sub])
    print(f"\n{name} (n={len(sub)}):")
    print(f"  lift, all hits      mean {L.mean():.2f}x  median {np.median(L):.2f}x  >1.5x: {(L>1.5).sum()}/{len(L)}")
    print(f"  lift, strong hits   mean {LS.mean():.2f}x  median {np.median(LS):.2f}x  >1.5x: {(LS>1.5).sum()}/{len(LS)}")
    print(f"  top-300 hit rate    {tk.mean():.2f}%   AF rounds 8-10  {af.mean():.2f}%   "
          f"pathway wins {(tk>af).sum()}/{len(tk)}")
    print(f"  total hits in 300   pathway {sum(r['top_k_hits'] for r in sub)}   "
          f"AF {sum(r['af_late_hits'] for r in sub)}")

agg(rows, "ALL")
agg([r for r in rows if not EASY(r["screen"])], "excluding the U_179x fitness family")
agg([r for r in rows if EASY(r["screen"])], "the U_179x fitness family only")
