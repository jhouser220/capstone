"""Compare the HSA arms with bank.mean kept vs zeroed, pick-for-pick and on EF."""
import json, statistics as st, sys

S = "output/runs/external"
import os
on  = {(r["arm"], r["screen"]): r for r in json.load(open(f"{S}/zeromean_meankept.json"))}
off = {(r["arm"], r["screen"]): r for r in json.load(open(f"{S}/zeromean_meanzeroed.json"))}
assert set(on) == set(off), "arm/screen sets differ"

for arm in sorted({k[0] for k in on}):
    keys = sorted(k for k in on if k[0] == arm)
    print(f"\n=== {arm} ===  ({len(keys)} screens)")
    print(f"{'screen':<15}{'hits(mean)':>11}{'hits(zero)':>11}{'d_hits':>8}"
          f"{'EF(mean)':>10}{'EF(zero)':>10}{'dEF':>9}{'picks same':>12}{'1st diff':>9}")
    dh, de, ident, ovl = [], [], 0, []
    for k in keys:
        a, b = on[k], off[k]
        pa, pb = a["pick_list"], b["pick_list"]
        same = pa == pb
        ident += same
        o = len(set(pa) & set(pb)); ovl.append(o)
        d = [i for i, (x, y) in enumerate(zip(pa, pb)) if x != y]
        dh.append(b["hits"] - a["hits"]); de.append(b["ef"] - a["ef"])
        print(f"{k[1]:<15}{a['hits']:>11}{b['hits']:>11}{b['hits']-a['hits']:>8}"
              f"{a['ef']:>10.4f}{b['ef']:>10.4f}{b['ef']-a['ef']:>+9.4f}"
              f"{('YES' if same else f'{o}/1000'):>12}{(d[0] if d else '-'):>9}")
    print(f"\n  identical pick sequences : {ident} of {len(keys)}")
    print(f"  median pick overlap      : {st.median(ovl):.0f} / 1000")
    print(f"  d_hits   total {sum(dh):+d}   mean {st.mean(dh):+.3f}   max |.| {max(map(abs,dh))}")
    print(f"  dEF      mean {st.mean(de):+.5f}   median {st.median(de):+.5f}   max |.| {max(map(abs,de)):.5f}")
    if len(de) > 1:
        print(f"  dEF      sd   {st.stdev(de):.5f}")
