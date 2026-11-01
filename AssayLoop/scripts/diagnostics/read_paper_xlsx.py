"""Read `paper_per_round_data.xlsx` — the AssayLoop paper's own per-round results.

Why this is hand-rolled. `openpyxl` is not in this project's environment, and an .xlsx is
a zip of XML, so reading it with the standard library costs less than a dependency.

What the workbook holds: 59 methods x 10 rounds on the paper's 20 held-out test screens,
which are EXACTLY this project's `paper_test` set — all 20 names match, and every screen's
total_hits and library_size agree with ours. So the paper's numbers and ours can be put on
one axis without re-normalising anything.

Sheets used here:
  Screens         per-screen total_hits, library_size (verified against ours)
  Mean_long       one row per method x round, averaged over screens
  By_screen_long  one row per method x screen x round — what the per-round plots need

A caveat that governs how EF is computed from this file. Our adjusted EF charges a budget
of `n_in_library + n_out_of_universe + n_unfilled`, forgiving picks that are real genes
outside the screen's own library. `By_screen_long` gives `n_acquired` and
`out_of_universe` per screen but NOT the forgiven count, which exists only as a
cross-screen mean in `Mean_long`. `ef_bounds()` therefore returns a bracket: the low end
forgives nothing, the high end forgives the method's mean. They differ by about 1.5%.
"""
from __future__ import annotations

import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
DEFAULT = "paper_per_round_data.xlsx"


def _sheets(path):
    z = zipfile.ZipFile(path)
    shared = ["".join(t.text or "" for t in si.iter(NS + "t"))
              for si in ET.fromstring(z.read("xl/sharedStrings.xml"))]
    names = [s.get("name") for s in ET.fromstring(z.read("xl/workbook.xml")).iter(NS + "sheet")]
    out = {}
    for i, nm in enumerate(names, start=1):
        member = f"xl/worksheets/sheet{i}.xml"
        if member not in z.namelist():
            continue
        rows = []
        for r in ET.fromstring(z.read(member)).iter(NS + "row"):
            cells = {}
            for c in r.iter(NS + "c"):
                v = c.find(NS + "v")
                if v is None:
                    continue
                key = re.match(r"[A-Z]+", c.get("r")).group(0)
                cells[key] = shared[int(v.text)] if c.get("t") == "s" else v.text
            rows.append(cells)
        if not rows:
            continue
        hdr = rows[0]
        out[nm] = [{hdr[c]: r.get(c) for c in hdr if c in r} for r in rows[1:]]
    return out


def load(path=DEFAULT):
    """{sheet name: [row dicts]}; values are strings, converted by the helpers below."""
    p = Path(path)
    if not p.is_file():
        raise SystemExit(f"{p} not found — the paper's per-round workbook is needed for "
                         f"these figures")
    return _sheets(p)


def per_round(sheets, display_names=None):
    """{display name: {screen: [cum_hits per round 1..10]}} plus the screen table.

    `By_screen_long` is keyed on the ORIGINAL method name, `Mean_long` carries the
    display name; this joins them so callers can ask for 'Gemini-3.1-pro' and get rows.
    """
    disp = {}
    for r in sheets["Mean_long"]:
        disp[r["method (original name)"]] = r["display name"]
    want = set(display_names) if display_names else None
    out, meta = {}, {}
    for r in sheets["By_screen_long"]:
        name = disp.get(r["method (original name)"], r["method (original name)"])
        if want and name not in want:
            continue
        rd = int(r["round"])
        d = out.setdefault(name, {}).setdefault(r["screen"], {})
        d[rd] = dict(cum_hits=float(r["cum_hits"]), n_acquired=float(r["n_acquired"]),
                     out_of_universe=float(r["out_of_universe"]),
                     total_hits=float(r["total_hits"]), library=float(r["library_size"]))
        meta[r["screen"]] = (float(r["total_hits"]), float(r["library_size"]))
    return out, meta


def forgiven_mean(sheets):
    """{display name: {round: mean in-universe-but-out-of-library picks}} — the only
    place the forgiven count exists, and only as a cross-screen mean."""
    out = {}
    for r in sheets["Mean_long"]:
        v = r.get("in_universe_not_lib")
        if v is None:
            continue
        out.setdefault(r["display name"], {})[int(r["round"])] = float(v)
    return out


def ef_bounds(entry, forgiven=0.0, budget_per_round=100):
    """(low, high) adjusted EF for one screen-round, bracketing the unknown forgiven count.

    low forgives nothing (charges the full budget); high forgives `forgiven` picks. Our
    own runs compute the exact value, which lies between the two.
    """
    b = budget_per_round * max(1, int(round(entry["n_acquired"] / budget_per_round)))
    H, D = entry["total_hits"], entry["library"]
    if H <= 0:
        return float("nan"), float("nan")
    lo = entry["cum_hits"] / (b * H / D)
    hi = entry["cum_hits"] / (max(b - forgiven, 1.0) * H / D)
    return lo, hi
