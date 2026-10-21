"""Absolute risks for the four main analyses (TriNetX exports).

Usage: python3 -I stats_absolute_risk.py <export root>

Prints a Markdown table; the main analyses are read from <export root>/<Main_Full|Female_Breast|Female_Endo_noHyst|Male>.

Per outcome and arm:
  - events;
  - person-years inside the outcome window (day 365 to the last recorded fact, death or day 3,650),
    interpolated from the 100-day follow-up distribution (lib_trinetx_export.window_person_years;
    follow-up is not truncated at the outcome), and incidence per 1,000 person-years;
  - cumulative incidence at 3 and 5 years from the index date = 1 - KM survival at day 1095 / 1825
    (the outcome window starts at day 365);
  - risk difference (SGLT2i - GLP-1RA) at 3 and 5 years, SE from the two KM CI half-widths.
Rounded once, ROUND_HALF_UP.
"""
import csv
import math
import re
import sys
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

from lib_trinetx_export import window_person_years

Z = 1.959963984540054
SOURCES = {
    "Main_Full": ["Kidney", "Bladder", "Colorectal", "Lung", "Liver", "Pancrea", "Esophag", "Thyroid", "Oral", "Death"],
    "Female_Breast": ["Breast", "Death"],
    "Female_Endo_noHyst": ["Endometrial", "Death"],
    "Male": ["Prostate", "Death"],
}


def q(x, p="0.01"):
    s = str(Decimal(str(x)).quantize(Decimal(p), rounding=ROUND_HALF_UP))
    return s[1:] if s.startswith("-") and set(s[1:]) <= set("0.") else s


def rows_after_header(path, startswith):
    rows = list(csv.reader(open(path, encoding="utf-8-sig", newline="")))
    i = next(k for k, r in enumerate(rows) if r and r[0].strip() == startswith or (r and startswith == "Time (Days)" and r[0].startswith(startswith)))
    return rows[i + 1:]


def km_at(path, day):
    data = [r for r in rows_after_header(path, "Time (Days)") if r and r[0].strip().isdigit()]
    last = None
    for r in data:
        if int(r[0]) <= day:
            last = r
        else:
            break
    s1, l1, u1, s2, l2, u2 = (float(x) for x in last[1:7])
    # cumulative incidence and CI (1 - survival; bounds swap)
    return (1 - s1, 1 - u1, 1 - l1), (1 - s2, 1 - u2, 1 - l2)


def main(root):
    root = Path(root)
    print("| Analysis | Outcome | Events S/G | Person-years S/G | Per 1,000 person-years S/G | 3-year cumulative incidence S/G (%) | 3-year risk difference (%) | 5-year cumulative incidence S/G (%) | 5-year risk difference (%) |")
    print("|---|---|---|---|---|---|---|---|---|")
    for sub, keys in SOURCES.items():
        d = next(p for p in sorted((root / sub).glob("IJC_*")) if p.is_dir())
        wpy = window_person_years(d)
        readme = (d / "_Read_Me.txt").read_text(encoding="utf-8-sig")
        names = {int(k): v.strip() for k, v in re.findall(r"Outcome (\d+) Name: (.+)", readme)}
        for key in keys:
            n = next((k for k, v in names.items() if key.lower() in v.lower()), None)
            if n is None or not (d / f"Outcome_{n}_Result_b_KM_table.csv").exists():
                sys.exit(f"missing outcome {key!r} in {d}")
            st = [r for r in rows_after_header(d / f"Outcome_{n}_Result_b_KM_table.csv", "Cohort") if r and r[0] in ("1", "2")]
            ev = {r[0]: int(r[3]) for r in st}
            py = {"1": wpy["S"], "2": wpy["G"]}
            ir = {c: ev[c] / py[c] * 1000 for c in ("1", "2")}
            out = []
            for day in (1095, 1825):
                a, b = km_at(d / f"Outcome_{n}_Result_b_KM_graph.csv", day)
                se1, se2 = (a[2] - a[1]) / (2 * Z), (b[2] - b[1]) / (2 * Z)
                rd = a[0] - b[0]
                se = math.sqrt(se1 ** 2 + se2 ** 2)
                out.append((f"{q(a[0] * 100)} / {q(b[0] * 100)}",
                            f"{q(rd * 100, '0.001')} ({q((rd - Z * se) * 100, '0.001')} to {q((rd + Z * se) * 100, '0.001')})"))
            print(f"| {sub} | {names[n]} | {ev['1']}/{ev['2']} | {q(py['1'], '1')}/{q(py['2'], '1')} | "
                  f"{q(ir['1'])} / {q(ir['2'])} | {out[0][0]} | {out[0][1]} | {out[1][0]} | {out[1][1]} |")


if __name__ == "__main__":
    main(sys.argv[1])
