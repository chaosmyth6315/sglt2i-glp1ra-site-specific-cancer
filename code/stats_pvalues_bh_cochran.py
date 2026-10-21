"""Main-analysis statistics from the four TriNetX main-analysis exports.

Usage: python3 -I stats_pvalues_bh_cochran.py <export root>

For the 12 site-specific cancers (9 from Main_Full, breast from Female_Breast,
endometrial from Female_Endo_noHyst, prostate from Male):
  - Wald p from the Cox HR and its 95% CI (SE = (ln U - ln L) / (2 * 1.959964)),
  - Benjamini-Hochberg q across the 12 Wald p values,
  - Cochran's Q and I^2 on the log HRs across the 12 sites.
Also prints deaths by arm in each of the four analyses.
HR/CI are rounded once, ROUND_HALF_UP, from the export values.
"""
import csv
import math
import re
import sys
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

Z = 1.959963984540054
SOURCES = {
    "Main_Full": ["Lung", "Colorectal", "Liver", "Kidney", "Pancrea", "Esophag", "Thyroid", "Oral", "Bladder"],
    "Female_Breast": ["Breast"],
    "Female_Endo_noHyst": ["Endometrial"],
    "Male": ["Prostate"],
}


def q(x, p="0.01"):
    return str(Decimal(str(x)).quantize(Decimal(p), rounding=ROUND_HALF_UP))


# HR (95% CI) display: a CI bound that would print as 1.00 while the CI excludes 1 is shown as 0.99 or
# 1.01, as in the figures.
def fmt_hr(hr, lo, hi):
    def bound(x):
        s = q(x)
        return ("0.99" if x < 1 else "1.01") if s == "1.00" and not (lo <= 1 <= hi) else s
    return f"{q(hr)} ({bound(lo)}–{bound(hi)})"


def sections(path):
    out, cur = {}, None
    with open(path, encoding="utf-8-sig", newline="") as f:
        for row in csv.reader(f):
            cells = [c.strip() for c in row]
            if not any(cells) or cells in (['" "'], [" "]):
                continue
            if len(cells) == 1:
                cur = cells[0]
                out[cur] = []
            elif cur is not None:
                out[cur].append(cells)
    return out


def norm_sf(z):
    return 0.5 * math.erfc(z / math.sqrt(2))


def read_export(d):
    readme = (d / "_Read_Me.txt").read_text(encoding="utf-8-sig")
    names = {int(k): v.strip() for k, v in re.findall(r"Outcome (\d+) Name: (.+)", readme)}
    res = {}
    for n, nm in names.items():
        kmf = d / f"Outcome_{n}_Result_b_KM_table.csv"
        if not kmf.exists():
            continue
        km = sections(kmf)
        st = km["Cohort Statistics"][1:]
        hr = km["Hazard Ratio"][1]
        res[nm] = {
            "n": (st[0][2], st[1][2]),
            "ev": (st[0][3], st[1][3]),
            "hr": float(hr[0]), "lo": float(hr[1]), "hi": float(hr[2]),
            "hr_s": hr[0], "lo_s": hr[1], "hi_s": hr[2],
        }
    return res


def main(root):
    root = Path(root)
    rows, deaths = [], []
    for sub, keys in SOURCES.items():
        d = next(p for p in sorted((root / sub).glob("IJC_*")) if p.is_dir())
        res = read_export(d)
        for nm, r in res.items():
            if nm.lower().startswith("death"):
                deaths.append((sub, r))
            elif any(k.lower() in nm.lower() for k in keys):
                b = math.log(r["hr"])
                se = (math.log(r["hi"]) - math.log(r["lo"])) / (2 * Z)
                rows.append({"src": sub, "name": nm, "b": b, "se": se, "p": 2 * norm_sf(abs(b / se)), **r})
        for k in keys:   # every site must be present exactly once, otherwise the BH family changes
            hits = [nm for nm in res if k.lower() in nm.lower()]
            if len(hits) != 1:
                sys.exit(f"{sub}: expected one outcome matching {k!r}, found {hits}")
    # Benjamini-Hochberg
    m = len(rows)
    assert m == sum(len(k) for k in SOURCES.values()), m
    order = sorted(range(m), key=lambda i: rows[i]["p"])
    qv = [0.0] * m
    prev = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        prev = min(prev, rows[i]["p"] * m / rank)
        qv[i] = prev
    # Cochran's Q
    w = [1 / r["se"] ** 2 for r in rows]
    bbar = sum(wi * r["b"] for wi, r in zip(w, rows)) / sum(w)
    Q = sum(wi * (r["b"] - bbar) ** 2 for wi, r in zip(w, rows))
    df = m - 1
    # chi-square upper tail via regularized gamma (series/continued fraction)
    pQ = chi2_sf(Q, df)
    I2 = max(0.0, (Q - df) / Q) * 100 if Q > 0 else 0.0

    print("| Cancer | Source | Events S/G | HR (95% CI) | Wald p | BH q |")
    print("|---|---|---|---|---|---|")
    for i in sorted(range(m), key=lambda i: rows[i]["p"]):
        r = rows[i]
        print(f"| {r['name']} | {r['src']} | {r['ev'][0]}/{r['ev'][1]} | "
              f"{fmt_hr(r['hr'], r['lo'], r['hi'])} | {fmt_p(r['p'])} | {fmt_p(qv[i])} |")
    print()
    print(f"Cochran's Q across {m} sites = {q(Q)}, df = {df}, p = {fmt_p(pQ)}, I² = {q(I2, '0.1')}%")
    print()
    print("| Analysis | Matched n per arm | Deaths SGLT2i | Deaths GLP-1RA | Death HR (95% CI) |")
    print("|---|---|---|---|---|")
    for sub, r in deaths:
        print(f"| {sub} | {r['n'][0]} | {r['ev'][0]} | {r['ev'][1]} | {fmt_hr(r['hr'], r['lo'], r['hi'])} |")


def fmt_p(p):
    if p < 0.001:
        return "<0.001"
    s = q(p, "0.001")
    if p < 0.05 and s == "0.050":
        return q(p, "0.0001")
    return s


def chi2_sf(x, k):
    a, xx = k / 2.0, x / 2.0
    if xx <= 0:
        return 1.0
    gln = math.lgamma(a)
    if xx < a + 1:
        ap, s, d = a, 1.0 / a, 1.0 / a
        for _ in range(500):
            ap += 1
            d *= xx / ap
            s += d
            if abs(d) < abs(s) * 1e-15:
                break
        return 1 - s * math.exp(-xx + a * math.log(xx) - gln)
    b = xx + 1 - a
    c, d = 1e300, 1 / b
    h = d
    for i in range(1, 500):
        an = -i * (i - a)
        b += 2
        d = an * d + b
        d = 1e-300 if abs(d) < 1e-300 else d
        c = b + an / c
        c = 1e-300 if abs(c) < 1e-300 else c
        d = 1 / d
        de = d * c
        h *= de
        if abs(de - 1) < 1e-15:
            break
    return math.exp(-xx + a * math.log(xx) - gln) * h


if __name__ == "__main__":
    main(sys.argv[1])
