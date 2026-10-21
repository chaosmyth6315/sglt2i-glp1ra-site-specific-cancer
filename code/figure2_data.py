"""Figure 2 data: one row per site-specific cancer, read from the four primary-analysis exports.

Usage: python3 -I figure2_data.py <export root> <out.csv>

Input layout: <export root>/<subfolder>/IJC_*/ for the subfolders Main_Full, Female_Breast,
Female_Endo_noHyst and Male; each export folder holds _Read_Me.txt (outcome names) and
Outcome_<k>_Result_b_KM_table.csv.

Columns: block, label, n_per_arm, ev1, ev2, cir1, cir2 (100 × (1 − KM survival at the end of the
window), %), hr, lo, hi (export values, 7 dp), wald_p, bh_q (Benjamini–Hochberg, across the 12 cancers),
evalue, evalue_ci, ph_p (Schoenfeld), dagger (ph_p < 0.05), src and outcome_no (export folder and
outcome number), and the display strings *_s.
Wald p from the HR and its 95% CI: SE = (ln hi − ln lo) / (2 × 1.96), p = erfc(|ln HR / SE| / √2).
E-value (VanderWeele & Ding 2017) from the unrounded HR and the CI bound nearer the null;
1.00 when the CI includes 1. Display strings are rounded once, here, with ROUND_HALF_UP; a CI
bound is never printed as 1.00 when the CI excludes 1 (fmt_bound). figure2_forest.R only lays them out.
"""
import csv
import math
import re
import sys
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

Z = 1.959963984540054
# (block label, export subfolder, [(outcome key, display label)]); the key is matched
# case-insensitively against the outcome names in _Read_Me.txt
BLOCKS = [
    ("All patients", "Main_Full", [
        ("Kidney", "Kidney cancer"), ("Bladder", "Bladder cancer"), ("Colorectal", "Colorectal cancer"),
        ("Lung", "Lung cancer"), ("Liver", "Liver and intrahepatic bile duct cancer"),
        ("Pancrea", "Pancreatic cancer"), ("Esophag", "Esophageal cancer"), ("Thyroid", "Thyroid cancer"),
        ("Oral", "Oral, oropharyngeal and hypopharyngeal cancer")]),
    ("Women", "Female_Breast", [("Breast", "Breast cancer")]),
    ("Women without prior hysterectomy", "Female_Endo_noHyst", [("Endometrial", "Endometrial cancer")]),
    ("Men", "Male", [("Prostate", "Prostate cancer")]),
]


def sections(path):
    out, cur = {}, None
    with open(path, encoding="utf-8-sig", newline="") as f:
        for row in csv.reader(f):
            cells = [c.strip() for c in row]
            if not any(cells) or cells in (['" "'], [" "]):
                continue
            if len(cells) == 1:
                cur = cells[0]; out[cur] = []
            elif cur is not None:
                out[cur].append(cells)
    return out


def q2(x, places="0.01"):
    return str(Decimal(str(x)).quantize(Decimal(places), rounding=ROUND_HALF_UP))


def fmt_bound(x, lo, hi):
    """2 dp. A bound that would print as 1.00 while the CI excludes 1 is moved one unit away from
    the null (0.9972 -> 0.99; 1.0035 -> 1.01) so that the display never contradicts the p value."""
    s = q2(x)
    if s == "1.00" and not (lo <= 1 <= hi):
        s = q2(Decimal("0.99")) if x < 1 else q2(Decimal("1.01"))
    return s


def fmt_p(p):
    if p < 0.001:
        return "<0.001"
    s = q2(p, "0.001")
    return q2(p, "0.0001") if (p < 0.05 and s == "0.050") else s


def evalue(rr):
    rr = 1 / rr if rr < 1 else rr
    return rr + math.sqrt(rr * (rr - 1))


def main(root, out):
    root = Path(root)
    rows = []
    for block, sub, keys in BLOCKS:
        d = next(p for p in sorted((root / sub).glob("IJC_*")) if p.is_dir())
        names = {int(k): v.strip() for k, v in re.findall(r"Outcome (\d+) Name: (.+)", (d / "_Read_Me.txt").read_text(encoding="utf-8-sig"))}
        for key, label in keys:
            k = next(i for i, v in names.items() if key.lower() in v.lower())
            km = sections(d / f"Outcome_{k}_Result_b_KM_table.csv")
            st = km["Cohort Statistics"][1:3]
            hr, lo, hi = (float(x) for x in km["Hazard Ratio"][1][:3])
            ph_p = float(km["Proportionality"][1][2])
            b, se = math.log(hr), (math.log(hi) - math.log(lo)) / (2 * Z)
            p = math.erfc(abs(b / se) / math.sqrt(2))
            ci_near = hi if hr < 1 else lo
            ev_ci = 1.0 if lo <= 1 <= hi else evalue(ci_near)
            rows.append(dict(block=block, label=label, n_per_arm=int(st[0][2]),
                             ev1=int(st[0][3]), ev2=int(st[1][3]),
                             cir1=100 * (1 - float(st[0][5])), cir2=100 * (1 - float(st[1][5])),
                             hr=hr, lo=lo, hi=hi, wald_p=p, evalue=evalue(hr), evalue_ci=ev_ci,
                             ph_p=ph_p, dagger=int(ph_p < 0.05), src=d.name, outcome_no=k))
    m = len(rows)
    order = sorted(range(m), key=lambda i: rows[i]["wald_p"])
    prev = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        prev = min(prev, rows[i]["wald_p"] * m / rank)
        rows[i]["bh_q"] = prev
    for r in rows:
        # en dash between the CI bounds; † (proportional-hazards violation, ph_p < 0.05) follows the HR
        r["hr_ci_s"] = (f"{q2(r['hr'])} ({fmt_bound(r['lo'], r['lo'], r['hi'])}–{fmt_bound(r['hi'], r['lo'], r['hi'])})"
                        + ("†" if r["dagger"] else ""))
        r["p_s"] = fmt_p(r["wald_p"])
        r["q_s"] = fmt_p(r["bh_q"])
        r["ev_s"] = f"{q2(r['evalue'])} ({q2(r['evalue_ci'])})"
        r["ev1_s"] = f"{r['ev1']:,} ({q2(r['cir1'])})"
        r["ev2_s"] = f"{r['ev2']:,} ({q2(r['cir2'])})"
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    for r in rows:
        print(f"{r['block'][:5]:5} {r['label'][:22]:22} n={r['n_per_arm']} ev={r['ev1']}/{r['ev2']} "
              f"CIR={r['cir1']:.2f}/{r['cir2']:.2f} HR={r['hr']:.4f} ({r['lo']:.4f}-{r['hi']:.4f}) "
              f"p={r['wald_p']:.4g} q={r['bh_q']:.4g} E={r['evalue']:.2f}({r['evalue_ci']:.2f}) PH={r['ph_p']:.3f}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
