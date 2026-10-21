"""Manifest and number-at-risk tables for Figure 3 and Figure S3.

Writes <outdir>/km_manifest.csv (figure, panel, label, cohort, km_csv path, log-rank p string, HR string)
and <outdir>/km_risk.csv (label, arm, year, n at risk; from lib_km_number_at_risk) for the 12 sites.
Figure 3: kidney, bladder (full cohort), breast (women), prostate (men).
Figure S3: colorectal, lung, liver, pancreatic, esophageal, thyroid, oral (full), endometrial (women
without prior hysterectomy).
Usage: python3 -P figure3_s3_km_data.py <export root> <outdir>
"""
import csv
import sys
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

from lib_trinetx_export import export_dir, find, read_outcomes
from lib_km_number_at_risk import risk_table


# Round half up to the given number of places.
def q(x, places="0.01"):
    return Decimal(str(x)).quantize(Decimal(places), rounding=ROUND_HALF_UP)


# P values: "<0.001" below 0.001, otherwise 3 decimals; a P value below 0.05 that would print
# as 0.050 is shown with 4 decimals.
def fmt_p(p):
    if p < 0.001:
        return "<0.001"
    s = str(q(p, "0.001"))
    return str(q(p, "0.0001")) if (p < 0.05 and s == "0.050") else s


def fmt_bound(x, lo, hi):
    """Format a CI bound to 2 decimals. A bound that would print as 1.00 while the CI (lo, hi)
    excludes 1 is printed as 0.99 or 1.01, so a CI that excludes 1 never shows 1.00."""
    s = str(q(x))
    if s == "1.00" and not (lo <= 1 <= hi):
        s = "0.99" if x < 1 else "1.01"
    return s


def fmt_hr(hr, lo, hi):
    return f"{q(hr)} ({fmt_bound(lo, lo, hi)}–{fmt_bound(hi, lo, hi)})"


PANELS = [
    ("3", "A", "Kidney cancer", "Main_Full", "Kidney"), ("3", "B", "Bladder cancer", "Main_Full", "Bladder"),
    ("3", "C", "Breast cancer (women)", "Female_Breast", "Breast"), ("3", "D", "Prostate cancer (men)", "Male", "Prostate"),
    ("S3", "A", "Colorectal cancer", "Main_Full", "Colorectal"), ("S3", "B", "Lung cancer", "Main_Full", "Lung"),
    ("S3", "C", "Liver and intrahepatic bile duct cancer", "Main_Full", "Liver"), ("S3", "D", "Pancreatic cancer", "Main_Full", "Pancrea"),
    ("S3", "E", "Esophageal cancer", "Main_Full", "Esophag"), ("S3", "F", "Thyroid cancer", "Main_Full", "Thyroid"),
    ("S3", "G", "Oral, oropharyngeal and hypopharyngeal cancer", "Main_Full", "Oral"),
    ("S3", "H", "Endometrial cancer (women without prior hysterectomy)", "Female_Endo_noHyst", "Endometrial"),
]
YEARS = list(range(1, 11))


def main(root, outdir):
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    man, risk = [], []
    for fig, panel, label, sub, key in PANELS:
        d = export_dir(root, sub)
        r = find(read_outcomes(d), key)
        km = d / f"Outcome_{r['no']}_Result_b_KM_graph.csv"
        p = r["logrank_p"]
        man.append(dict(figure=fig, panel=panel, label=label, cohort=sub, km_csv=str(km),
                        logrank=("Log-rank p < 0.001" if p < 0.001 else f"Log-rank p = {fmt_p(p)}"),
                        hr="HR " + fmt_hr(r['hr'], r['lo'], r['hi']).replace("(", "(95% CI ", 1), n_per_arm=r["n"][0], ev1=r["ev"][0], ev2=r["ev"][1]))
        rt = risk_table(km, [y * 365 for y in YEARS])
        for arm in ("S", "G"):
            for y in YEARS:
                risk.append(dict(label=label, arm=arm, year=y, n=rt[arm][y * 365]))
        print(f"{fig:2} {panel} {label[:30]:30} S:{[rt['S'][y*365] for y in YEARS]}")
    with open(outdir / "km_manifest.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(man[0].keys())); w.writeheader(); w.writerows(man)
    with open(outdir / "km_risk.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(risk[0].keys())); w.writeheader(); w.writerows(risk)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
