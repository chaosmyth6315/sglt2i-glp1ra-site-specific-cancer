"""Pooled sensitivity estimates quoted in the text, and the landmark column of Figure 2.

Usage: python3 -I stats_sensitivity_pooled.py <sensitivity export root> <out landmark csv>

Analyses that exceeded the platform's complexity limit in the full cohort were run in the separately matched
cohorts of women and men, and the two log hazard ratios were combined by fixed-effect inverse-variance pooling.
The same pooling combines the two calendar periods (calendar-stratified estimate). Expected folders under
<sensitivity export root>, one TriNetX "Export All" per analysis:

  Landmark_Female_Breast, Landmark_Male                       outcome window from day 365, patients with the
                                                              outcome before the window excluded
  Cal2013-2018_Full, Cal2013-2018_Female_Breast,              index dates 2013-2018
  Cal2013-2018_Female_Endo_noHyst, Cal2013-2018_Male
  Cal2019-2024_Female_Breast, Cal2019-2024_Female_Endo_noHyst, index dates 2019-2024
  Cal2019-2024_Male
  PreIdx_Female, PreIdx_Male                                  at least one encounter and one type 2 diabetes
                                                              diagnosis on or before the index prescription
  PSA_Male                                                    PSA testing, day 1 to day 3,650 (men)

Prints the estimates; writes <out landmark csv> (label, lm_s) for figure2_forest.R.
Rounded once, ROUND_HALF_UP; a confidence limit that excludes 1 is never printed as 1.00.
"""
import csv
import math
import sys
from decimal import Decimal, ROUND_HALF_UP, ROUND_DOWN, ROUND_UP

from lib_trinetx_export import Z, export_dir, find, read_outcomes

ROOT, OUT = sys.argv[1], sys.argv[2]
_cache = {}


def outcome(sub, key):
    if sub not in _cache:
        _cache[sub] = read_outcomes(export_dir(ROOT, sub))
    r = find(_cache[sub], key)
    assert r is not None and r["b"] is not None, (sub, key)
    return r


def pool(items):
    """Fixed-effect inverse-variance pooling of (log HR, SE) pairs."""
    w = [1 / se ** 2 for _, se in items]
    b = sum(wi * bi for wi, (bi, _) in zip(w, items)) / sum(w)
    return b, 1 / math.sqrt(sum(w))


def q(x, places="0.01", mode=ROUND_HALF_UP):
    return Decimal(repr(x)).quantize(Decimal(places), rounding=mode)


def hr_ci(b, se):
    h, lo, hi = math.exp(b), math.exp(b - Z * se), math.exp(b + Z * se)
    los, his = q(lo), q(hi)
    if str(his) == "1.00" and hi < 1:
        his = q(hi, mode=ROUND_DOWN)
    if str(los) == "1.00" and lo > 1:
        los = q(lo, mode=ROUND_UP)
    return f"{q(h)} ({los}–{his})"


def p_txt(b, se):
    p = math.erfc(abs(b / se) / math.sqrt(2))
    if p < 0.001:
        return "<0.001"
    s = q(p, "0.001")
    return str(q(p, "0.0001")) if (p < 0.05 and str(s) == "0.050") else str(s)


def estimate(subs, key):
    rs = [outcome(s, key) for s in subs]
    return pool([(r["b"], r["se"]) for r in rs]), rs


def line(label, b, se, rs):
    n = tuple(sum(r["n"][i] for r in rs) for i in (0, 1))
    ev = tuple(sum(r["ev"][i] for r in rs) for i in (0, 1))
    print(f"  {label:42s} {hr_ci(b, se):18s} p {p_txt(b, se):7s} analysed n S/G {n[0]:,}/{n[1]:,}  events S/G {ev[0]:,}/{ev[1]:,}")


SEX = ("Female_Breast", "Male")
FULL_SITES = [("Kidney cancer", "Kidney"), ("Bladder cancer", "Bladder"), ("Colorectal cancer", "Colorectal")]

# landmark (Figure 2 column)
print("Landmark analysis (day 365-3,650, patients with the outcome before the window excluded)")
rows = []
for label, key in FULL_SITES + [("Breast cancer", "Breast"), ("Prostate cancer", "Prostat")]:
    subs = [f"Landmark_{s}" for s in SEX] if label in dict(FULL_SITES) else \
        ["Landmark_Female_Breast" if key == "Breast" else "Landmark_Male"]
    (b, se), rs = estimate(subs, key)
    line(label, b, se, rs)
    rows.append((label, hr_ci(b, se)))
with open(OUT, "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f); w.writerow(["label", "lm_s"]); w.writerows(rows)

# calendar-stratified
print("\nCalendar-stratified (index dates 2013-2018 and 2019-2024 combined)")
CAL = FULL_SITES + [("Breast cancer", "Breast"), ("Endometrial cancer", "Endometri"), ("Prostate cancer", "Prostat")]
for label, key in CAL:
    if label in dict(FULL_SITES):
        s18, s24 = ["Cal2013-2018_Full"], [f"Cal2019-2024_{s}" for s in SEX]
    else:
        cohort = {"Breast": "Female_Breast", "Endometri": "Female_Endo_noHyst", "Prostat": "Male"}[key]
        s18, s24 = [f"Cal2013-2018_{cohort}"], [f"Cal2019-2024_{cohort}"]
    (b18, se18), r18 = estimate(s18, key)
    (b24, se24), r24 = estimate(s24, key)
    b, se = pool([(b18, se18), (b24, se24)])
    zd = (b18 - b24) / math.sqrt(se18 ** 2 + se24 ** 2)
    line(label, b, se, r18 + r24)
    print(f"  {'':42s} 2013-2018 {hr_ci(b18, se18)}; 2019-2024 {hr_ci(b24, se24)}; "
          f"P for period difference {p_txt(zd, 1.0)}")

# encounter and diagnosis before the index prescription
print("\nAt least one encounter and one type 2 diabetes diagnosis on or before the index prescription")
for label, key in FULL_SITES + [("Breast cancer", "Breast"), ("Prostate cancer", "Prostat")]:
    subs = ["PreIdx_Female", "PreIdx_Male"] if label in dict(FULL_SITES) else \
        ["PreIdx_Female" if key == "Breast" else "PreIdx_Male"]
    (b, se), rs = estimate(subs, key)
    line(label, b, se, rs)

# PSA testing
print("\nPSA testing, matched cohort of men (day 1-3,650)")
r = outcome("PSA_Male", "PSA")
n = r["n"][0]
pct = [q(e * 100 / n, "0.01") for e in r["ev"]]
print(f"  {n:,} pairs; SGLT2i {r['ev'][0]:,} ({pct[0]}%) vs GLP-1 RA {r['ev'][1]:,} ({pct[1]}%); "
      f"HR {hr_ci(r['b'], r['se'])}, p {p_txt(r['b'], r['se'])}")
