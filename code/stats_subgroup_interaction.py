"""Subgroup table with P for interaction (TriNetX subgroup exports).

Usage: python3 -I stats_subgroup_interaction.py <subgroup export root> <main export root> [Full|F|M]

Cohorts: Full (kidney, bladder, colorectal), F = women (breast), M = men (prostate).
Each subgroup is a separately matched analysis (Sub_<cohort>_<stratum>).
Sex (full cohort only) uses the Female_Breast and Male main-analysis exports.
P for interaction: two strata -> z test on the difference of log HRs; three or more -> Cochran's Q
(SE of log HR from the 95% CI). Cells with no estimate (masked counts) are reported as such.
HR/CI rounded once, ROUND_HALF_UP; P to 3 decimals (4 when it would round to 0.050).
"""
import csv
import math
import re
import sys
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

Z = 1.959963984540054
CANCERS = {"Full": ["Kidney", "Bladder", "Colorectal"], "F": ["Breast"], "M": ["Prostate"]}
GROUPS = [
    ("Age at index", [("18–64 y", "IdxAge18-64"), ("≥65 y", "IdxAge65plus")]),
    ("Sex", [("Female", "SEX_F"), ("Male", "SEX_M")]),
    ("Race", [("White", "White"), ("Black", "Black"), ("Asian", "Asian")]),
    ("BMI", [("<30 kg/m²", "BMIlt30"), ("≥30 kg/m²", "BMIge30")]),
    ("Smoking", [("Yes", "SmokeY"), ("No", "SmokeN")]),
    ("eGFR", [("<60", "eGFRlt60"), ("≥60", "eGFRge60")]),
    ("HbA1c", [("<7%", "HbA1clt7"), ("≥7%", "HbA1cge7")]),
    ("UACR", [("<30 mg/g", "UACRlt30"), ("≥30 mg/g", "UACRge30")]),
    ("CKD", [("Yes", "CKDY"), ("No", "CKDN")]),
    ("Dyslipidaemia", [("Yes", "DyslipY"), ("No", "DyslipN")]),
    ("Heart failure", [("Yes", "HFY"), ("No", "HFN")]),
    ("Hypertension", [("Yes", "HTNY"), ("No", "HTNN")]),
    ("Ischaemic heart disease", [("Yes", "IHDY"), ("No", "IHDN")]),
    ("Insulin", [("Yes", "InsulinY"), ("No", "InsulinN")]),
    ("Metformin", [("Yes", "MetforminY"), ("No", "MetforminN")]),
]


def q(x, p="0.01"):
    return str(Decimal(str(x)).quantize(Decimal(p), rounding=ROUND_HALF_UP))


# HR (95% CI) display: a CI bound that would print as 1.00 while the CI excludes 1 is shown as 0.99 or
# 1.01, as in the figures.
def fmt_hr(hr, lo, hi):
    def bound(x):
        s = q(x)
        return ("0.99" if x < 1 else "1.01") if s == "1.00" and not (lo <= 1 <= hi) else s
    return f"{q(hr)} ({bound(lo)}–{bound(hi)})"


def fmt_p(p):
    if p < 0.001:
        return "<0.001"
    s = q(p, "0.001")
    return q(p, "0.0001") if (p < 0.05 and s == "0.050") else s


def sections(path):
    out, cur = {}, None
    with open(path, encoding="utf-8-sig", newline="") as f:
        for row in csv.reader(f):
            cells = [c.strip() for c in row]
            if not any(cells):
                continue
            if len(cells) == 1:
                cur = cells[0]
                out[cur] = []
            elif cur is not None:
                out[cur].append(cells)
    return out


def read(folder, keys):
    dirs = [p for p in sorted(folder.glob("IJC_*")) if p.is_dir()] if folder.exists() else []
    if not dirs:
        return None
    d = dirs[0]
    names = {int(k): v.strip() for k, v in re.findall(r"Outcome (\d+) Name: (.+)", (d / "_Read_Me.txt").read_text(encoding="utf-8-sig"))}
    res = {}
    for key in keys:
        k = next((i for i, v in names.items() if key.lower() in v.lower()), None)
        f = d / f"Outcome_{k}_Result_b_KM_table.csv" if k else None
        if not f or not f.exists():
            res[key] = None
            continue
        km = sections(f)
        st = km["Cohort Statistics"][1:3]
        n, ev = (st[0][2], st[1][2]), (st[0][3], st[1][3])
        surv = (st[0][5], st[1][5])   # KM survival at the end of the window (for the cumulative incidence, CIR)
        try:
            h, lo, hi = (float(x) for x in km["Hazard Ratio"][1][:3])
            res[key] = {"n": n, "ev": ev, "surv": surv, "hr": h, "lo": lo, "hi": hi, "b": math.log(h), "se": (math.log(hi) - math.log(lo)) / (2 * Z)}
        except (KeyError, ValueError, IndexError):
            res[key] = {"n": n, "ev": ev, "surv": surv, "hr": None}
    return res


def p_int(ests):
    ests = [e for e in ests if e and e.get("hr")]
    if len(ests) < 2:
        return None
    if len(ests) == 2:
        a, b = ests
        z = (a["b"] - b["b"]) / math.sqrt(a["se"] ** 2 + b["se"] ** 2)
        return math.erfc(abs(z) / math.sqrt(2))
    w = [1 / e["se"] ** 2 for e in ests]
    bp = sum(wi * e["b"] for wi, e in zip(w, ests)) / sum(w)
    Q = sum(wi * (e["b"] - bp) ** 2 for wi, e in zip(w, ests))
    df = len(ests) - 1
    # chi-square upper tail (df small): regularized gamma via series
    return chi2_sf(Q, df)


def chi2_sf(x, k):
    a, x2 = k / 2, x / 2
    if x2 <= 0:
        return 1.0
    term = s = 1 / a
    n = 1
    while term > 1e-15 * s:
        term *= x2 / (a + n)
        s += term
        n += 1
    lower = math.exp(-x2 + a * math.log(x2) - math.lgamma(a)) * s
    return max(0.0, 1 - lower)


def main(root, root_main, cohort):
    root, root_main = Path(root), Path(root_main)
    keys = CANCERS[cohort]
    print(f"\n## {cohort}\n")
    print("| Subgroup | " + " | ".join(f"{k}: events S/G | {k}: HR (95% CI) | {k}: P-int" for k in keys) + " |")
    print("|" + "---|" * (1 + 3 * len(keys)))
    for gname, strata in GROUPS:
        if gname == "Sex" and cohort != "Full":
            continue
        rows = []
        for lab, code in strata:
            if code == "SEX_F":
                folder = root_main / "Female_Breast"
            elif code == "SEX_M":
                folder = root_main / "Male"
            else:
                folder = root / f"Sub_{cohort}_{code}"
            rows.append((lab, read(folder, keys)))
        pints = {k: p_int([r[1][k] if r[1] else None for r in rows]) for k in keys}
        print(f"| **{gname}** | " + " | ".join(f" | | {fmt_p(pints[k]) if pints[k] is not None else ''}" for k in keys) + " |")
        for lab, r in rows:
            cells = []
            for k in keys:
                e = r[k] if r else None
                if r is None:
                    cells += ["not run", "", ""]
                elif e is None:
                    cells += ["–", "–", ""]
                elif e.get("hr") is None:
                    cells += [f"{e['ev'][0]}/{e['ev'][1]}", "not estimable", ""]
                else:
                    cells += [f"{e['ev'][0]}/{e['ev'][1]}", fmt_hr(e['hr'], e['lo'], e['hi']), ""]
            print(f"| {lab} | " + " | ".join(cells) + " |")


if __name__ == "__main__":
    for c in (sys.argv[3:] or ["Full", "F", "M"]):
        main(sys.argv[1], sys.argv[2], c)
