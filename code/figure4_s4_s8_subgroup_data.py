"""Subgroup data for Figure 4 and Supplementary Figures S4–S8: long CSV, one row per cancer × stratum.

Usage: python3 -I figure4_s4_s8_subgroup_data.py <subgroup export root> <primary export root> <out.csv>

Input layout: <subgroup export root>/Sub_<cohort>_<stratum code>/IJC_*/ (cohort Full, F or M);
<primary export root>/Female_Breast/IJC_*/ and <primary export root>/Male/IJC_*/.
Cancers: Kidney, Bladder, Colorectal (full cohort, Sub_Full_*; the sex strata are the Female_Breast
and Male primary exports), Breast (Sub_F_*), Prostate (Sub_M_*).
Reads the KM tables through stats_subgroup_interaction.read(); P for interaction via p_int()
(two strata: z test on the difference of the log HRs; three strata (race): Cochran's Q). Display
strings are rounded once (ROUND_HALF_UP; a CI bound that would print 1.00 while the CI excludes 1 is
moved to 0.99 / 1.01).
Strata with no export (Sub_Full_SmokeN, blocked by the platform's complexity limit) get
status "not_run"; strata that do not apply (sex in the cohorts of women and men) are absent.
"""
import csv
import sys
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from stats_subgroup_interaction import GROUPS, read, p_int, fmt_p  # noqa: E402

LABELS = {
    "Age at index": "Age at index (yrs)", "BMI": "BMI (kg/m²)", "eGFR": "eGFR (mL/min/1.73 m²)",
    "HbA1c": "HbA1c (%)", "UACR": "UACR (mg/g)", "Dyslipidaemia": "Dyslipidemia",
    "Ischaemic heart disease": "Ischemic heart disease",
}
STRATA = {"18–64 y": "18–64", "≥65 y": "≥65", "<30 kg/m²": "<30", "≥30 kg/m²": "≥30",
          "<7%": "<7", "≥7%": "≥7", "<30 mg/g": "<30", "≥30 mg/g": "≥30", "Black": "African American"}
CANCERS = [("Kidney", "Full"), ("Bladder", "Full"), ("Colorectal", "Full"), ("Breast", "F"), ("Prostate", "M")]


def q2(x, places="0.01"):
    return str(Decimal(str(x)).quantize(Decimal(places), rounding=ROUND_HALF_UP))


def fmt_bound(x, lo, hi):
    """2 dp. A bound that would print as 1.00 while the CI excludes 1 is moved one unit away from
    the null (0.9972 -> 0.99; 1.0035 -> 1.01) so that the display never contradicts the p value."""
    s = q2(x)
    if s == "1.00" and not (lo <= 1 <= hi):
        s = q2(Decimal("0.99")) if x < 1 else q2(Decimal("1.01"))
    return s


def cir(e):
    """10-year cumulative incidence (%) = 100 × (1 − KM survival at day 3,650), rounded once."""
    return {f"cir{i + 1}_s": str(q2(100 * (1 - Decimal(e["surv"][i])))) for i in range(2)}


def main(root, root_main, out):
    root, root_main = Path(root), Path(root_main)
    rows = []
    for cancer, cohort in CANCERS:
        for gname, strata in GROUPS:
            if gname == "Sex" and cohort != "Full":
                continue
            ests = []
            for lab, code in strata:
                if code == "SEX_F":
                    folder = root_main / "Female_Breast"
                elif code == "SEX_M":
                    folder = root_main / "Male"
                else:
                    folder = root / f"Sub_{cohort}_{code}"
                r = read(folder, [cancer])
                e = r[cancer] if r else None
                ests.append((lab, code, folder.name, e))
            pi = p_int([e for _, _, _, e in ests])
            for lab, code, src, e in ests:
                row = dict(cancer=cancer, cohort=cohort, group=LABELS.get(gname, gname),
                           stratum=STRATA.get(lab, lab), code=code, src=src,
                           p_int=fmt_p(pi) if pi is not None else "")
                if e is None:
                    row.update(status="not_run", n_per_arm="", ev1="", ev2="", hr="", lo="", hi="", hr_ci_s="", cir1_s="", cir2_s="")
                elif e.get("hr") is None:
                    row.update(status="not_estimable", n_per_arm=e["n"][0], ev1=e["ev"][0], ev2=e["ev"][1],
                               hr="", lo="", hi="", hr_ci_s="", **cir(e))
                else:
                    row.update(status="ok", n_per_arm=e["n"][0], ev1=e["ev"][0], ev2=e["ev"][1],
                               hr=e["hr"], lo=e["lo"], hi=e["hi"],
                               hr_ci_s=f"{q2(e['hr'])} ({fmt_bound(e['lo'], e['lo'], e['hi'])}-{fmt_bound(e['hi'], e['lo'], e['hi'])})", **cir(e))
                rows.append(row)
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    for r in rows:
        print(f"{r['cancer']:10} {r['group']:24} {r['stratum']:17} {r['status']:13} n={r['n_per_arm']:>7} "
              f"ev={r['ev1']}/{r['ev2']} {r['hr_ci_s']:20} p-int={r['p_int']}")


if __name__ == "__main__":
    main(*sys.argv[1:4])
