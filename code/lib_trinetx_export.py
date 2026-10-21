"""Reader for the TriNetX "Export All" folders of the analyses.

export_dir(root, sub)        -> the IJC_* analysis folder inside <export root>/<sub>
read_outcomes(d)             -> {outcome name: {n, ev, hr, lo, hi, se, b, wald_p, ph_p, logrank_p,
                                 surv_end, risk}} (n/ev/risk as (S, G) tuples; None HR when masked)
matched_n(d)                 -> (S, G) matched counts from the baseline table row "AI"
follow_up(d)                 -> {"S"|"G": row of Follow_up_time_metrics_matched_results_table.csv as a dict
                                 of strings keyed by the platform's column names}
follow_up_quartiles(d)       -> 25th, 50th and 75th percentiles of follow-up (days) per arm
window_person_years(d)       -> person-years per arm inside the outcome window (day 365 to day 3,650)
find(res, key)               -> the outcome whose name contains key (case-insensitive)
"""
import csv
import math
import re
from pathlib import Path

Z = 1.959963984540054


def export_dir(root, sub):
    return next(p for p in sorted((Path(root) / sub).glob("IJC_*")) if p.is_dir())


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


def outcome_names(d):
    txt = (Path(d) / "_Read_Me.txt").read_text(encoding="utf-8-sig")
    return {int(k): v.strip() for k, v in re.findall(r"Outcome (\d+) Name: (.+)", txt)}


def read_outcomes(d):
    d = Path(d)
    res = {}
    for k, nm in outcome_names(d).items():
        kmf = d / f"Outcome_{k}_Result_b_KM_table.csv"
        if not kmf.exists():
            continue
        km = sections(kmf)
        st = km["Cohort Statistics"][1:3]
        r = {"no": k, "n": (int(st[0][2]), int(st[1][2])), "ev": (int(st[0][3]), int(st[1][3])),
             "surv_end": (float(st[0][5]), float(st[1][5])) if len(st[0]) > 5 and st[0][5] else None,
             "hr": None, "lo": None, "hi": None, "se": None, "b": None, "wald_p": None, "ph_p": None,
             "logrank_p": None, "risk": None}
        try:
            h, lo, hi = (float(x) for x in km["Hazard Ratio"][1][:3])
            r.update(hr=h, lo=lo, hi=hi, b=math.log(h), se=(math.log(hi) - math.log(lo)) / (2 * Z))
            r["wald_p"] = math.erfc(abs(r["b"] / r["se"]) / math.sqrt(2))
        except (KeyError, ValueError, IndexError, ZeroDivisionError):
            pass
        try:
            r["ph_p"] = float(km["Proportionality"][1][2])
        except (KeyError, ValueError, IndexError):
            pass
        try:
            r["logrank_p"] = float(km["Log-Rank Test"][1][2])
        except (KeyError, ValueError, IndexError):
            pass
        moa = d / f"Outcome_{k}_Result_a_MOA_table.csv"
        if moa.exists():
            try:
                ms = sections(moa)["Cohort Statistics"][1:3]
                r["risk"] = (float(ms[0][4]), float(ms[1][4]))
            except (KeyError, ValueError, IndexError):
                pass
        res[nm] = r
    return res


def find(res, key):
    k = next((nm for nm in res if key.lower() in nm.lower()), None)
    return res[k] if k else None


def matched_n(d):
    rows = list(csv.reader(open(Path(d) / "Baseline_Patient_Characteristics.csv", encoding="utf-8-sig", newline="")))
    ai = next(r for r in rows if r and r[0] == "AI")
    return int(ai[21]), int(ai[29])


def before_n(d):
    rows = list(csv.reader(open(Path(d) / "Baseline_Patient_Characteristics.csv", encoding="utf-8-sig", newline="")))
    ai = next(r for r in rows if r and r[0] == "AI")
    return int(ai[3]), int(ai[11])


def follow_up(d):
    rows = list(csv.reader(open(Path(d) / "Follow_up_time_metrics_matched_results_table.csv", encoding="utf-8-sig", newline="")))
    hdr = next(i for i, r in enumerate(rows) if r and r[0].strip() == "Cohort")
    out = {}
    for r in rows[hdr + 1:]:
        if r and r[0] in ("1", "2"):
            out["S" if r[0] == "1" else "G"] = dict(zip(rows[hdr], r))
    return out


def follow_up_quartiles(d, tol_days=5):
    """25th, 50th and 75th percentiles of follow-up (days) per arm, {"S": (q1, q2, q3), "G": ...}.

    The platform gives the median and only the WIDTH of the IQR, so the quartiles are interpolated from the
    100-day follow-up distribution (..._matched_results_graph.csv). Its bins are labelled by their midpoints
    (label 400 = 350–450 days; first bin 0–50, last 3,550–3,650). Checked against the platform: the
    interpolated median and IQR width must agree within `tol_days`.
    """
    rows = list(csv.reader(open(Path(d) / "Follow_up_time_metrics_matched_results_graph.csv", encoding="utf-8-sig", newline="")))
    i = next(k for k, r in enumerate(rows) if r and r[0].startswith("Time (Days)"))
    data = [(int(r[0]), int(r[1]), int(r[2])) for r in rows[i + 1:] if r and r[0].strip().lstrip("-").isdigit()]
    edges = [(max(t - 50, 0), min(t + 50, 3650)) for t, _, _ in data]

    def q(cnt, p):
        target, c = p * sum(cnt), 0
        for (a, b), n in zip(edges, cnt):
            if c + n >= target:
                return a + (target - c) / n * (b - a)
            c += n

    plat, out = follow_up(d), {}
    for arm, col in (("S", 1), ("G", 2)):
        cnt = [r[col] for r in data]
        q1, q2, q3 = (q(cnt, p) for p in (0.25, 0.5, 0.75))
        med, iqr = float(plat[arm]["Median Follow-up (Days)"]), float(plat[arm]["Interquartile Range"])
        assert abs(q2 - med) < tol_days and abs((q3 - q1) - iqr) < tol_days, (d, arm, q1, q2, q3, med, iqr)
        out[arm] = (q1, q2, q3)
    return out


def window_person_years(d, start=365, end=3650, rel_tol=0.01):
    """Person-years inside the outcome window (day `start` to day `end`) per arm, {"S": py, "G": py}.

    Interpolated from the 100-day follow-up distribution (bins labelled by midpoint; first bin 0–50, last
    3,550–3,650), assuming follow-up is spread evenly within each bin. Check: the same calculation from day 0
    must reproduce matched n × the platform's mean follow-up within `rel_tol`.
    """
    rows = list(csv.reader(open(Path(d) / "Follow_up_time_metrics_matched_results_graph.csv", encoding="utf-8-sig", newline="")))
    i = next(k for k, r in enumerate(rows) if r and r[0].startswith("Time (Days)"))
    data = [(int(r[0]), int(r[1]), int(r[2])) for r in rows[i + 1:] if r and r[0].strip().lstrip("-").isdigit()]
    edges = [(max(t - 50, 0), min(t + 50, end)) for t, _, _ in data]

    def mean_excess(a, b, s):   # E[max(0, x − s)] for x uniform on [a, b]
        if b <= s:
            return 0.0
        if a >= s:
            return (a + b) / 2 - s
        return (b - s) ** 2 / (2 * (b - a))

    plat, out = follow_up(d), {}
    for arm, col in (("S", 1), ("G", 2)):
        n = sum(r[col] for r in data)
        total = sum(r[col] * mean_excess(a, b, 0) for r, (a, b) in zip(data, edges))
        ref = n * float(plat[arm]["Mean Follow-up (Days)"])
        assert abs(total - ref) / ref < rel_tol, (d, arm, total, ref)
        out[arm] = sum(r[col] * mean_excess(a, b, start) for r, (a, b) in zip(data, edges)) / 365.25
    return out
