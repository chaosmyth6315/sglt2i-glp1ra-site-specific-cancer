"""Figure S2: distribution of follow-up time in the full matched cohort.

Source: <export root>/Main_Full/IJC_*/Follow_up_time_metrics_matched_results_graph.csv (patients per
100-day bin of follow-up, by arm) and Follow_up_time_metrics_matched_results_table.csv (mean, SD,
median, IQR width), both in the same export folder.
The platform labels each bin by its midpoint (bin "400" = 350–450 days; first bin 0–50, last
3,550–3,650): treating the labels as midpoints reproduces the platform's mean, median and IQR width
within 2 days, whereas treating them as lower edges misses the median by 50 days. Bars are drawn over
the true intervals. Summary statistics printed in the plot: mean, SD and median from the platform table;
the 25th and 75th percentiles are interpolated from the bins (asserted to give the platform's median and
IQR width within 5 days). x ticks every year; y "Number of patients". Arial, monochrome, 600 dpi.
Usage: python3 -P figure_s2_followup.py <..._graph.csv> <..._table.csv> <out.png>
"""
import csv
import sys
from decimal import Decimal, ROUND_HALF_UP

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

from lib_trinetx_export import follow_up_quartiles

YEAR = Decimal("365.25")


def yrs(days):
    return str((Decimal(str(days)) / YEAR).quantize(Decimal("0.01"), rounding=ROUND_HALF_UP))


def main(graph, table, out):
    rows = list(csv.reader(open(graph, encoding="utf-8-sig", newline="")))
    i = next(k for k, r in enumerate(rows) if r and r[0].startswith("Time (Days)"))
    data = [(int(r[0]), int(r[1]), int(r[2])) for r in rows[i + 1:] if r and r[0].strip().lstrip("-").isdigit()]
    edges = [(max(d - 50, 0), min(d + 50, 3650)) for d, _, _ in data]      # labels are bin midpoints
    arms = {1: [a for _, a, _ in data], 2: [b for _, _, b in data]}
    trow = list(csv.reader(open(table, encoding="utf-8-sig", newline="")))
    plat = {int(r[0]): dict(mean=float(r[2]), sd=float(r[3]), med=float(r[4]), iqr=float(r[5]))
            for r in trow if r and r[0] in ("1", "2")}
    txt = []
    qs = follow_up_quartiles(Path(graph).parent)          # checked against the platform
    for k, name, arm in ((1, "SGLT2i", "S"), (2, "GLP-1 RA", "G")):
        q1, q2, q3 = qs[arm]
        p = plat[k]
        txt.append(f"{name}\n   Mean {yrs(p['mean'])} (SD {yrs(p['sd'])}) years\n"
                   f"   Median {yrs(p['med'])} (IQR {yrs(q1)}\u2013{yrs(q3)}) years")
        print(name, "Q1 %.0f median %.0f/%.0f Q3 %.0f IQR %.0f/%.0f" % (q1, q2, p["med"], q3, q3 - q1, p["iqr"]))
    plt.rcParams["font.family"] = "Arial"
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    lw = [(b - a) / 365.25 for a, b in edges]
    left = [a / 365.25 for a, _ in edges]
    ax.bar([l + 0.05 * w + 0.225 * w for l, w in zip(left, lw)], arms[1], width=[0.45 * w for w in lw],
           facecolor="white", edgecolor="black", linewidth=0.6, label="SGLT2i")
    ax.bar([l + 0.5 * w + 0.225 * w for l, w in zip(left, lw)], arms[2], width=[0.45 * w for w in lw],
           facecolor="white", edgecolor="black", linewidth=0.6, hatch="////", label="GLP-1 RA")
    ax.set_xlabel("Follow-up time, years (100-day bins)", fontsize=9)
    ax.set_ylabel("Number of patients", fontsize=9)
    ax.tick_params(labelsize=8)
    ax.set_xlim(0, 10.05); ax.set_xticks(range(0, 11))
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax.legend(fontsize=8, frameon=False, loc="upper left")
    ax.text(6.4, 18000, "\n\n".join(txt), ha="left", va="top", fontsize=7.5, linespacing=1.4)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(out, dpi=600)
    print("bins", len(data), "totals", sum(arms[1]), sum(arms[2]))
    print("\n".join(txt))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
