"""Figure S1: Love plot of standardized mean differences (SMDs) before and after matching.

Two side-by-side panels so that the figure fits one page at 7 pt (178 mm wide); the x axis starts at
−0.01 so that markers at 0 are not clipped. Male and female are separate rows: both are matching
covariates, and their SMDs differ because some patients have unknown sex.

Source: <export root>/Main_Full/IJC_*/Baseline_Patient_Characteristics.csv. Of the 85 matched
covariates, 81 have an SMD in the platform's balance table; four were constant in both arms (any
antidiabetic drug, ATC A10, 100%; dependence on renal dialysis Z99.2, immunosuppressants L04A and
immunostimulants L03, 0%, because they are exclusion criteria) and have no SMD. Statins and fibrates are
in the platform's table but were not matching covariates; they are plotted with the label "(not matched)".
Laboratory values and BMI are plotted with the SMD of their mean row. The platform matched them as one
present/absent indicator per prespecified range (Supplementary Method S1); the SMDs of those range rows are
not plotted.
Absolute SMDs sorted by the before-matching value; dashed vertical line at 0.1. Arial, monochrome, 600 dpi.
Usage: python3 -P figure_s1_love_plot.py <Baseline_Patient_Characteristics.csv> <out.png>
"""
import csv
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# Read the baseline export; lines before the header row (first cell "Characteristic ID") are skipped.
def load(path):
    with open(path, encoding="utf-8-sig", newline="") as f:
        rows = list(csv.reader(f))
    hdr_i = next(i for i, r in enumerate(rows) if r and r[0] == "Characteristic ID")
    hdr = rows[hdr_i]
    return [dict(zip(hdr, r)) for r in rows[hdr_i + 1:] if len(r) >= len(hdr)]


NAMES = {"AI": "Age at index", "2106-3": "White", "M": "Male", "F": "Female", "2054-5": "African American",
         "2131-1": "Other race", "2028-9": "Asian", "I10-I1A": "Hypertension", "E78": "Dyslipidemia",
         "E66": "Overweight and obesity", "Z12": "Cancer screening encounter", "Z00.0": "General medical examination",
         "I20-I25": "Ischemic heart disease", "E11.4": "T2DM, neurological complications",
         "E11.2": "T2DM, kidney complications", "N18": "Chronic kidney disease", "F41": "Anxiety disorders",
         "I50": "Heart failure", "Z87.891": "History of nicotine dependence", "F17": "Nicotine dependence",
         "I48": "Atrial fibrillation and flutter", "E11.5": "T2DM, circulatory complications",
         "I60-I69": "Cerebrovascular diseases", "J44": "COPD", "E11.3": "T2DM, ophthalmic complications",
         "K76.0": "Fatty liver", "I73": "Peripheral vascular diseases", "I42": "Cardiomyopathy",
         "I63": "Cerebral infarction", "Z72.0": "Tobacco use", "F10": "Alcohol-related disorders",
         "I70.2": "Atherosclerosis of extremity arteries", "K74": "Liver cirrhosis",
         "K74.6": "Other and unspecified cirrhosis", "G45": "Transient cerebral ischemic attacks",
         "I71": "Aortic aneurysm", "K75.81": "Nonalcoholic steatohepatitis", "K85": "Acute pancreatitis",
         "K22.7": "Barrett's esophagus", "K76.9": "Liver disease, unspecified", "B18.2": "Chronic hepatitis C",
         "I74": "Arterial embolism and thrombosis", "K86.1": "Other chronic pancreatitis",
         "B18.1": "Chronic hepatitis B", "K70.0": "Alcoholic fatty liver",
         "K21.01": "GERD with esophagitis, with bleeding", "C10AA": "Statins", "A10BA": "Biguanides",
         "C09": "Renin-angiotensin system agents", "A10A": "Insulins and analogues", "C07": "Beta blockers",
         "B01": "Antithrombotic agents", "N06A": "Antidepressants", "C08": "Calcium channel blockers",
         "M01A": "NSAIDs", "A10BB": "Sulfonylureas", "B01AC": "Platelet aggregation inhibitors",
         "A02BC": "Proton pump inhibitors", "A10BH": "DPP-4 inhibitors", "B01AF": "Direct factor Xa inhibitors",
         "C10AB": "Fibrates", "A10BG": "Thiazolidinediones", "C04": "Peripheral vasodilators",
         "9083": "BMI", "9029": "Sodium", "8001": "eGFR", "9028": "Potassium", "9022": "Calcium",
         "9030": "Blood urea nitrogen", "9037": "Hemoglobin A1c", "9044": "Alanine aminotransferase",
         "9047": "Aspartate aminotransferase", "9046": "Alkaline phosphatase", "9045": "Albumin",
         "9014": "Hemoglobin", "9001": "HDL cholesterol", "9004": "Triglyceride", "9015": "Leukocytes",
         "9002": "LDL cholesterol", "9016": "Lymphocytes", "9063": "C-reactive protein",
         "1013712": "Emergency department visit", "1013659": "Inpatient or observation care"}
DISPLAY_ONLY = {"A10", "L04A", "L03", "Z99.2"}
NOT_MATCHED = {"C10AA", "C10AB"}


def main(csv_path, out):
    recs = load(csv_path)
    rows = []
    seen = set()
    for r in recs:
        cid = r["Characteristic ID"]
        if cid in DISPLAY_ONLY or r["Category"] or cid in seen:
            continue
        if r["Cohort 1 Before: Mean"] == "" and cid.isdigit() and cid != "1013712" and cid != "1013659":
            continue  # laboratory "any value" count rows (the mean row is the covariate)
        if not r["Before: Standardized Mean Difference"]:
            continue
        seen.add(cid)
        name = NAMES.get(cid, r["Characteristic Name"]) + (" (not matched)" if cid in NOT_MATCHED else "")
        rows.append((name, abs(float(r["Before: Standardized Mean Difference"])),
                     abs(float(r["After: Standardized Mean Difference"]))))
    assert len(rows) == 83 and sum("not matched" in x[0] for x in rows) == 2, (len(rows), [x[0] for x in rows])
    rows.sort(key=lambda x: -x[1])                      # largest imbalance first
    half = (len(rows) + 1) // 2
    cols = [rows[:half], rows[half:]]
    plt.rcParams["font.family"] = "Arial"
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 6.9), constrained_layout=True)
    xmax = max(r[1] for r in rows) * 1.05
    for k, (ax, part) in enumerate(zip(axes, cols)):
        y = [half - 1 - i for i in range(len(part))]   # both panels aligned at the top
        for yi, r in zip(y, part):
            ax.plot([r[2], r[1]], [yi, yi], color="0.75", linewidth=0.6, zorder=2)
        ax.scatter([r[1] for r in part], y, marker="o", facecolors="white", edgecolors="black", s=14,
                   linewidths=0.7, label="Before matching", zorder=3)
        ax.scatter([r[2] for r in part], y, marker="o", color="black", s=14, label="After matching", zorder=4)
        ax.axvline(0.1, color="black", linestyle="--", linewidth=0.7)
        ax.set_yticks(y); ax.set_yticklabels([r[0] for r in part], fontsize=7)
        ax.set_ylim(-0.8, half - 0.2); ax.set_xlim(-0.01, xmax)
        ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5])
        ax.set_xlabel("Absolute standardized mean difference", fontsize=7.5)
        ax.tick_params(axis="x", labelsize=7); ax.tick_params(axis="y", length=2)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_linewidth(0.6)
    axes[1].legend(loc="lower right", fontsize=7.5, frameon=False)
    fig.savefig(out, dpi=600)
    print("max after-matching SMD:", max(r[2] for r in rows), "covariates:", len(rows))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
