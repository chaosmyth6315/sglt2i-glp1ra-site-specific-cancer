# SGLT2 inhibitors versus GLP-1 receptor agonists and site-specific cancer risk

Analysis and figure code for the paper *Comparative Risk of Site-Specific Cancers Among SGLT2 Inhibitor and GLP-1 Receptor Agonist Initiators With Type 2 Diabetes* (International Journal of Cancer, under review).

The study emulates a target trial in the TriNetX US Collaborative Network. Adult new users of an SGLT2 inhibitor or a GLP-1 receptor agonist (2013 to 2024) were matched 1:1 on propensity scores. Incident cancer at 12 sites was followed from year 1 to year 10 after the index date. The cohorts, outcomes and Cox models were run on the TriNetX platform. The code in this repository reads the platform's "Export All" output and builds the figures of the paper and its supplement, and the derived statistics quoted in the text (multiplicity-adjusted and heterogeneity p values, P for interaction and absolute risks).

The code documents how each figure and derived statistic was computed. Running it requires TriNetX exports, which cannot be shared.

## Data

The data are de-identified electronic health records licensed from TriNetX. Patient-level data and the platform exports cannot be redistributed, so none are included here. Researchers with TriNetX access can rebuild the cohorts and outcomes from the code lists in Supplementary Method S1 of the paper and export the analyses with "Export All". Aggregate data are available from the corresponding author on reasonable request.

Each TriNetX query rebuilds and re-matches its cohorts on the network as it stands that day, so a new run returns cohorts of a different size from those in the paper.

Expected layout of the exports (one folder per analysis, each holding the platform's `IJC_*` analysis folder):

```
<primary export root>/
  Main_Full/  Female_Breast/  Female_Endo_noHyst/  Male/        primary analyses
  Figure1_AnalyzeCriteria/                                         Criteria Analysis exports (Figure 1)
<subgroup export root>/
  Sub_<cohort>_<stratum>/                                          separately matched subgroup analyses
```

`<cohort>` is `Full`, `F` (women) or `M` (men).

## Requirements

- Python 3 with the packages in `requirements.txt` (`pip install -r requirements.txt`)
- R with forestploter, ggplot2, patchwork, ragg and zoo
- draw.io desktop (command-line export) to turn the Figure 1 `.drawio` file into PDF
- Arial installed for the figures

Tested with Python 3.14, matplotlib 3.10.8, R 4.4.2, forestploter 1.1.2, ggplot2 4.0.2, patchwork 1.3.2, ragg 1.5.0 and zoo 1.8.15.

## Running

```bash
./run_all.sh <primary export root> <subgroup export root> <out dir>
```

This writes the figure data (`.csv`), the Figure 1 draw.io file, the figures (`.png`, 600 dpi) and text logs of the statistics quoted in the text. The published figures are these images with the white margin trimmed. Each script can also be run on its own from the `code/` folder with `PYTHONPATH=.`; its usage line is at the top of the file.

## Scripts

Each script is named after what it produces. `figure*` scripts build the numbered figures, `stats_*` scripts print the statistics quoted in the text, and `lib_*` files are shared helpers that are imported, not run. The tables of the paper were built from the same exports and are not part of this repository.

| Script | Produces |
|---|---|
| `figure1_flowchart.py` | Figure 1, cohort flow chart (draw.io) |
| `figure2_data.py`, `figure2_forest.R` | Figure 2, site-specific hazard ratios |
| `figure3_s3_km_data.py`, `figure3_s3_km.R` | Figure 3 and Figure S3, Kaplan–Meier curves with numbers at risk |
| `figure4_s4_s8_subgroup_data.py`, `figure4_subgroup_forest.R`, `figure_s4_s8_subgroup_forest.R` | Figure 4 and Figures S4 to S8, subgroup analyses with P for interaction |
| `figure_s1_love_plot.py` | Figure S1, standardized mean differences before and after matching |
| `figure_s2_followup.py` | Figure S2, distribution of follow-up time |
| `stats_pvalues_bh_cochran.py` | Wald p values, Benjamini–Hochberg q values and Cochran's Q across the 12 sites; deaths by arm |
| `stats_subgroup_interaction.py` | Subgroup hazard ratios with P for interaction (also used for Figure 4) |
| `stats_absolute_risk.py` | Absolute risks in the four primary analyses |
| `lib_trinetx_export.py` | Reader for the TriNetX "Export All" folders |
| `lib_km_number_at_risk.py` | Numbers at risk recovered from the Kaplan–Meier exports |

## Conventions

- Every number is rounded once, half up, from the most precise value in the export.
- Hazard ratios and confidence limits have two decimals. A confidence limit that excludes 1 is never printed as 1.00 (0.9972 prints as 0.99).
- p values are Wald p values computed from the exported hazard ratio and its 95% CI. Benjamini–Hochberg q values are computed across the 12 sites.
- Cumulative incidence is 100 × (1 − Kaplan–Meier survival) at day 3,650. The outcome window runs from day 365 to day 3,650.
- TriNetX masks small counts. Strata the platform could not estimate are shown as NA.

## License

MIT; see `LICENSE`.
