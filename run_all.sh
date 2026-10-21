#!/usr/bin/env bash
# Rebuild the figure data, figures and the statistics quoted in the text from TriNetX "Export All" folders.
#
# usage: ./run_all.sh <primary export root> <subgroup export root> <out dir>
#   primary export root : Main_Full, Female_Breast, Female_Endo_noHyst, Male, Figure1_AnalyzeCriteria
#   subgroup export root: Sub_<cohort>_<stratum> analyses (Figure 4, Figures S4 to S8)
set -euo pipefail
if [ $# -ne 3 ]; then
  sed -n '4,6p' "$0"; exit 1
fi
P=$(cd "$1" && pwd); SUB=$(cd "$2" && pwd)
mkdir -p "$3"; O=$(cd "$3" && pwd)
mkdir -p "$O/figure_data/km" "$O/figures" "$O/logs"
cd "$(dirname "$0")/code"
export PYTHONPATH=.
export LC_ALL=en_US.UTF-8

analysis() { ls -d "$P/$1"/IJC_* | head -n 1; }      # the analysis folder inside one export
BASE_FULL="$(analysis Main_Full)/Baseline_Patient_Characteristics.csv"

echo "Figure data and statistics"
python3 figure2_data.py "$P" "$O/figure_data/fig2_data.csv" > "$O/logs/fig2_data.txt"
python3 figure4_s4_s8_subgroup_data.py "$SUB" "$P" "$O/figure_data/fig4_data.csv" > "$O/logs/fig4_data.txt"
python3 figure3_s3_km_data.py "$P" "$O/figure_data/km" > "$O/logs/km_data.txt"
python3 stats_pvalues_bh_cochran.py "$P" > "$O/logs/main_stats.txt"
python3 stats_absolute_risk.py "$P" > "$O/logs/absolute_risks.txt"

echo "Figures"
python3 figure1_flowchart.py "$P/Figure1_AnalyzeCriteria" "$P" "$O/figures/Figure 1_flow chart.drawio" >/dev/null
FIGOUT="$O/figures" FIG2DATA="$O/figure_data/fig2_data.csv" Rscript figure2_forest.R
FIGOUT="$O/figures" KMDIR="$O/figure_data/km" Rscript figure3_s3_km.R
FIGOUT="$O/figures" FIG4DATA="$O/figure_data/fig4_data.csv" Rscript figure4_subgroup_forest.R
FIGOUT="$O/figures" FIG4DATA="$O/figure_data/fig4_data.csv" Rscript figure_s4_s8_subgroup_forest.R
python3 figure_s1_love_plot.py "$BASE_FULL" "$O/figures/Figure S1_Love plot.png" >/dev/null
FU="$(analysis Main_Full)"
python3 figure_s2_followup.py "$FU/Follow_up_time_metrics_matched_results_graph.csv" \
  "$FU/Follow_up_time_metrics_matched_results_table.csv" "$O/figures/Figure S2_follow-up time.png" >/dev/null

echo "Done: $O"
