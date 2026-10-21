#################################################################
# Supplementary Figures S4–S8: detailed subgroup forest plots, one per site:
#   colorectal (S4), kidney (S5) and bladder (S6) cancer in the full cohort, breast cancer
#   in women (S7) and prostate cancer in men (S8). They give the numbers behind main-text
#   Figure 4, which shows the forest plots only. Output: FIGOUT/Figure S<n>_subgroup <site>.png.
# Columns: subgroup | matched pairs | SGLT2i events (CIR, %) | GLP-1 RA events (CIR, %) | forest |
#   HR (95% CI) | P for interaction (printed on the group row). CIR as in Figure 2:
#   100 × (1 − KM survival at day 3,650), rounded once in figure4_s4_s8_subgroup_data.py.
# Data: fig4_data.csv from figure4_s4_s8_subgroup_data.py; strings are rounded once there
#   (ROUND_HALF_UP; a CI bound is never printed as 1.00 when the CI excludes 1);
#   the hyphen in the CI is shown as an en dash.
# "NA" = analysis blocked by the platform's complexity limit (full cohort, non-smokers);
#   sex strata are absent in the cohorts of women and men.
# Three-line table borders; print size <=178 mm wide, Arial 8 pt; log x axis 0.6–1.6 as in
#   Figure 4 (forestploter draws an arrow where a CI runs past the axis).
# SITES (comma-separated, e.g. "Kidney,Breast") selects the sites; default all five.
# Usage: FIGOUT=<dir> FIG4DATA=<csv> Rscript figure_s4_s8_subgroup_forest.R (en_US.UTF-8)
#################################################################
library(forestploter); library(grid); library(ragg)
pdf(NULL)   # null device for text measurement, so no Rplots.pdf is written
outdir <- Sys.getenv("FIGOUT")
dat <- read.csv(Sys.getenv("FIG4DATA"), stringsAsFactors = FALSE, check.names = FALSE,
                encoding = "UTF-8", colClasses = "character")
sites <- c(Kidney = "kidney cancer", Bladder = "bladder cancer",
           Breast = "breast cancer", Prostate = "prostate cancer",
           Colorectal = "colorectal cancer")   # names = site keys in fig4_data.csv and in the file names
fig_no <- c(Colorectal = "S4", Kidney = "S5", Bladder = "S6", Breast = "S7", Prostate = "S8")

build <- function(nm) {
  d <- dat[dat$cancer == nm, ]
  rows <- list()
  for (g in unique(d$group)) {
    sub <- d[d$group == g, ]
    pint <- unique(sub$p_int); stopifnot(length(pint) == 1)
    rows[[length(rows) + 1]] <- data.frame(Subgroup = g, pairs = "", e1 = "", e2 = "", hr = NA, lo = NA,
                                           hi = NA, ci = "", p = pint, hdr = TRUE, stringsAsFactors = FALSE)
    for (i in seq_len(nrow(sub))) {
      r <- sub[i, ]
      run <- r$status == "ok"
      if (run) stopifnot(nzchar(r$cir1_s), nzchar(r$cir2_s))
      rows[[length(rows) + 1]] <- data.frame(
        Subgroup = paste0("   ", r$stratum),
        pairs = if (run) format(as.integer(r$n_per_arm), big.mark = ",") else "",
        e1 = if (run) sprintf("%s (%s)", format(as.integer(r$ev1), big.mark = ","), r$cir1_s) else "",
        e2 = if (run) sprintf("%s (%s)", format(as.integer(r$ev2), big.mark = ","), r$cir2_s) else "",
        hr = if (run) as.numeric(r$hr) else NA, lo = if (run) as.numeric(r$lo) else NA,
        hi = if (run) as.numeric(r$hi) else NA,
        ci = if (run) gsub("-", "–", r$hr_ci_s, fixed = TRUE) else "NA",
        p = "", hdr = FALSE, stringsAsFactors = FALSE)
    }
  }
  df <- do.call(rbind, rows)
  stopifnot(sum(!df$hdr) == nrow(d))
  ok <- !is.na(df$hr)
  stopifnot(all(df$hr[ok] > 0.6 & df$hr[ok] < 1.6))
  out <- data.frame(Subgroup = df$Subgroup, check.names = FALSE, stringsAsFactors = FALSE)
  out$`Matched\npairs` <- df$pairs
  out$`SGLT2i\nEvents (CIR, %)` <- df$e1
  out$`GLP-1 RA\nEvents (CIR, %)` <- df$e2
  out$` ` <- paste(rep(" ", 30), collapse = "")
  out$`HR (95% CI)` <- df$ci
  out$`P for\ninteraction` <- df$p

  tm <- forest_theme(base_size = 8, base_family = "Arial",
    refline_gp = gpar(lwd = 0.75, lty = "dashed", col = "grey35"),
    ci_pch = 15, ci_col = "black", ci_fill = "black", ci_lwd = 0.9, ci_Theight = 0.25,
    arrow_type = "closed", arrow_length = 0.06, arrow_gp = gpar(fontsize = 7.5),
    xaxis_gp = gpar(fontsize = 7.5, lwd = 0.75),
    core = list(padding = unit(c(1.2, 2.2), "mm"), bg_params = list(fill = "white")),
    colhead = list(fg_params = list(hjust = 0, x = 0.04, fontface = "bold"),
                   padding = unit(c(1.2, 2.2), "mm")))
  p <- forest(out, est = df$hr, lower = df$lo, upper = df$hi, sizes = 0.45, ci_column = 5,
              ref_line = 1, x_trans = "log", xlim = c(0.6, 1.6), ticks_at = c(0.7, 1, 1.4),
              arrow_lab = c("Favors SGLT2i", "Favors GLP-1 RA"), theme = tm)
  TH <- 1.33; TN <- 0.67
  p <- edit_plot(p, row = which(df$hdr), gp = gpar(fontface = "bold"))
  p <- add_border(p, part = "header", row = 1, where = "top", gp = gpar(lwd = TH, col = "black"))
  p <- add_border(p, part = "header", row = 1, where = "bottom", gp = gpar(lwd = TN, col = "black"))
  p <- add_border(p, row = nrow(out), where = "bottom", gp = gpar(lwd = TH, col = "black"))
  wh <- get_wh(p)
  f <- file.path(outdir, sprintf("Figure %s_subgroup %s.png", fig_no[[nm]], sites[[nm]]))
  agg_png(f, width = wh[1] + 0.08, height = wh[2] + 0.08, units = "in", res = 600, background = "white")
  plot(p); invisible(dev.off())
  cat(sprintf("%-9s rows %d (strata %d, run %d, NA %d) size %.2f x %.2f in\n", nm, nrow(df), sum(!df$hdr),
              sum(ok), sum(df$ci == "NA"), wh[1], wh[2]))
}
which_sites <- strsplit(Sys.getenv("SITES", paste(names(sites), collapse = ",")), ",")[[1]]
for (nm in which_sites) build(nm)
