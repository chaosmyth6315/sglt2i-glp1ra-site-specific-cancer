#################################################################
# Figure 2: forest plot of site-specific cancer HRs, SGLT2i versus GLP-1 RA
# Input: fig2_data.csv written by figure2_data.py from the four TriNetX
#   exports under <export root> (Main_Full, Female_Breast,
#   Female_Endo_noHyst, Male). Every displayed string is rounded once
#   there (ROUND_HALF_UP); this script only lays the table out.
#   p = Cox Wald p (from HR and 95% CI); q = Benjamini–Hochberg across the
#   12 sites; † = proportional-hazards violation (Schoenfeld p < 0.05);
#   CIR % = 100 × (1 − KM survival) at day 3,650; E-value from unrounded HR/CI.
# Layout: print size (at most 178 mm wide, 8 pt Arial), log x axis 0.7–1.4,
#   dashed reference line at 1, 'SGLT2i versus GLP-1 RA' title row, cohort
#   rows as '<cohort> (<n> matched pairs)', three-line rules, † set as a
#   superscript after the HR (95% CI).
# Output: 'Figure 2_Clinical Outcomes.png' (600 dpi) in FIGOUT.
# Usage: FIGOUT=<dir> FIG2DATA=<csv> Rscript figure2_forest.R
#   (LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8; ragg needs an ASCII output path)
#################################################################
library(forestploter); library(grid); library(ragg)
pdf(NULL)   # null device for text measurement, so no Rplots.pdf is written
DAG <- "†"
outdir <- Sys.getenv("FIGOUT")
dat <- read.csv(Sys.getenv("FIG2DATA"), stringsAsFactors = FALSE, check.names = FALSE,
                encoding = "UTF-8")
stopifnot(nrow(dat) == 12)

blocks <- unique(dat$block)
blank <- function(lab) data.frame(Outcomes = lab, ev1 = "", ev2 = "", hr = NA, lo = NA, hi = NA,
                                  hrci = "", p = "", q = "", ev = "", hdr = TRUE, stringsAsFactors = FALSE)
rows <- list(blank("SGLT2i versus GLP-1 RA"))
for (b in blocks) {
  sub <- dat[dat$block == b, ]
  n <- format(sub$n_per_arm[1], big.mark = ",")
  name <- if (b == "All patients") "Full cohort" else b
  rows[[length(rows) + 1]] <- blank(sprintf("%s (%s matched pairs)", name, n))
  for (i in seq_len(nrow(sub))) {
    r <- sub[i, ]
    rows[[length(rows) + 1]] <- data.frame(Outcomes = r$label, ev1 = r$ev1_s, ev2 = r$ev2_s,
                                           hr = r$hr, lo = r$lo, hi = r$hi, hrci = r$hr_ci_s,
                                           p = r$p_s, q = r$q_s, ev = r$ev_s, hdr = FALSE,
                                           stringsAsFactors = FALSE)
  }
}
df <- do.call(rbind, rows)
hdr_rows <- which(df$hdr)

wrap_outcome <- function(x, is_sub, width_main = 36, width_sub = 30, indent = "   ") {
  mapply(function(txt, sub) {
    w <- strwrap(txt, width = if (sub) width_sub else width_main)
    if (sub) w <- paste0(indent, w)
    paste(w, collapse = "\n")
  }, x, is_sub, USE.NAMES = FALSE)
}
df$Outcomes <- wrap_outcome(df$Outcomes, !df$hdr)

df$`SGLT2i\nEvents (CIR, %)` <- df$ev1
df$`GLP-1 RA\nEvents (CIR, %)` <- df$ev2
df$`HR (95% CI)` <- df$hrci
df$`p value` <- df$p
df$`q value` <- df$q
df$`E-value` <- df$ev
df$` ` <- paste(rep(" ", 28), collapse = "")

show_cols <- c("Outcomes", "SGLT2i\nEvents (CIR, %)", "GLP-1 RA\nEvents (CIR, %)", " ",
               "HR (95% CI)", "p value", "q value", "E-value")

tm <- forest_theme(base_size = 8, base_family = "Arial",
  refline_gp = gpar(lwd = 0.75, lty = "dashed", col = "grey35"),
  ci_pch = 15, ci_col = "black", ci_fill = "black", ci_lwd = 0.9, ci_Theight = 0.25,
  arrow_type = "closed", arrow_length = 0.06, arrow_gp = gpar(fontsize = 7.5),
  xaxis_gp = gpar(fontsize = 7.5, lwd = 0.75),
  core = list(padding = unit(c(1.3, 2.6), "mm"), bg_params = list(fill = "white")),
  colhead = list(fg_params = list(hjust = 0, x = 0.02, fontface = "bold"),
                 padding = unit(c(1.6, 2.4), "mm")))

p <- forest(df[, show_cols],
  est = df$hr, lower = df$lo, upper = df$hi,
  sizes = 0.45, ci_column = 4,
  ref_line = 1, x_trans = "log", arrow_lab = c("Favors SGLT2i", "Favors GLP-1 RA"),
  xlim = c(0.7, 1.4), ticks_at = c(0.8, 1, 1.25), theme = tm)

TH <- 1.33; TN <- 0.67   # 1 pt and 0.5 pt rules
p <- edit_plot(p, row = hdr_rows, gp = gpar(fontface = "bold"))
p <- add_border(p, part = "header", row = 1, where = "top", gp = gpar(lwd = TH, col = "black"))
p <- add_border(p, part = "header", row = 1, where = "bottom", gp = gpar(lwd = TN, col = "black"))
p <- add_border(p, row = nrow(df), where = "bottom", gp = gpar(lwd = TH, col = "black"))

wh <- get_wh(p)
agg_png(file.path(outdir, "Figure 2_Clinical Outcomes.png"),
        width = wh[1] + 0.08, height = wh[2] + 0.08, units = "in", res = 600, background = "white")
plot(p)
# superscript trailing daggers
invisible(grid.force())
esc <- function(s) { s <- gsub("\\\\", "\\\\\\\\", s); gsub("'", "\\\\'", s, fixed = TRUE) }
make_expr <- function(s) {
  if (is.na(s) || s == "") return(parse(text = "''")[[1]])
  if (endsWith(s, DAG)) {
    n <- esc(substr(s, 1, nchar(s) - nchar(DAG)))
    parse(text = paste0("'", n, "'^'", DAG, "'"))[[1]]
  } else parse(text = paste0("'", esc(s), "'"))[[1]]
}
for (n in grid.ls(print = FALSE)$name) {
  obj <- try(grid.get(n), silent = TRUE)
  if (inherits(obj, "try-error") || !inherits(obj, "text") || is.null(obj$label)) next
  lab <- as.character(obj$label)
  if (!any(endsWith(lab, DAG))) next
  grid.edit(n, label = as.expression(lapply(lab, make_expr)))
}
dev.off()
cat("fig2 OK", wh, "\n")
