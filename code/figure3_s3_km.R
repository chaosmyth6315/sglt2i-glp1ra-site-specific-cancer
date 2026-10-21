#################################################################
# Figure 3 (4-panel Kaplan–Meier plot) and Supplementary Figure S3 (8-panel Kaplan–Meier plot).
# Inputs (written by figure3_s3_km_data.py): KMDIR/km_manifest.csv (panel -> exported KM csv, log-rank p, HR)
#   and KMDIR/km_risk.csv (number at risk per arm at years 1–10, recovered from the KM
#   confidence intervals in the export; see lib_km_number_at_risk.py).
# Style: exposure (SGLT2i) red solid line with a pink 95% CI band; comparator (GLP-1 RA)
#   black dashed line with a light-blue band; bands are drawn under the lines. The same y axis is
#   used on every panel. x axis from year 1 (start of the outcome window) to year 10.
#   Number-at-risk table under each panel. Arial, 600 dpi.
# Layout: 2 columns at print size (178 mm wide); numbers at risk every 2 years (2–10) in black
#   at 7 pt; one shared legend; x title on the bottom row, y title on the left column.
#   The BUILD environment variable selects '3', 'S3' or both (default both).
# Usage: FIGOUT=<dir> KMDIR=<dir> Rscript figure3_s3_km.R   (en_US.UTF-8 locale)
#################################################################
library(ggplot2); library(zoo); library(ragg); library(patchwork)
pdf(NULL)   # null device for text measurement, so no Rplots.pdf is written
outdir <- Sys.getenv("FIGOUT"); kmdir <- Sys.getenv("KMDIR")
man  <- read.csv(file.path(kmdir, "km_manifest.csv"), stringsAsFactors = FALSE, encoding = "UTF-8")
risk <- read.csv(file.path(kmdir, "km_risk.csv"), stringsAsFactors = FALSE, encoding = "UTF-8")

read_km <- function(path){
  lines <- readLines(path, encoding = "UTF-8")
  hi <- grep("^Time \\(Days\\)", lines)[1]
  d <- read.csv(text = paste(lines[hi:length(lines)], collapse = "\n"), check.names = FALSE)
  names(d) <- c("Time", "case", "case_lower", "case_upper", "control", "control_lower", "control_upper")
  cols <- names(d)[-1]
  d[cols] <- lapply(d[cols], function(x) na.locf(as.numeric(x), na.rm = FALSE))
  d$years <- d$Time / 365.25
  d
}
YMIN <- 0.92   # common y axis for all 12 panels (lowest 10-year survival: prostate cancer ≈ 0.94)
labels <- c(case = "SGLT2i", control = "GLP-1 RA")

COL <- c(case = "red", control = "black")
RISK_YEARS <- c(2, 4, 6, 8, 10)
pt <- function(x) x / ggplot2::.pt   # font size in pt -> geom_text size

km_panel <- function(m, last = TRUE, left = TRUE){
  d <- read_km(m$km_csv)
  stopifnot(min(c(d$case_lower, d$control_lower), na.rm = TRUE) > YMIN)
  title <- paste0("(", m$panel, ") ", m$label)
  if (nchar(title) > 44) title <- paste(strwrap(title, width = 44), collapse = "\n")
  p <- ggplot(d, aes(x = years)) +
    geom_ribbon(aes(ymin = control_lower, ymax = control_upper, fill = "control"), alpha = 0.5) +
    geom_ribbon(aes(ymin = case_lower, ymax = case_upper, fill = "case"), alpha = 0.5) +
    geom_line(aes(y = control, color = "control", linetype = "control"), linewidth = 0.4) +
    geom_line(aes(y = case, color = "case", linetype = "case"), linewidth = 0.4) +
    scale_color_manual(name = NULL, values = COL, labels = labels, breaks = c("case", "control")) +
    scale_fill_manual(name = NULL, values = c(case = "pink", control = "lightblue"), labels = labels,
                      breaks = c("case", "control")) +
    scale_linetype_manual(name = NULL, values = c(case = "solid", control = "dashed"), labels = labels,
                          breaks = c("case", "control")) +
    scale_x_continuous(limits = c(0.55, 10.45), breaks = RISK_YEARS, expand = c(0, 0)) +
    scale_y_continuous(limits = c(YMIN, 1), breaks = seq(YMIN, 1, by = 0.02)) +
    labs(x = if (last) "Years since the index date" else NULL, y = if (left) "Cancer-free probability" else NULL, title = title) +
    annotate("text", x = 0.75, y = YMIN + 0.003, label = paste0(m$hr, "\n", m$logrank), hjust = 0, vjust = 0,
             size = pt(7.5), family = "Arial", lineheight = 1.05) +
    theme_classic(base_family = "Arial", base_size = 8) +
    theme(plot.title = element_text(hjust = 0, size = 8.5, face = "bold"),
          axis.title = element_text(size = 8), axis.text = element_text(size = 7.5, colour = "black"),
          axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35),
          legend.position = "top", legend.text = element_text(size = 8),
          legend.key.width = unit(1.0, "cm"), legend.key.height = unit(0.3, "cm"),
          plot.margin = margin(4, 16, 0, 2))
  r <- risk[risk$label == m$label & risk$year %in% RISK_YEARS, ]
  stopifnot(nrow(r) == 2 * length(RISK_YEARS))
  r$arm <- factor(r$arm, levels = c("G", "S"), labels = c("GLP-1 RA", "SGLT2i"))
  tbl <- ggplot(r, aes(x = year, y = arm, label = format(n, big.mark = ","))) +
    geom_text(size = pt(7), family = "Arial") +
    scale_x_continuous(limits = c(0.55, 10.45), breaks = RISK_YEARS, expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL, title = "Number at risk") +
    theme_void(base_family = "Arial", base_size = 8) +
    theme(axis.text.y = element_text(size = 7.5, hjust = 1, colour = "black", margin = margin(r = 4)),
          plot.title = element_text(size = 7.5, hjust = 0, margin = margin(b = 2)),
          plot.margin = margin(2, 16, 8, 2))
  list(p, tbl)
}

compose3 <- function(ms, file, width, height = 6.6, kmtbl = c(4, 1)){
  n <- nrow(ms)
  panels <- lapply(seq_len(n), function(i) {
    pr <- km_panel(ms[i, ], last = i > n - 2, left = i %% 2 == 1)
    pr[[1]] / pr[[2]] + plot_layout(heights = kmtbl)
  })
  p <- wrap_plots(c(list(guide_area()), list(wrap_plots(panels, ncol = 2))), ncol = 1,
                  heights = c(0.25, 10), guides = "collect") &
       theme(legend.position = "top", legend.direction = "horizontal")
  agg_png(file.path(outdir, file), width = width, height = height, units = "in", res = 600, background = "white")
  print(p); dev.off()
}

BUILD <- Sys.getenv("BUILD", "3,S3")
if (grepl("(^|,)3($|,)", BUILD)) { compose3(man[man$figure == "3", ], "Figure 3_KM.png", 7.0); cat("fig3 OK\n") }
if (grepl("S3", BUILD)) {
  # Figure S3 uses the Figure 3 layout (2 columns, 178 mm, at risk every 2 years, shared legend)
  compose3(man[man$figure == "S3", ], "Figure S3_KM other outcomes.png", 7.0, 8.8, kmtbl = c(3, 1))
  cat("s3 OK\n")
}
