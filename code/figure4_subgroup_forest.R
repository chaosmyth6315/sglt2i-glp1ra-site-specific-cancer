# ---------------------------------------------------------------
# Figure 4 (main text): subgroup forest plots.
# Four panels: kidney and bladder cancer (full cohort), breast cancer (women),
#   prostate cancer (men). Forest plots only, without numbers: HR, CI, N, events
#   and P for interaction are shown in Supplementary Figures S5–S8.
#   "*" on a group row = P for interaction < 0.05 (two strata: z test on the difference
#   of the log HRs; race: Cochran's Q), defined in the legend.
#   "Not applicable" = sex strata in the cohorts of women and men;
#   "NA" = analysis blocked by the platform's complexity limit (full cohort, non-smokers).
# Laid out at print size: 178 mm wide, Arial 7–8 pt; log x axis 0.6–1.6 (a CI that runs
#   beyond the axis ends in an arrow); dashed reference line at 1 (grey35).
# Data: fig4_data.csv from figure4_s4_s8_subgroup_data.py.
# Usage: FIGOUT=<dir> FIG4DATA=<csv> Rscript figure4_subgroup_forest.R (en_US.UTF-8 locale)
# ---------------------------------------------------------------
library(ggplot2); library(patchwork); library(ragg)
pdf(NULL)   # null device for text measurement, so no Rplots.pdf is written
outdir <- Sys.getenv("FIGOUT")
dat <- read.csv(Sys.getenv("FIG4DATA"), stringsAsFactors = FALSE, check.names = FALSE,
                encoding = "UTF-8", colClasses = "character")
cancers <- c("Kidney", "Bladder", "Breast", "Prostate")
titles  <- c(Kidney = "Kidney cancer", Bladder = "Bladder cancer",
             Breast = "Breast cancer (women)", Prostate = "Prostate cancer (men)")
pt <- function(x) x / ggplot2::.pt
XL <- c(0.6, 1.6); XB <- c(0.7, 1, 1.4)

# row skeleton in the order of the kidney (full-cohort) rows
k <- dat[dat$cancer == "Kidney", ]
skel <- data.frame(group = character(), stratum = character(), hdr = logical(), stringsAsFactors = FALSE)
for (g in unique(k$group)) {
  skel <- rbind(skel, data.frame(group = g, stratum = "", hdr = TRUE, stringsAsFactors = FALSE))
  for (s in k$stratum[k$group == g])
    skel <- rbind(skel, data.frame(group = g, stratum = s, hdr = FALSE, stringsAsFactors = FALSE))
}
nr <- nrow(skel)
skel$y <- rev(seq_len(nr))

pnum <- function(s) ifelse(s == "<0.001", 0.0005, suppressWarnings(as.numeric(s)))

panel <- function(nm, first) {
  d <- dat[dat$cancer == nm, ]
  idx <- match(paste(skel$group, skel$stratum), paste(d$group, d$stratum))
  r <- skel
  r$hr <- as.numeric(d$hr[idx]); r$lo <- as.numeric(d$lo[idx]); r$hi <- as.numeric(d$hi[idx])
  r$status <- d$status[idx]
  r$txt <- ifelse(!r$hdr & is.na(idx), "Not applicable", ifelse(!r$hdr & r$status %in% "not_run", "NA", ""))
  r$hr[r$hdr | r$txt != ""] <- NA
  r$sig <- FALSE
  for (i in which(r$hdr)) {
    v <- unique(d$p_int[d$group == r$group[i]])
    stopifnot(length(v) <= 1)
    r$sig[i] <- length(v) == 1 && v != "" && pnum(v) < 0.05
  }
  ok <- !is.na(r$hr)
  stopifnot(all(r$hr[ok] > XL[1] & r$hr[ok] < XL[2]))
  r$lo_c <- pmax(r$lo, XL[1] * 1.001); r$hi_c <- pmin(r$hi, XL[2] / 1.001)
  r$cut_lo <- ok & r$lo <= XL[1]; r$cut_hi <- ok & r$hi >= XL[2]
  cat(sprintf("%-9s arrows: %s\n", nm, paste(c(paste0(r$stratum[r$cut_lo], "<"), paste0(r$stratum[r$cut_hi], ">"))[
      c(r$cut_lo[r$cut_lo], r$cut_hi[r$cut_hi])], collapse = ", ")))
  cap <- 0.17
  cat(sprintf("%-9s plotted %d strata; not applicable %d; NA %d; * %s\n", nm, sum(ok), sum(r$txt == "Not applicable"),
              sum(r$txt == "NA"), paste(r$group[r$sig], collapse = ", ")))
  ggplot(r, aes(y = y)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey35", linewidth = 0.3) +
    geom_segment(data = r[ok, ], aes(x = lo_c, xend = hi_c, yend = y), linewidth = 0.3) +
    geom_segment(data = r[ok & !r$cut_lo, ], aes(x = lo, xend = lo, y = y - cap, yend = y + cap), linewidth = 0.3) +
    geom_segment(data = r[ok & !r$cut_hi, ], aes(x = hi, xend = hi, y = y - cap, yend = y + cap), linewidth = 0.3) +
    geom_segment(data = r[r$cut_lo, ], aes(x = hr, xend = lo_c, yend = y), linewidth = 0.3,
                 arrow = arrow(length = unit(1.2, "mm"), type = "closed", angle = 25)) +
    geom_segment(data = r[r$cut_hi, ], aes(x = hr, xend = hi_c, yend = y), linewidth = 0.3,
                 arrow = arrow(length = unit(1.2, "mm"), type = "closed", angle = 25)) +
    geom_point(data = r[ok, ], aes(x = hr), shape = 15, size = 1.1) +
    geom_label(data = r[r$txt != "", ], aes(x = 1, label = txt), size = pt(7), family = "Arial",
               fill = "white", label.size = 0, label.padding = unit(0.6, "mm")) +
    geom_text(data = r[r$sig, ], aes(x = XL[2] * 0.95, label = "*"), size = pt(10), family = "Arial",
              vjust = 0.75) +
    annotate("text", x = 0.95, y = -2.1, label = "Favors\nSGLT2i", hjust = 1, size = pt(7),
             family = "Arial", lineheight = 0.95) +
    annotate("text", x = 1.05, y = -2.1, label = "Favors\nGLP-1 RA", hjust = 0, size = pt(7),
             family = "Arial", lineheight = 0.95) +
    scale_x_log10(limits = XL, breaks = XB, labels = c("0.7", "1", "1.4"), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim = c(0.4, nr + 0.6), clip = "off") +
    labs(title = titles[nm], x = NULL, y = NULL) +
    theme_classic(base_family = "Arial", base_size = 8) +
    theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5, margin = margin(b = 4)),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
          axis.text.x = element_text(size = 7, colour = "black"),
          axis.line.x = element_line(linewidth = 0.3), axis.ticks.x = element_line(linewidth = 0.3),
          plot.margin = margin(2, 9, 28, 9))
}

lab <- ggplot(skel, aes(y = y)) +
  geom_text(aes(x = ifelse(hdr, 0, 0.06), label = ifelse(hdr, group, stratum),
                fontface = ifelse(hdr, "bold", "plain")), hjust = 0, size = pt(7), family = "Arial") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(ylim = c(0.4, nr + 0.6), clip = "off") +
  labs(title = "Subgroup", x = NULL, y = NULL) +
  theme_void(base_family = "Arial", base_size = 8) +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0, margin = margin(b = 4)),
        plot.margin = margin(2, 2, 28, 2))

ps <- c(list(lab), lapply(seq_along(cancers), function(i) panel(cancers[i], i == 1)))
p <- wrap_plots(ps, nrow = 1, widths = c(1.15, 1, 1, 1, 1))
agg_png(file.path(outdir, "Figure 4_subgroup.png"), width = 7.0, height = 6.4, units = "in",
        res = 600, background = "white")
print(p); dev.off()
cat("fig4 OK\n")
