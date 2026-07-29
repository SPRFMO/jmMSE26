#!/usr/bin/env Rscript

# SC14 trade-off figures for the 500-iteration annual-change candidate set.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(ggthemes)
})

jm_root <- Sys.getenv("JMMSE_ROOT", "/Users/jim/_mymods/sprfmo/jmMSE")
out_dir <- file.path("doc", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source_files <- list(
  reference = file.path(
    "output", "candidate-performance", "reference", "performance_with_vb.rds"
  ),
  robustness = file.path(
    "output", "candidate-performance", "robustness", "performance_with_vb.rds"
  )
)

missing_files <- names(source_files)[!file.exists(unlist(source_files))]
if (length(missing_files)) {
  stop("Missing candidate performance files: ",
    paste(missing_files, collapse = ", "))
}

read_results <- function(keys) {
  rbindlist(lapply(keys, function(key) {
    ans <- as.data.table(readRDS(source_files[[key]]))
    ans[, source_set := key]
    ans
  }), use.names = TRUE)
}

labels <- data.table(
  mp_code = paste0("tun", c(29, 43, 45, 47, 32, 44, 46, 48)),
  mp = paste0("MP", c(29, 43, 45, 47, 32, 44, 46, 48)),
  short_label = c(
    "HS+20", "HS-20", "HSsym", "HS-30",
    "PR+20", "PR-20", "PRsym", "PR-30"
  ),
  family = rep(c("Hockey-stick", "Power ramp"), each = 4)
)
labels[, display := paste0(short_label, " (", mp, ")")]
labels[, display := factor(display, levels = display)]

om_labels <- c(
  h1_0.16 = "Reference (om11)",
  h1_0.16_selex = "Alt. selectivity (om11_1)",
  h1_0.16_lowrec = "Recruitment crash (om11_2)",
  h1_0.16_cycle = "Cyclic recruitment (om11_3)",
  h1_0.16h = "Steepness 0.8 (om12)",
  h1_1.14 = "Model 1.14 (om13)",
  h2_0.16 = "Two stock (om21)",
  h2_0.16_mov = "Movement (om21_1)",
  h2_0.16h = "Two stock, h=0.8 (om22)",
  h2_1.14 = "Two stock, model 1.14 (om23)"
)

summarize_tradeoff <- function(dat) {
  dat[, mp_code := sub("_om.*$", "", mp)]
  dat <- dat[mp_code %in% labels$mp_code &
    statistic %in% c("C", "SBMSY", "IACC") &
    ((statistic == "C" & year %in% 2026:2030) |
      (statistic %in% c("SBMSY", "IACC") & year %in% 2036:2040))]

  by_iter <- dat[, .(value = mean(data, na.rm = TRUE)),
    by = .(om, biol, mp_code, iter, statistic)]
  wide <- dcast(
    by_iter, om + biol + mp_code + iter ~ statistic,
    value.var = "value"
  )
  setnames(
    wide,
    c("C", "IACC", "SBMSY"),
    c("catch_near", "iacc_long", "sb_long")
  )

  ans <- wide[, .(
    catch_mean = mean(catch_near, na.rm = TRUE),
    catch_q10 = quantile(catch_near, 0.10, na.rm = TRUE),
    catch_q90 = quantile(catch_near, 0.90, na.rm = TRUE),
    sb_mean = mean(sb_long, na.rm = TRUE),
    sb_q10 = quantile(sb_long, 0.10, na.rm = TRUE),
    sb_q90 = quantile(sb_long, 0.90, na.rm = TRUE),
    iacc_mean = mean(iacc_long, na.rm = TRUE),
    iacc_q10 = quantile(iacc_long, 0.10, na.rm = TRUE),
    iacc_q90 = quantile(iacc_long, 0.90, na.rm = TRUE),
    n_replicates = uniqueN(iter)
  ), by = .(om, biol, mp_code)]

  ans <- labels[ans, on = "mp_code"]
  ans[, om_panel := factor(
    unname(om_labels[om]),
    levels = unname(om_labels)
  )]
  ans[]
}

reference <- summarize_tradeoff(read_results("reference"))
robustness <- summarize_tradeoff(read_results("robustness"))
tradeoff <- rbind(reference, robustness, use.names = TRUE)

if (any(tradeoff$n_replicates != 500L)) {
  stop("Expected 500 iterations in every SC14 trade-off summary")
}

fwrite(
  tradeoff[order(om_panel, biol, family, display)],
  file.path("doc", "data", "sc14-preliminary-tradeoffs.csv")
)

palette <- c(`Hockey-stick` = "#0f4c8a", `Power ramp` = "#ea9b1f")

tradeoff_plot <- function(dat, facet = NULL, x_label = NULL) {
  p <- ggplot(
    dat,
    aes(catch_mean, sb_mean, colour = family, shape = family)
  ) +
    geom_hline(yintercept = 1, linewidth = 0.35, linetype = 2,
      colour = "grey55") +
    geom_errorbar(
      aes(ymin = sb_q10, ymax = sb_q90),
      width = 0, linewidth = 0.35, alpha = 0.55
    ) +
    geom_errorbar(
      aes(xmin = catch_q10, xmax = catch_q90),
      orientation = "y", width = 0, linewidth = 0.35, alpha = 0.55
    ) +
    geom_point(size = 2.5) +
    geom_text_repel(
      aes(label = short_label),
      size = 2.7, seed = 20260724, box.padding = 0.25,
      point.padding = 0.15, min.segment.length = 0,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    scale_colour_manual(values = palette) +
    labs(
      x = if (is.null(x_label)) {
        "Mean catch, 2026-2030 (kt)"
      } else {
        x_label
      },
      y = expression("Mean " * SB/SB[MSY] * ", 2036-2040"),
      colour = "HCR family",
      shape = "HCR family",
      caption = paste(
        "Points are means; bars show 10th-90th percentiles among 500",
        "posterior iterations. Dashed line: SB/SBMSY = 1."
      )
    ) +
    ggthemes::theme_few(base_size = 10) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      plot.caption = element_text(size = 8, colour = "grey35")
    )
  if (!is.null(facet)) {
    p <- p + facet
  }
  p
}

ggsave(
  file.path(out_dir, "sc14-tradeoff-reference.png"),
  tradeoff_plot(reference),
  width = 7.2, height = 4.7, dpi = 180, bg = "white"
)

biomass_iacc_plot <- ggplot(
  reference,
  aes(sb_mean, iacc_mean, colour = family, shape = family)
) +
  geom_vline(xintercept = 1, linewidth = 0.35, linetype = 2,
    colour = "grey55") +
  geom_errorbar(
    aes(ymin = iacc_q10, ymax = iacc_q90),
    width = 0, linewidth = 0.35, alpha = 0.55
  ) +
  geom_errorbar(
    aes(xmin = sb_q10, xmax = sb_q90),
    orientation = "y", width = 0, linewidth = 0.35, alpha = 0.55
  ) +
  geom_point(size = 2.5) +
  geom_text_repel(
    aes(label = short_label),
    size = 2.7, seed = 20260724, box.padding = 0.25,
    point.padding = 0.15, min.segment.length = 0,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = palette) +
  labs(
    x = expression("Mean " * SB/SB[MSY] * ", 2036-2040"),
    y = "Mean IACC, 2036-2040 (%)",
    colour = "HCR family",
    shape = "HCR family",
    caption = paste(
      "IACC is the absolute year-to-year percentage change in catch.",
      "\nPoints are means; bars show 10th-90th percentiles among 500",
      "posterior iterations. Dashed line: SB/SBMSY = 1."
    )
  ) +
  ggthemes::theme_few(base_size = 10) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.caption = element_text(size = 8, colour = "grey35")
  )

ggsave(
  file.path(out_dir, "sc14-tradeoff-reference-biomass-iacc.png"),
  biomass_iacc_plot,
  width = 7.2, height = 4.7, dpi = 180, bg = "white"
)

single_stock <- robustness[biol == "CJM"]
ggsave(
  file.path(out_dir, "sc14-tradeoff-single-stock-robustness.png"),
  tradeoff_plot(single_stock, facet_wrap(~om_panel, scales = "free", ncol = 2)),
  width = 8.2, height = 8.2, dpi = 180, bg = "white"
)

two_stock <- robustness[biol %in% c("Southern", "North")]
ggsave(
  file.path(out_dir, "sc14-tradeoff-two-stock-robustness.png"),
  tradeoff_plot(
    two_stock,
    facet_grid(biol ~ om_panel, scales = "free"),
    x_label = "Mean component catch, 2026-2030 (kt)"
  ),
  width = 11.5, height = 6.5, dpi = 180, bg = "white"
)

message("Wrote ", nrow(tradeoff), " OM-stock-CMP summaries and four figures.")
