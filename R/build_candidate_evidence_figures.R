suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(ggthemes)
  library(ggh4x)
})

output_dir <- file.path("doc", "figures", "candidates")
summary_dir <- file.path("doc", "data", "candidates")
performance_root <- Sys.getenv("JMMSE_PERFORMANCE_ROOT",
  unset = file.path("output", "candidate-performance-500"))
reference_runs_file <- Sys.getenv("JMMSE_REFERENCE_RUNS",
  unset = file.path("..", "jmMSE-500-refine", "model", "tune",
    "refine_500_from_100", "runs.rds"))
reference_om_file <- Sys.getenv("JMMSE_REFERENCE_OM",
  unset = file.path("..", "jmMSE-500-refine", "data",
    "om11_h1_0.16_065.rds"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

candidate_labels <- c(
  tun29 = "HS+20", tun43 = "HS-20", tun45 = "HSsym", tun47 = "HS-30",
  tun32 = "PR+20", tun44 = "PR-20", tun46 = "PRsym", tun48 = "PR-30"
)
candidate_label_levels <- unname(candidate_labels[c(
  "tun29", "tun43", "tun45", "tun47",
  "tun32", "tun44", "tun46", "tun48"
)])

candidate_vb_plot <- function(mode) {
  input_file <- file.path(performance_root, mode,
    "performance_with_vb.rds")
  if (!file.exists(input_file)) stop("Missing candidate results: ", input_file)

  perf <- as.data.table(readRDS(input_file))
  needed <- c("VB2025", "VBMSY")
  if (!all(needed %in% unique(perf$statistic))) {
    stop("Missing vulnerable-biomass statistics in ", input_file)
  }

  # The robustness vulnerable-biomass boxplot is a like-for-like summary of
  # the five single-stock CJM robustness OMs. North and Southern components
  # from the four two-stock OMs are examined in the stock-specific figures.
  vb_perf <- if (mode == "robustness") perf[biol == "CJM"] else perf

  periods <- rbindlist(list(
    vb_perf[statistic %in% needed & year %in% 2026:2035][,
      .(data = mean(data)), by = .(om, biol, mp, statistic, iter)][,
        period := "Short term (2026-2035)"],
    vb_perf[statistic %in% needed & year %in% 2041:2050][,
      .(data = mean(data)), by = .(om, biol, mp, statistic, iter)][,
        period := "Long term (2041-2050)"]
  ))
  periods[, metric := factor(statistic, levels = needed,
    labels = c("VB / VB[2025]", "VB / VB[MSY]"))]
  periods[, mp := factor(candidate_labels[sub("_.*$", "", mp)],
    levels = candidate_label_levels)]

  summary <- periods[, .(
    median = median(data),
    lower = quantile(data, 0.05),
    upper = quantile(data, 0.95)
  ), by = .(om, biol, mp, period, metric)]
  fwrite(summary, file.path(summary_dir, paste0("vb_", mode, "_summary.csv")))

  y_caps <- periods[, .(
    upper = max(
      quantile(data, 0.75, na.rm = TRUE) +
        1.5 * IQR(data, na.rm = TRUE),
      na.rm = TRUE
    )
  ), by = metric]
  vb2025_cap <- y_caps[metric == "VB / VB[2025]", upper] * 1.03
  vbmsy_cap <- y_caps[metric == "VB / VB[MSY]", upper] * 1.03

  plot <- ggplot(periods, aes(mp, data, fill = mp)) +
    geom_boxplot(outlier.shape = NA, linewidth = 0.25) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey45") +
    facet_grid(metric ~ period, scales = "free_y") +
    labs(
      title = paste("Candidate vulnerable-biomass performance:", mode),
      x = NULL, y = "Ratio", fill = NULL
    ) +
    ggh4x::facetted_pos_scales(y = list(
      metric == "VB / VB[2025]" ~
        scale_y_continuous(limits = c(0, vb2025_cap)),
      metric == "VB / VB[MSY]" ~
        scale_y_continuous(limits = c(0, vbmsy_cap))
    )) +
    ggthemes::theme_few(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  ggsave(file.path(output_dir, paste0("vb_", mode, ".png")), plot,
    width = 9, height = 6, dpi = 160)

  near_term <- perf[
    statistic %in% c("C", "VBMSY") & year %in% 2026:2035,
    .(data = mean(data)),
    by = .(om, biol, mp, statistic, iter)
  ]
  near_term[, scenario := fifelse(
    grepl("^tun[0-9]+_", mp),
    sub("^tun[0-9]+_", "", mp),
    "reference"
  )]
  near_term[, mp := factor(candidate_labels[sub("_.*$", "", mp)],
    levels = candidate_label_levels)]
  paired <- dcast(near_term, scenario + om + biol + mp + iter ~ statistic,
    value.var = "data")
  tradeoff <- paired[, .(
    catch = median(C),
    catch_lower = quantile(C, 0.05),
    catch_upper = quantile(C, 0.95),
    vbmsy = median(VBMSY),
    vbmsy_lower = quantile(VBMSY, 0.05),
    vbmsy_upper = quantile(VBMSY, 0.95),
    n = .N
  ), by = .(biol, mp)]
  fwrite(tradeoff,
    file.path(summary_dir, paste0("catch_vbmsy_", mode, "_summary.csv")))

  tradeoff_plot <- ggplot(tradeoff, aes(catch, vbmsy, colour = mp)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey45") +
    geom_errorbar(aes(ymin = vbmsy_lower, ymax = vbmsy_upper),
      width = 0, linewidth = 0.35) +
    geom_errorbar(aes(xmin = catch_lower, xmax = catch_upper),
      orientation = "y", width = 0, linewidth = 0.35) +
    geom_point(size = 2.5) +
    geom_text_repel(aes(label = mp), size = 2.7, show.legend = FALSE,
      seed = 12000, min.segment.length = 0, max.overlaps = Inf) +
    facet_wrap(~biol, scales = "free") +
    labs(
      title = paste("Near-term catch and vulnerable-biomass status:", mode),
      subtitle = paste(
        "Points are MP medians; bars are marginal 5th-95th percentiles",
        "for 2026-2035 means"
      ),
      x = "Mean annual catch, 2026-2035 (model catch units)",
      y = "Mean VB / VB[MSY], 2026-2035",
      colour = NULL
    ) +
    ggthemes::theme_few(base_size = 10) +
    theme(legend.position = "none")

  ggsave(file.path(output_dir, paste0("catch_vbmsy_", mode, ".png")),
    tradeoff_plot, width = 9, height = 5.5, dpi = 160)

  long_term <- perf[
    statistic %in% c("C", "VB2025") & year %in% 2041:2050,
    .(data = mean(data)),
    by = .(om, biol, mp, statistic, iter)
  ]
  long_term[, scenario := fifelse(
    grepl("^tun[0-9]+_", mp),
    sub("^tun[0-9]+_", "", mp),
    "reference"
  )]
  long_term[, mp := factor(candidate_labels[sub("_.*$", "", mp)],
    levels = candidate_label_levels)]
  paired_long <- dcast(
    long_term,
    scenario + om + biol + mp + iter ~ statistic,
    value.var = "data")
  long_tradeoff <- paired_long[, .(
    catch = median(C),
    catch_lower = quantile(C, 0.05),
    catch_upper = quantile(C, 0.95),
    vb2025 = median(VB2025),
    vb2025_lower = quantile(VB2025, 0.05),
    vb2025_upper = quantile(VB2025, 0.95),
    n = .N
  ), by = .(biol, mp)]
  fwrite(long_tradeoff,
    file.path(summary_dir, paste0("catch_vb2025_", mode, "_summary.csv")))

  long_tradeoff_plot <- ggplot(long_tradeoff,
    aes(catch, vb2025, colour = mp)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey45") +
    geom_errorbar(aes(ymin = vb2025_lower, ymax = vb2025_upper),
      width = 0, linewidth = 0.35) +
    geom_errorbar(aes(xmin = catch_lower, xmax = catch_upper),
      orientation = "y", width = 0, linewidth = 0.35) +
    geom_point(size = 2.5) +
    geom_text_repel(aes(label = mp), size = 2.7, show.legend = FALSE,
      seed = 12000, min.segment.length = 0, max.overlaps = Inf) +
    facet_wrap(~biol, scales = "free") +
    labs(
      title = paste("Long-term catch and vulnerable biomass vs 2025:", mode),
      subtitle = paste(
        "Points are MP medians; bars are marginal 5th-95th percentiles",
        "for 2041-2050 means"
      ),
      x = "Mean annual catch, 2041-2050 (model catch units)",
      y = "Mean VB / VB[2025], 2041-2050",
      colour = NULL
    ) +
    ggthemes::theme_few(base_size = 10) +
    theme(legend.position = "none")

  ggsave(file.path(output_dir, paste0("catch_vb2025_", mode, ".png")),
    long_tradeoff_plot, width = 9, height = 5.5, dpi = 160)
}

candidate_vb_plot("reference")
candidate_vb_plot("robustness")

# Reference-set catch trajectories for the two asymmetric hockey-stick MPs.
# The shorthand refers to the binding side of each TAC-change limit:
# tun43 is -20%/+15%, while tun29 is -15%/+20%.
reference_file <- file.path(performance_root, "reference",
  "performance_with_vb.rds")
reference_perf <- as.data.table(readRDS(reference_file))

# Recover the dynamic biomass denominator used in CMP tuning. For the complete
# 500-draw workflow all eight runs are stored together; the fallback preserves
# compatibility with the archived split 100-draw artifacts.
if (nzchar(reference_runs_file)) {
  all_candidate_runs <- readRDS(reference_runs_file)
  candidate_runs <- all_candidate_runs@.Data[
    match(names(candidate_labels), names(all_candidate_runs))]
  names(candidate_runs) <- names(candidate_labels)
} else {
  candidate_runs <- c(
    readRDS(file.path("..", "jmMSE", "model", "candidates", "reference",
      "runs.rds"))[c("tun29", "tun43", "tun32", "tun44")]@.Data,
    readRDS(file.path("..", "jmMSE", "model", "candidates",
      "additional-change-limits", "reference", "runs.rds"))[
        c("tun45", "tun47", "tun46", "tun48")]@.Data
  )
  names(candidate_runs) <- c(
    "tun29", "tun43", "tun32", "tun44",
    "tun45", "tun47", "tun46", "tun48"
  )
}
unfished_ssb <- readRDS(reference_om_file)$unfishedSSB
flq_year_iter_table <- function(x, run_code, statistic_code) {
  years <- as.integer(dimnames(x)$year)
  iterations <- as.integer(dimnames(x)$iter)
  data.table(
    mp = run_code,
    iter = rep(iterations, each = length(years)),
    year = rep(years, times = length(iterations)),
    statistic = statistic_code,
    data = as.numeric(x)
  )
}
dynamic_kobe_annual <- rbindlist(lapply(names(candidate_runs), function(run_code) {
  projected_om <- mse::om(candidate_runs[[run_code]])
  metric <- FLCore::metrics(projected_om)$CJM
  rp <- FLCore::refpts(projected_om)
  years <- as.integer(dimnames(metric$SB)$year)
  sb0msy <- metric$SB /
    ((rp["SBMSY", ] / rp["SB0", ]) * unfished_ssb$CJM[, as.character(years)])
  ffmsy <- metric$F / rp["FMSY", ]
  rbindlist(list(
    flq_year_iter_table(sb0msy, run_code, "SB0MSY"),
    flq_year_iter_table(ffmsy, run_code, "FMSY")
  ))
}))

# Reference-OM Kobe plots for the near- and long-term reporting periods.
# Kobe convention: biomass status is horizontal and fishing pressure vertical.
kobe_mps <- c(
  "tun29", "tun43", "tun45", "tun47",
  "tun32", "tun44", "tun46", "tun48"
)
kobe_labels <- c(
  tun29 = "HS+20", tun43 = "HS-20", tun45 = "HSsym", tun47 = "HS-30",
  tun32 = "PR+20", tun44 = "PR-20", tun46 = "PRsym", tun48 = "PR-30"
)
kobe_periods <- rbindlist(list(
  dynamic_kobe_annual[
    mp %in% kobe_mps & statistic %in% c("SB0MSY", "FMSY") &
      year %in% 2026:2035,
    .(data = mean(data, na.rm = TRUE)),
    by = .(mp, iter, statistic)
  ][, period := "Near term (2026-2035)"],
  dynamic_kobe_annual[
    mp %in% kobe_mps & statistic %in% c("SB0MSY", "FMSY") &
      year %in% 2041:2050,
    .(data = mean(data, na.rm = TRUE)),
    by = .(mp, iter, statistic)
  ][, period := "Long term (2041-2050)"]
))
kobe_wide <- dcast(kobe_periods, mp + iter + period ~ statistic,
  value.var = "data")
kobe_summary <- kobe_wide[, .(
  f_fmsy = median(FMSY, na.rm = TRUE),
  sb_sbmsy = median(SB0MSY, na.rm = TRUE),
  n = sum(is.finite(FMSY) & is.finite(SB0MSY))
), by = .(mp, period)]
kobe_summary[, cmp := factor(kobe_labels[mp],
  levels = unname(kobe_labels[kobe_mps]))]
kobe_summary[, period := factor(period,
  levels = c("Near term (2026-2035)", "Long term (2041-2050)"))]
fwrite(kobe_summary,
  file.path(summary_dir, "kobe_reference_periods_summary.csv"))

kobe_quadrants <- rbindlist(lapply(
  levels(kobe_summary$period),
  function(period_name) data.table(
    period = period_name,
    quadrant = c("Green", "Yellow", "Orange", "Red"),
    xmin = c(1, 1, 0, 0),
    xmax = c(Inf, Inf, 1, 1),
    ymin = c(0, 1, 0, 1),
    ymax = c(1, Inf, 1, Inf),
    fill = c("#8BCB88", "#F2D46F", "#E7A35B", "#D97B72")
  )
))
kobe_quadrants[, period := factor(period, levels = levels(kobe_summary$period))]

kobe_plot <- ggplot() +
  geom_rect(
    data = kobe_quadrants,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
    alpha = 0.16, inherit.aes = FALSE
  ) +
  scale_fill_identity() +
  geom_vline(xintercept = 1, colour = "grey35", linewidth = 0.45) +
  geom_hline(yintercept = 1, colour = "grey35", linewidth = 0.45) +
  geom_point(
    data = kobe_summary,
    aes(sb_sbmsy, f_fmsy, colour = cmp),
    size = 2.8
  ) +
  geom_text_repel(
    data = kobe_summary,
    aes(sb_sbmsy, f_fmsy, label = cmp, colour = cmp),
    seed = 12000, size = 3, min.segment.length = 0,
    max.overlaps = Inf, show.legend = FALSE
  ) +
  facet_wrap(~period, nrow = 1) +
  coord_cartesian(xlim = c(0, 1.9), ylim = c(0, 1.35), expand = FALSE) +
  labs(
    title = "Reference-OM Kobe status by candidate management procedure",
    subtitle = paste(
      "Points are medians across 500 posterior iterations of each",
      "iteration's period mean"
    ),
    x = "Spawning biomass relative to dynamic SB[MSY]",
    y = "Fishing mortality relative to F[MSY]",
    colour = "CMP"
  ) +
  ggthemes::theme_few(base_size = 10) +
  theme(legend.position = "none")

ggsave(file.path(output_dir, "kobe_reference_periods.png"),
  kobe_plot, width = 9, height = 4.8, dpi = 160)

catch_trajectories <- reference_perf[
  statistic == "C" & mp %in% c("tun43", "tun29") & year %in% 2025:2050
]

available_iters <- sort(unique(catch_trajectories$iter))
if (length(available_iters) < 15L) {
  stop("Fewer than 15 reference-set iterations are available")
}
set.seed(12000)
selected_iters <- sort(sample(available_iters, 15L))
catch_trajectories <- catch_trajectories[iter %in% selected_iters]

ids_by_mp <- catch_trajectories[, .(ids = list(sort(unique(iter)))), by = mp]
if (nrow(ids_by_mp) != 2L ||
    !setequal(ids_by_mp[mp == "tun43", ids][[1]],
      ids_by_mp[mp == "tun29", ids][[1]]) ||
    !setequal(ids_by_mp[mp == "tun43", ids][[1]], selected_iters)) {
  stop("Spaghetti panels do not contain the same 15 iteration IDs")
}

mp_labels <- c(
  tun43 = "HS -20% (MP43; limits -20%/+15%)",
  tun29 = "HS +20% (MP29; limits -15%/+20%)"
)
catch_trajectories[, panel := factor(mp_labels[mp],
  levels = unname(mp_labels))]

fwrite(data.table(
  iter = selected_iters,
  hs_minus_20_run = "tun43",
  hs_minus_20_present = selected_iters %in%
    ids_by_mp[mp == "tun43", ids][[1]],
  hs_plus_20_run = "tun29",
  hs_plus_20_present = selected_iters %in%
    ids_by_mp[mp == "tun29", ids][[1]],
  projection_years_per_panel = 26L
),
  file.path(summary_dir, "catch_spaghetti_reference_iterations.csv"))

catch_spaghetti_plot <- ggplot(catch_trajectories,
  aes(year, data, group = interaction(mp, iter), colour = factor(iter))) +
  geom_line(linewidth = 0.70, alpha = 0.84) +
  facet_wrap(~panel, ncol = 1, scales = "fixed") +
  scale_colour_viridis_d(option = "turbo", guide = "none") +
  labs(
    title = "Projected catch trajectories: reference operating model",
    subtitle = paste(
      "The same 15 reproducibly selected iterations and matching colors",
      "are shown for each MP, 2025-2050"
    ),
    x = "Projection year",
    y = "Catch (model catch units)",
    colour = "Iteration"
  ) +
  ggthemes::theme_few(base_size = 10)

ggsave(file.path(output_dir, "catch_spaghetti_hs_reference.png"),
  catch_spaghetti_plot, width = 9, height = 8, dpi = 160)

# Long-term reference-set quilt. Scores are normalized independently within
# each metric and show relative performance among the eight annual-change
# variants listed in the report's naming-convention table.
additional_reference_file <- file.path(performance_root,
  "additional-change-limits", "reference", "performance_with_vb.rds")
if (file.exists(additional_reference_file)) {
  additional_reference_perf <- as.data.table(readRDS(additional_reference_file))
  quilt_perf <- rbindlist(list(
    reference_perf[mp %in% c("tun29", "tun43", "tun32", "tun44")],
    additional_reference_perf[mp %in% c("tun45", "tun47", "tun46", "tun48")]
  ), use.names = TRUE)
} else {
  # The 500-draw reference file already contains all eight candidates.
  quilt_perf <- reference_perf[mp %in% names(candidate_labels)]
}

dynamic_kobe_wide <- dcast(
  dynamic_kobe_annual[year %in% 2041:2050],
  mp + iter + year ~ statistic, value.var = "data")
dynamic_kobe <- dynamic_kobe_wide[, .(
  SB0green = mean(SB0MSY >= 1 & FMSY <= 1, na.rm = TRUE),
  SB0red = mean(SB0MSY < 1 & FMSY > 1, na.rm = TRUE)
), by = mp][, melt(.SD, id.vars = "mp", variable.name = "statistic",
  value.name = "value")]

catch_diagnostics <- quilt_perf[
  statistic == "C" & year %in% 2026:2050,
  .(catch = mean(data, na.rm = TRUE)),
  by = .(mp, iter, year)
][order(mp, iter, year)]
catch_diagnostics[, previous_catch := shift(catch), by = .(mp, iter)]
catch_diagnostics[, catch_reduction := fifelse(
  is.finite(previous_catch) & previous_catch > 0,
  pmax((previous_catch - catch) / previous_catch, 0) * 100,
  NA_real_
)]
catch_summary <- rbindlist(list(
  # Probability over all valid iteration-year outcomes in the long-term
  # window. For example, if 10% of iterations are below 270 kt in five of
  # ten years, the pooled probability is 0.10 * 0.50 = 0.05.
  quilt_perf[statistic == "C" & year %in% 2041:2050,
    .(value = mean(data < 270, na.rm = TRUE)), by = mp][
      , statistic := "PC270"],
  catch_diagnostics[year %in% 2026:2050,
    .(iteration_mean = mean(catch_reduction, na.rm = TRUE)), by = .(mp, iter)][
      , .(value = median(iteration_mean, na.rm = TRUE)), by = mp][
        , statistic := "Creduction"]
), use.names = TRUE)

metric_definitions <- data.table(
  statistic = c(
    "SBMSY", "FMSY", "C", "IACC", "VB2025", "VBMSY",
    "SB0red", "PC270", "Creduction"
  ),
  metric = c(
    "SB / SB[MSY]", "F / F[MSY]", "Catch", "IACC",
    "VB / VB[2025]", "VB / VB[MSY]", "P(Kobe red)",
    "P(catch < 270 kt)", "Mean catch reduction"
  ),
  direction = c(
    "higher is better", "lower is better", "higher is better",
    "lower is better", "higher is better", "higher is better",
    "lower is better", "lower is better", "lower is better"
  ),
  higher_is_better = c(TRUE, FALSE, TRUE, FALSE, TRUE, TRUE,
    FALSE, FALSE, FALSE)
)

base_quilt <- quilt_perf[
  year %in% 2041:2050,
  .(iteration_mean = mean(data, na.rm = TRUE)),
  by = .(mp, statistic, iter)
][
  , .(value = median(iteration_mean, na.rm = TRUE)),
  by = .(mp, statistic)
]
quilt <- rbindlist(list(base_quilt, dynamic_kobe, catch_summary),
  use.names = TRUE, fill = TRUE)[metric_definitions, on = "statistic"]

quilt[, `:=`(
  metric_min = min(value),
  metric_max = max(value)
), by = statistic]
quilt[, relative_badness := fifelse(
  metric_max == metric_min,
  0,
  fifelse(
    higher_is_better,
    (metric_max - value) / (metric_max - metric_min),
    (value - metric_min) / (metric_max - metric_min)
  )
)]
quilt[, `:=`(
  mp = factor(
    c(
      tun29 = "HS+20 (MP29)", tun43 = "HS-20 (MP43)",
      tun45 = "HSsym (MP45)", tun47 = "HS-30 (MP47)",
      tun32 = "PR+20 (MP32)", tun44 = "PR-20 (MP44)",
      tun46 = "PRsym (MP46)", tun48 = "PR-30 (MP48)"
    )[mp],
    levels = rev(c(
      "HS+20 (MP29)", "HS-20 (MP43)", "HSsym (MP45)", "HS-30 (MP47)",
      "PR+20 (MP32)", "PR-20 (MP44)", "PRsym (MP46)", "PR-30 (MP48)"
    ))
  ),
  metric = factor(metric, levels = metric_definitions$metric),
  years = "2041-2050",
  summary = "Median across iteration-specific 2041-2050 means"
)]
quilt[statistic == "SB0red",
  summary := "Proportion across all 2041-2050 iteration-years"]
quilt[statistic == "PC270",
  summary := "Proportion of 2041-2050 iteration-years with catch below 270 kt"]
quilt[statistic == "Creduction",
  summary := "Median across iteration-specific 2026-2050 mean reductions"]
quilt[, display_value := fifelse(
  statistic == "C", sprintf("%.0f", value),
  fifelse(statistic %in% c("IACC", "Creduction"), sprintf("%.1f", value),
  fifelse(statistic %in% c("SB0green", "SB0red", "PC270"),
    sprintf("%.1f%%", 100 * value),
    sprintf("%.2f", value))
))]
quilt[, relative_preference := 1 - relative_badness]

fwrite(quilt[, .(
  mp, statistic, metric, value, display_value, direction,
  relative_badness, relative_preference, metric_min, metric_max, years, summary
)], file.path(summary_dir, "candidate_quilt_reference_summary.csv"))

# Editable input layer for the report's scorecard method. A score of 100 is
# relatively best within the displayed candidate set and 0 relatively worst.
# Weights are defaults only and are intentionally not aggregated here.
scorecard_input <- quilt[, .(
  mp,
  statistic,
  metric,
  raw_value = value,
  preferred_direction = direction,
  normalized_score = 100 * relative_preference,
  # Include Kobe-red risk and mean catch reduction in addition to the six
  # established metrics. Low-catch risk remains visible and editable but is
  # not silently added to the worked composite without an agreed weight.
  include = statistic %in% c(
    "SBMSY", "FMSY", "C", "IACC", "VB2025", "VBMSY", "SB0red",
    "Creduction"
  ),
  weight = 1,
  years,
  summary
)]
fwrite(scorecard_input,
  file.path(summary_dir, "candidate_scorecard_input_reference.csv"))
jsonlite::write_json(
  scorecard_input,
  file.path(summary_dir, "candidate_scorecard_input_reference.json"),
  dataframe = "rows", pretty = TRUE, auto_unbox = TRUE, na = "null"
)

quilt_plot <- ggplot(quilt, aes(metric, mp, fill = relative_preference)) +
  geom_tile(colour = "grey80", linewidth = 0.35) +
  geom_text(aes(label = display_value), size = 4.2) +
  scale_fill_gradient(
    low = "#7952A8", high = "#F4EFF8", limits = c(0, 1),
    name = "Relative\npreference"
  ) +
  labs(
    title = "Long-term candidate performance metrics: reference OM",
    subtitle = paste(
      "Labels are long-term summaries; probabilities pool 2041-2050",
      "iteration-years; color is normalized within each metric"
    ),
    x = NULL, y = NULL,
    caption = paste(
      "Purple = lower and light lavender = higher relative preference.",
      "This is not an absolute acceptability score."
    )
  ) +
  ggthemes::theme_few(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    plot.title = element_text(size = 17),
    plot.subtitle = element_text(size = 12),
    plot.caption = element_text(size = 10),
    panel.grid = element_blank(),
    legend.position = "right"
  )

ggsave(file.path(output_dir, "candidate_quilt_reference.png"),
  quilt_plot, width = 15, height = 7.2, dpi = 160)
