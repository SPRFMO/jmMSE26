#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggthemes)
})

args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args)) args[[1]] else
  file.path("output", "candidate-performance-500", "reference", "babs.rds")
performance_input <- if (length(args) >= 2) args[[2]] else
  file.path(dirname(input), "performance_babs.rds")
if (!file.exists(input)) stop("Banking-and-borrowing runs not found: ", input)
if (!file.exists(performance_input))
  stop("Banking-and-borrowing performance not found: ", performance_input)

runs <- readRDS(input)
run_meta <- data.table(
  run = names(runs),
  scenario = rep(c("Reference OM", "Low recruitment", "Cyclic recruitment"),
    each = 2),
  cmp = rep(c("HS+20", "PR+20"), 3)
)

tracking <- rbindlist(lapply(seq_along(runs), function(i) {
  x <- as.data.table(runs[[i]]@tracking)
  x[metric %chin% c("banking.isys", "borrowing.isys"), .(
    run = names(runs)[[i]],
    year = as.integer(year),
    iter = as.integer(iter),
    event = fifelse(metric == "banking.isys", "Banking applied",
      "Borrowing applied"),
    amount = data
  )]
}))
tracking <- merge(tracking, run_meta, by = "run")
n_iter <- max(tracking$iter, na.rm = TRUE)
event_counts <- tracking[amount > 0, .(count = .N),
  by = .(scenario, cmp, year, event)]
grid <- CJ(
  scenario = unique(run_meta$scenario),
  cmp = unique(run_meta$cmp),
  year = 2026:2049,
  event = c("Banking applied", "Borrowing applied"),
  unique = TRUE
)
event_counts <- merge(grid, event_counts,
  by = c("scenario", "cmp", "year", "event"), all.x = TRUE)
event_counts[is.na(count), count := 0L]
event_counts[, scenario := factor(scenario,
  levels = c("Reference OM", "Low recruitment", "Cyclic recruitment"))]
event_counts[, cmp := factor(cmp, levels = c("HS+20", "PR+20"))]
y_upper <- max(50, ceiling(max(event_counts$count) / 50) * 50)

out_data <- file.path("doc", "data", "banking_borrowing_eligibility_counts.csv")
fwrite(event_counts, out_data)

plot <- ggplot(event_counts,
  aes(year, count, shape = event, colour = event)) +
  geom_point(size = 2.1, alpha = 0.85) +
  facet_grid(scenario ~ cmp) +
  scale_shape_manual(values = c("Banking applied" = 16,
    "Borrowing applied" = 17)) +
  scale_colour_manual(values = c("Banking applied" = "#4C78A8",
    "Borrowing applied" = "#7A5195")) +
  scale_x_continuous(breaks = seq(2026, 2049, 4)) +
  scale_y_continuous(limits = c(0, y_upper),
    breaks = scales::breaks_pretty(5)) +
  labs(
    title = "Recorded banking and borrowing applications",
    subtitle = paste0("Counts among ", n_iter,
      " iterations in the updated run archive"),
    x = "Projection year",
    y = "Number of iterations",
    shape = NULL,
    colour = NULL
  ) +
  theme_few(base_size = 11) +
  theme(legend.position = "bottom")

out_plot <- file.path("doc", "figures", "banking-borrowing-eligibility.png")
ggsave(out_plot, plot, width = 9, height = 7, dpi = 180)

# Use the supplied performance table after confirming that it matches the six
# runs in the updated archive.
perf <- as.data.table(readRDS(performance_input))
if (!"run" %in% names(perf))
  stop("performance_babs.rds has no 'run' field")
if (!setequal(as.character(unique(perf$run)), names(runs)))
  stop("Run IDs differ between babs.rds and performance_babs.rds")

# Save the long-term dynamic Kobe-green results calculated from the same run
# archive used for the figures. This is the reproducible source for the report
# table.
sb0green <- as.data.table(perf)[
  statistic == "SB0green" & year %in% 2041:2050,
  .(SB0green = mean(data, na.rm = TRUE)), by = mp
]
sb0green[, `:=`(
  scenario = fcase(
    grepl("om11_2", as.character(mp)), "Low recruitment",
    grepl("om11_3", as.character(mp)), "Cyclic recruitment",
    default = "Reference OM"
  ),
  cmp = fifelse(grepl("tun29", as.character(mp)), "HS+20", "PR+20")
)]
sb0green[, `:=`(
  scenario_order = match(scenario,
    c("Reference OM", "Low recruitment", "Cyclic recruitment")),
  cmp_order = match(cmp, c("HS+20", "PR+20"))
)]
setorder(sb0green, scenario_order, cmp_order)
out_sb0green <- file.path("doc", "data",
  "banking_borrowing_sb0green_summary.csv")
fwrite(sb0green[, .(scenario, cmp, SB0green)], out_sb0green)

sb_period <- as.data.table(perf)[statistic == "SBMSY" & year %in% 2026:2050]
sb_period[, `:=`(
  scenario = fcase(
    grepl("om11_2", as.character(mp)), "Low recruitment",
    grepl("om11_3", as.character(mp)), "Cyclic recruitment",
    default = "Reference OM"
  ),
  cmp = fifelse(grepl("tun29", as.character(mp)), "HS+20", "PR+20"),
  period = fifelse(year %in% 2026:2035, "Near term (2026–2035)",
    fifelse(year %in% 2041:2050, "Long term (2041–2050)", NA_character_))
)]
sb_period <- sb_period[!is.na(period), .(
  mean_sb_sbmsy = mean(data, na.rm = TRUE)
), by = .(scenario, cmp, period, iter)]
sb_period[, scenario := factor(scenario,
  levels = c("Reference OM", "Low recruitment", "Cyclic recruitment"))]
sb_period[, cmp := factor(cmp, levels = c("HS+20", "PR+20"))]
sb_period[, period := factor(period,
  levels = c("Near term (2026–2035)", "Long term (2041–2050)"))]

out_sb_data <- file.path("doc", "data",
  "banking_borrowing_sbmsy_period_summary.csv")
fwrite(sb_period, out_sb_data)

sb_plot <- ggplot(sb_period, aes(scenario, mean_sb_sbmsy, fill = cmp)) +
  geom_hline(yintercept = 1, linewidth = 0.5, linetype = 2,
    colour = "grey40") +
  geom_boxplot(position = position_dodge(width = 0.75), width = 0.62,
    outlier.alpha = 0.35) +
  facet_wrap(~period, ncol = 2) +
  scale_fill_manual(values = c("HS+20" = "#4C78A8", "PR+20" = "#7A5195")) +
  labs(
    x = NULL,
    y = "Mean spawning biomass / SBMSY",
    fill = "CMP"
  ) +
  theme_few(base_size = 11) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 20, hjust = 1)
  )

out_sb_plot <- file.path("doc", "figures",
  "banking-borrowing-sbmsy-periods.png")
ggsave(out_sb_plot, sb_plot, width = 9, height = 4.8, dpi = 180)

# Confirm that the archive contains recorded applications.
recorded <- rbindlist(lapply(seq_along(runs), function(i) {
  x <- as.data.table(runs[[i]]@tracking)
  x[metric %chin% c("banking.isys", "borrowing.isys"), .(
    run = names(runs)[[i]], metric, nonzero = sum(data != 0, na.rm = TRUE)
  ), by = metric]
}))
if (!all(recorded[, any(nonzero > 0), by = metric]$V1))
  warning("At least one transaction type has no non-zero records")

message("Wrote ", paste(c(out_data, out_plot, out_sb0green, out_sb_data,
  out_sb_plot),
  collapse = ", "))
