#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggthemes)
})

args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args)) args[[1]] else
  path.expand("~/Downloads/babs/babs.rds")
if (!file.exists(input)) stop("Banking-and-borrowing runs not found: ", input)

runs <- readRDS(input)
run_meta <- data.table(
  run = names(runs),
  scenario = rep(c("Reference OM", "Low recruitment", "Cyclic recruitment"),
    each = 2),
  cmp = rep(c("HS+20", "PR+20"), 3),
  healthy_rule = rep(c(3, 2), 3)
)

tracking <- rbindlist(lapply(seq_along(runs), function(i) {
  x <- as.data.table(runs[[i]]@tracking)
  wide <- dcast(
    x[metric %chin% c("hcr", "rule.hcr")],
    year + iter ~ metric,
    value.var = "data"
  )
  setorder(wide, iter, year)
  wide[, previous_hcr := shift(hcr), by = iter]
  wide[, run := names(runs)[[i]]]
  wide
}))
tracking <- merge(tracking, run_meta, by = "run")
tracking[, event := fcase(
  !is.na(previous_hcr) & rule.hcr >= healthy_rule &
    hcr < previous_hcr * 0.85, "Borrowing eligible",
  !is.na(previous_hcr) & rule.hcr >= healthy_rule &
    hcr > previous_hcr * 1.15, "Banking eligible",
  default = "Neither"
)]

event_counts <- tracking[event != "Neither", .(count = .N),
  by = .(scenario, cmp, year = as.integer(year), event)]
grid <- CJ(
  scenario = unique(run_meta$scenario),
  cmp = unique(run_meta$cmp),
  year = 2026:2049,
  event = c("Banking eligible", "Borrowing eligible"),
  unique = TRUE
)
event_counts <- merge(grid, event_counts,
  by = c("scenario", "cmp", "year", "event"), all.x = TRUE)
event_counts[is.na(count), count := 0L]
event_counts[, scenario := factor(scenario,
  levels = c("Reference OM", "Low recruitment", "Cyclic recruitment"))]
event_counts[, cmp := factor(cmp, levels = c("HS+20", "PR+20"))]

out_data <- file.path("doc", "data", "banking_borrowing_eligibility_counts.csv")
fwrite(event_counts, out_data)

plot <- ggplot(event_counts,
  aes(year, count, shape = event, colour = event)) +
  geom_point(size = 2.1, alpha = 0.85) +
  facet_grid(scenario ~ cmp) +
  scale_shape_manual(values = c("Banking eligible" = 16,
    "Borrowing eligible" = 17)) +
  scale_colour_manual(values = c("Banking eligible" = "#4C78A8",
    "Borrowing eligible" = "#7A5195")) +
  scale_x_continuous(breaks = seq(2026, 2049, 4)) +
  scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, 10)) +
  labs(
    title = "Years eligible for banking or borrowing",
    subtitle = paste0(
      "Counts among 100 iterations, reconstructed from >15% changes in ",
      "annual HCR advice\nand the CMP-specific healthy-rule condition"
    ),
    x = "Projection year",
    y = "Number of eligible iterations",
    shape = NULL,
    colour = NULL
  ) +
  theme_few(base_size = 11) +
  theme(legend.position = "bottom")

out_plot <- file.path("doc", "figures", "banking-borrowing-eligibility.png")
ggsave(out_plot, plot, width = 9, height = 7, dpi = 180)

# Summarize spawning biomass for the same six B&B runs. Each plotted value is
# an iteration-level mean over a reporting period, rather than an annual value.
perf <- mse::performance(runs)
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

# The transferred run archive currently contains no non-zero recorded
# transactions. Keep this check explicit so eligibility is not mislabeled as
# implemented banking or borrowing.
recorded <- rbindlist(lapply(seq_along(runs), function(i) {
  x <- as.data.table(runs[[i]]@tracking)
  x[metric %chin% c("banking.isys", "borrowing.isys"), .(
    run = names(runs)[[i]], metric, nonzero = sum(data != 0, na.rm = TRUE)
  ), by = metric]
}))
if (any(recorded$nonzero != 0L))
  warning("The archive contains non-zero recorded transactions; review caption")

message("Wrote ", paste(c(out_data, out_plot, out_sb_data, out_sb_plot),
  collapse = ", "))
