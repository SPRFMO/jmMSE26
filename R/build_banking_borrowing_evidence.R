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

message("Wrote ", out_data, " and ", out_plot)
