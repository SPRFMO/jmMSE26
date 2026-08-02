#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggthemes)
})

input_file <- "output/candidate-performance-500/reference/performance_with_vb.rds"
figure_file <- "doc/figures/candidates/reference_catch_four_high_variability_iterations.png"
summary_file <- "doc/data/candidates/reference_catch_four_high_variability_iterations.csv"
trajectory_file <- "doc/data/candidates/reference_catch_four_high_variability_trajectories.csv"

candidate_labels <- c(
  tun29 = "HS+20 (MP29)", tun43 = "HS-20 (MP43)",
  tun45 = "HSsym (MP45)", tun47 = "HS-30 (MP47)",
  tun32 = "PR+20 (MP32)", tun44 = "PR-20 (MP44)",
  tun46 = "PRsym (MP46)", tun48 = "PR-30 (MP48)"
)
candidate_colors <- c(
  "HS+20 (MP29)" = "#C87A8A", "HS-20 (MP43)" = "#B6875B",
  "HSsym (MP45)" = "#909646", "HS-30 (MP47)" = "#55A067",
  "PR+20 (MP32)" = "#00A396", "PR-20 (MP44)" = "#409BBB",
  "PRsym (MP46)" = "#9189C7", "PR-30 (MP48)" = "#BE7AB4"
)

performance <- as.data.table(readRDS(input_file))
catch <- performance[
  statistic == "C" & mp %in% names(candidate_labels) & year %in% 2025:2050,
  .(catch = mean(data, na.rm = TRUE)),
  by = .(mp, iter, year)
]
setorder(catch, mp, iter, year)
catch[, annual_change := 100 * abs(catch / shift(catch) - 1),
  by = .(mp, iter)]

cmp_variability <- catch[year %in% 2026:2050,
  .(mean_iacc = mean(annual_change, na.rm = TRUE)), by = .(iter, mp)]
iteration_variability <- cmp_variability[, .(
  q25_cmp_iacc = as.numeric(quantile(mean_iacc, 0.25)),
  median_cmp_iacc = median(mean_iacc),
  minimum_cmp_iacc = min(mean_iacc),
  maximum_cmp_iacc = max(mean_iacc)
), by = iter]

# Require high variability across most CMPs by retaining the top decile of the
# 25th percentile across the eight CMP-specific IACC values. Within that set,
# use farthest-point sampling on standardized mean catch trajectories to select
# four contrasting temporal shapes reproducibly.
eligibility_cutoff <- as.numeric(quantile(
  iteration_variability$q25_cmp_iacc, 0.90))
eligible <- iteration_variability[q25_cmp_iacc >= eligibility_cutoff]

mean_shapes <- catch[iter %in% eligible$iter & year %in% 2026:2050,
  .(catch = mean(catch)), by = .(iter, year)]
mean_shapes[, standardized_catch := (catch - mean(catch)) / sd(catch), by = iter]
shape_wide <- dcast(mean_shapes, iter ~ year,
  value.var = "standardized_catch")
shape_ids <- shape_wide$iter
shape_matrix <- as.matrix(shape_wide[, -"iter"])

selected <- eligible[order(-q25_cmp_iacc, -median_cmp_iacc), iter][1L]
while (length(selected) < 4L) {
  remaining <- setdiff(shape_ids, selected)
  minimum_distance <- vapply(remaining, function(candidate_iter) {
    min(vapply(selected, function(selected_iter) {
      sqrt(sum((
        shape_matrix[match(candidate_iter, shape_ids), ] -
        shape_matrix[match(selected_iter, shape_ids), ]
      )^2))
    }, numeric(1)))
  }, numeric(1))
  selected <- c(selected, remaining[which.max(minimum_distance)])
}

selection <- iteration_variability[iter %in% selected]
selection[, selection_order := match(iter, selected)]
setorder(selection, selection_order)
selection[, eligibility_q25_iacc_cutoff := eligibility_cutoff]

plot_data <- catch[iter %in% selected & year %in% 2026:2050]
plot_data[, `:=`(
  cmp = factor(candidate_labels[mp], levels = unname(candidate_labels)),
  panel = factor(
    sprintf("Iteration %d", iter),
    levels = sprintf("Iteration %d", selected)
  )
)]

plot <- ggplot(plot_data, aes(year, catch, colour = cmp, group = cmp)) +
  geom_line(linewidth = 0.85, alpha = 0.92) +
  facet_wrap(~panel, ncol = 2, scales = "fixed") +
  scale_colour_manual(values = candidate_colors, drop = FALSE) +
  labs(
    title = "Contrasting high-variability catch projections: reference OM",
    subtitle = paste(
      "Four contrasting shapes selected from iterations in the top 10%",
      "of across-CMP catch variability"
    ),
    x = "Projection year",
    y = "Catch (model catch units)",
    colour = "CMP",
    caption = paste(
      "Eligibility uses the 25th percentile across eight CMP-specific mean",
      "IACC values; panels are then selected for contrasting temporal shapes."
    )
  ) +
  ggthemes::theme_few(base_size = 12) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.text = element_text(size = 9),
    strip.text = element_text(face = "bold")
  ) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))

dir.create(dirname(figure_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(summary_file), recursive = TRUE, showWarnings = FALSE)
ggsave(figure_file, plot, width = 11, height = 8.2, dpi = 180)
fwrite(selection, summary_file)
fwrite(plot_data[, .(selection_order = match(iter, selected), iter, mp,
  cmp = as.character(cmp), year, catch)], trajectory_file)

message("Selected iterations: ", paste(selected, collapse = ", "))
message("Wrote figure to ", figure_file)
message("Wrote selection summary to ", summary_file)
