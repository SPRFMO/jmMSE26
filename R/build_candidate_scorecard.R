#!/usr/bin/env Rscript

# Build the worked reference-OM scorecard used in Section 5.3 of the SC14 MSE
# report. The calculation deliberately starts from the editable input layer
# behind Figure 4 so that metric inclusion and weights remain explicit.

args <- commandArgs(trailingOnly = TRUE)
input_path <- if (length(args) >= 1L) {
  args[[1L]]
} else {
  "doc/data/candidates/candidate_scorecard_input_reference.csv"
}
output_path <- if (length(args) >= 2L) {
  args[[2L]]
} else {
  "doc/data/candidates/candidate_scorecard_result_reference.csv"
}
contribution_path <- if (length(args) >= 3L) {
  args[[3L]]
} else {
  "doc/data/candidates/candidate_scorecard_contributions_reference.csv"
}
table_path <- if (length(args) >= 4L) {
  args[[4L]]
} else {
  "doc/data/candidates/candidate_scorecard_result_reference.md"
}

scorecard <- read.csv(input_path, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("mp", "metric", "raw_value", "preferred_direction",
              "normalized_score", "include", "weight", "years", "summary")
missing_columns <- setdiff(required, names(scorecard))
if (length(missing_columns)) {
  stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
}

scorecard$include <- as.logical(scorecard$include)
scorecard$weight <- as.numeric(scorecard$weight)
input_normalized_score <- as.numeric(scorecard$normalized_score)
scorecard$raw_value <- as.numeric(scorecard$raw_value)

# Recalculate the direction-aware min--max score from the raw Figure 4 values.
# The existing normalized_score column is retained as an independent check.
scorecard$normalized_score <- NA_real_
for (metric_name in unique(scorecard$metric)) {
  metric_rows <- scorecard$metric == metric_name
  values <- scorecard$raw_value[metric_rows]
  directions <- unique(scorecard$preferred_direction[metric_rows])
  if (length(directions) != 1L) {
    stop("Metric has inconsistent preferred directions: ", metric_name)
  }
  metric_range <- max(values) - min(values)
  if (metric_range == 0) {
    scorecard$normalized_score[metric_rows] <- 100
  } else if (directions == "higher is better") {
    scorecard$normalized_score[metric_rows] <-
      100 * (values - min(values)) / metric_range
  } else if (directions == "lower is better") {
    scorecard$normalized_score[metric_rows] <-
      100 * (max(values) - values) / metric_range
  } else {
    stop("Unknown preferred direction for ", metric_name, ": ", directions)
  }
}
if (any(abs(scorecard$normalized_score - input_normalized_score) > 1e-8)) {
  stop("Recalculated normalized scores do not match the Figure 4 input layer.")
}

included <- scorecard[scorecard$include & scorecard$weight > 0, , drop = FALSE]
if (!nrow(included)) {
  stop("No metrics have include = TRUE and weight > 0.")
}
if (any(!is.finite(included$normalized_score)) ||
    any(!is.finite(included$weight))) {
  stop("Included normalized scores and weights must be finite.")
}

metric_settings <- unique(included[c("metric", "weight")])
if (anyDuplicated(metric_settings$metric)) {
  stop("Each included metric must have one common weight across CMPs.")
}
total_weight <- sum(metric_settings$weight)
if (total_weight <= 0) {
  stop("The total included metric weight must be positive.")
}

cmp_order <- unique(included$mp)
expected_metrics <- c(
  "SB / SB[MSY]", "F / F[MSY]", "Catch", "IACC",
  "VB / VB[2025]", "VB / VB[MSY]", "P(Kobe red)",
  "Mean catch reduction"
)
if (!setequal(metric_settings$metric, expected_metrics)) {
  stop(
    "The balanced example requires exactly these six metrics: ",
    paste(expected_metrics, collapse = ", ")
  )
}

weight_schemes <- list(
  equal = setNames(rep(1 / 8, 8), expected_metrics),
  balanced = setNames(
    c(0.10, 0.10, 1 / 6, 1 / 6, 0.10, 0.10, 0.10, 1 / 6),
    expected_metrics
  )
)

# A candidate-set-dependent sensitivity weighting based on how strongly each
# metric distinguishes the included CMP point estimates. The square root
# moderates the otherwise overwhelming influence of IACC under direct-CV
# weighting. These weights describe observed contrast, not management value.
metric_cv <- vapply(
  expected_metrics,
  function(metric_name) {
    values <- included$raw_value[included$metric == metric_name]
    stats::sd(values) / abs(mean(values))
  },
  numeric(1)
)
if (any(!is.finite(metric_cv)) || sum(metric_cv) <= 0) {
  stop("Dispersion weights require finite, positive metric CVs.")
}
weight_schemes$dispersion <- sqrt(metric_cv) / sum(sqrt(metric_cv))

calculate_scheme <- function(scheme_name, scheme_weights) {
  rows <- included
  rows$weight_scheme <- scheme_name
  rows$scheme_weight <- unname(scheme_weights[rows$metric])
  rows$weighted_contribution <-
    rows$normalized_score * rows$scheme_weight
  scores <- vapply(
    cmp_order,
    function(candidate) {
      candidate_rows <- rows[rows$mp == candidate, , drop = FALSE]
      if (nrow(candidate_rows) != length(scheme_weights) ||
          !setequal(candidate_rows$metric, names(scheme_weights))) {
        stop(
          "Every CMP must contain exactly one row for every included metric: ",
          candidate
        )
      }
      sum(candidate_rows$weighted_contribution)
    },
    numeric(1)
  )
  list(rows = rows, scores = scores)
}

equal <- calculate_scheme("Equal weights", weight_schemes$equal)
balanced <- calculate_scheme(
  "Balanced: 50% fishing performance / 50% stock condition",
  weight_schemes$balanced
)
dispersion <- calculate_scheme(
  "Dispersion sensitivity: square-root-CV weights",
  weight_schemes$dispersion
)

result <- data.frame(
  mp = cmp_order,
  equal_weight_rank = rank(-equal$scores, ties.method = "min"),
  equal_weight_score = unname(equal$scores),
  balanced_rank = rank(-balanced$scores, ties.method = "min"),
  balanced_score = unname(balanced$scores),
  dispersion_rank = rank(-dispersion$scores, ties.method = "min"),
  dispersion_score = unname(dispersion$scores),
  included_metrics = length(expected_metrics),
  years = paste(sort(unique(included$years)), collapse = "; "),
  stringsAsFactors = FALSE
)
result <- result[order(result$equal_weight_rank, result$mp), , drop = FALSE]
result[c(
  "equal_weight_score", "balanced_score", "dispersion_score"
)] <-
  round(result[c(
    "equal_weight_score", "balanced_score", "dispersion_score"
  )], 2)

contributions <- rbind(
  equal$rows, balanced$rows, dispersion$rows
)
contributions <- contributions[
  order(
    contributions$weight_scheme,
    match(contributions$mp, result$mp),
    match(contributions$metric, expected_metrics)
  ),
  c("weight_scheme", "mp", "metric", "raw_value",
    "preferred_direction", "normalized_score", "include", "scheme_weight",
    "weighted_contribution", "years", "summary")
]
contributions$normalized_score <- round(contributions$normalized_score, 4)
contributions$weighted_contribution <-
  round(contributions$weighted_contribution, 4)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write.csv(result, output_path, row.names = FALSE, na = "")
write.csv(contributions, contribution_path, row.names = FALSE, na = "")

table_lines <- c(
  paste(
    "| CMP | Equal weight | Balanced | Dispersion weighted |"
  ),
  "|:---|---:|---:|---:|",
  sprintf(
    "| %s | %.2f (%d) | %.2f (%d) | %.2f (%d) |",
    result$mp,
    result$equal_weight_score,
    result$equal_weight_rank,
    result$balanced_score,
    result$balanced_rank,
    result$dispersion_score,
    result$dispersion_rank
  ),
  "",
  paste(
    ": Relative scorecard sensitivity results for the eight",
    "performance-metric-quilt CMPs under the reference OM, 2041--2050.",
    "The balanced score assigns one-half across Catch, IACC, and mean catch",
    "reduction, and one-half across SB/SBMSY, F/FMSY, VB/VB[2025],",
    "VB/VB[MSY], and P(Kobe red). The",
    "dispersion-weighted sensitivity derives weights from the square root of",
    "each metric's across-CMP coefficient of variation. Each cell reports",
    "score (rank).",
    "{#tbl-scorecard-results}"
  )
)
writeLines(table_lines, table_path)

message("Wrote ", nrow(result), " CMP scores to ", output_path)
message("Wrote ", nrow(contributions), " metric contributions to ",
        contribution_path)
message("Wrote report table to ", table_path)
