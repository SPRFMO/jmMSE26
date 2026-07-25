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

included$weighted_contribution <-
  included$normalized_score * included$weight / total_weight

cmp_order <- unique(included$mp)
score <- vapply(
  cmp_order,
  function(candidate) {
    rows <- included[included$mp == candidate, , drop = FALSE]
    if (nrow(rows) != nrow(metric_settings) ||
        !setequal(rows$metric, metric_settings$metric)) {
      stop("Every CMP must contain exactly one row for every included metric: ",
           candidate)
    }
    sum(rows$weighted_contribution)
  },
  numeric(1)
)

result <- data.frame(
  rank = rank(-score, ties.method = "min"),
  mp = cmp_order,
  weighted_score = unname(score),
  included_metrics = nrow(metric_settings),
  total_weight = total_weight,
  weight_scheme = if (length(unique(metric_settings$weight)) == 1L) {
    "Equal weights"
  } else {
    "User-specified weights"
  },
  years = paste(sort(unique(included$years)), collapse = "; "),
  stringsAsFactors = FALSE
)
result <- result[order(result$rank, result$mp), , drop = FALSE]
result$weighted_score <- round(result$weighted_score, 2)

contributions <- included[
  order(match(included$mp, result$mp), included$metric),
  c("mp", "metric", "raw_value", "preferred_direction", "normalized_score",
    "include", "weight", "weighted_contribution", "years", "summary")
]
contributions$normalized_score <- round(contributions$normalized_score, 4)
contributions$weighted_contribution <-
  round(contributions$weighted_contribution, 4)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
write.csv(result, output_path, row.names = FALSE, na = "")
write.csv(contributions, contribution_path, row.names = FALSE, na = "")

table_lines <- c(
  "| Relative order | CMP | Equal-weight score (0--100) |",
  "|---:|:---|---:|",
  sprintf(
    "| %d | %s | %.2f |",
    result$rank,
    result$mp,
    result$weighted_score
  )
)
writeLines(table_lines, table_path)

message("Wrote ", nrow(result), " CMP scores to ", output_path)
message("Wrote ", nrow(contributions), " metric contributions to ",
        contribution_path)
message("Wrote report table to ", table_path)
