#!/usr/bin/env Rscript

# Check agreement among the files in the named SC14 MSE release candidate.

suppressPackageStartupMessages({
  library(yaml)
  library(jsonlite)
  library(data.table)
  library(Slick)
})

args <- commandArgs(trailingOnly = TRUE)
manifest_path <- if (length(args)) args[[1L]] else
  file.path("release", "sc14-mse-2026-rc1.yml")
manifest <- yaml::read_yaml(manifest_path)
tolerance <- as.numeric(manifest$expected$numeric_tolerance)
checks <- list()

add_check <- function(area, check, outcome, expected = "", observed = "",
                      source = "", note = "") {
  checks[[length(checks) + 1L]] <<- data.frame(
    area = area,
    check = check,
    outcome = outcome,
    expected = paste(expected, collapse = "; "),
    observed = paste(observed, collapse = "; "),
    source = source,
    note = note,
    stringsAsFactors = FALSE
  )
}

same_set <- function(x, y) setequal(as.character(x), as.character(y))
fmt_set <- function(x) paste(sort(unique(as.character(x))), collapse = ", ")

sha256 <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  out <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
  sub("[[:space:]].*$", "", out[[1L]])
}

git_value <- function(path, args) {
  out <- suppressWarnings(system2("git", c("-C", shQuote(path), args),
    stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0L) return(NA_character_)
  paste(out, collapse = "\n")
}

for (repo_name in names(manifest$repositories)) {
  repo <- manifest$repositories[[repo_name]]
  current <- git_value(repo$path, "rev-parse HEAD")
  expected_revision <- repo$revision
  resolved_revision <- git_value(repo$path,
    paste("rev-parse", paste0(expected_revision, "^{commit}")))
  revision_ok <- !is.na(current) && !is.na(resolved_revision) &&
    identical(current, resolved_revision)
  observed_revision <- if (revision_ok &&
      !grepl("^[0-9a-f]{40}$", expected_revision)) {
    paste(expected_revision, "resolves to current commit")
  } else {
    current
  }
  add_check("release identity", paste(repo$name, "revision"),
    if (revision_ok) "Pass" else "Not ready",
    expected_revision, observed_revision, repo$path,
    if (is.na(resolved_revision)) {
      paste("Create the", expected_revision, "tag on the approved release commit.")
    } else if (!revision_ok) {
      paste("The tag resolves to", resolved_revision, "rather than the current commit.")
    } else "")
  dirty <- git_value(repo$path, "status --porcelain --untracked-files=no")
  clean <- !is.na(dirty) && !nzchar(dirty)
  add_check("release identity", paste(repo$name, "tracked files"),
    if (clean) "Pass" else "Not ready", "no uncommitted tracked changes",
    if (clean) "none" else dirty, repo$path)
}

candidate_ids <- unlist(manifest$expected$candidates, use.names = FALSE)
candidate_codes <- unlist(manifest$expected$candidate_codes, use.names = FALSE)
registry_path <- manifest$authoritative_files$cmp_registry
registry <- read.csv(registry_path, check.names = FALSE,
  stringsAsFactors = FALSE)
release_registry <- registry[tolower(registry$sc14_release) == "yes", , drop = FALSE]
add_check("CMP definitions", "SC14 CMP membership",
  if (same_set(release_registry$cmp_id, candidate_ids)) "Pass" else "Fail",
  fmt_set(candidate_ids), fmt_set(release_registry$cmp_id), registry_path)

controls_path <- manifest$authoritative_files$simulation_controls
controls <- read.csv(controls_path, check.names = FALSE,
  stringsAsFactors = FALSE)
controls <- controls[controls$mp %in% candidate_codes, , drop = FALSE]
add_check("CMP definitions", "simulation-control membership",
  if (same_set(controls$mp, candidate_codes)) "Pass" else "Fail",
  fmt_set(candidate_codes), fmt_set(controls$mp), controls_path)

control_fields <- c("target", "min", "lim", "trigger", "dlow", "dupp")
definition_errors <- character()
for (code in candidate_codes) {
  cmp_id <- sub("^tun", "MP", code)
  a <- release_registry[release_registry$cmp_id == cmp_id, , drop = FALSE]
  b <- controls[controls$mp == code, , drop = FALSE]
  if (nrow(a) != 1L || nrow(b) != 1L) next
  if (sub("\\.hcr$", "", a$hcr) != b$hcr) {
    definition_errors <- c(definition_errors,
      paste(cmp_id, "HCR", a$hcr, "!=", b$hcr))
  }
  for (field in control_fields) {
    av <- as.numeric(a[[field]])
    bv <- as.numeric(b[[field]])
    if (!isTRUE(all.equal(av, bv, tolerance = tolerance))) {
      definition_errors <- c(definition_errors,
        paste(cmp_id, field, format(av, digits = 12), "!=",
          format(bv, digits = 12)))
    }
  }
}
add_check("CMP definitions", "registry agrees with simulation controls",
  if (!length(definition_errors)) "Pass" else "Fail",
  "all eight HCRs and numerical settings agree",
  if (length(definition_errors)) definition_errors else "all agree",
  paste(registry_path, controls_path, sep = "; "))

tuning_path <- manifest$authoritative_files$tuning_summary
tuning <- read.csv(tuning_path, stringsAsFactors = FALSE)
tuning <- tuning[tuning$mp %in% candidate_codes, , drop = FALSE]
tuning_errors <- character()
for (code in candidate_codes) {
  cmp_id <- sub("^tun", "MP", code)
  a <- release_registry[release_registry$cmp_id == cmp_id, , drop = FALSE]
  b <- tuning[tuning$mp == code, , drop = FALSE]
  if (nrow(a) != 1L || nrow(b) != 1L ||
      abs(as.numeric(a$trigger) - as.numeric(b$trigger)) > tolerance ||
      a$label != b$short_label) {
    tuning_errors <- c(tuning_errors, cmp_id)
  }
}
add_check("tuning", "registry agrees with 500-draw tuning summary",
  if (!length(tuning_errors) && nrow(tuning) == length(candidate_codes)) "Pass" else "Fail",
  "eight matching labels and triggers",
  if (length(tuning_errors)) tuning_errors else "all agree", tuning_path)

read_performance <- function(path) {
  x <- readRDS(path)
  if (isS4(x) && "performance" %in% methods::slotNames(x)) x <- x@performance
  x <- data.table::as.data.table(x)
  x$run <- as.character(x$run)
  x$mp <- as.character(x$mp)
  x$om <- as.character(x$om)
  x$biol <- as.character(x$biol)
  x$statistic <- as.character(x$statistic)
  x$year <- as.integer(as.character(x$year))
  x$iter <- as.integer(as.character(x$iter))
  x$data <- as.numeric(x$data)
  x$cmp_code <- sub("_.*$", "", x$run)
  x
}

reference_path <- manifest$authoritative_files$reference_results
robustness_path <- manifest$authoritative_files$robustness_results
reference <- read_performance(reference_path)
robustness <- read_performance(robustness_path)
data.table::setkey(reference, run, biol, statistic, year, iter)
data.table::setkey(robustness, run, biol, statistic, year, iter)

check_performance <- function(x, set_name, path, expected_oms) {
  expected_stats <- unlist(manifest$expected$statistics, use.names = FALSE)
  add_check("simulation results", paste(set_name, "CMP membership"),
    if (same_set(x$cmp_code, candidate_codes)) "Pass" else "Fail",
    fmt_set(candidate_codes), fmt_set(x$cmp_code), path)
  add_check("simulation results", paste(set_name, "draws"),
    if (identical(sort(unique(x$iter)),
      seq_len(as.integer(manifest$expected$draws)))) "Pass" else "Fail",
    manifest$expected$draws, length(unique(x$iter)), path)
  add_check("simulation results", paste(set_name, "years"),
    if (min(x$year) == manifest$expected$result_first_year &&
      max(x$year) == manifest$expected$result_last_year) "Pass" else "Fail",
    paste(manifest$expected$result_first_year,
      manifest$expected$result_last_year, sep = "-") ,
    paste(min(x$year), max(x$year), sep = "-"), path)
  add_check("simulation results", paste(set_name, "statistics"),
    if (same_set(x$statistic, expected_stats)) "Pass" else "Fail",
    fmt_set(expected_stats), fmt_set(x$statistic), path)
  if (!is.null(expected_oms)) {
    run_oms <- unique(sub(paste0("^(?:", paste(candidate_codes,
      collapse = "|"), ")_"), "", x$run))
    add_check("simulation results", paste(set_name, "operating models"),
      if (same_set(run_oms, expected_oms)) "Pass" else "Fail",
      fmt_set(expected_oms), fmt_set(run_oms), path)
  }
}

check_performance(reference, "reference", reference_path, NULL)
check_performance(robustness, "robustness", robustness_path,
  unlist(manifest$expected$robustness_om_codes, use.names = FALSE))

run_index_path <- file.path("release", "sc14-run-index.csv")
run_index <- read.csv(run_index_path, stringsAsFactors = FALSE)
add_check("simulation results", "published run list candidate membership",
  if (same_set(run_index$cmp_code, candidate_codes)) "Pass" else "Fail",
  fmt_set(candidate_codes), fmt_set(run_index$cmp_code), run_index_path)
add_check("simulation results", "published run list draw count",
  if (all(run_index$draws == manifest$expected$draws)) "Pass" else "Fail",
  manifest$expected$draws, fmt_set(run_index$draws), run_index_path)
expected_design_rows <- length(candidate_codes) * (1L + 5L + 4L * 2L)
add_check("simulation results", "published run and stock combinations",
  if (nrow(run_index) == expected_design_rows) "Pass" else "Fail",
  expected_design_rows, nrow(run_index), run_index_path,
  "Two-stock operating models have one row for each stock component.")

slick_path <- manifest$authoritative_files$slick_file
slick <- readRDS(slick_path)
slick_check <- try(capture.output(Slick::Check(slick)), silent = TRUE)
slick_valid <- !inherits(slick_check, "try-error")
add_check("Slick file", "Slick structural check",
  if (slick_valid) "Pass" else "Fail", "Slick::Check completes",
  if (slick_valid) "complete" else as.character(slick_check), slick_path)

mp_meta <- Slick::Metadata(Slick::MPs(slick))
expected_labels <- release_registry$label[match(sub("^tun", "MP",
  candidate_codes), release_registry$cmp_id)]
add_check("Slick file", "CMP codes and labels",
  if (identical(as.character(mp_meta$Code), candidate_codes) &&
    identical(as.character(mp_meta$Label), expected_labels)) "Pass" else "Fail",
  paste(candidate_codes, expected_labels, sep = "=", collapse = ", "),
  paste(mp_meta$Code, mp_meta$Label, sep = "=", collapse = ", "), slick_path)

design <- Slick::Design(Slick::OMs(slick))
add_check("Slick file", "operating-model and stock rows",
  if (nrow(design) == 14L && all(c("om23", "North", "Southern") %in%
    unlist(design, use.names = FALSE))) "Pass" else "Fail",
  "14 rows including om23 North and Southern",
  paste(nrow(design), "rows"), slick_path)

# Compare all future numerical values stored in the Slick time-series array
# with the named reference and robustness performance files.
ts <- Slick::Value(Slick::Timeseries(slick))
ts_codes <- Slick::Code(Slick::Timeseries(slick))
ts_years <- Slick::Time(Slick::Timeseries(slick))
compare_codes <- intersect(unlist(manifest$expected$statistics,
  use.names = FALSE), ts_codes)
max_difference <- 0
mismatch_count <- 0L
for (om_i in seq_len(nrow(design))) {
  for (mp_i in seq_along(candidate_codes)) {
    code <- candidate_codes[mp_i]
    if (design$Set[om_i] == "Reference") {
      dat <- reference[run == code & biol == design$Stock[om_i]]
    } else {
      run_name <- paste(code, design$OM[om_i], sep = "_")
      dat <- robustness[run == run_name & biol == design$Stock[om_i]]
    }
    for (stat in compare_codes) {
      for (yr in intersect(2026:manifest$expected$result_last_year,
        ts_years)) {
        z <- dat[statistic == stat & year == yr, .(iter, data)]
        data.table::setorder(z, iter)
        expected_values <- z$data
        if (stat == "FMSY") expected_values <- pmin(expected_values, 4)
        observed_values <- ts[, om_i, mp_i, match(stat, ts_codes),
          match(yr, ts_years)]
        both_na <- is.na(expected_values) & is.na(observed_values)
        differences <- abs(expected_values - observed_values)
        differences[both_na] <- 0
        bad <- is.na(differences) | differences > tolerance
        mismatch_count <- mismatch_count + sum(bad)
        if (any(is.finite(differences))) {
          max_difference <- max(max_difference, differences, na.rm = TRUE)
        }
      }
    }
  }
}
add_check("Slick file", "future numerical values agree with result files",
  if (mismatch_count == 0L) "Pass" else "Fail", "no differences",
  paste(mismatch_count, "differences; maximum absolute difference",
    format(max_difference, scientific = TRUE)),
  paste(reference_path, robustness_path, slick_path, sep = "; "))

kobe <- Slick::Value(Slick::Kobe(slick))
kobe_years <- Slick::Time(Slick::Kobe(slick))
tune_years <- manifest$expected$tuning_first_year:
  manifest$expected$tuning_last_year
tune_i <- which(kobe_years %in% tune_years)
recalculated <- vapply(seq_along(candidate_codes), function(i) {
  mean(kobe[, 1L, i, 1L, tune_i] > 1 &
    kobe[, 1L, i, 2L, tune_i] < 1, na.rm = TRUE)
}, numeric(1))
tuning_ordered <- tuning[match(candidate_codes, tuning$mp), ]
tuning_difference <- max(abs(recalculated - tuning_ordered$tuning_probability))
add_check("tuning", "reported tuning probabilities recalculate from Slick",
  if (tuning_difference <= tolerance) "Pass" else "Fail",
  paste(format(tuning_ordered$tuning_probability, digits = 6), collapse = ", "),
  paste(format(recalculated, digits = 6), collapse = ", "),
  paste(tuning_path, slick_path, sep = "; "),
  paste("maximum absolute difference", format(tuning_difference,
    scientific = TRUE)))

# Check the public summary tables that feed figures and the scorecard.
display_names <- paste0(expected_labels, " (",
  sub("^tun", "MP", candidate_codes), ")")
summary_files <- list(
  quilt = "doc/data/candidates/candidate_quilt_reference_summary.csv",
  scorecard_input = "doc/data/candidates/candidate_scorecard_input_reference.csv",
  scorecard_result = "doc/data/candidates/candidate_scorecard_result_reference.csv",
  kobe = "doc/data/candidates/kobe_reference_periods_summary.csv",
  vb_reference = "doc/data/candidates/vb_reference_summary.csv",
  vb_robustness = "doc/data/candidates/vb_robustness_summary.csv"
)
summary_membership <- list(
  quilt = display_names,
  scorecard_input = display_names,
  scorecard_result = display_names,
  kobe = expected_labels,
  vb_reference = expected_labels,
  vb_robustness = expected_labels
)
summary_errors <- character()
summary_data <- list()
for (nm in names(summary_files)) {
  x <- read.csv(summary_files[[nm]], stringsAsFactors = FALSE,
    check.names = FALSE)
  summary_data[[nm]] <- x
  column <- if (nm == "kobe") "cmp" else "mp"
  if (!same_set(x[[column]], summary_membership[[nm]])) {
    summary_errors <- c(summary_errors, nm)
  }
}
add_check("public summaries", "CMP membership in figure and scorecard tables",
  if (!length(summary_errors)) "Pass" else "Fail",
  "the same eight CMPs in all six tables",
  if (length(summary_errors)) summary_errors else "all agree",
  paste(unlist(summary_files), collapse = "; "))

quilt_values <- summary_data$quilt[c("mp", "statistic", "metric", "value")]
score_values <- summary_data$scorecard_input[
  c("mp", "statistic", "metric", "raw_value")]
joined <- merge(quilt_values, score_values,
  by = c("mp", "statistic", "metric"), all = TRUE)
score_difference <- max(abs(joined$value - joined$raw_value), na.rm = TRUE)
score_missing <- anyNA(joined[c("value", "raw_value")])
add_check("public summaries", "scorecard inputs agree with the figure summary",
  if (!score_missing && score_difference <= tolerance) "Pass" else "Fail",
  "identical numerical inputs", paste("maximum absolute difference",
    format(score_difference, scientific = TRUE)),
  paste(summary_files$quilt, summary_files$scorecard_input, sep = "; "))

kobe_summary <- summary_data$kobe
kobe_difference <- 0
for (i in seq_len(nrow(kobe_summary))) {
  mp_i <- match(kobe_summary$cmp[i], expected_labels)
  period_years <- if (grepl("2041", kobe_summary$period[i]))
    2041:2050 else 2026:2035
  year_i <- which(kobe_years %in% period_years)
  expected_f <- median(rowMeans(kobe[, 1L, mp_i, 2L, year_i], na.rm = TRUE))
  expected_sb <- median(rowMeans(kobe[, 1L, mp_i, 1L, year_i], na.rm = TRUE))
  kobe_difference <- max(kobe_difference,
    abs(expected_f - kobe_summary$f_fmsy[i]),
    abs(expected_sb - kobe_summary$sb_sbmsy[i]))
}
add_check("public summaries", "Kobe figure table agrees with the Slick file",
  if (kobe_difference <= tolerance) "Pass" else "Fail",
  "identical period summaries", paste("maximum absolute difference",
    format(kobe_difference, scientific = TRUE)),
  paste(summary_files$kobe, slick_path, sep = "; "))

page_expectations <- list(
  MP29 = list(path = "doc/cmp/mp29.qmd", trigger = "2.125", old = character()),
  MP32 = list(path = "doc/cmp/mp32.qmd", trigger = "1.640625", old = "1.421875"),
  MP43 = list(path = "doc/cmp/mp43.qmd", trigger = "1.9375", old = "3.0625"),
  MP44 = list(path = "doc/cmp/mp44.qmd", trigger = "1.6015625", old = "1.9140625")
)
page_errors <- character()
for (cmp_id in names(page_expectations)) {
  item <- page_expectations[[cmp_id]]
  text <- paste(readLines(item$path, warn = FALSE), collapse = "\n")
  if (!grepl(item$trigger, text, fixed = TRUE) ||
      any(vapply(item$old, grepl, logical(1), x = text, fixed = TRUE))) {
    page_errors <- c(page_errors, cmp_id)
  }
}
add_check("public pages", "individual CMP pages use the refined triggers",
  if (!length(page_errors)) "Pass" else "Fail",
  "current triggers and no earlier trigger values",
  if (length(page_errors)) page_errors else "all agree",
  paste(vapply(page_expectations, `[[`, character(1), "path"),
    collapse = "; "))

file_list_path <- file.path("release", "sc14-file-list.csv")
file_list <- read.csv(file_list_path, stringsAsFactors = FALSE)
current_hash <- vapply(file_list$path, sha256, character(1))
missing_files <- file_list$path[is.na(current_hash)]
changed_files <- file_list$path[!is.na(current_hash) &
  current_hash != file_list$sha256]
add_check("published files", "all listed files exist",
  if (!length(missing_files)) "Pass" else "Fail", "all present",
  if (length(missing_files)) missing_files else "all present", file_list_path)
add_check("published files", "listed files match the accepted baseline",
  if (!length(changed_files)) "Pass" else "Not ready", "all fingerprints agree",
  if (length(changed_files)) changed_files else "all agree", file_list_path,
  if (length(changed_files)) paste("Review the changes, then intentionally run",
    "R/prepare_sc14_release.R to accept a new baseline.") else "")

for (path in unlist(manifest$published_outputs, use.names = FALSE)) {
  add_check("published files", paste("published output", path),
    if (file.exists(path)) "Pass" else "Not ready", "file exists",
    if (file.exists(path)) "present" else "missing", path)
}

superseded <- unlist(manifest$superseded_current_products, use.names = FALSE)
source_text <- paste(vapply(c("README.md", list.files("doc",
  pattern = "\\.(qmd|md)$", recursive = TRUE, full.names = TRUE)),
  function(path) paste(readLines(path, warn = FALSE), collapse = "\n"),
  character(1)), collapse = "\n")
stale_mentions <- superseded[vapply(superseded,
  function(path) grepl(path, source_text, fixed = TRUE), logical(1))]
add_check("published files", "superseded filenames are not presented as current",
  if (!length(stale_mentions)) "Pass" else "Not ready", "none",
  if (length(stale_mentions)) stale_mentions else "none", manifest_path,
  "Historical references may remain, but they must be labelled as replaced.")

results <- do.call(rbind, checks)
dir.create("release", showWarnings = FALSE)
write.csv(results, file.path("release", "sc14-check-results.csv"),
  row.names = FALSE, na = "")
public_release_dir <- file.path("doc", "data", "release")
dir.create(public_release_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(results, file.path(public_release_dir, "sc14-check-results.csv"),
  row.names = FALSE, na = "")
json_result <- list(
  release_id = manifest$release_id,
  # Keep the published verification byte-stable. The release manifest date
  # identifies this accepted baseline; rerunning the same checks must not
  # alter the release merely because the clock or commit SHA changed.
  checked = as.character(manifest$prepared),
  summary = as.list(table(factor(results$outcome,
    levels = c("Pass", "Not ready", "Fail")))),
  checks = results
)
jsonlite::write_json(json_result, file.path("release", "sc14-check-results.json"),
  pretty = TRUE, auto_unbox = TRUE, na = "null")
jsonlite::write_json(json_result,
  file.path(public_release_dir, "sc14-check-results.json"),
  pretty = TRUE, auto_unbox = TRUE, na = "null")

summary <- table(factor(results$outcome,
  levels = c("Pass", "Not ready", "Fail")))
cat(manifest$release_id, "\n")
cat("  Pass: ", summary[["Pass"]], "\n", sep = "")
cat("  Not ready: ", summary[["Not ready"]], "\n", sep = "")
cat("  Fail: ", summary[["Fail"]], "\n", sep = "")
cat("Wrote release/sc14-check-results.csv and .json\n")

if (summary[["Fail"]] > 0L) quit(status = 1L)
if (summary[["Not ready"]] > 0L) quit(status = 2L)
