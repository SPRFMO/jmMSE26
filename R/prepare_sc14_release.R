#!/usr/bin/env Rscript

# Prepare the file and run lists for the named SC14 MSE release candidate.
# Run this only when intentionally accepting a new release baseline. Routine
# checks are performed by R/check_sc14_release.R and do not update the baseline.

suppressPackageStartupMessages({
  library(yaml)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
manifest_path <- if (length(args)) args[[1L]] else
  file.path("release", "sc14-mse-2026-rc1.yml")
manifest <- yaml::read_yaml(manifest_path)

sha256 <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  unname(tools::md5sum(path)) # replaced below by the platform SHA-256 command
}

sha256 <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  out <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
  sub("[[:space:]].*$", "", out[[1L]])
}

read_performance <- function(path) {
  x <- readRDS(path)
  if (isS4(x) && "performance" %in% methods::slotNames(x)) x <- x@performance
  data.table::as.data.table(x)
}

candidate_codes <- unlist(manifest$expected$candidate_codes, use.names = FALSE)
candidate_ids <- unlist(manifest$expected$candidates, use.names = FALSE)
registry <- read.csv(manifest$authoritative_files$cmp_registry,
  check.names = FALSE, stringsAsFactors = FALSE)
labels <- setNames(registry$label[match(candidate_ids, registry$cmp_id)],
  sub("^MP", "tun", candidate_ids))

make_index <- function(path, set_name) {
  x <- read_performance(path)
  x$run <- as.character(x$run)
  x$om <- as.character(x$om)
  x$biol <- as.character(x$biol)
  x$statistic <- as.character(x$statistic)
  x$year <- as.integer(as.character(x$year))
  x$iter <- as.integer(as.character(x$iter))
  x$cmp_code <- sub("_.*$", "", x$run)
  keep <- x$cmp_code %in% candidate_codes
  x <- x[keep, , drop = FALSE]
  index <- x[, .(
    draws = data.table::uniqueN(iter),
    first_year = min(year),
    last_year = max(year),
    statistics = paste(sort(unique(statistic)), collapse = ";")
  ), by = .(run, cmp_code, source_om = om, stock = biol)]
  index[, `:=`(
    set = set_name,
    cmp_id = sub("^tun", "MP", cmp_code),
    label = unname(labels[cmp_code]),
    om_code = if (set_name == "reference") {
      manifest$expected$reference_om_code
    } else {
      substring(run, nchar(cmp_code) + 2L)
    },
    source_file = path
  )]
  as.data.frame(index[, .(set, run, cmp_code, cmp_id, label, om_code,
    source_om, stock, draws, first_year, last_year, statistics,
    source_file)])
}

run_index <- rbind(
  make_index(manifest$authoritative_files$reference_results, "reference"),
  make_index(manifest$authoritative_files$robustness_results, "robustness")
)
run_index <- run_index[order(match(run_index$cmp_code, candidate_codes),
  run_index$set, run_index$om_code, run_index$stock), , drop = FALSE]
dir.create("release", showWarnings = FALSE)
write.csv(run_index, file.path("release", "sc14-run-index.csv"),
  row.names = FALSE, na = "")
public_release_dir <- file.path("doc", "data", "release")
dir.create(public_release_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(run_index, file.path(public_release_dir, "sc14-run-index.csv"),
  row.names = FALSE, na = "")

entries <- manifest$release_files
file_list <- do.call(rbind, lapply(entries, function(entry) {
  path <- entry$path
  info <- if (file.exists(path)) file.info(path) else NULL
  data.frame(
    category = entry$category,
    path = path,
    bytes = if (is.null(info)) NA_real_ else unname(info$size),
    sha256 = sha256(path),
    stringsAsFactors = FALSE
  )
}))
write.csv(file_list, file.path("release", "sc14-file-list.csv"),
  row.names = FALSE, na = "")
write.csv(file_list, file.path(public_release_dir, "sc14-file-list.csv"),
  row.names = FALSE, na = "")

cat("Prepared ", manifest$release_id, "\n", sep = "")
cat("  run records: ", nrow(run_index), "\n", sep = "")
cat("  listed files: ", nrow(file_list), "\n", sep = "")
cat("Now review the two CSV files and run R/check_sc14_release.R.\n")
