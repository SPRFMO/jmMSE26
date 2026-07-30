#!/usr/bin/env Rscript

# Restore the eight-CMP evidence layer from the completed 100-iteration
# reference and nine-OM robustness runs retained in the adjacent jmMSE
# repository. The base comparison and the symmetric/30%-decrease variants
# were saved separately, so this script augments each set with vulnerable-
# biomass metrics and then combines only the eight report CMPs.

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
jmMSE_dir <- if (length(args)) args[[1L]] else "../jmMSE"
jmMSE_dir <- normalizePath(jmMSE_dir, mustWork = TRUE)

source_dirs <- list(
  reference_base = file.path(jmMSE_dir, "model", "candidates",
    "reference"),
  reference_variants = file.path(jmMSE_dir, "model", "candidates",
    "additional-change-limits", "reference"),
  robustness_base = file.path(jmMSE_dir, "model", "candidates",
    "robustness"),
  robustness_variants = file.path(jmMSE_dir, "model", "candidates",
    "additional-change-limits", "robustness")
)

missing_dirs <- names(source_dirs)[!vapply(
  source_dirs,
  function(path) all(file.exists(file.path(path,
    c("runs.rds", "performance.rds")))),
  logical(1)
)]
if (length(missing_dirs)) {
  stop("Missing saved 100-iteration inputs: ",
    paste(missing_dirs, collapse = ", "))
}

work_dir <- tempfile("candidate-performance-100-")
dir.create(work_dir, recursive = TRUE)
augment_script <- file.path("R", "augment_candidate_vb_performance.R")

augmented <- setNames(vector("list", length(source_dirs)), names(source_dirs))
for (source_name in names(source_dirs)) {
  output_file <- file.path(work_dir, paste0(source_name, ".rds"))
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(augment_script, source_dirs[[source_name]], output_file)
  )
  if (status != 0L || !file.exists(output_file)) {
    stop("VB augmentation failed for ", source_name)
  }
  augmented[[source_name]] <- as.data.table(readRDS(output_file))
}

candidate_codes <- c(
  "tun29", "tun43", "tun45", "tun47",
  "tun32", "tun44", "tun46", "tun48"
)
base_codes <- c("tun29", "tun43", "tun32", "tun44")
variant_codes <- setdiff(candidate_codes, base_codes)

select_candidates <- function(x, codes) {
  normalized_mp <- sub("_.*$", "", x$mp)
  x[normalized_mp %in% codes]
}

reference <- rbindlist(list(
  select_candidates(augmented$reference_base, base_codes),
  select_candidates(augmented$reference_variants, variant_codes)
), use.names = TRUE, fill = TRUE)
robustness <- rbindlist(list(
  select_candidates(augmented$robustness_base, base_codes),
  select_candidates(augmented$robustness_variants, variant_codes)
), use.names = TRUE, fill = TRUE)

validate <- function(x, mode, expected_oms) {
  normalized_mp <- sub("_.*$", "", x$mp)
  if (!identical(sort(unique(x$iter)), 1:100)) {
    stop(mode, " iterations are not exactly 1:100")
  }
  if (!setequal(unique(normalized_mp), candidate_codes)) {
    stop(mode, " does not contain the expected eight CMPs")
  }
  if (uniqueN(x$om) != expected_oms) {
    stop(mode, " contains ", uniqueN(x$om), " OMs; expected ", expected_oms)
  }
  if (!setequal(unique(x$statistic),
      c("SBMSY", "FMSY", "C", "IACC", "VB2025", "VBMSY"))) {
    stop(mode, " does not contain the expected six performance metrics")
  }
  if (anyNA(x$data)) stop(mode, " contains missing performance values")
  invisible(TRUE)
}

validate(reference, "reference", 1L)
validate(robustness, "robustness", 9L)

output_files <- c(
  reference = file.path("output", "candidate-performance", "reference",
    "performance_with_vb.rds"),
  robustness = file.path("output", "candidate-performance", "robustness",
    "performance_with_vb.rds")
)
dir.create(dirname(output_files[["reference"]]), recursive = TRUE,
  showWarnings = FALSE)
dir.create(dirname(output_files[["robustness"]]), recursive = TRUE,
  showWarnings = FALSE)
saveRDS(reference, output_files[["reference"]], compress = "xz")
saveRDS(robustness, output_files[["robustness"]], compress = "xz")

message("Restored 100-iteration reference evidence: ", nrow(reference),
  " rows")
message("Restored 100-iteration robustness evidence: ", nrow(robustness),
  " rows")
