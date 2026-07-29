#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Slick)
})

candidate_codes <- c(
  "tun29", "tun43", "tun45", "tun47",
  "tun32", "tun44", "tun46", "tun48"
)
statistics <- c("SBMSY", "FMSY", "C", "IACC", "VB2025", "VBMSY")

reference <- as.data.table(readRDS(
  "output/candidate-performance/reference/performance_with_vb.rds"
))
robustness <- as.data.table(readRDS(
  "output/candidate-performance/robustness/performance_with_vb.rds"
))
slick <- readRDS("output/jm_candidates.slick")
Check(slick)

source_codes <- unique(c(
  sub("_.*$", "", reference$run),
  sub("_om.*$", "", robustness$run)
))
if (!setequal(source_codes, candidate_codes))
  stop("Performance inputs and SC14 candidate list differ")
if (!identical(Code(MPs(slick)), candidate_codes))
  stop("Slick MP order and SC14 candidate list differ")
if (!setequal(Metadata(Timeseries(slick))$Code, statistics))
  stop("Slick and report performance-measure sets differ")

ts_dim <- dim(Value(Timeseries(slick)))
if (!identical(ts_dim, c(500L, 14L, 8L, 6L, 81L)))
  stop("Unexpected Slick dimensions: ", paste(ts_dim, collapse = " x "))
if (!identical(sort(unique(reference$iter)), 1:500) ||
    !identical(sort(unique(robustness$iter)), 1:500))
  stop("Performance inputs do not both contain iterations 1--500")
if (anyNA(reference$data) || anyNA(robustness$data))
  stop("Performance inputs contain missing numerical values")

registry <- fread("doc/data/cmp-registry.csv")
if (!setequal(registry[cmp_id %in% sub("^tun", "MP", candidate_codes), cmp_id],
    sub("^tun", "MP", candidate_codes)))
  stop("Registry and SC14 candidate list differ")

cat("SC14/Slick parity validation passed\n")
cat("Candidates:", paste(candidate_codes, collapse = ", "), "\n")
cat("Slick dimensions:", paste(ts_dim, collapse = " x "), "\n")
cat("Measures:", paste(statistics, collapse = ", "), "\n")
