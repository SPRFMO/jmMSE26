#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(FLCore)
  library(FLBRP)
  library(FLFishery)
  library(mse)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop(paste(
    "Usage: Rscript R/augment_candidate_vb_performance.R",
    "<source-candidate-results-dir> <jmMSE26-output-file>"
  ))
}

results_dir <- normalizePath(args[[1]], mustWork = TRUE)
runs_file <- file.path(results_dir, "runs.rds")
performance_file <- file.path(results_dir, "performance.rds")
output_file <- args[[2]]
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

if (!file.exists(runs_file) || !file.exists(performance_file)) {
  stop("Expected runs.rds and performance.rds in ", results_dir)
}
if (file.exists(output_file)) {
  stop("Refusing to overwrite existing output: ", output_file)
}

bevholtss3_to_bevholt <- function(stock_recruit) {
  steepness <- as.numeric(params(stock_recruit)["s", ])
  r0 <- as.numeric(params(stock_recruit)["R0", ])
  sb0 <- as.numeric(params(stock_recruit)["v", ])
  out <- stock_recruit
  model(out) <- bevholt()$model
  params(out) <- FLPar(
    a = (4 * steepness * r0) / (5 * steepness - 1),
    b = sb0 * (1 - steepness) / (5 * steepness - 1)
  )
  out
}

vulnerable_biomass_at_msy <- function(projected_om, biol_code,
  start_year = 2025L) {
  stock_i <- match(biol_code, names(biols(projected_om)))
  if (is.na(stock_i)) stop("Biological component not found: ", biol_code)

  equilibrium_stock <- window(stock(projected_om)[[stock_i]],
    start = start_year)
  stock_recruit <- bevholtss3_to_bevholt(sr(biols(projected_om)[[stock_i]]))
  equilibrium <- brp(FLBRP(equilibrium_stock, stock_recruit))

  f_grid <- fbar(equilibrium)
  vb_grid <- vb(equilibrium)
  f_msy <- as.numeric(refpts(equilibrium)["msy", "harvest", ])
  n_iter <- dim(vb_grid)[6]

  vapply(seq_len(n_iter), function(iter_i) {
    approx(
      x = as.numeric(f_grid[, , , , , iter_i]),
      y = as.numeric(vb_grid[, , , , , iter_i]),
      xout = f_msy[[iter_i]],
      rule = 2,
      ties = mean
    )$y
  }, numeric(1))
}

vb_performance <- function(run, run_code, type_code) {
  projected_om <- om(run)
  stocks <- stock(projected_om)

  rbindlist(lapply(names(stocks), function(biol_code) {
    vb_annual <- vb(stocks[[biol_code]])
    years <- as.integer(dimnames(vb_annual)$year)
    baseline_i <- match(2025L, years)
    if (is.na(baseline_i)) stop("Year 2025 not found for ", run_code)

    vb_2025 <- as.numeric(vb_annual[, baseline_i, , , , ])
    vb_msy <- vulnerable_biomass_at_msy(projected_om, biol_code)
    values <- list(
      VB2025 = sweep(vb_annual, 6, vb_2025, "/"),
      VBMSY = sweep(vb_annual, 6, vb_msy, "/")
    )
    labels <- c(
      VB2025 = "VB/VB[2025]",
      VBMSY = "VB/VB[MSY]"
    )
    descriptions <- c(
      VB2025 = "Vulnerable biomass relative to vulnerable biomass in 2025.",
      VBMSY = paste(
        "Vulnerable biomass relative to equilibrium vulnerable biomass",
        "at the FLBRP MSY fishing mortality."
      )
    )

    rbindlist(lapply(names(values), function(stat_code) {
      dat <- as.data.table(as.data.frame(values[[stat_code]]))
      dat[, .(
        om = as.character(name(projected_om)),
        biol = biol_code,
        statistic = stat_code,
        name = labels[[stat_code]],
        desc = descriptions[[stat_code]],
        year = as.integer(year),
        iter = as.integer(iter),
        data = as.numeric(data),
        type = type_code,
        run = run_code,
        mp = run_code
      )]
    }))
  }))
}

runs <- readRDS(runs_file)
original <- as.data.table(readRDS(performance_file))
run_types <- original[, .(type = unique(type)[1]), by = run]
type_map <- setNames(run_types$type, run_types$run)
added <- rbindlist(Map(vb_performance, runs, names(runs),
  unname(type_map[names(runs)])), use.names = TRUE)
combined <- rbindlist(list(original, added), fill = TRUE, use.names = TRUE)

saveRDS(combined, output_file, compress = "xz")
roundtrip <- readRDS(output_file)
expected <- c(sort(unique(original$statistic)), "VB2025", "VBMSY")
if (!setequal(unique(roundtrip$statistic), expected)) {
  stop("Saved statistic codes do not match expected values")
}
if (roundtrip[statistic == "VB2025" & year == 2025L,
    max(abs(data - 1), na.rm = TRUE)] > 1e-10) {
  stop("VB2025 is not exactly one in the baseline year")
}

message("Wrote ", output_file)
message("Statistics: ", paste(sort(unique(roundtrip$statistic)), collapse = ", "))
