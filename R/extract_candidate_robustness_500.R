#!/usr/bin/env Rscript

# Extract compact performance data from the checkpoint-only 500-draw
# robustness run.  Each checkpoint is processed independently so that the
# complete set of large FLmse objects is never held in memory or duplicated.

suppressPackageStartupMessages({
  library(data.table)
  library(FLCore)
  library(FLBRP)
  library(FLFishery)
  library(mse)
})

args <- commandArgs(trailingOnly = TRUE)
checkpoint_dir <- if (length(args) >= 1L) args[[1]] else
  "../jmMSE-500-refine/model/candidates/robustness_500/checkpoints"
output_file <- if (length(args) >= 2L) args[[2]] else
  "output/candidate-performance-500/robustness/performance_with_vb.rds"

checkpoint_dir <- normalizePath(checkpoint_dir, mustWork = TRUE)
checkpoint_files <- sort(list.files(checkpoint_dir, "\\.rds$", full.names = TRUE))
expected_oms <- c("om11_1", "om11_2", "om11_3", "om12", "om13",
  "om21", "om21_1", "om22", "om23")
if (!setequal(tools::file_path_sans_ext(basename(checkpoint_files)), expected_oms))
  stop("Robustness checkpoints do not match the expected nine OMs")
if (file.exists(output_file)) stop("Refusing to overwrite: ", output_file)
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

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
  equilibrium_stock <- window(stock(projected_om)[[stock_i]], start = start_year)
  stock_recruit <- bevholtss3_to_bevholt(sr(biols(projected_om)[[stock_i]]))
  equilibrium <- brp(FLBRP(equilibrium_stock, stock_recruit))
  f_grid <- fbar(equilibrium)
  vb_grid <- vb(equilibrium)
  f_msy <- as.numeric(refpts(equilibrium)["msy", "harvest", ])
  vapply(seq_len(dim(vb_grid)[6]), function(iter_i) {
    approx(
      x = as.numeric(f_grid[, , , , , iter_i]),
      y = as.numeric(vb_grid[, , , , , iter_i]),
      xout = f_msy[[iter_i]], rule = 2, ties = mean
    )$y
  }, numeric(1))
}

performance_one <- function(run, run_code, type_code = "rob") {
  projected_om <- om(run)
  rp <- refpts(projected_om)
  metric_groups <- metrics(projected_om)
  rbindlist(lapply(names(metric_groups), function(biol_code) {
    met <- metric_groups[[biol_code]]
    biol_rp <- if (is(rp, "FLPars")) rp[[biol_code]] else rp
    ref_quant <- function(code, q) {
      vals <- as.numeric(biol_rp[code, ])
      FLQuant(array(rep(vals, each = prod(dim(q)[1:5])), dim = dim(q),
        dimnames = dimnames(q)))
    }
    values <- list(
      SBMSY = met$SB / ref_quant("SBMSY", met$SB),
      FMSY = met$F / ref_quant("FMSY", met$F), C = met$C,
      IACC = 100 * abs(met$C[, -1] / met$C[, -dim(met$C)[2]] - 1)
    )
    labels <- c(SBMSY = "SB/SBMSY", FMSY = "F/FMSY", C = "Catch",
      IACC = "Interannual catch change")
    core <- rbindlist(lapply(names(values), function(stat_code) {
      dat <- as.data.table(as.data.frame(values[[stat_code]]))
      dat[, .(om = as.character(name(projected_om)), biol = biol_code,
        statistic = stat_code, name = labels[[stat_code]],
        desc = labels[[stat_code]], year = as.integer(year),
        iter = as.integer(iter), data = as.numeric(data), type = type_code,
        run = run_code, mp = run_code)]
    }))

    vb_annual <- vb(stock(projected_om)[[biol_code]])
    years <- as.integer(dimnames(vb_annual)$year)
    baseline_i <- match(2025L, years)
    if (is.na(baseline_i)) stop("Year 2025 not found for ", run_code)
    vb_values <- list(
      VB2025 = sweep(vb_annual, 6, as.numeric(vb_annual[, baseline_i, , , , ]), "/"),
      VBMSY = sweep(vb_annual, 6,
        vulnerable_biomass_at_msy(projected_om, biol_code), "/")
    )
    vb_labels <- c(VB2025 = "VB/VB[2025]", VBMSY = "VB/VB[MSY]")
    vb_desc <- c(
      VB2025 = "Vulnerable biomass relative to vulnerable biomass in 2025.",
      VBMSY = paste("Vulnerable biomass relative to equilibrium vulnerable",
        "biomass at the FLBRP MSY fishing mortality."))
    added <- rbindlist(lapply(names(vb_values), function(stat_code) {
      dat <- as.data.table(as.data.frame(vb_values[[stat_code]]))
      dat[, .(om = as.character(name(projected_om)), biol = biol_code,
        statistic = stat_code, name = vb_labels[[stat_code]],
        desc = vb_desc[[stat_code]], year = as.integer(year),
        iter = as.integer(iter), data = as.numeric(data), type = type_code,
        run = run_code, mp = run_code)]
    }))
    rbindlist(list(core, added), use.names = TRUE)
  }))
}

pieces <- vector("list", length(checkpoint_files))
for (i in seq_along(checkpoint_files)) {
  om_code <- tools::file_path_sans_ext(basename(checkpoint_files[[i]]))
  message(sprintf("[%d/%d] Extracting %s", i, length(checkpoint_files), om_code))
  runs <- readRDS(checkpoint_files[[i]])
  if (length(runs) != 8L) stop(om_code, " does not contain eight CMPs")
  pieces[[i]] <- rbindlist(Map(performance_one, runs, names(runs)),
    use.names = TRUE, fill = TRUE)
  rm(runs)
  gc(verbose = FALSE)
}

combined <- rbindlist(pieces, use.names = TRUE, fill = TRUE)
combined[, iter := as.integer(as.character(iter))]
if (!identical(sort(unique(combined$iter)), 1:500))
  stop("Extracted performance does not contain iterations 1--500")
if (uniqueN(combined$run) != 72L) stop("Expected 72 CMP--OM runs")
if (!setequal(unique(combined$statistic),
    c("SBMSY", "FMSY", "C", "IACC", "VB2025", "VBMSY")))
  stop("Unexpected performance metric set")
if (combined[statistic == "VB2025" & year == 2025L,
    max(abs(data - 1), na.rm = TRUE)] > 1e-10)
  stop("VB2025 is not exactly one in 2025")
if (anyNA(combined$data)) stop("Extracted performance contains missing values")

saveRDS(combined, output_file, compress = "xz")
message("Wrote ", output_file)
message("Rows: ", format(nrow(combined), big.mark = ","),
  "; OMs: ", uniqueN(combined$om), "; OM-stock series: ",
  uniqueN(combined[, paste(om, biol)]), "; runs: ", uniqueN(combined$run))
