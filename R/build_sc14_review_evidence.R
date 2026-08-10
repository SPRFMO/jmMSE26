#!/usr/bin/env Rscript

# Build the machine-readable evidence tables requested by the independent
# SC14 review. Run from the jmMSE26 repository root.

suppressPackageStartupMessages({
  library(Slick)
})

slick_path <- "output/jm_candidates_500.slick"
output_path <- "doc/data/candidates/candidate_biological_decision_table.csv"

if (!file.exists(slick_path)) stop("Missing checked Slick file: ", slick_path)

slick <- readRDS(slick_path)
Slick::Check(slick)

design <- slick@OMs@Design
mp_codes <- slick@MPs@Code
mp_labels <- slick@MPs@Label
kobe <- slick@Kobe@Value
quilt <- slick@Quilt@Value

long_years <- 2041:2050
long_index <- match(long_years, slick@Kobe@Time)
if (anyNA(long_index)) stop("Kobe array does not contain 2041--2050")

kobe_codes <- slick@Kobe@Code
sb_index <- match("SB/dynamic SBMSY", kobe_codes)
f_index <- match("F/FMSY", kobe_codes)
if (anyNA(c(sb_index, f_index))) stop("Required Kobe measures are missing")

quilt_codes <- slick@Quilt@Code
catch_index <- match("Catch 2041-2050", quilt_codes)
iacc_index <- match("IACC", quilt_codes)
if (anyNA(c(catch_index, iacc_index))) stop("Required quilt measures are missing")

rows <- vector("list", nrow(design) * length(mp_codes))
row_id <- 0L

for (om_i in seq_len(nrow(design))) {
  for (mp_i in seq_along(mp_codes)) {
    row_id <- row_id + 1L
    sb <- kobe[, om_i, mp_i, sb_index, long_index, drop = FALSE]
    f <- kobe[, om_i, mp_i, f_index, long_index, drop = FALSE]
    dim(sb) <- c(dim(kobe)[1L], length(long_index))
    dim(f) <- c(dim(kobe)[1L], length(long_index))

    valid <- is.finite(sb) & is.finite(f)
    green <- valid & sb >= 1 & f <= 1
    below <- is.finite(sb) & sb < 1
    overfishing <- is.finite(f) & f > 1

    annual_below <- colMeans(below, na.rm = TRUE)
    any_below <- apply(below, 1L, any, na.rm = TRUE)

    rows[[row_id]] <- data.frame(
      set = design$Set[om_i],
      om = design$OM[om_i],
      stock = design$Stock[om_i],
      cmp_code = mp_codes[mp_i],
      cmp = mp_labels[mp_i],
      p_dynamic_kobe_green = mean(green[valid]),
      p_sb_below_dynamic_sbmsy = mean(below[is.finite(sb)]),
      max_annual_p_sb_below_dynamic_sbmsy = max(annual_below, na.rm = TRUE),
      p_any_year_sb_below_dynamic_sbmsy = mean(any_below),
      p_f_above_fmsy = mean(overfishing[is.finite(f)]),
      median_long_term_catch = quilt[om_i, mp_i, catch_index],
      median_long_term_iacc = quilt[om_i, mp_i, iacc_index],
      annual_use_ready = "No -- operational specification remains",
      stringsAsFactors = FALSE
    )
  }
}

result <- do.call(rbind, rows)
numeric_columns <- vapply(result, is.numeric, logical(1))
result[numeric_columns] <- lapply(result[numeric_columns], function(x) round(x, 6))
write.csv(result, output_path, row.names = FALSE, na = "")

message("Wrote ", nrow(result), " CMP-by-OM/stock rows to ", output_path)

# Index-error parameters and historical coverage are extracted from the exact
# analysis repository named by the release manifest. The projection method is
# the same oem_index_dev_pars()/rmvlnormAR1() path used by constructJJOM().
manifest <- yaml::read_yaml("release/sc14-mse-2026-rc1.yml")
analysis_path <- manifest$repositories$analysis$path
parameter_source <- file.path(
  analysis_path, "report/data/index_residuals/om11_index_dev_pars.csv")
correlation_source <- file.path(
  analysis_path, "report/data/index_residuals/om11_index_dev_corr.csv")

if (!file.exists(parameter_source) || !file.exists(correlation_source)) {
  stop("Missing index-error exports in the named analysis repository")
}

parameters <- read.csv(parameter_source, stringsAsFactors = FALSE)
correlations <- read.csv(correlation_source, stringsAsFactors = FALSE,
  check.names = FALSE)

old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(analysis_path)
source("config.R")
source("utilities.R")
om11 <- readRDS("data/om11_h1_0.16_065.rds")
index_errors <- oem_index_dev_pars(om11, iy = 2025)
setwd(old_wd)

coverage <- lapply(names(index_errors$idev), function(index_name) {
  values <- as.array(window(index_errors$idev[[index_name]], end = 2025))
  years <- as.integer(dimnames(values)$year)
  has_value <- apply(is.finite(values), 2L, any)
  data.frame(
    index = index_name,
    first_year = min(years[has_value]),
    last_year = max(years[has_value]),
    years_with_data = sum(has_value),
    missing_years = paste(years[!has_value], collapse = ";"),
    stringsAsFactors = FALSE
  )
})
coverage <- do.call(rbind, coverage)
parameters <- merge(coverage, parameters, by = "index", sort = FALSE)
parameters$estimation_end_year <- 2025L
parameters$sd_scale <- "log"
parameters$method <- "iteration-specific multivariate lognormal AR(1)"

write.csv(parameters, "doc/data/index-error-parameters.csv",
  row.names = FALSE, na = "")
write.csv(correlations, "doc/data/index-error-correlations.csv",
  row.names = FALSE, na = "")
message("Wrote index-error parameter and correlation tables")
