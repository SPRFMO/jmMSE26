# Build a validated Slick object for the current jmMSE candidate set.

suppressPackageStartupMessages({
  library(data.table)
  library(Slick)
})

read_performance <- function(path) {
  if (!file.exists(path)) stop("Performance input not found: ", path)
  x <- readRDS(path)
  if (inherits(x, "data.frame")) return(as.data.table(x))
  if (isS4(x) && "performance" %in% methods::slotNames(x)) {
    out <- as.data.table(x@performance)
    if (!nrow(out)) stop("The performance slot is empty: ", path)
    return(out)
  }
  stop("Expected a performance table or an S4 object with a performance slot")
}

candidate_run_map <- function(perf, candidate_codes) {
  needed <- c("run", "om", "biol", "statistic", "year", "iter", "data")
  missing <- setdiff(needed, names(perf))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))

  runs <- unique(perf[, .(run, source_om = as.character(om),
    biol = as.character(biol))])
  pattern <- paste0("^(", paste(candidate_codes, collapse = "|"), ")(?:_(.+))?$")
  runs[, mp := sub(pattern, "\\1", run)]
  runs <- runs[mp %in% candidate_codes]
  if (!nrow(runs)) {
    stop("No candidate runs found. Expected: ",
      paste(candidate_codes, collapse = ", "))
  }
  runs[, suffix := sub(pattern, "\\2", run)]
  runs[suffix == run | suffix == "", suffix := NA_character_]
  runs[, om_code := fifelse(is.na(suffix), source_om, suffix)]
  runs[, OM := paste(om_code, biol, sep = " / ")]
  unique(runs[, .(run, mp, OM, om_code, biol, source_om)])
}

performance_matrix <- function(perf, statistic_code, run_code, biol_code,
                               years, iters) {
  x <- perf[run == run_code & biol == biol_code &
    statistic == statistic_code & year %in% years & iter %in% iters,
    .(data = mean(data, na.rm = TRUE)), by = .(iter, year)]
  wide <- dcast(x, iter ~ year, value.var = "data")
  setorder(wide, iter)
  wide <- merge(data.table(iter = iters), wide, by = "iter", all.x = TRUE,
    sort = FALSE)
  for (yr in as.character(years)) {
    if (!yr %in% names(wide)) wide[, (yr) := NA_real_]
  }
  as.matrix(wide[, as.character(years), with = FALSE])
}

finite_mean <- function(x) {
  x[!is.finite(x)] <- NA_real_
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

timeseries_matrix <- function(ts, om_i, mp_i, pi_i, year_i) {
  out <- Value(ts)[, om_i, mp_i, pi_i, year_i, drop = FALSE]
  matrix(as.numeric(out), nrow = dim(Value(ts))[1], ncol = length(year_i))
}

mean_catch_reduction <- function(catch_matrix) {
  if (ncol(catch_matrix) < 2) return(NA_real_)
  previous <- catch_matrix[, -ncol(catch_matrix), drop = FALSE]
  current <- catch_matrix[, -1, drop = FALSE]
  reduction <- ifelse(previous > 0,
    pmax((previous - current) / previous, 0) * 100,
    NA_real_)
  finite_mean(reduction)
}

make_om_metadata <- function(run_map, source_file) {
  design <- unique(run_map[, .(
    OM, Model = om_code, Stock = biol, SourceOM = source_om,
    Source = basename(source_file)
  )])
  setorder(design, Model, Stock)
  factors <- rbindlist(lapply(names(design), function(nm) data.table(
    Factor = nm,
    Level = unique(as.character(design[[nm]])),
    Description = unique(as.character(design[[nm]]))
  )))
  OMs(Factors = as.data.frame(factors), Design = as.data.frame(design))
}

make_mp_metadata <- function(candidate_codes, registry_file) {
  registry <- fread(registry_file)
  registry <- registry[cmp_id %in% sub("^tun", "MP", candidate_codes)]
  meta <- data.table(code = candidate_codes)
  meta[, cmp_id := sub("^tun", "MP", code)]
  meta <- merge(meta, registry, by = "cmp_id", all.x = TRUE, sort = FALSE)
  meta[, Label := fifelse(is.na(label), cmp_id, label)]
  meta[, Description := fifelse(is.na(hcr),
    paste("jmMSE candidate", cmp_id),
    sprintf("%s with %s; target %s, trigger %s, min %s, lim %s.",
      estimator, hcr, target, trigger, min, lim))]
  MPs(Code = meta$code, Label = meta$Label,
    Description = meta$Description)
}

build_candidate_slick <- function(
  performance_file,
  out_file = file.path("output", "jm_candidates.slick"),
  registry_file = file.path("doc", "data", "cmp-registry.csv"),
  candidate_codes = c("tun29", "tun32"),
  time_now = 2025L
) {
  perf <- read_performance(performance_file)
  setDT(perf)
  perf[, `:=`(
    run = as.character(run),
    om = as.character(om),
    biol = as.character(biol),
    statistic = as.character(statistic),
    year = as.integer(as.character(year)),
    iter = as.integer(as.character(iter)),
    data = as.numeric(data)
  )]

  run_map <- candidate_run_map(perf, candidate_codes)
  available <- intersect(c("SBMSY", "FMSY", "C", "IACC"),
    unique(perf$statistic))
  if (!all(c("SBMSY", "FMSY", "C") %in% available)) {
    stop("Slick export requires SBMSY, FMSY, and C performance statistics")
  }

  coverage <- perf[run %in% run_map$run & statistic %in% available,
    .N, by = .(run, biol, statistic, year)]
  expected <- nrow(run_map) * length(available)
  years <- sort(coverage[, .N, by = year][N == expected, year])
  years <- years[years >= time_now]
  if (!length(years)) stop("No common years across candidate runs")
  iters <- sort(unique(perf[run %in% run_map$run & statistic == "C", iter]))

  om_codes <- unique(run_map$OM)
  mp_codes <- candidate_codes[candidate_codes %in% run_map$mp]
  pi_meta <- data.frame(
    Code = available,
    Label = c(SBMSY = "SB/SBMSY", FMSY = "F/FMSY", C = "Catch",
      IACC = "IACC")[available],
    Description = c(
      SBMSY = "Spawning biomass relative to SBMSY.",
      FMSY = "Fishing mortality relative to FMSY.",
      C = "Catch.",
      IACC = "Interannual catch change."
    )[available],
    stringsAsFactors = FALSE
  )

  ts <- Timeseries()
  Metadata(ts) <- pi_meta
  Time(ts) <- years
  TimeNow(ts) <- time_now
  TimeLab(ts) <- "Year"
  Value(ts) <- array(NA_real_, dim = c(length(iters), length(om_codes),
    length(mp_codes), nrow(pi_meta), length(years)))

  for (i in seq_len(nrow(run_map))) {
    om_i <- match(run_map$OM[i], om_codes)
    mp_i <- match(run_map$mp[i], mp_codes)
    if (is.na(mp_i)) next
    for (pi_i in seq_along(available)) {
      Value(ts)[, om_i, mp_i, pi_i, ] <- performance_matrix(
        perf, available[pi_i], run_map$run[i], run_map$biol[i], years, iters)
    }
  }
  Target(ts) <- ifelse(available %in% c("SBMSY", "FMSY"), 1, NA_real_)
  Limit(ts) <- ifelse(available == "SBMSY", 0.5, NA_real_)
  Check(ts)

  kobe <- Kobe(
    Code = c("SB/SBMSY", "F/FMSY"),
    Label = c("SB/SBMSY", "F/FMSY"),
    Description = c("Spawning biomass relative to SBMSY.",
      "Fishing mortality relative to FMSY.")
  )
  projection <- years > time_now
  Time(kobe) <- years[projection]
  Value(kobe) <- Value(ts)[, , , match(c("SBMSY", "FMSY"), available),
    projection, drop = FALSE]
  Target(kobe) <- c(1, 1)
  Limit(kobe) <- c(0.5, NA_real_)
  Check(kobe)

  terminal_year <- max(years)
  terminal_i <- match(terminal_year, years)
  box_meta <- data.frame(
    Code = available,
    Label = pi_meta$Label,
    Description = c(
      SBMSY = paste("Terminal-year", terminal_year,
        "spawning biomass relative to SBMSY."),
      FMSY = paste("Terminal-year", terminal_year,
        "fishing mortality relative to FMSY."),
      C = paste("Terminal-year", terminal_year, "catch."),
      IACC = paste("Terminal-year", terminal_year,
        "interannual percentage change in catch.")
    )[available],
    stringsAsFactors = FALSE
  )
  boxplot <- Boxplot()
  Metadata(boxplot) <- box_meta
  box_value <- array(
    as.numeric(Value(ts)[, , , , terminal_i, drop = FALSE]),
    dim = c(length(iters), length(om_codes), length(mp_codes),
      length(available)))
  box_value[!is.finite(box_value)] <- NA_real_
  Value(boxplot) <- box_value
  Check(boxplot)

  short_years <- years[years %in% 2026:2030]
  long_years <- years[years %in% 2036:2040]
  projection_years <- years[years > time_now]
  if (!length(short_years) || !length(long_years)) {
    stop("Boxplot, quilt, and trade-off summaries require 2026-2030 and 2036-2040")
  }
  short_i <- match(short_years, years)
  long_i <- match(long_years, years)
  projection_i <- match(projection_years, years)
  sb_i <- match("SBMSY", available)
  f_i <- match("FMSY", available)
  catch_i <- match("C", available)
  iacc_i <- match("IACC", available)

  summary_codes <- c("P(green)", "P(yellow)", "P(orange)", "P(red)",
    "P(<270)", "Catch 2026-2030", "Catch 2036-2040",
    "Catch reduction", "SB/SBMSY", "F/FMSY", "IACC")
  summary_labels <- c("P(green)", "P(yellow)", "P(orange)", "P(red)",
    "P(catch below 270)", "Short-term catch", "Mean catch",
    "Mean catch reduction", "Mean SB/SBMSY", "Mean F/FMSY", "Mean IACC")
  summary_descriptions <- c(
    "Probability of Kobe green status over 2036-2040.",
    "Probability of Kobe yellow status over 2036-2040.",
    "Probability of Kobe orange status over 2036-2040.",
    "Probability of Kobe red status over 2036-2040.",
    "Probability that catch is below 270 over 2036-2040.",
    "Mean catch over 2026-2030.",
    "Mean catch over 2036-2040.",
    "Mean annual percentage catch reduction over 2026-2050.",
    "Mean spawning biomass relative to SBMSY over 2036-2040.",
    "Mean fishing mortality relative to FMSY over 2036-2040.",
    "Mean interannual percentage change in catch over 2036-2040."
  )
  summaries <- array(NA_real_,
    dim = c(length(om_codes), length(mp_codes), length(summary_codes)))

  for (om_i in seq_along(om_codes)) {
    for (mp_i in seq_along(mp_codes)) {
      sb_long <- timeseries_matrix(ts, om_i, mp_i, sb_i, long_i)
      f_long <- timeseries_matrix(ts, om_i, mp_i, f_i, long_i)
      catch_short <- timeseries_matrix(ts, om_i, mp_i, catch_i, short_i)
      catch_long <- timeseries_matrix(ts, om_i, mp_i, catch_i, long_i)
      catch_projection <- timeseries_matrix(ts, om_i, mp_i, catch_i,
        projection_i)
      iacc_long <- if (is.na(iacc_i)) matrix(NA_real_, nrow = length(iters),
        ncol = length(long_i)) else
        timeseries_matrix(ts, om_i, mp_i, iacc_i, long_i)
      valid_kobe <- is.finite(sb_long) & is.finite(f_long)
      summaries[om_i, mp_i, ] <- c(
        finite_mean(ifelse(valid_kobe, sb_long >= 1 & f_long <= 1, NA)),
        finite_mean(ifelse(valid_kobe, sb_long < 1 & f_long <= 1, NA)),
        finite_mean(ifelse(valid_kobe, sb_long >= 1 & f_long > 1, NA)),
        finite_mean(ifelse(valid_kobe, sb_long < 1 & f_long > 1, NA)),
        finite_mean(catch_long < 270),
        finite_mean(catch_short),
        finite_mean(catch_long),
        mean_catch_reduction(catch_projection),
        finite_mean(sb_long),
        finite_mean(f_long),
        finite_mean(iacc_long)
      )
    }
  }

  quilt <- Quilt(Code = summary_codes, Label = summary_labels,
    Description = summary_descriptions)
  Value(quilt) <- summaries
  Check(quilt)

  tradeoff_codes <- c("P(green)", "P(red)", "Catch 2026-2030",
    "Catch 2036-2040", "Catch reduction", "SB/SBMSY", "F/FMSY", "IACC")
  tradeoff_i <- match(tradeoff_codes, summary_codes)
  tradeoff <- Tradeoff(Code = tradeoff_codes,
    Label = summary_labels[tradeoff_i],
    Description = summary_descriptions[tradeoff_i])
  Value(tradeoff) <- summaries[, , tradeoff_i, drop = FALSE]
  Check(tradeoff)

  slick <- Slick()
  Title(slick) <- "SPRFMO Jack Mackerel Candidate MPs"
  Subtitle(slick) <- paste(sub("^tun", "MP", mp_codes), collapse = " and ")
  Date(slick) <- Sys.Date()
  Author(slick) <- "jmMSE"
  Institution(slick) <- "SPRFMO"
  Introduction(slick) <- paste(
    "Candidate MP performance exported from", basename(performance_file),
    "using explicit iteration x OM x MP x indicator x year dimensions."
  )
  MPs(slick) <- make_mp_metadata(mp_codes, registry_file)
  OMs(slick) <- make_om_metadata(run_map, performance_file)
  Boxplot(slick) <- boxplot
  Quilt(slick) <- quilt
  Timeseries(slick) <- ts
  Kobe(slick) <- kobe
  Tradeoff(slick) <- tradeoff

  Check(slick)
  methods::validObject(slick)
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(slick, out_file)
  message("Wrote validated Slick object: ", out_file)
  invisible(slick)
}
