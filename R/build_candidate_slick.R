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
  Timeseries(slick) <- ts
  Kobe(slick) <- kobe

  Check(slick)
  methods::validObject(slick)
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  saveRDS(slick, out_file)
  message("Wrote validated Slick object: ", out_file)
  invisible(slick)
}

