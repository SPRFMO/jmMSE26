# Profile the completed 500-draw CMP controls immediately below and above
# each selected trigger, using the same OM, OEM, random stream, and tuning
# statistic as the refinement run.

argv <- commandArgs(trailingOnly = TRUE)
option <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), argv, value = TRUE)
  if (!length(hit)) return(default)
  sub(paste0("^--", name, "="), "", hit[[length(hit)]])
}

project_root <- normalizePath(getwd())
jm_root <- normalizePath(option("jm-root", file.path("..", "jmMSE-500-refine")),
  mustWork = TRUE)
output_dir <- file.path(project_root, option(
  "output-dir", file.path("output", "candidate-trigger-tolerance-500")
))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cores <- as.integer(option("cores", "5"))
if (!is.finite(cores) || cores < 1L) stop("--cores must be positive")

old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(jm_root)
source("config.R")

library(data.table)

refinement_dir <- file.path("model", "tune", "refine_500_from_100")
summary <- fread(file.path(refinement_dir, "tuning_summary.csv"))
trace <- fread(file.path(refinement_dir, "tuning_trace.csv"))
tuned <- readRDS(file.path(refinement_dir, "runs.rds"))
omc <- readRDS(file.path("data", "om11_h1_0.16_065.rds"))

stopifnot(
  dims(omc$om)$iter == 500L,
  oem_iters(omc$oem) == 500L,
  setequal(summary$mp, names(tuned))
)

tuning_years <- 2041:2050
objective <- 0.60

kobe_green_probability <- function(run, dSB0, years) {
  projected_om <- om(run)
  metric_groups <- metrics(projected_om)
  reference_points <- refpts(projected_om)
  if (length(metric_groups) != 1L) {
    stop("Reference tolerance profile expects one managed stock")
  }
  stock_code <- names(metric_groups)[1L]
  stock_refpts <- if (is(reference_points, "FLPars")) {
    reference_points[[stock_code]]
  } else {
    reference_points
  }
  dynamic_bmsy <- stock_refpts["SBMSY", ] / stock_refpts["SB0", ] *
    dSB0[[stock_code]][, ac(years)]
  stock_metrics <- metric_groups[[stock_code]]
  green <- stock_metrics$SB[, ac(years)] / dynamic_bmsy >= 1 &
    stock_metrics$F[, ac(years)] / stock_refpts["FMSY", ] <= 1
  mean(as.numeric(green), na.rm = TRUE)
}

# Use the final bisection resolution for each rule family: one sixteenth for
# hockeystick triggers and one thirty-second for power-ramp triggers.
grid <- summary[, .(
  trigger = trigger + c(-1, 0, 1) *
    if (hcr == "hockeystick") 0.0625 else 0.03125,
  position = c("lower", "selected", "upper"),
  step_size = if (hcr == "hockeystick") 0.0625 else 0.03125
), by = .(mp, label, short_label, hcr, selected_trigger = trigger,
  selected_probability = tuning_probability)]

existing <- trace[, .SD[.N], by = .(mp, trigger)]
grid[, `:=`(tuning_probability = NA_real_, source = NA_character_)]

for (i in seq_len(nrow(grid))) {
  row <- grid[i]
  if (row$position == "selected") {
    grid[i, `:=`(
      tuning_probability = selected_probability,
      source = "selected control"
    )]
    next
  }

  hit <- existing[mp == row$mp & abs(trigger - row$trigger) < 1e-12]
  if (nrow(hit)) {
    grid[i, `:=`(
      tuning_probability = hit$tuning_probability[[1]],
      source = "existing refinement trace"
    )]
    next
  }

  message(
    "Profiling ", row$mp, " ", row$position,
    ": trigger=", format(row$trigger, digits = 8)
  )
  candidate_control <- control(tuned[[row$mp]])
  args(candidate_control$hcr)$trigger <- row$trigger
  set.seed(43000L + sum(utf8ToInt(row$mp)))
  run <- mp(
    omc$om,
    oem = omc$oem,
    ctrl = candidate_control,
    args = list(iy = 2025L, fy = 2050L),
    parallel = cores > 1L,
    verbose = FALSE
  )
  probability <- kobe_green_probability(
    run, omc$unfishedSSB, tuning_years
  )
  grid[i, `:=`(
    tuning_probability = probability,
    source = "new local profile run"
  )]
  fwrite(grid, file.path(output_dir, "trigger_tolerance_profile.csv"))
  rm(run)
  gc()
}

grid[, position_order := match(position, c("lower", "selected", "upper"))]
setorder(grid, mp, position_order)
grid[, position_order := NULL]
grid[, difference_from_objective := tuning_probability - objective]
grid[, difference_from_selected :=
  tuning_probability - selected_probability]
fwrite(grid, file.path(output_dir, "trigger_tolerance_profile.csv"))

profile_summary <- dcast(
  grid,
  mp + label + short_label + hcr + selected_trigger +
    selected_probability + step_size ~ position,
  value.var = c("trigger", "tuning_probability")
)
profile_summary[, `:=`(
  probability_change_lower =
    tuning_probability_lower - tuning_probability_selected,
  probability_change_upper =
    tuning_probability_upper - tuning_probability_selected,
  local_slope =
    (tuning_probability_upper - tuning_probability_lower) /
    (trigger_upper - trigger_lower),
  monotonic = tuning_probability_lower <= tuning_probability_selected &
    tuning_probability_selected <= tuning_probability_upper
)]
fwrite(profile_summary,
  file.path(output_dir, "trigger_tolerance_summary.csv"))
writeLines(capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo.txt"))

if (anyNA(grid$tuning_probability)) stop("Tolerance profile is incomplete")
if (!all(profile_summary$monotonic)) {
  warning("At least one local profile is not monotonic")
}

print(profile_summary)
message("Saved complete trigger-tolerance profile to ", output_dir)
