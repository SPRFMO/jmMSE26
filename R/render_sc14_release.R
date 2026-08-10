#!/usr/bin/env Rscript

# Render the SC14 package with a fixed creation time so repeat PDF builds from
# the same source have the same file fingerprint.

suppressPackageStartupMessages(library(yaml))

args <- commandArgs(trailingOnly = TRUE)
manifest_path <- file.path("release", "sc14-mse-2026-rc1.yml")
manifest <- yaml::read_yaml(manifest_path)
release_date <- as.Date(manifest$prepared)
release_epoch <- as.numeric(as.POSIXct(release_date, tz = "UTC"))
Sys.setenv(SOURCE_DATE_EPOCH = format(release_epoch, scientific = FALSE))

render_target <- if (length(args)) args[[1L]] else character()
old_working_directory <- getwd()
on.exit(setwd(old_working_directory), add = TRUE)
setwd("doc")

command_args <- c("render", render_target)
status <- system2("quarto", command_args)
if (!identical(status, 0L)) {
  stop("Quarto release render failed with status ", status)
}

cat("Rendered ", manifest$release_id,
  " with SOURCE_DATE_EPOCH=", format(release_epoch, scientific = FALSE),
  "\n", sep = "")
