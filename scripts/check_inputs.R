#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(yaml)
})

config_file <- "config/config.yaml"

fail <- function(msg) {
  cat(msg, "\n")
  quit(status = 1)
}

ok <- function(msg) cat(sprintf("✅ %s\n", msg))
warn <- function(msg) cat(sprintf("⚠️ %s\n", msg))
info <- function(msg) cat(sprintf("ℹ️ %s\n", msg))

if (!file.exists(config_file)) {
  fail("❌ Missing config/config.yaml")
}

cfg <- yaml::read_yaml(config_file)

if (is.null(cfg$paths)) {
  fail("❌ 'paths' section not found in config/config.yaml")
}

required_items <- list(
  ortholog_map = cfg$paths$ortholog_map,
  human_series_matrix = cfg$paths$human_series_matrix,
  human_gct_matrix = cfg$paths$human_gct_matrix
)

optional_items <- list(
  human_series_matrix_external = cfg$paths$human_series_matrix_external,
  human_gct_matrix_external = cfg$paths$human_gct_matrix_external
)

cat("========================================\n")
cat("Input File Check\n")
cat("========================================\n")

info("Checking required input files...")
required_ok <- TRUE

for (nm in names(required_items)) {
  p <- required_items[[nm]]

  if (is.null(p) || !nzchar(p)) {
    cat(sprintf("❌ %s is empty in config/config.yaml\n", nm))
    required_ok <- FALSE
    next
  }

  if (file.exists(p)) {
    ok(sprintf("%s -> %s", nm, p))
  } else {
    cat(sprintf("❌ %s not found -> %s\n", nm, p))
    required_ok <- FALSE
  }
}

info("Checking optional external validation files...")

for (nm in names(optional_items)) {
  p <- optional_items[[nm]]

  if (is.null(p) || !nzchar(p)) {
    warn(sprintf("%s is empty; external validation may be skipped", nm))
    next
  }

  if (file.exists(p)) {
    ok(sprintf("%s -> %s", nm, p))
  } else {
    warn(sprintf("%s not found -> %s", nm, p))
  }
}

if (!required_ok) {
  fail("❌ Required input file check failed.")
}

cat("========================================\n")
ok("All required inputs are available.")
cat("========================================\n")
