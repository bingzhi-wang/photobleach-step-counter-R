#!/usr/bin/env Rscript
# Batch photobleaching step counting from the command line.
#
# Author: Bingzhi Wang (2023), ported to R (2026). Built on the step-finding code of
# Konstantinos Tsekouras and Sina Jazani (Copyright (C) 2015, GNU GPL v3 or later).
#
# Usage
# -----
#   # 1st run: set the priors by clicking on each trace, then count steps
#   Rscript count_photobleaching.R traces.txt -o results
#
#   # re-run later without clicking (e.g. with another window size)
#   Rscript count_photobleaching.R traces.txt -o results_w50 --priors results/priors.csv -w 50
#
# Options
#   -o, --outdir DIR     output directory (default: results)
#   -w, --window N       window size in frames (default: 100)
#   --priors FILE        priors.csv from an earlier run (skips the clicking)
#   --only I [J ...]     analyse only these 0-based trajectory indices
#   -h, --help           show this help

usage <- function() {
  cat("Usage: Rscript count_photobleaching.R TRACES [-o DIR] [-w N] [--priors FILE] [--only I J ...]\n",
      "\n  TRACES             text file, one trajectory per column",
      "\n  -o, --outdir DIR   output directory (default: results)",
      "\n  -w, --window N     window size in frames (default: 100)",
      "\n  --priors FILE      priors.csv from an earlier run (skips the clicking)",
      "\n  --only I [J ...]   analyse only these 0-based trajectory indices\n", sep = "")
}

# --- load the package: prefer the source next to this script (a cloned repo) ---------
script_dir <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else getwd()
})
src <- file.path(script_dir, "R")
if (dir.exists(src) && file.exists(file.path(script_dir, "DESCRIPTION"))) {
  for (f in list.files(src, pattern = "\\.R$", full.names = TRUE)) source(f)
} else {
  suppressPackageStartupMessages(library(photobleach))
}

# --- parse arguments -----------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || any(args %in% c("-h", "--help"))) { usage(); quit(status = if (length(args)) 0 else 1) }

opt <- list(traces = NULL, outdir = "results", window = 100, priors = NULL, only = NULL)
i <- 1
while (i <= length(args)) {
  a <- args[i]
  if (a %in% c("-o", "--outdir"))      { opt$outdir <- args[i + 1]; i <- i + 2 }
  else if (a %in% c("-w", "--window")) { opt$window <- as.integer(args[i + 1]); i <- i + 2 }
  else if (a == "--priors")            { opt$priors <- args[i + 1]; i <- i + 2 }
  else if (a == "--only") {
    i <- i + 1
    while (i <= length(args) && grepl("^[0-9]+$", args[i])) {
      opt$only <- c(opt$only, as.integer(args[i])); i <- i + 1
    }
  }
  else if (is.null(opt$traces) && !startsWith(a, "-")) { opt$traces <- a; i <- i + 1 }
  else { usage(); stop("unknown argument: ", a, call. = FALSE) }
}
if (is.null(opt$traces)) { usage(); quit(status = 1) }

invisible(run_batch(opt$traces, outdir = opt$outdir, window = opt$window,
                    priors = opt$priors, only = opt$only))
