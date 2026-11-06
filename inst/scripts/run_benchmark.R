#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
bench <- if (length(args) >= 1L) args[[1L]] else "B1_shared_strength"
out <- if (length(args) >= 2L) args[[2L]] else paste0("benchmark_", bench)
reps <- if (length(args) >= 3L) as.integer(args[[3L]]) else 10L

suppressPackageStartupMessages(library(HeteroSpaSim))
generate_benchmark(bench, n_replicates = reps, out_dir = out)
cat("Benchmark written to", normalizePath(out), "
")
