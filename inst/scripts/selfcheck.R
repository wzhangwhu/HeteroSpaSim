#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(HeteroSpaSim))
cat("HeteroSpaSim version:", as.character(packageVersion("HeteroSpaSim")), "\n")
cat("R version:", R.version.string, "\n")
cat("Matrix version:", as.character(packageVersion("Matrix")), "\n")

cat("1. Small multimodal simulation...\n")
cfg <- heterospa_config(
    n_cells = 150,
    n_domains = 3,
    n_cell_types = 4,
    modalities = c("RNA", "ATAC", "Protein"),
    n_features = c(RNA = 30, ATAC = 40, Protein = 8),
    n_spatial_basis = 8,
    chunk_size = 20,
    resolution = c(RNA = 1, ATAC = 2, Protein = 3),
    pairing = 0.5,
    capture_fraction = 0.4,
    overlap = 0.7,
    misalignment = "mild",
    missing_rate = 0.1,
    seed = 101
)
sim <- simulate_spatial_multiomics(config = cfg)
stopifnot(isTRUE(validate_simulation(sim)))
print(summary(sim))

cat("1a. Metadata interface...\n")
stopifnot(
    length(sim$metadata$domain) == sim$parameters$n_cells,
    identical(sim$metadata$domain, sim$metadata$cells$domain),
    identical(sim$metadata$cell_type, sim$metadata$cells$cell_type),
    identical(sim$metadata$batch, sim$metadata$cells$batch)
)

cat("1b. Sparse-correspondence compatibility...\n")
A0 <- Matrix::sparseMatrix(i = 1:2, j = 1:2, x = 1, dims = c(2, 4))
B0 <- Matrix::sparseMatrix(i = 1:2, j = 3:4, x = 1, dims = c(2, 4))
C0 <- HeteroSpaSim:::.correspondence_sparse(A0, B0)
stopifnot(inherits(C0, "Matrix"), identical(dim(C0), c(2L, 2L)), Matrix::nnzero(C0) == 0L)

cat("2. Metrics...\n")
z <- sim$truth$measurement_annotations$RNA$observed_domain
stopifnot(abs(score_domains(z, z)["ARI"] - 1) < 1e-12)
xy <- sim$spatial$observed_true_coordinates$RNA
stopifnot(abs(score_alignment(xy, xy)["RMSE"]) < 1e-12)
stopifnot(nrow(sim$correspondence$RNA__ATAC) == nrow(sim$assays$RNA))
stopifnot(ncol(sim$correspondence$RNA__ATAC) == nrow(sim$assays$ATAC))
lat_self <- score_latent_recovery(sim$latent$shared, sim$latent$shared)
stopifnot(unname(lat_self["MeanCanonicalCorrelation"]) == 1)
stopifnot(unname(lat_self["NormalizedProjectionDistance"]) == 0)
# Rotation-equivalent representations span the same subspace and should also
# have zero projection distance up to numerical precision.
set.seed(991)
Qrot <- qr.Q(qr(matrix(stats::rnorm(ncol(sim$latent$shared)^2),
                            ncol(sim$latent$shared), ncol(sim$latent$shared))))
lat_rot <- score_latent_recovery(sim$latent$shared, sim$latent$shared %*% Qrot)
stopifnot(unname(lat_rot["NormalizedProjectionDistance"]) < 1e-7)

cat("3. Generic export...\n")
out <- tempfile("heterospa_selfcheck_")
export_benchmark(sim, out)
stopifnot(file.exists(file.path(out, "RNA", "counts.mtx")))
stopifnot(file.exists(file.path(out, "truth", "latent_cells.csv")))
unlink(out, recursive = TRUE)

cat("4. Binary ATAC at spot resolution...\n")
at <- simulate_spatial_multiomics(
    n_cells = 120, modalities = "ATAC", n_features = c(ATAC = 20),
    resolution = c(ATAC = 4), n_domains = 3, n_cell_types = 3,
    n_spatial_basis = 6, chunk_size = 20, seed = 102
)
stopifnot(all(unique(as.numeric(at$assays$ATAC)) %in% c(0, 1)))

cat("5. Three-modality mosaic overlap...\n")
mo <- simulate_spatial_multiomics(
    n_cells = 600, modalities = c("RNA", "ATAC", "Protein"),
    n_features = c(RNA = 12, ATAC = 15, Protein = 6), overlap = 0.5,
    n_domains = 3, n_cell_types = 4, n_spatial_basis = 6, chunk_size = 12, seed = 103
)
ov <- unlist(mo$truth$realized$spatial_overlap)
stopifnot(all(abs(ov - 0.5) < 0.05))

cat("6. Benchmark definitions...\n")
stopifnot(length(benchmark_presets()) == 6L)

cat("HeteroSpaSim self-check completed successfully.\n")
