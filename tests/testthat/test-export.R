test_that("generic export creates expected files", {
    sim <- simulate_spatial_multiomics(
        n_cells = 60, modalities = c("RNA", "ATAC"),
        n_features = c(RNA = 12, ATAC = 15), n_domains = 2,
        n_cell_types = 3, n_spatial_basis = 5, chunk_size = 12, seed = 3
    )
    d <- tempfile("heterospa_export_")
    export_benchmark(sim, d)
    expect_true(file.exists(file.path(d, "public_inputs", "RNA", "counts.mtx")))
    expect_true(file.exists(file.path(d, "evaluator_truth", "truth", "latent_cells.csv")))
    expect_true(file.exists(file.path(d, "evaluator_truth", "config.json")) || file.exists(file.path(d, "evaluator_truth", "config.dput")))
})
