test_that("small simulation produces coherent dimensions", {
    cfg <- heterospa_config(
        n_cells = 120, n_domains = 3, n_cell_types = 4,
        modalities = c("RNA", "ATAC", "Protein"),
        n_features = c(RNA = 30, ATAC = 40, Protein = 8),
        n_spatial_basis = 8, chunk_size = 20,
        resolution = c(RNA = 1, ATAC = 2, Protein = 3),
        seed = 7
    )
    sim <- simulate_spatial_multiomics(config = cfg)
    expect_s3_class(sim, "HeteroSpaSim")
    expect_true(validate_simulation(sim))
    expect_equal(ncol(sim$assays$RNA), 30)
    expect_equal(ncol(sim$measurement$operators$ATAC), 120)
    expect_equal(nrow(sim$assays$Protein), nrow(sim$spatial$observed_coordinates$Protein))
})

test_that("complete truth is retained by default", {
    sim <- simulate_spatial_multiomics(
        n_cells = 80, modalities = c("RNA", "Protein"),
        n_features = c(RNA = 20, Protein = 6),
        n_domains = 3, n_cell_types = 3,
        n_spatial_basis = 6, chunk_size = 20, seed = 2
    )
    expect_true(all(c("RNA", "Protein") %in% names(sim$truth$complete_assays)))
})
