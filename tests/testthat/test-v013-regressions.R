test_that("observed correspondence follows modality missingness", {
    cfg <- heterospa_config(
        n_cells = 160, n_domains = 3, n_cell_types = 4,
        modalities = c("RNA", "ATAC", "Protein"),
        n_features = c(RNA = 20, ATAC = 25, Protein = 8),
        missing_rate = c(RNA = 0.15, ATAC = 0.25, Protein = 0.10),
        n_spatial_basis = 6, chunk_size = 20, seed = 21
    )
    sim <- simulate_spatial_multiomics(config = cfg)
    expect_equal(nrow(sim$correspondence$RNA__ATAC), nrow(sim$assays$RNA))
    expect_equal(ncol(sim$correspondence$RNA__ATAC), nrow(sim$assays$ATAC))
    expect_equal(nrow(sim$spatial$observed_true_coordinates$Protein), nrow(sim$assays$Protein))
    expect_gt(nrow(sim$truth$complete_correspondence$RNA__ATAC), nrow(sim$correspondence$RNA__ATAC))
})

test_that("multi-modal mosaic overlap is close to the target for every modality pair", {
    cfg <- heterospa_config(
        n_cells = 900, modalities = c("RNA", "ATAC", "Protein"),
        n_features = c(RNA = 10, ATAC = 12, Protein = 6),
        overlap = 0.50, n_domains = 3, n_cell_types = 4,
        n_spatial_basis = 6, chunk_size = 10, seed = 5
    )
    expect_warning(sim <- simulate_spatial_multiomics(config = cfg), "pairing is constrained")
    vals <- unlist(sim$truth$realized$spatial_overlap)
    expect_true(all(abs(vals - 0.50) < 0.03))
})

test_that("binary ATAC remains binary after spot aggregation", {
    cfg <- heterospa_config(
        n_cells = 180, modalities = c("ATAC"), n_features = c(ATAC = 30),
        resolution = c(ATAC = 5), atac_binary = TRUE,
        n_domains = 3, n_cell_types = 4, n_spatial_basis = 6,
        chunk_size = 15, seed = 12
    )
    sim <- simulate_spatial_multiomics(config = cfg)
    vals <- unique(as.numeric(sim$assays$ATAC))
    expect_true(all(vals %in% c(0, 1)))
})

test_that("invalid positive-scale parameters fail early", {
    expect_error(heterospa_config(rna_dispersion = -1))
    expect_error(heterospa_config(shared_length_scale = 0))
    expect_error(heterospa_config(protein_ambient = -0.1))
})
