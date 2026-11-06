test_that("sparse correspondence works for multiple modalities", {
  cfg <- heterospa_config(
    n_cells = 120,
    n_domains = 3,
    n_cell_types = 4,
    modalities = c("RNA", "ATAC", "Protein"),
    n_features = c(RNA = 20, ATAC = 25, Protein = 10),
    resolution = c(RNA = 1, ATAC = 2, Protein = 3),
    pairing = 0.5,
    overlap = 0.7,
    seed = 11
  )
  sim <- simulate_spatial_multiomics(config = cfg)
  expect_s3_class(sim, "HeteroSpaSim")
  expect_true(length(sim$correspondence) == 3L)
  expect_true(all(vapply(sim$correspondence, inherits, logical(1), what = "Matrix")))
})
