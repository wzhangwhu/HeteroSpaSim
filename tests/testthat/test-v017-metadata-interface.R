test_that("metadata aliases are present and consistent", {
  sim <- simulate_spatial_multiomics(
    n_cells = 120,
    modalities = c("RNA", "ATAC"),
    n_features = c(RNA = 20, ATAC = 30),
    n_domains = 3,
    n_cell_types = 4,
    seed = 17
  )
  expect_length(sim$metadata$domain, 120)
  expect_length(sim$metadata$cell_type, 120)
  expect_length(sim$metadata$batch, 120)
  expect_identical(sim$metadata$domain, sim$metadata$cells$domain)
  expect_identical(sim$metadata$cell_type, sim$metadata$cells$cell_type)
  expect_identical(sim$metadata$batch, sim$metadata$cells$batch)
  expect_equal(sim$metadata$cells$x, sim$spatial$latent_coordinates[,1])
  expect_equal(sim$metadata$cells$y, sim$spatial$latent_coordinates[,2])
})
