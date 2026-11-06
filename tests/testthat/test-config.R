test_that("configuration is valid and named", {
    cfg <- heterospa_config(n_cells = 100, modalities = c("RNA", "ATAC"), n_features = c(40, 50))
    expect_equal(cfg$n_cells, 100L)
    expect_equal(names(cfg$n_features), c("RNA", "ATAC"))
    expect_true(all(cfg$shared_strength >= 0 & cfg$shared_strength <= 1))
})

test_that("invalid parameters fail", {
    expect_error(heterospa_config(shared_strength = 1.2))
    expect_error(heterospa_config(overlap = -0.1))
})
