test_that("benchmark definitions are complete", {
    b <- benchmark_presets()
    expect_equal(length(b), 6)
    expect_true(all(paste0("B", 1:6) %in% substr(names(b), 1, 2)))
})
