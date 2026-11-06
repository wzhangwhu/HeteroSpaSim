test_that("domain scores are perfect for identical labels", {
    z <- c(1, 1, 2, 2, 3, 3)
    sc <- score_domains(z, z)
    expect_equal(unname(sc["ARI"]), 1)
    expect_equal(unname(sc["NMI"]), 1)
    expect_equal(unname(sc["Purity"]), 1)
})

test_that("alignment score is zero for identical coordinates", {
    xy <- cbind(c(0, 1), c(0, 1))
    expect_equal(unname(score_alignment(xy, xy)["RMSE"]), 0)
})

test_that("latent recovery recognizes identical subspaces", {
    set.seed(1)
    h <- matrix(rnorm(100), 20, 5)
    sc <- score_latent_recovery(h, h)
    expect_equal(unname(sc["NormalizedProjectionDistance"]), 0)
    expect_equal(unname(sc["MeanCanonicalCorrelation"]), 1)
})


test_that("latent recovery is invariant to full-rank rotation", {
    set.seed(2)
    h <- matrix(rnorm(120), 24, 5)
    R <- qr.Q(qr(matrix(rnorm(25), 5, 5)))
    sc <- score_latent_recovery(h, h %*% R)
    expect_lt(unname(sc["NormalizedProjectionDistance"]), 1e-7)
    expect_gt(unname(sc["MeanCanonicalCorrelation"]), 1 - 1e-10)
})

test_that("binary imputation rejects mismatched lengths", {
    expect_error(score_binary_imputation(c(0, 1, 0), c(0.1, 0.9)), "same number")
})
