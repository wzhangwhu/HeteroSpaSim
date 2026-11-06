test_that("sparse correspondence is compatible with Matrix sparseSummary behaviour", {
    cfg <- heterospa_config(
        n_cells = 120,
        n_domains = 3,
        n_cell_types = 4,
        modalities = c("RNA", "ATAC"),
        n_features = c(RNA = 15, ATAC = 20),
        resolution = c(RNA = 1, ATAC = 2),
        pairing = 0.5,
        overlap = 0.7,
        n_spatial_basis = 6,
        chunk_size = 15,
        seed = 1401
    )
    expect_no_error(sim <- simulate_spatial_multiomics(config = cfg))
    expect_s3_class(sim, "HeteroSpaSim")
    expect_equal(nrow(sim$correspondence$RNA__ATAC), nrow(sim$assays$RNA))
    expect_equal(ncol(sim$correspondence$RNA__ATAC), nrow(sim$assays$ATAC))
})

test_that("empty sparse correspondence returns a valid Matrix with dimnames", {
    A <- Matrix::sparseMatrix(i = 1:2, j = 1:2, x = 1, dims = c(2, 4),
                              dimnames = list(c("a1", "a2"), NULL))
    B <- Matrix::sparseMatrix(i = 1:2, j = 3:4, x = 1, dims = c(2, 4),
                              dimnames = list(c("b1", "b2"), NULL))
    C <- HeteroSpaSim:::.correspondence_sparse(A, B)
    expect_true(inherits(C, "Matrix"))
    expect_equal(dim(C), c(2, 2))
    expect_equal(Matrix::nnzero(C), 0)
    expect_equal(rownames(C), c("a1", "a2"))
    expect_equal(colnames(C), c("b1", "b2"))
})


test_that("correspondence scoring accepts tied maximum-overlap targets", {
    C <- Matrix::Matrix(matrix(c(0.5, 0.5, 0, 0,
                                 0,   0,   1, 0), nrow = 2, byrow = TRUE), sparse = TRUE)
    sc <- score_correspondence(C, c(2L, 3L))
    expect_equal(unname(sc["Top1Accuracy"]), 1)
})

test_that("benchmark replicates use matched seeds across parameter levels", {
    b <- generate_benchmark("B1_shared_strength", n_replicates = 2,
                            n_cells = 40, modalities = c("RNA", "ATAC"),
                            base_config = heterospa_config(
                                n_cells = 40, modalities = c("RNA", "ATAC"),
                                n_features = c(RNA = 10, ATAC = 12),
                                n_spatial_basis = 5, chunk_size = 10,
                                seed = 1
                            ), seed_start = 17)
    man <- b$manifest
    expect_equal(unique(man$seed[man$replicate == 1]), 17L)
    expect_equal(unique(man$seed[man$replicate == 2]), 18L)
})
