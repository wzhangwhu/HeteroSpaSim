small_cfg <- function(...) heterospa_config(n_cells = 60, n_domains = 3, n_cell_types = 4,
    modalities = c("RNA", "ATAC", "Protein"), n_features = c(RNA = 20, ATAC = 24, Protein = 8), n_spatial_basis = 5, chunk_size = 10, seed = 9, ...)

test_that("singleton sampling cannot leave the eligible pool", {
    f <- HeteroSpaSim:::.sample_without
    set.seed(1);expect_identical(replicate(20, f(7L, 1L)), rep(7L, 20))
})
test_that("isolated streams fix unrelated biology/loadings/counts", {
    cfg <- small_cfg();x <- simulate_spatial_multiomics(cfg)
    for (nm in c("rna_mean", "rna_dispersion", "rna_library_sd")) {
        c <- cfg;c[[nm]] <- cfg[[nm]] * 1.5;y <- simulate_spatial_multiomics(c)
        expect_identical(x$latent, y$latent)
        expect_identical(x$truth$feature$RNA$shared_loadings, y$truth$feature$RNA$shared_loadings)
        expect_identical(x$assays$ATAC, y$assays$ATAC);expect_identical(x$assays$Protein, y$assays$Protein)
    }
    c <- cfg;c$misalignment <- "severe";y <- simulate_spatial_multiomics(c)
    expect_identical(x$assays, y$assays)
    c <- cfg;c$missing_rate["RNA"] <- .3;y <- simulate_spatial_multiomics(c)
    expect_identical(x$assays$ATAC, y$assays$ATAC);expect_identical(x$assays$Protein, y$assays$Protein)
    c <- cfg;c$chunk_size <- 17L;y <- simulate_spatial_multiomics(c)
    expect_identical(x$assays, y$assays)
    set.seed(99);before <- .Random.seed;simulate_spatial_multiomics(cfg);expect_identical(.Random.seed, before)
})
test_that("validator rejects corruption and nonfinite inputs early", {
    x <- simulate_spatial_multiomics(small_cfg());y <- x;y$assays$RNA@x[1] <- -1
    expect_error(validate_simulation(y, TRUE), "nonnegative")
    y <- x;rownames(y$assays$RNA) <- rev(rownames(y$assays$RNA));expect_error(validate_simulation(y, TRUE), "ID")
    y <- x;y$truth$measurement_annotations$RNA$observed_domain[1] <- 99;expect_error(validate_simulation(y, TRUE), "truth")
    expect_error(heterospa_config(n_cells = 60.5), "integers")
    expect_error(heterospa_config(atac_binary = NA), "Boolean")
    expect_error(heterospa_config(misalignment = "invalid"), "misalignment")
    c <- small_cfg();c$shared_strength[1] <- NA;expect_error(simulate_spatial_multiomics(c), "finite")
})
test_that("domain/type reference labels are independently preserved and IDs align", {
    x <- simulate_spatial_multiomics(small_cfg());ref <- list()
    for (m in c("RNA", "ATAC")) {
        counts <- x$assays[[m]];rownames(counts) <- sprintf("s%d", 1:60)
        ref[[m]] <- list(counts = counts, coordinates = x$spatial$latent_coordinates, domain = x$metadata$domain)
    }
    y <- simulate_from_reference(ref, seed = 9);expect_identical(y$metadata$domain, as.integer(factor(ref$RNA$domain)))
    a <- calibrate_reference(ref);ref$ATAC$counts <- ref$ATAC$counts[60:1, ];b <- calibrate_reference(ref)
    expect_equal(a$shared_strength, b$shared_strength)
    ref$RNA$cell_type <- x$metadata$cell_type;y <- simulate_from_reference(ref, seed = 9);expect_identical(y$metadata$cell_type, as.integer(factor(x$metadata$cell_type)))
    ref$RNA$domain <- NULL;y <- simulate_from_reference(ref, seed = 9);expect_identical(y$metadata$cell_type, as.integer(factor(x$metadata$cell_type)))
    ref$RNA$counts@x[1] <- -2;expect_error(calibrate_reference(ref), "raw counts")
})
test_that("metric ties and zero-rank truth are handled", {
    expect_equal(unname(score_binary_imputation(c(1, 0), c(.5, .5))["AUPRC"]), .5)
    expect_equal(score_binary_imputation(c(1, 0), c(.5, .5)), score_binary_imputation(c(0, 1), c(.5, .5)))
    expect_error(score_latent_recovery(matrix(0, 20, 2), matrix(0, 20, 2)), "non-zero rank")
})
test_that("V1 visible pairs and V2 hidden maps have distinct legal scopes", {
    x <- simulate_spatial_multiomics(small_cfg());v1 <- tempfile();v2 <- tempfile()
    export_benchmark(x, v1, design = "V1");export_benchmark(x, v2, design = "V2")
    expect_true(file.exists(file.path(v1, "public_inputs/visible_pairs/RNA__ATAC.csv")))
    expect_false(dir.exists(file.path(v2, "public_inputs/visible_pairs")))
    expect_false(file.exists(file.path(v2, "public_inputs/RNA/coordinates_true.csv")))
    map <- read.csv(file.path(v2, "evaluator_truth/public_unit_map_RNA.csv"))
    counts <- Matrix::readMM(file.path(v2, "public_inputs/RNA/counts.mtx"))
    expect_equal(as.matrix(counts), unname(as.matrix(x$assays$RNA[map$original_row, ])), ignore_attr = TRUE)
    expect_error(export_benchmark(x, v2), "nonempty")
})
test_that("constant-coordinate spatial batches and pairing constraints are explicit", {
    c <- small_cfg();c$batch_count <- 3L;c$batch_strength <- .2
    xy <- cbind(rep(1, 60), seq_len(60));x <- simulate_spatial_multiomics(c, coordinates = xy);expect_true(validate_simulation(x))
    c <- small_cfg();c$pairing <- 0;c$capture_fraction <- 1
    expect_warning(simulate_spatial_multiomics(c), "constrained")
})
