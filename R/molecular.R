.feature_prefix <- function(modality) {
    switch(modality, RNA = "Gene", ATAC = "Peak", Protein = "Protein", "Feature")
}

.generate_counts_block <- function(modality, eta, idx, cfg, batch, batch_effects, baseline_all, dispersion_all) {
    n <- nrow(eta)
    if (!is.null(batch_effects)) eta <- eta + batch_effects[batch, idx, drop = FALSE]
    eta <- .clip(eta, -6, 6)
    baseline <- baseline_all[idx]

    if (modality == "RNA") {
        lib <- attr(baseline_all, "library")
        mu <- exp(sweep(eta, 2L, baseline, "+"))
        mu <- mu * lib
        mu <- .clip(mu, 1e-8, 1e5)
        th <- dispersion_all[idx]
        x <- if (identical(cfg$rng_mode, "isolated")) as.vector(vapply(seq_along(idx), function(j) .with_rng_stream(cfg, paste0("counts/", modality, "/", idx[j]), stats::rnbinom(n, size = th[j], mu = mu[, j])), numeric(n))) else stats::rnbinom(length(mu), size = rep(th, each = n), mu = as.vector(mu))
        matrix(x, nrow = n, ncol = length(idx))
    } else if (modality == "ATAC") {
        prob <- stats::plogis(sweep(eta, 2L, baseline, "+"))
        if (cfg$atac_binary) {
            x <- if (identical(cfg$rng_mode, "isolated")) as.vector(vapply(seq_along(idx), function(j) .with_rng_stream(cfg, paste0("counts/", modality, "/", idx[j]), stats::rbinom(n, size = 1, prob = prob[, j])), numeric(n))) else stats::rbinom(length(prob), size = 1, prob = as.vector(prob))
        } else {
            lambda <- pmax(as.vector(prob) * cfg$atac_depth, 1e-8)
            x <- if (identical(cfg$rng_mode, "isolated")) as.vector(vapply(seq_along(idx), function(j) .with_rng_stream(cfg, paste0("counts/", modality, "/", idx[j]), stats::rpois(n, lambda = lambda[(j - 1L) * n + seq_len(n)])), numeric(n))) else stats::rpois(length(lambda), lambda = lambda)
        }
        matrix(x, nrow = n, ncol = length(idx))
    } else {
        lib <- attr(baseline_all, "library")
        mu <- exp(sweep(eta, 2L, baseline, "+"))
        mu <- .clip(mu * lib, 1e-8, 1e5)
        th <- dispersion_all[idx]
        x <- if (identical(cfg$rng_mode, "isolated")) as.vector(vapply(seq_along(idx), function(j) .with_rng_stream(cfg, paste0("counts/", modality, "/", idx[j]), stats::rnbinom(n, size = th[j], mu = mu[, j])), numeric(n))) else stats::rnbinom(length(mu), size = rep(th, each = n), mu = as.vector(mu))
        ambient <- if (identical(cfg$rng_mode, "isolated")) as.vector(vapply(idx, function(j) .with_rng_stream(cfg, paste0("ambient/", modality, "/", j), stats::rpois(n, lambda = cfg$protein_ambient)), numeric(n))) else stats::rpois(length(mu), lambda = cfg$protein_ambient)
        matrix(x + ambient, nrow = n, ncol = length(idx))
    }
}

.simulate_cell_assay <- function(modality, cfg, H0, Hs, batch, A = NULL, keep_cell = FALSE) {
    if (!is.null(cfg$reference_distribution_profile)) {
        return(.simulate_reference_assay(modality, cfg, H0, Hs, batch, A, keep_cell))
    }
    n <- nrow(H0)
    p <- cfg$n_features[[modality]]
    q0 <- ncol(H0); qs <- ncol(Hs)
    W0 <- .simulate_loadings(p, q0, density = cfg$loading_density, cfg = cfg, key = paste0("loadings/", modality, "/shared"))
    Ws <- .simulate_loadings(p, qs, density = cfg$loading_density, cfg = cfg, key = paste0("loadings/", modality, "/specific"))
    rho <- cfg$shared_strength[[modality]]

    if (modality == "RNA") {
        baseline <- .with_rng_stream(cfg, paste0("baseline/", modality), stats::rnorm(p, log(cfg$rna_mean), 0.8))
        dispersion <- .with_rng_stream(cfg, paste0("dispersion/", modality), exp(stats::rnorm(p, log(cfg$rna_dispersion), 0.25)))
        attr(baseline, "library") <- .with_rng_stream(cfg, paste0("library/", modality), stats::rlnorm(n, 0, cfg$rna_library_sd))
    } else if (modality == "ATAC") {
        baseline <- .with_rng_stream(cfg, paste0("baseline/", modality), stats::rnorm(p, stats::qlogis(.clip(cfg$atac_accessibility, 1e-4, 1 - 1e-4)), 0.9))
        dispersion <- rep(NA_real_, p)
    } else {
        baseline <- .with_rng_stream(cfg, paste0("baseline/", modality), stats::rnorm(p, log(cfg$protein_mean), 0.55))
        dispersion <- .with_rng_stream(cfg, paste0("dispersion/", modality), exp(stats::rnorm(p, log(cfg$protein_dispersion), 0.20)))
        attr(baseline, "library") <- .with_rng_stream(cfg, paste0("library/", modality), stats::rlnorm(n, 0, 0.20))
    }

    batch_effects <- NULL
    if (cfg$batch_count > 1L && cfg$batch_strength > 0) {
        batch_effects <- .with_rng_stream(cfg, paste0("batch_effects/", modality), matrix(stats::rnorm(cfg$batch_count * p, 0, cfg$batch_strength), cfg$batch_count, p))
    }

    starts <- seq.int(1L, p, by = cfg$chunk_size)
    measured_blocks <- vector("list", length(starts))
    cell_blocks <- if (keep_cell) vector("list", length(starts)) else NULL
    for (b in seq_along(starts)) {
        idx <- starts[[b]]:min(p, starts[[b]] + cfg$chunk_size - 1L)
        u0 <- H0 %*% t(W0[idx, , drop = FALSE])
        us <- Hs %*% t(Ws[idx, , drop = FALSE])
        u0 <- .standardize_cols(u0)
        us <- .standardize_cols(us)
        eta <- cfg$signal_strength * (sqrt(rho) * u0 + sqrt(1 - rho) * us)
        x <- .generate_counts_block(modality, eta, idx, cfg, batch, batch_effects, baseline, dispersion)
        xb <- Matrix::Matrix(x, sparse = TRUE)
        if (is.null(A)) {
            measured_blocks[[b]] <- xb
        } else {
            agg <- A %*% xb
            # Preserve the declared ATAC binary observation mode after
            # cell-to-spot aggregation.  A spot is accessible if at least one
            # contributing latent cell is accessible at that peak.
            if (modality == "ATAC" && isTRUE(cfg$atac_binary)) {
                agg <- Matrix::Matrix(1 * (agg > 0), sparse = TRUE)
            }
            measured_blocks[[b]] <- agg
        }
        if (keep_cell) cell_blocks[[b]] <- xb
        rm(u0, us, eta, x, xb)
    }
    measured <- .sparse_cbind(measured_blocks)
    cell_counts <- if (keep_cell) .sparse_cbind(cell_blocks) else NULL
    feat <- sprintf("%s_%06d", .feature_prefix(modality), seq_len(p))
    colnames(measured) <- feat
    if (!is.null(A)) rownames(measured) <- rownames(A) else rownames(measured) <- sprintf("Cell_%06d", seq_len(n))
    if (!is.null(cell_counts)) {
        colnames(cell_counts) <- feat
        rownames(cell_counts) <- sprintf("Cell_%06d", seq_len(n))
    }
    rownames(W0) <- feat; rownames(Ws) <- feat
    list(counts = measured, cell_counts = cell_counts, shared_loadings = W0, specific_loadings = Ws,
         baseline = as.numeric(baseline), dispersion = dispersion,
         feature_names = feat, shared_strength = rho,
         library_factors = if (identical(cfg$rng_mode, "isolated")) attr(baseline, "library") else NULL,
         batch_effects = if (identical(cfg$rng_mode, "isolated")) batch_effects else NULL)
}

.aggregate_assay <- function(cell_counts, A) {
    out <- A %*% cell_counts
    if (!inherits(out, "Matrix")) out <- Matrix::Matrix(out, sparse = TRUE)
    rownames(out) <- rownames(A)
    colnames(out) <- colnames(cell_counts)
    out
}

.missing_unit_mask <- function(rate, mechanism, A, coords, domains) {
    n <- nrow(A)
    if (rate <= 0) return(rep(TRUE, n))
    if (mechanism == "mcar") {
        miss <- stats::runif(n) < rate
    } else if (mechanism == "domain") {
        dz <- .row_mode_from_operator(A, domains)
        score <- ifelse(dz == max(domains), 1.8, 0.6)
        prob <- .clip(rate * score / mean(score), 0, 0.95)
        miss <- stats::runif(n) < prob
    } else {
        x <- coords[, 1L]
        sx <- stats::sd(x)
        z <- if (is.finite(sx) && sx > 1e-8) (x - mean(x)) / sx else rep(0, length(x))
        score <- exp(z)
        prob <- .clip(rate * score / mean(score), 0, 0.95)
        miss <- stats::runif(n) < prob
    }
    if (all(miss)) miss[sample.int(n, 1L)] <- FALSE
    !miss
}
