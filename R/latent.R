.factor_component <- function(domains, celltypes, coords, q, n_domains, n_cell_types,
                              n_basis, length_scale, factor_noise,
                              domain_sd = 1, celltype_sd = 0.8, spatial_sd = 0.8,
                              cfg = NULL, key = "latent") {
    dom_eff <- .with_rng_stream(cfg, paste0(key, "/domain"), matrix(stats::rnorm(n_domains * q, 0, domain_sd), n_domains, q))
    ct_eff <- .with_rng_stream(cfg, paste0(key, "/celltype"), matrix(stats::rnorm(n_cell_types * q, 0, celltype_sd), n_cell_types, q))
    basis <- .rbf_basis(coords, n_basis = n_basis, length_scale = length_scale)
    coef <- .with_rng_stream(cfg, paste0(key, "/spatial"), matrix(stats::rnorm(ncol(basis$phi) * q, 0, spatial_sd), ncol(basis$phi), q))
    spatial <- basis$phi %*% coef
    H <- dom_eff[domains, , drop = FALSE] + ct_eff[celltypes, , drop = FALSE] + spatial +
        .with_rng_stream(cfg, paste0(key, "/noise"), matrix(stats::rnorm(nrow(coords) * q, 0, factor_noise), nrow(coords), q))
    list(H = .standardize_cols(H), domain_effect = dom_eff, celltype_effect = ct_eff,
         spatial_coefficients = coef, knots = basis$knots)
}

.orthogonalize_specific <- function(Hs, H0) {
    H0 <- .standardize_cols(H0)
    qr0 <- qr(H0)
    Q <- qr.Q(qr0)[, seq_len(qr0$rank), drop = FALSE]
    Hs2 <- Hs - Q %*% crossprod(Q, Hs)
    .standardize_cols(Hs2)
}

.simulate_latent_biology <- function(cfg, coords, domains, celltypes) {
    shared <- .factor_component(
        domains, celltypes, coords, q = cfg$shared_dim,
        n_domains = cfg$n_domains, n_cell_types = cfg$n_cell_types,
        n_basis = cfg$n_spatial_basis, length_scale = cfg$shared_length_scale,
        factor_noise = cfg$factor_noise, cfg = cfg, key = "latent/shared"
    )
    specific <- list()
    specific_meta <- list()
    for (m in cfg$modalities) {
        sp <- .factor_component(
            domains, celltypes, coords, q = cfg$specific_dim[[m]],
            n_domains = cfg$n_domains, n_cell_types = cfg$n_cell_types,
            n_basis = cfg$n_spatial_basis,
            length_scale = cfg$specific_length_scale,
            factor_noise = cfg$factor_noise,
            domain_sd = 0.8, celltype_sd = 0.7, spatial_sd = 0.7,
            cfg = cfg, key = paste0("latent/specific/", m)
        )
        specific[[m]] <- .orthogonalize_specific(sp$H, shared$H)
        specific_meta[[m]] <- sp[names(sp) != "H"]
    }
    list(shared = shared$H, specific = specific,
         shared_meta = shared[names(shared) != "H"], specific_meta = specific_meta)
}

.simulate_loadings <- function(p, q, density = 0.15, sd = 1, cfg = NULL, key = "loadings") {
    W <- .with_rng_stream(cfg, paste0(key, "/normal"), matrix(stats::rnorm(p * q, 0, sd), p, q))
    mask <- .with_rng_stream(cfg, paste0(key, "/mask"), matrix(stats::runif(p * q) < density, p, q))
    W <- W * mask
    empty <- which(rowSums(abs(W)) == 0)
    if (length(empty)) {
        for (i in empty) {
            if (is.null(cfg) || identical(cfg$rng_mode %||% "legacy", "legacy")) {
                W[i, sample.int(q, 1L)] <- stats::rnorm(1L, 0, sd)
                next
            }
            rescue <- .with_rng_stream(cfg, paste0(key, "/rescue/", i), c(sample.int(q, 1L), stats::rnorm(1L, 0, sd)))
            W[i, as.integer(rescue[1L])] <- rescue[2L]
        }
    }
    W
}
