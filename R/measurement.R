.coverage_masks <- function(coords, modalities, overlap) {
    M <- length(modalities)
    n <- nrow(coords)
    overlap <- .clip(as.numeric(overlap), 0, 1)
    if (M == 1L || overlap >= 0.999999) {
        out <- replicate(M, rep(TRUE, n), simplify = FALSE)
        names(out) <- modalities
        return(out)
    }

    # Construct a shared spatial core plus modality-specific, non-overlapping
    # peripheral sectors.  If each modality contains the common core C and one
    # unique sector U_m, pairwise Jaccard overlap is |C|/(|C|+|U_m|+|U_l|).
    # Choosing the core fraction below makes this equal to the requested
    # overlap (up to integer rounding), for any number of modalities.
    core_fraction <- (2 * overlap / M) / (1 - overlap + 2 * overlap / M)
    n_core <- max(0L, min(n, round(core_fraction * n)))

    cx <- stats::median(coords[, 1L]); cy <- stats::median(coords[, 2L])
    r2 <- (coords[, 1L] - cx)^2 + (coords[, 2L] - cy)^2
    core_idx <- if (n_core > 0L) order(r2)[seq_len(n_core)] else integer()
    remaining <- setdiff(seq_len(n), core_idx)

    # Angular ordering gives contiguous peripheral sectors rather than random
    # missing cells, preserving a spatially interpretable mosaic geometry.
    if (length(remaining)) {
        ang <- atan2(coords[remaining, 2L] - cy, coords[remaining, 1L] - cx)
        remaining <- remaining[order(ang)]
        grp <- cut(seq_along(remaining), breaks = M, labels = FALSE, include.lowest = TRUE)
    } else {
        grp <- integer()
    }

    out <- vector('list', M)
    for (i in seq_len(M)) {
        unique_idx <- if (length(remaining)) remaining[grp == i] else integer()
        mask <- rep(FALSE, n)
        mask[c(core_idx, unique_idx)] <- TRUE
        out[[i]] <- mask
    }
    names(out) <- modalities
    out
}

.paired_cell_subsets <- function(masks, modalities, pairing, capture_fraction = NULL) {
    base_pool <- which(masks[[1L]])
    auto_cap <- if (is.null(capture_fraction)) {
        if (pairing >= 0.999) 1 else min(0.95, 0.95 / (2 - pairing))
    } else {
        capture_fraction
    }
    nbase <- max(1L, floor(length(base_pool) * auto_cap))
    base_sel <- if (nbase >= length(base_pool)) base_pool else sample(base_pool, nbase, replace = FALSE)
    out <- list(); out[[modalities[1L]]] <- sort(base_sel)
    if (length(modalities) == 1L) return(out)

    for (m in modalities[-1L]) {
        pool <- which(masks[[m]])
        target_n <- max(1L, min(length(pool), floor(length(pool) * auto_cap)))
        common_pool <- intersect(base_sel, pool)
        npair_target <- min(round(pairing * min(length(base_sel), target_n)), length(common_pool))
        paired <- if (npair_target > 0L) sample(common_pool, npair_target, replace = FALSE) else integer()
        need_unique <- target_n - length(paired)
        unique <- .sample_without(pool, need_unique, exclude = base_sel)
        if (length(unique) < need_unique) {
            # Pairing target may be infeasible under strong mosaic non-overlap or high capture.
            extra <- .sample_without(pool, need_unique - length(unique), exclude = c(paired, unique))
            unique <- c(unique, extra)
        }
        sel <- sort(unique(c(paired, unique)))
        if (length(sel) > target_n) sel <- sel[seq_len(target_n)]
        out[[m]] <- sel
    }
    out
}

.build_operator <- function(coords, selected_cells, resolution, modality) {
    n <- nrow(coords)
    selected_cells <- sort(unique(selected_cells))
    if (!length(selected_cells)) stop("A modality has no cells after spatial coverage selection.", call. = FALSE)
    if (resolution <= 1.05) {
        nr <- length(selected_cells)
        A <- Matrix::sparseMatrix(i = seq_len(nr), j = selected_cells, x = 1,
                                  dims = c(nr, n))
        unit_id <- sprintf("%s_%06d", modality, seq_len(nr))
    } else {
        sub <- coords[selected_cells, , drop = FALSE]
        target_units <- max(1L, round(length(selected_cells) / resolution))
        g <- max(1L, ceiling(sqrt(target_units)))
        ix <- pmin(floor(sub[, 1L] * g) + 1L, g)
        iy <- pmin(floor(sub[, 2L] * g) + 1L, g)
        grp0 <- (iy - 1L) * g + ix
        lev <- sort(unique(grp0))
        grp <- match(grp0, lev)
        A <- Matrix::sparseMatrix(i = grp, j = selected_cells, x = 1,
                                  dims = c(length(lev), n))
        unit_id <- sprintf("%s_%06d", modality, seq_len(nrow(A)))
    }
    rownames(A) <- unit_id
    list(A = A, unit_id = unit_id)
}

.operator_centers <- function(A, coords) {
    rs <- as.numeric(Matrix::rowSums(A))
    rs[rs == 0] <- 1
    xy <- cbind(
        as.numeric(A %*% coords[, 1L]) / rs,
        as.numeric(A %*% coords[, 2L]) / rs
    )
    colnames(xy) <- c("x", "y")
    rownames(xy) <- rownames(A)
    xy
}

.resolve_misalignment <- function(x) {
    if (is.list(x)) {
        defaults <- list(rotation_deg = 0, translation = c(0, 0), scale = 1,
                         warp = 0, jitter = 0)
        return(utils::modifyList(defaults, x))
    }
    x <- match.arg(as.character(x), c("none", "mild", "moderate", "severe"))
    switch(x,
        none = list(rotation_deg = 0, translation = c(0, 0), scale = 1, warp = 0, jitter = 0),
        mild = list(rotation_deg = 5, translation = c(0.015, -0.01), scale = 1.01, warp = 0.01, jitter = 0.002),
        moderate = list(rotation_deg = 15, translation = c(0.04, -0.025), scale = 1.04, warp = 0.03, jitter = 0.004),
        severe = list(rotation_deg = 35, translation = c(0.10, -0.07), scale = 1.10, warp = 0.08, jitter = 0.008)
    )
}

.apply_transform <- function(xy, spec, modality_index = 1L) {
    if (modality_index == 1L) {
        truth <- list(rotation_deg = 0, translation = c(0, 0), scale = 1, warp = 0, jitter = 0)
        return(list(coords = xy, truth = truth))
    }
    p <- .resolve_misalignment(spec)
    sign_flip <- if (modality_index %% 2L == 0L) 1 else -1
    angle <- sign_flip * p$rotation_deg * pi / 180
    R <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), 2L, 2L, byrow = TRUE)
    centered <- sweep(xy, 2L, c(0.5, 0.5), "-")
    out <- p$scale * centered %*% t(R)
    out <- sweep(out, 2L, c(0.5, 0.5) + sign_flip * p$translation, "+")
    if (p$warp > 0) {
        warp <- cbind(sin(2 * pi * xy[, 2L]), sin(2 * pi * xy[, 1L])) * p$warp * sign_flip
        out <- out + warp
    }
    if (p$jitter > 0) out <- out + matrix(stats::rnorm(length(out), 0, p$jitter), nrow(out), 2L)
    colnames(out) <- c("x", "y"); rownames(out) <- rownames(xy)
    list(coords = out, truth = c(p, list(signed_rotation_deg = sign_flip * p$rotation_deg,
                                          signed_translation = sign_flip * p$translation)))
}

.correspondence_sparse <- function(A, B) {
    # A and B are measurement operators, usually sparse Matrix S4 objects.
    # Do not use base::t(B) here: on some R 4.4 / Matrix combinations the
    # transpose can fall through to t.default() during package execution.
    # tcrossprod() computes A %*% t(B) directly and preserves sparsity.
    if (!inherits(A, "Matrix")) A <- Matrix::Matrix(A, sparse = TRUE)
    if (!inherits(B, "Matrix")) B <- Matrix::Matrix(B, sparse = TRUE)
    if (ncol(A) != ncol(B)) {
        stop("Measurement operators must have the same latent-cell dimension.", call. = FALSE)
    }
    O <- Matrix::tcrossprod(A, B)
    # Matrix::summary() may return a sparseSummary object whose dim/nrow
    # behaviour varies across Matrix versions.  Test structural emptiness on
    # the sparse matrix itself instead of applying logical negation to nrow().
    if (Matrix::nnzero(O) == 0L) {
        return(Matrix::sparseMatrix(
            i = integer(), j = integer(), x = numeric(),
            dims = c(nrow(A), nrow(B)),
            dimnames = list(rownames(A), rownames(B))
        ))
    }
    # Avoid summary(O) entirely. Some Matrix builds on Windows/R 4.4 can
    # fail inside sparse summary methods. Convert to triplet representation
    # and read the structural slots directly.
    Ot <- methods::as(O, "TsparseMatrix")
    ii <- as.integer(Ot@i) + 1L
    jj <- as.integer(Ot@j) + 1L
    xx <- as.numeric(Ot@x)
    ra <- as.numeric(Matrix::rowSums(A))
    rb <- as.numeric(Matrix::rowSums(B))
    den <- ra[ii] + rb[jj] - xx
    val <- xx / pmax(den, 1e-12)
    Matrix::sparseMatrix(i = ii, j = jj, x = val, dims = c(nrow(A), nrow(B)),
                         dimnames = list(rownames(A), rownames(B)))
}

.build_measurements <- function(cfg, coords) {
    masks <- .coverage_masks(coords, cfg$modalities, cfg$overlap)
    subsets <- .with_rng_stream(cfg, "measurement/selection", .paired_cell_subsets(masks, cfg$modalities, cfg$pairing, cfg$capture_fraction))
    ops <- list(); complete_xy <- list(); observed_xy <- list(); transforms <- list()
    for (i in seq_along(cfg$modalities)) {
        m <- cfg$modalities[[i]]
        op <- .build_operator(coords, subsets[[m]], cfg$resolution[[m]], m)
        xy <- .operator_centers(op$A, coords)
        tr <- .with_rng_stream(cfg, paste0("measurement/transform/", m), .apply_transform(xy, cfg$misalignment, i))
        ops[[m]] <- op$A
        complete_xy[[m]] <- xy
        observed_xy[[m]] <- tr$coords
        transforms[[m]] <- tr$truth
    }
    corr <- list()
    if (length(cfg$modalities) >= 2L) {
        for (i in seq_len(length(cfg$modalities) - 1L)) {
            for (j in (i + 1L):length(cfg$modalities)) {
                a <- cfg$modalities[[i]]; b <- cfg$modalities[[j]]
                corr[[paste(a, b, sep = "__")]] <- .correspondence_sparse(ops[[a]], ops[[b]])
            }
        }
    }
    list(operators = ops, coordinates_complete = complete_xy,
         coordinates_transformed = observed_xy, coverage_masks = masks,
         selected_cells = subsets, transformations = transforms,
         correspondence = corr)
}

.measurement_annotations <- function(A, domains, celltypes, n_domains, n_cell_types, batches = NULL) {
    n <- ncol(A)
    dom_hot <- Matrix::sparseMatrix(i = seq_len(n), j = domains, x = 1, dims = c(n, n_domains))
    ct_hot <- Matrix::sparseMatrix(i = seq_len(n), j = celltypes, x = 1, dims = c(n, n_cell_types))
    dc <- as.matrix(A %*% dom_hot)
    cc <- as.matrix(A %*% ct_hot)
    rsd <- pmax(rowSums(dc), 1)
    rsc <- pmax(rowSums(cc), 1)
    dp <- dc / rsd
    cp <- cc / rsc
    colnames(dp) <- paste0("Domain", seq_len(n_domains))
    colnames(cp) <- paste0("CellType", seq_len(n_cell_types))
    rownames(dp) <- rownames(A); rownames(cp) <- rownames(A)
    out <- list(
        domain = max.col(dp, ties.method = "first"),
        cell_type = max.col(cp, ties.method = "first"),
        domain_composition = dp,
        celltype_composition = cp
    )
    if (!is.null(batches)) out$batch <- .row_mode_from_operator(A, batches)
    out
}
