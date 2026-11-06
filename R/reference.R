.reference_counts <- function(obj) {
    if (is.list(obj) && !inherits(obj, "Matrix") && !is.data.frame(obj) && !is.null(obj$counts)) {
        obj$counts
    } else {
        obj
    }
}

.reference_coordinates <- function(obj) {
    if (is.list(obj) && !inherits(obj, "Matrix") && !is.data.frame(obj) && !is.null(obj$coordinates)) {
        as.matrix(obj$coordinates)
    } else {
        NULL
    }
}

.reference_labels <- function(obj) {
    if (!is.list(obj) || inherits(obj, "Matrix") || is.data.frame(obj)) return(NULL)
    if (is.null(obj$domain) && is.null(obj$cell_type)) return(NULL)
    list(domain = obj$domain, cell_type = obj$cell_type)
}

.reference_spatial_scale <- function(obj, max_units = 500L) {
    xy <- .reference_coordinates(obj)
    x <- .reference_counts(obj)
    if (is.null(xy) || nrow(xy) != nrow(x) || nrow(x) < 10L) return(NA_real_)
    idx <- unique(round(seq(1, nrow(x), length.out = min(max_units, nrow(x)))))
    xy <- .normalize_coordinates(xy[idx, , drop = FALSE])
    emb <- .reference_embedding(x[idx, , drop = FALSE], q = 2L, max_cells = length(idx), max_features = 200L)
    if (is.null(emb) || ncol(emb) < 1L) return(NA_real_)
    d <- as.numeric(stats::dist(xy))
    g <- as.numeric(stats::dist(matrix(emb[, 1L], ncol = 1L)))^2 / 2
    ok <- is.finite(d) & is.finite(g) & d > 0
    d <- d[ok]; g <- g[ok]
    if (length(d) < 20L) return(NA_real_)
    br <- unique(stats::quantile(d, probs = seq(0, 1, length.out = 11L), names = FALSE))
    if (length(br) < 4L) return(NA_real_)
    bin <- cut(d, breaks = br, include.lowest = TRUE, labels = FALSE)
    dm <- tapply(d, bin, stats::median)
    gm <- tapply(g, bin, stats::median)
    sill <- stats::median(g[d >= stats::quantile(d, 0.8)], na.rm = TRUE)
    if (!is.finite(sill) || sill <= 0) return(NA_real_)
    hit <- which(gm >= 0.63 * sill)
    ell <- if (length(hit)) dm[min(hit)] else stats::median(d)
    .clip(as.numeric(ell), 0.03, 0.60)
}

.matrix_reference_stats <- function(x, modality) {
    if (!inherits(x, "Matrix")) x <- Matrix::Matrix(x, sparse = TRUE)
    n <- nrow(x); p <- ncol(x)
    lib <- as.numeric(Matrix::rowSums(x))
    mu <- as.numeric(Matrix::colMeans(x))
    ex2 <- as.numeric(Matrix::colMeans(x^2))
    va <- pmax(ex2 - mu^2, 0)
    nz <- as.numeric(Matrix::colSums(x != 0))
    zero <- 1 - nz / max(n, 1)
    disp <- mu^2 / pmax(va - mu, 1e-8)
    disp[!is.finite(disp) | disp <= 0] <- NA_real_
    list(
        modality = modality, n_units = n, n_features = p,
        library_median = stats::median(lib),
        library_log_sd = stats::sd(log1p(lib)),
        feature_mean_median = if (any(mu > 0)) stats::median(mu[mu > 0]) else 0,
        dispersion_median = if (any(is.finite(disp))) stats::median(disp[is.finite(disp)]) else NA_real_,
        zero_rate_median = stats::median(zero)
    )
}

.reference_embedding <- function(x, q = 5L, max_cells = 2000L, max_features = 500L) {
    if (!inherits(x, "Matrix")) x <- Matrix::Matrix(x, sparse = TRUE)
    nr <- nrow(x); nc <- ncol(x)
    ridx <- if (nr > max_cells) sample.int(nr, max_cells) else seq_len(nr)
    vars <- as.numeric(Matrix::colMeans(x^2)) - as.numeric(Matrix::colMeans(x))^2
    fidx <- order(vars, decreasing = TRUE)[seq_len(min(max_features, nc))]
    z <- as.matrix(x[ridx, fidx, drop = FALSE])
    lib <- rowSums(z); lib[lib <= 0] <- 1
    z <- log1p(z / lib * 1e4)
    rank_use <- min(q, nrow(z) - 1L, ncol(z))
    if (rank_use < 1L) return(NULL)
    stats::prcomp(z, center = TRUE, scale. = FALSE, rank. = rank_use)$x
}

.estimate_reference_rho <- function(reference, modalities) {
    if (length(modalities) < 2L) return(setNames(0.7, modalities))
    pair_scores <- setNames(vector("list", length(modalities)), modalities)
    for (i in seq_len(length(modalities) - 1L)) {
        for (j in (i + 1L):length(modalities)) {
            a <- modalities[[i]]; b <- modalities[[j]]
            xa <- .reference_counts(reference[[a]])
            xb <- .reference_counts(reference[[b]])
            # Operational calibration requires row-wise paired units.
            if (!is.null(rownames(xa)) || !is.null(rownames(xb))) {
                if (is.null(rownames(xa)) || is.null(rownames(xb))) next
                common <- intersect(rownames(xa), rownames(xb))
                if (length(common) < 3L) next
                xa <- xa[match(common, rownames(xa)), , drop = FALSE]
                xb <- xb[match(common, rownames(xb)), , drop = FALSE]
            }
            if (nrow(xa) != nrow(xb)) next
            nuse <- min(nrow(xa), 2000L)
            idx <- seq_len(nuse)
            ea <- .reference_embedding(xa[idx, , drop = FALSE], max_cells = nuse)
            eb <- .reference_embedding(xb[idx, , drop = FALSE], max_cells = nuse)
            if (is.null(ea) || is.null(eb)) next
            cc <- try(stats::cancor(ea, eb)$cor, silent = TRUE)
            if (inherits(cc, "try-error") || !length(cc)) next
            sc <- .clip(mean(cc^2), 0.05, 0.95)
            pair_scores[[a]] <- c(pair_scores[[a]], sc)
            pair_scores[[b]] <- c(pair_scores[[b]], sc)
        }
    }
    all_scores <- unlist(pair_scores, use.names = FALSE)
    fallback <- if (length(all_scores)) stats::median(all_scores, na.rm = TRUE) else 0.7
    if (!is.finite(fallback)) fallback <- 0.7
    out <- setNames(vapply(modalities, function(m) {
        x <- pair_scores[[m]]
        if (length(x)) stats::median(x, na.rm = TRUE) else fallback
    }, numeric(1)), modalities)
    .clip(out, 0.05, 0.95)
}

#' Calibrate operational simulation parameters from reference data
#'
#' @param reference Named list. Each modality is either a units-by-features matrix or a list with `counts` and optional `coordinates`.
#' @param estimate_shared Logical; estimate an operational shared-signal parameter when paired units are available.
#' @return A calibration list.
#' @export
calibrate_reference <- function(reference, estimate_shared = TRUE, atac_binary = TRUE,
                                calibration = c("moment", "legacy", "distribution_v2")) {
    calibration <- match.arg(calibration)
    if (calibration == "distribution_v2") return(fit_reference_profile(reference, atac_binary, estimate_shared))
    if (!is.list(reference) || is.null(names(reference))) stop("reference must be a named list of modalities.", call. = FALSE)
    reference <- .prepare_reference(reference, atac_binary = atac_binary)
    modalities <- intersect(names(reference), c("RNA", "ATAC", "Protein"))
    if (!length(modalities)) stop("Reference must contain RNA, ATAC and/or Protein.", call. = FALSE)
    modalities <- names(reference)
    stats <- lapply(modalities, function(m) {
        x <- .reference_counts(reference[[m]])
        out <- .matrix_reference_stats(x, m)
        out$spatial_length_scale <- .reference_spatial_scale(reference[[m]])
        out
    })
    names(stats) <- modalities
    rho <- if (estimate_shared) .estimate_reference_rho(reference, modalities) else setNames(rep(0.7, length(modalities)), modalities)
    pairing_diagnostics <- list()
    if (length(modalities) >= 2L) for (i in seq_len(length(modalities) - 1L)) for (j in (i + 1L):length(modalities)) {
        a <- .reference_counts(reference[[modalities[i]]]);b <- .reference_counts(reference[[modalities[j]]])
        named <- !is.null(rownames(a)) && !is.null(rownames(b))
        common <- if (named) length(intersect(rownames(a), rownames(b))) else if (is.null(rownames(a)) && is.null(rownames(b)) && nrow(a) == nrow(b)) nrow(a) else 0L
        mode <- if (named && common >= 3L) "ID_matched" else if (!named && common >= 3L) "row_order_assumed" else "unestablished"
        pairing_diagnostics[[paste(modalities[c(i, j)], collapse = "__")]] <- list(mode = mode, paired_units = common)
        if (estimate_shared && mode == "unestablished") warning("Reference pairing not established for ", paste(modalities[c(i, j)], collapse = "/"), "; operational shared_strength uses fallback.", call. = FALSE)
    }
    n_units <- vapply(stats, function(s) s$n_units, integer(1))
    spatial_scales <- vapply(stats, function(s) s$spatial_length_scale, numeric(1))
    spatial_scale <- stats::median(spatial_scales[is.finite(spatial_scales)], na.rm = TRUE)
    if (!is.finite(spatial_scale)) spatial_scale <- 0.22

    labels <- .reference_labels(reference[[modalities[1L]]])
    composition <- NULL; nk <- NULL; nc <- NULL
    if (!is.null(labels) && length(labels$domain) == stats[[modalities[1L]]]$n_units &&
        length(labels$cell_type) == stats[[modalities[1L]]]$n_units) {
        dz <- as.integer(as.factor(labels$domain))
        cz <- as.integer(as.factor(labels$cell_type))
        nk <- max(dz); nc <- max(cz)
        tab <- table(factor(dz, levels = seq_len(nk)), factor(cz, levels = seq_len(nc)))
        composition <- prop.table(tab + 1e-3, margin = 1L)
    }

    list(modalities = modalities, per_modality = stats, shared_strength = rho,
         shared_length_scale = spatial_scale, domain_celltype_composition = composition,
         n_domains = nk, n_cell_types = nc, suggested_n_cells = max(n_units),
         calibrated_at = as.character(Sys.time()), pairing_diagnostics = pairing_diagnostics,
         note = "Shared strength and spatial length scale are operational calibration parameters, not uniquely identifiable biological truth.")
}

.prepare_reference <- function(reference, atac_binary = TRUE) {
    if (!is.list(reference) || is.null(names(reference)) || anyDuplicated(names(reference))) stop("Reference must be a uniquely named list.", call. = FALSE)
    if ("ADT" %in% names(reference)) {
        if ("Protein" %in% names(reference)) stop("Use either ADT or Protein, not both.", call. = FALSE)
        names(reference)[names(reference) == "ADT"] <- "Protein"
    }
    if (any(!names(reference) %in% c("RNA", "ATAC", "Protein"))) stop("Unknown reference modality.", call. = FALSE)
    for (m in names(reference)) {
        v <- reference[[m]];x <- .reference_counts(v)
        if (length(dim(x)) != 2L || nrow(x) < 2L || ncol(x) < 2L) stop(m, " reference needs at least two units and features.", call. = FALSE)
        x <- Matrix::Matrix(x, sparse = TRUE)
        vals <- if (inherits(x, "sparseMatrix")) x@x else as.numeric(x)
        if (any(!is.finite(vals)) || any(vals < 0) || any(vals != round(vals))) stop(m, " reference must contain finite nonnegative integer counts; use the raw counts layer, not normalized X.", call. = FALSE)
        for (ids in list(rownames(x), colnames(x))) if (!is.null(ids) && (anyNA(ids) || anyDuplicated(ids) || any(!nzchar(ids)))) stop(m, " reference IDs must be unique and nonmissing.", call. = FALSE)
        if (!is.list(v)) v <- list(counts = x) else v$counts <- x
        if (!is.null(v$coordinates)) {
            if (nrow(as.matrix(v$coordinates)) != nrow(x)) stop(m, " coordinates must align with count rows.", call. = FALSE)
            .normalize_coordinates(v$coordinates)
        }
        for (nm in c("domain", "cell_type")) if (!is.null(v[[nm]]) && (length(v[[nm]]) != nrow(x) || anyNA(v[[nm]]))) stop(m, " ", nm, " must align with count rows.", call. = FALSE)
        if (m == "ATAC" && isTRUE(atac_binary)) v$counts <- Matrix::Matrix(1 * (x > 0), sparse = TRUE)
        reference[[m]] <- v
    }
    reference
}

#' Simulate from a real-data calibration
#'
#' @param reference Reference list accepted by [calibrate_reference()].
#' @param preserve_geometry Use coordinates from the first reference modality when available and dimensionally compatible.
#' @param preserve_labels Reuse domain and cell-type labels from the first reference modality when available.
#' @param ... Overrides passed to [heterospa_config()].
#' @return A `HeteroSpaSim` object.
#' @export
simulate_from_reference <- function(reference, preserve_geometry = TRUE, preserve_labels = TRUE,
                                    ..., calibration = c("moment", "legacy", "distribution_v2"),
                                    reference_profile = NULL, depth_model = c("empirical", "lognormal"),
                                    distribution_control = list()) {
    calibration <- match.arg(calibration)
    args <- list(...)
    if (calibration == "distribution_v2") {
        return(.simulate_distribution_reference(reference, preserve_geometry, preserve_labels,
            args, reference_profile, match.arg(depth_model), distribution_control))
    }
    binary <- args$atac_binary %||% TRUE
    reference <- .prepare_reference(reference, atac_binary = binary)
    fit <- calibrate_reference(reference, atac_binary = binary)
    mods <- fit$modalities
    nf <- setNames(vapply(fit$per_modality, function(x) x$n_features, integer(1)), mods)
    base <- list(
        n_cells = fit$suggested_n_cells,
        modalities = mods,
        n_features = nf,
        shared_strength = fit$shared_strength,
        shared_length_scale = fit$shared_length_scale,
        specific_length_scale = max(0.03, 0.75 * fit$shared_length_scale)
    )
    if (!is.null(fit$domain_celltype_composition)) {
        base$n_domains <- fit$n_domains
        base$n_cell_types <- fit$n_cell_types
        base$celltype_composition <- fit$domain_celltype_composition
    }
    if (("n_domains" %in% names(args) || "n_cell_types" %in% names(args)) &&
        !"celltype_composition" %in% names(args)) {
        base$celltype_composition <- NULL
    }
    if ("RNA" %in% mods) {
        s <- fit$per_modality$RNA
        base$rna_mean <- max(s$feature_mean_median, if (calibration == "legacy") 0.1 else 0.001)
        base$rna_dispersion <- ifelse(is.finite(s$dispersion_median), s$dispersion_median, 8)
        if (calibration == "moment") base$rna_dispersion <- .clip(base$rna_dispersion, 0.05, 100)
        base$rna_library_sd <- max(s$library_log_sd, 0.05)
    }
    if ("ATAC" %in% mods) {
        s <- fit$per_modality$ATAC
        base$atac_accessibility <- .clip(1 - s$zero_rate_median, 0.005, 0.5)
    }
    if ("Protein" %in% mods) {
        s <- fit$per_modality$Protein
        base$protein_mean <- max(s$feature_mean_median, if (calibration == "legacy") 0.1 else 0.001)
        base$protein_dispersion <- ifelse(is.finite(s$dispersion_median), s$dispersion_median, 12)
        if (calibration == "moment") base$protein_dispersion <- .clip(base$protein_dispersion, 0.05, 100)
    }
    first_coords <- .reference_coordinates(reference[[mods[1L]]])
    if (isTRUE(preserve_geometry) && !is.null(first_coords) && !"n_cells" %in% names(args)) {
        base$n_cells <- nrow(first_coords)
    }
    cfg <- do.call(heterospa_config, utils::modifyList(base, args))
    coords_use <- NULL
    if (isTRUE(preserve_geometry) && !is.null(first_coords)) {
        if (nrow(first_coords) == cfg$n_cells) coords_use <- first_coords
        else warning("Reference coordinates were not used because their row count does not match n_cells.", call. = FALSE)
    }
    labels_use <- if (isTRUE(preserve_labels)) .reference_labels(reference[[mods[1L]]]) else NULL
    if (!is.null(labels_use) && !is.null(labels_use$domain) && length(labels_use$domain) != cfg$n_cells) labels_use$domain <- NULL
    if (!is.null(labels_use) && !is.null(labels_use$cell_type) && length(labels_use$cell_type) != cfg$n_cells) labels_use$cell_type <- NULL
    run <- function(cfg) simulate_spatial_multiomics(
        config = cfg, coordinates = coords_use,
        domain_labels = if (is.null(labels_use)) NULL else labels_use$domain,
        cell_type_labels = if (is.null(labels_use)) NULL else labels_use$cell_type
    )
    sim <- run(cfg)
    history <- list()
    safe_design <- all(cfg$resolution <= 1.05) && cfg$pairing == 1 && cfg$overlap == 1 && all(cfg$missing_rate == 0) && (is.null(cfg$capture_fraction) || cfg$capture_fraction == 1)
    if (calibration == "moment" && !safe_design) warning("Moment correction skipped for heterogeneous measurement overrides; calibrate a matched baseline first, then vary measurement parameters.", call. = FALSE)
    for (attempt in 0:2) {
        history[[attempt + 1L]] <- list(config = sim$parameters, median_counts = vapply(sim$assays, function(x) stats::median(Matrix::rowSums(x)), numeric(1)))
        if (attempt >= 2L || calibration != "moment" || !safe_design) break
        for (m in mods) {
            param <- switch(m, RNA = "rna_mean", Protein = "protein_mean", ATAC = if (binary) "atac_accessibility" else "atac_depth")
            if (param %in% names(args)) next
            target <- fit$per_modality[[m]]$library_median
            if (target <= 0) next
            observed <- stats::median(Matrix::rowSums(sim$assays[[m]]))
            ratio <- .clip(target / max(observed, 1e-8), 0.1, 10)
            cfg[[param]] <- max(1e-5, cfg[[param]] * ratio)
            if (param == "atac_accessibility") cfg[[param]] <- .clip(cfg[[param]], 0.0001, 0.8)
        }
        sim <- run(cfg)
    }
    sim$reference_calibration <- fit
    sim$reference_calibration$mode <- calibration
    sim$reference_calibration$moment_history <- history
    sim$reference_calibration$input_preprocessing <- if (binary && "ATAC" %in% mods) "ATAC counts binarized before fitting" else "Count-scale reference"
    sim
}
