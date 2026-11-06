#' Simulate heterogeneous spatial multi-omics data
#'
#' @param config Optional configuration from [heterospa_config()].
#' @param scenario Optional easy, moderate or hard preset used when config is NULL.
#' @param modalities Modalities used by a scenario preset.
#' @param n_cells Number of latent cells used by a scenario preset.
#' @param seed Random seed. Overrides the seed in config when supplied.
#' @param coordinates Optional n-by-2 tissue template coordinates. They are rescaled to [0,1]^2 internally.
#' @param domain_labels Optional template spatial-domain labels for the latent cells.
#' @param cell_type_labels Optional template cell-type labels for the latent cells.
#' @param ... Additional arguments forwarded to [heterospa_config()] when config is NULL and scenario is NULL.
#' @return An object of class `HeteroSpaSim`.
#' @export
simulate_spatial_multiomics <- function(config = NULL, scenario = NULL,
                                        modalities = c("RNA", "ATAC", "Protein"),
                                        n_cells = 3000L, seed = NULL, coordinates = NULL,
                                        domain_labels = NULL, cell_type_labels = NULL, ...) {
    if (is.null(config)) {
        if (!is.null(scenario)) {
            config <- scenario_preset(scenario, modalities = modalities, n_cells = n_cells,
                                      seed = seed %||% 1L)
        } else {
            config <- heterospa_config(n_cells = n_cells, modalities = modalities,
                                       seed = seed %||% 1L, ...)
        }
    } else {
        config <- validate_config(config)
        if (!is.null(seed)) {
            if (length(seed) != 1L || !is.numeric(seed) || !is.finite(seed) || seed != round(seed) || abs(seed) > .Machine$integer.max) stop("seed must be one finite R integer.", call. = FALSE)
            config$seed <- as.integer(seed)
        }
    }
    if (identical(config$rng_mode, "isolated")) {
        had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
        on.exit({if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv) else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)}, add = TRUE)
    }
    set.seed(config$seed)

    template_original <- NULL
    if (is.null(coordinates)) {
        coords <- .with_rng_stream(config, "tissue/coordinates", .simulate_coordinates(config$n_cells, config$geometry))
    } else {
        template_original <- as.matrix(coordinates)
        if (nrow(template_original) != config$n_cells) {
            stop("coordinates must have n_cells rows.", call. = FALSE)
        }
        coords <- .normalize_coordinates(template_original)
    }
    domain_levels <- NULL; celltype_levels <- NULL
    if (is.null(domain_labels)) {
        domains <- .with_rng_stream(config, "tissue/domains", .simulate_domains(coords, config$n_domains, config$domain_pattern, config$spatial_strength))
    } else {
        if (length(domain_labels) != config$n_cells || anyNA(domain_labels)) stop("domain_labels must contain one non-missing label per latent cell.", call. = FALSE)
        fd <- factor(domain_labels)
        domains <- as.integer(fd); domain_levels <- levels(fd)
        Knew <- nlevels(fd)
        if (!is.null(config$celltype_composition) && nrow(config$celltype_composition) != Knew) {
            stop("Template domain_labels are incompatible with celltype_composition rows.", call. = FALSE)
        }
        config$n_domains <- Knew
    }

    if (is.null(cell_type_labels)) {
        ct <- .with_rng_stream(config, "tissue/celltypes", .simulate_celltypes(domains, config$n_domains, config$n_cell_types, config$domain_purity,
                                  composition = config$celltype_composition))
        celltypes <- ct$cell_type
    } else {
        if (length(cell_type_labels) != config$n_cells || anyNA(cell_type_labels)) stop("cell_type_labels must contain one non-missing label per latent cell.", call. = FALSE)
        fc <- factor(cell_type_labels)
        celltypes <- as.integer(fc); celltype_levels <- levels(fc)
        Cnew <- nlevels(fc)
        if (!is.null(config$celltype_composition) && ncol(config$celltype_composition) != Cnew) {
            # When both labels are supplied, empirical composition supersedes the fitted matrix.
            if (!is.null(domain_labels)) config$celltype_composition <- NULL else stop("Template cell_type_labels are incompatible with celltype_composition columns.", call. = FALSE)
        }
        config$n_cell_types <- Cnew
        comp <- .composition_from_labels(domains, celltypes, config$n_domains, config$n_cell_types)
        ct <- list(cell_type = celltypes, composition = comp)
    }
    batch <- .with_rng_stream(config, "measurement/batches", .assign_batches(coords, config$batch_count, config$batch_mode))
    latent <- .simulate_latent_biology(config, coords, domains, celltypes)
    meas <- .build_measurements(config, coords)

    assays <- list(); truth_complete <- list(); feature_truth <- list(); masks <- list();
    spatial_observed <- list(); spatial_true_observed <- list(); spatial_complete <- meas$coordinates_transformed
    cell_truth <- list(); measurement_annotations <- list()

    for (m in config$modalities) {
        mol <- .simulate_cell_assay(
            m, config, latent$shared, latent$specific[[m]], batch,
            A = meas$operators[[m]], keep_cell = config$keep_cell_truth
        )
        complete <- mol$counts
        keep <- .with_rng_stream(config, paste0("measurement/missing/", m), .missing_unit_mask(config$missing_rate[[m]], config$missing_mechanism,
                                   meas$operators[[m]], meas$coordinates_transformed[[m]], domains))
        assays[[m]] <- complete[keep, , drop = FALSE]
        spatial_observed[[m]] <- meas$coordinates_transformed[[m]][keep, , drop = FALSE]
        spatial_true_observed[[m]] <- meas$coordinates_complete[[m]][keep, , drop = FALSE]
        masks[[m]] <- list(observed = keep, missing = !keep,
                           complete_unit_ids = rownames(complete), observed_unit_ids = rownames(complete)[keep])
        ann <- .measurement_annotations(
            meas$operators[[m]], domains, celltypes, config$n_domains, config$n_cell_types, batches = batch
        )
        ann$observed_domain <- ann$domain[keep]
        ann$observed_cell_type <- ann$cell_type[keep]
        if (!is.null(ann$batch)) ann$observed_batch <- ann$batch[keep]
        ann$observed_domain_composition <- ann$domain_composition[keep, , drop = FALSE]
        ann$observed_celltype_composition <- ann$celltype_composition[keep, , drop = FALSE]
        measurement_annotations[[m]] <- ann
        if (config$keep_complete_assays) truth_complete[[m]] <- complete
        if (config$keep_cell_truth) cell_truth[[m]] <- mol$cell_counts
        feature_truth[[m]] <- mol[!names(mol) %in% c("counts", "cell_counts")]
    }

    # Correspondence matrices exposed at the top level must match the actually
    # observed assays after modality-level missingness.  The complete
    # pre-missingness matrices are retained separately as ground truth.
    correspondence_complete <- meas$correspondence
    correspondence_observed <- list()
    if (length(config$modalities) >= 2L) {
        for (i in seq_len(length(config$modalities) - 1L)) {
            for (j in (i + 1L):length(config$modalities)) {
                a <- config$modalities[[i]]; b <- config$modalities[[j]]
                nm <- paste(a, b, sep = "__")
                C <- correspondence_complete[[nm]]
                correspondence_observed[[nm]] <- C[masks[[a]]$observed, masks[[b]]$observed, drop = FALSE]
            }
        }
    }

    cell_metadata <- data.frame(
        cell_id = sprintf("Cell_%06d", seq_len(config$n_cells)),
        domain = domains,
        cell_type = celltypes,
        batch = batch,
        x = coords[, 1L], y = coords[, 2L],
        stringsAsFactors = FALSE
    )
    rownames(cell_metadata) <- cell_metadata$cell_id

    obj <- list(
        assays = assays,
        spatial = list(latent_coordinates = coords,
                       template_coordinates_original = template_original,
                       true_measurement_coordinates = meas$coordinates_complete,
                       observed_true_coordinates = spatial_true_observed,
                       complete_coordinates = spatial_complete,
                       observed_coordinates = spatial_observed),
        metadata = list(
            cells = cell_metadata,
            domain = domains,
            cell_type = celltypes,
            batch = batch,
            domain_celltype_composition = ct$composition,
            domain_levels = domain_levels,
            cell_type_levels = celltype_levels
        ),
        latent = list(shared = latent$shared, specific = latent$specific,
                      shared_meta = latent$shared_meta, specific_meta = latent$specific_meta),
        measurement = list(operators = meas$operators, selected_cells = meas$selected_cells,
                           coverage_masks = meas$coverage_masks),
        correspondence = correspondence_observed,
        transformations = meas$transformations,
        masks = masks,
        truth = list(feature = feature_truth,
                     complete_assays = truth_complete,
                     cell_assays = cell_truth,
                     complete_correspondence = correspondence_complete,
                     measurement_annotations = measurement_annotations),
        parameters = config,
        provenance = list(package = "HeteroSpaSim", version = as.character(utils::packageVersion("HeteroSpaSim")),
                          seed = config$seed, generated_at = as.character(Sys.time()))
    )
    class(obj) <- "HeteroSpaSim"
    obj$truth$realized <- .realized_diagnostics(obj)
    obj$provenance$rng_mode <- config$rng_mode
    obj$provenance$rng_protocol <- if (config$rng_mode == "isolated") "isolated-v1" else "legacy-v0.1.7"
    obj$provenance$rng_kind <- RNGkind()
    if (length(config$modalities) > 1L && all(config$resolution <= 1.05)) {
        vals <- unlist(obj$truth$realized$exact_pairing)
        # Pairing is defined against the first modality, not between all secondary pairs.
        vals <- vals[startsWith(names(vals), paste0(config$modalities[1L], "__"))]
        if (any(abs(vals - config$pairing) > max(0.1, 2 / config$n_cells))) warning("Requested pairing is constrained by capture_fraction/overlap; inspect truth$realized$exact_pairing.", call. = FALSE)
    }
    validate_simulation(obj, stop_on_error = TRUE)
    obj
}

.realized_diagnostics <- function(obj) {
    mods <- obj$parameters$modalities
    resolution <- setNames(vapply(mods, function(m) {
        mean(as.numeric(Matrix::rowSums(obj$measurement$operators[[m]])))
    }, numeric(1)), mods)
    missing <- setNames(vapply(mods, function(m) mean(obj$masks[[m]]$missing), numeric(1)), mods)
    alignment_rmse <- setNames(vapply(mods, function(m) {
        a <- obj$spatial$true_measurement_coordinates[[m]]
        b <- obj$spatial$complete_coordinates[[m]]
        sqrt(mean(rowSums((a - b)^2)))
    }, numeric(1)), mods)
    overlap <- list(); pairing <- list()
    if (length(mods) >= 2L) {
        for (i in seq_len(length(mods) - 1L)) {
            for (j in (i + 1L):length(mods)) {
                a <- mods[[i]]; b <- mods[[j]]; nm <- paste(a, b, sep = "__")
                ma <- obj$measurement$coverage_masks[[a]]
                mb <- obj$measurement$coverage_masks[[b]]
                overlap[[nm]] <- sum(ma & mb) / max(sum(ma | mb), 1)
                C <- obj$truth$complete_correspondence[[nm]]
                # Avoid summary() on sparse Matrix objects for compatibility
                # with Windows/R 4.4 Matrix builds.
                Ct <- methods::as(C, "TsparseMatrix")
                cx <- as.numeric(Ct@x)
                nexact <- sum(cx >= 1 - 1e-10, na.rm = TRUE)
                pairing[[nm]] <- nexact / max(min(nrow(C), ncol(C)), 1)
            }
        }
    }
    list(resolution = resolution, missing_rate = missing, alignment_rmse = alignment_rmse,
         spatial_overlap = overlap, exact_pairing = pairing)
}

#' Validate a HeteroSpaSim object
#'
#' @param x A simulation object.
#' @param stop_on_error Stop if validation fails.
#' @return Invisibly TRUE on success, or a character vector of problems.
#' @export
validate_simulation <- function(x, stop_on_error = FALSE) {
    problems <- character()
    add <- function(ok, message) if (!isTRUE(ok)) problems <<- c(problems, message)
    tryCatch({
        add(inherits(x, "HeteroSpaSim"), "Object does not inherit from HeteroSpaSim.")
        req <- c("assays", "spatial", "metadata", "latent", "measurement", "parameters", "truth", "masks", "correspondence")
        add(all(req %in% names(x)), "Missing simulation fields.")
        if (length(problems)) stop("Incomplete object.")
        cfg <- validate_config(x$parameters);mods <- cfg$modalities
        cells <- x$metadata$cells
        add(nrow(cells) == cfg$n_cells && !anyNA(cells) && !anyDuplicated(cells$cell_id), "Invalid latent cell metadata.")
        add(all(cells$domain %in% seq_len(cfg$n_domains)) && all(cells$cell_type %in% seq_len(cfg$n_cell_types)), "Invalid latent labels.")
        add(identical(x$metadata$domain, cells$domain) && identical(x$metadata$cell_type, cells$cell_type) && identical(x$metadata$batch, cells$batch), "Metadata aliases disagree.")
        add(all(is.finite(x$spatial$latent_coordinates)) && nrow(x$spatial$latent_coordinates) == cfg$n_cells, "Invalid latent coordinates.")
        add(all(is.finite(x$latent$shared)) && nrow(x$latent$shared) == cfg$n_cells && ncol(x$latent$shared) == cfg$shared_dim, "Invalid shared latent matrix.")
        for (m in mods) {
            X <- x$assays[[m]];A <- x$measurement$operators[[m]];mask <- x$masks[[m]]$observed
            values <- if (inherits(X, "sparseMatrix")) X@x else as.numeric(X)
            ids <- rownames(X);feat <- colnames(X);xy <- x$spatial$observed_coordinates[[m]]
            add(length(dim(X)) == 2L && ncol(X) == cfg$n_features[[m]] && nrow(X) > 0L, paste(m, "invalid dimensions."))
            add(all(is.finite(values)) && all(values >= 0) && all(values == round(values)), paste(m, "counts must be finite nonnegative integers."))
            if (m == "ATAC" && cfg$atac_binary) add(all(values %in% c(0, 1)), "Binary ATAC contains non-binary entries.")
            add(length(ids) == nrow(X) && !anyNA(ids) && !anyDuplicated(ids) && all(nzchar(ids)), paste(m, "invalid spot IDs."))
            add(length(feat) == ncol(X) && !anyNA(feat) && !anyDuplicated(feat) && all(nzchar(feat)), paste(m, "invalid feature IDs."))
            add(identical(dim(xy), c(nrow(X), 2L)) && all(is.finite(xy)) && identical(rownames(xy), ids), paste(m, "coordinate IDs/dimensions mismatch."))
            add(identical(rownames(x$spatial$observed_true_coordinates[[m]]), ids) && all(is.finite(x$spatial$observed_true_coordinates[[m]])), paste(m, "true coordinate IDs mismatch."))
            av <- if (inherits(A, "sparseMatrix")) A@x else as.numeric(A)
            add(ncol(A) == cfg$n_cells && all(is.finite(av)) && all(av %in% c(0, 1)) && all(Matrix::rowSums(A) > 0) && all(Matrix::colSums(A) <= 1), paste(m, "invalid measurement operator."))
            add(is.logical(mask) && length(mask) == nrow(A) && !anyNA(mask) && sum(mask) == nrow(X), paste(m, "mask dimensions mismatch."))
            add(identical(ids, rownames(A)[mask]) && identical(ids, x$masks[[m]]$observed_unit_ids), paste(m, "observed ID/mask mismatch."))
            add(identical(x$masks[[m]]$complete_unit_ids, rownames(A)) && identical(x$masks[[m]]$missing, !mask), paste(m, "complete ID/missing mask mismatch."))
            ann <- x$truth$measurement_annotations[[m]]
            add(length(ann$domain) == nrow(A) && identical(ann$observed_domain, ann$domain[mask]) && identical(ann$observed_cell_type, ann$cell_type[mask]), paste(m, "observed truth alignment mismatch."))
            expected <- .measurement_annotations(A, cells$domain, cells$cell_type, cfg$n_domains, cfg$n_cell_types, batches = cells$batch)
            add(identical(ann$domain, expected$domain) && identical(ann$cell_type, expected$cell_type) && isTRUE(all.equal(ann$domain_composition, expected$domain_composition, tolerance = 1e-12)) && isTRUE(all.equal(ann$celltype_composition, expected$celltype_composition, tolerance = 1e-12)), paste(m, "operator-derived truth mismatch."))
            if (!is.null(x$truth$complete_assays[[m]])) add(identical(dim(x$truth$complete_assays[[m]]), c(nrow(A), ncol(X))) && identical(rownames(x$truth$complete_assays[[m]]), rownames(A)) && isTRUE(all.equal(x$truth$complete_assays[[m]][mask, , drop = FALSE], X, tolerance = 0)), paste(m, "complete/observed assays disagree."))
            if (!is.null(x$truth$cell_assays[[m]])) {
                cell <- x$truth$cell_assays[[m]];full <- A %*% cell
                if (m == "ATAC" && cfg$atac_binary) full <- Matrix::Matrix(1 * (full > 0), sparse = TRUE)
                add(isTRUE(all.equal(full[mask, , drop = FALSE], X, tolerance = 0, check.attributes = FALSE)), paste(m, "cell-to-observed aggregation mismatch."))
            }
            add(nrow(x$latent$specific[[m]]) == cfg$n_cells && ncol(x$latent$specific[[m]]) == cfg$specific_dim[[m]] && all(is.finite(x$latent$specific[[m]])), paste(m, "invalid specific latent matrix."))
        }
        if (length(mods) >= 2L) for (i in seq_len(length(mods) - 1L)) for (j in (i + 1L):length(mods)) {
            a <- mods[i];b <- mods[j];nm <- paste(a, b, sep = "__")
            full <- .correspondence_sparse(x$measurement$operators[[a]], x$measurement$operators[[b]])
            observed <- full[x$masks[[a]]$observed, x$masks[[b]]$observed, drop = FALSE]
            add(isTRUE(all.equal(full, x$truth$complete_correspondence[[nm]], tolerance = 1e-12)), paste(nm, "complete correspondence inconsistent with operator."))
            add(isTRUE(all.equal(observed, x$correspondence[[nm]], tolerance = 1e-12)), paste(nm, "observed correspondence/IDs mismatch."))
        }
    }, error = function(e) problems <<- c(problems, paste("Validation could not complete:", conditionMessage(e))))
    if (length(problems)) {
        problems <- unique(problems)
        if (stop_on_error) stop(paste(problems, collapse = "\n"), call. = FALSE)
        return(problems)
    }
    invisible(TRUE)
}

#' @export
print.HeteroSpaSim <- function(x, ...) {
    cat("HeteroSpaSim object
")
    cat("  latent cells: ", x$parameters$n_cells, "
", sep = "")
    cat("  domains / cell types: ", x$parameters$n_domains, " / ", x$parameters$n_cell_types, "
", sep = "")
    cat("  modalities: ", paste(x$parameters$modalities, collapse = ", "), "
", sep = "")
    for (m in x$parameters$modalities) {
        cat("    ", m, ": ", nrow(x$assays[[m]]), " units x ", ncol(x$assays[[m]]), " features",
            "; rho=", format(x$parameters$shared_strength[[m]], digits = 3),
            "; resolution=", format(x$parameters$resolution[[m]], digits = 3), "
", sep = "")
    }
    invisible(x)
}

#' @export
summary.HeteroSpaSim <- function(object, ...) {
    mods <- object$parameters$modalities
    data.frame(
        modality = mods,
        n_observed_units = vapply(mods, function(m) nrow(object$assays[[m]]), integer(1)),
        n_complete_units = vapply(mods, function(m) nrow(object$measurement$operators[[m]]), integer(1)),
        n_features = vapply(mods, function(m) ncol(object$assays[[m]]), integer(1)),
        shared_strength = as.numeric(object$parameters$shared_strength[mods]),
        resolution = as.numeric(object$parameters$resolution[mods]),
        missing_rate_target = as.numeric(object$parameters$missing_rate[mods]),
        missing_rate_realized = vapply(mods, function(m) mean(object$masks[[m]]$missing), numeric(1)),
        stringsAsFactors = FALSE
    )
}

#' @export
plot.HeteroSpaSim <- function(x, what = c("domains", "celltypes", "modalities"), ...) {
    what <- match.arg(what)
    if (what == "domains" || what == "celltypes") {
        md <- x$metadata$cells
        lab <- if (what == "domains") md$domain else md$cell_type
        graphics::plot(md$x, md$y, col = lab, pch = 16, asp = 1,
                       xlab = "x", ylab = "y", main = if (what == "domains") "Spatial domains" else "Cell types", ...)
    } else {
        mods <- x$parameters$modalities
        xyall <- do.call(rbind, x$spatial$observed_coordinates)
        xr <- range(xyall[, 1L]); yr <- range(xyall[, 2L])
        graphics::plot(NA, xlim = xr, ylim = yr, asp = 1, xlab = "x", ylab = "y", main = "Observed modality coordinates")
        for (i in seq_along(mods)) {
            xy <- x$spatial$observed_coordinates[[mods[[i]]]]
            graphics::points(xy[, 1L], xy[, 2L], pch = 16, cex = 0.45, col = i)
        }
        graphics::legend("topright", legend = mods, col = seq_along(mods), pch = 16, bty = "n")
    }
    invisible(x)
}
