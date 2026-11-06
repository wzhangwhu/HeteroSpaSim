# Write JSON when jsonlite is available; otherwise write an R-readable
# text fallback so the core simulator and benchmark export do not depend
# on a compiled optional package.
.write_metadata_sidecar <- function(object, json_file, fallback_file) {
    if (requireNamespace("jsonlite", quietly = TRUE)) {
        jsonlite::write_json(object, json_file, pretty = TRUE, auto_unbox = TRUE, digits = NA)
        return(invisible(json_file))
    }
    dput(object, file = fallback_file)
    warning(
        "Optional package 'jsonlite' is unavailable; wrote an R dput sidecar instead of JSON: ",
        fallback_file,
        call. = FALSE
    )
    invisible(fallback_file)
}

.write_matrix_bundle <- function(mat, dir_path, prefix = "counts") {
    .safe_dir(dir_path)
    file_mtx <- file.path(dir_path, paste0(prefix, ".mtx"))
    Matrix::writeMM(mat, file_mtx)
    utils::write.table(data.frame(feature = colnames(mat)), file.path(dir_path, "features.tsv"),
                       sep = "	", quote = FALSE, row.names = FALSE, col.names = FALSE)
    utils::write.table(data.frame(unit = rownames(mat)), file.path(dir_path, "units.tsv"),
                       sep = "	", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

#' Export a language-neutral benchmark directory
#'
#' @param x A HeteroSpaSim object.
#' @param path Destination directory.
#' @param include_complete Write complete pre-missingness assays when available.
#' @return Invisibly, the normalized destination path.
#' @export
export_benchmark <- function(x, path, include_complete = TRUE,
                             layout = c("isolated", "legacy"), design = c("V1", "V2")) {
    validate_simulation(x, stop_on_error = TRUE)
    layout <- match.arg(layout);design <- match.arg(design)
    if (dir.exists(path) && length(list.files(path, all.files = TRUE, no.. = TRUE))) stop("Destination is nonempty; use a new export directory.", call. = FALSE)
    path <- .safe_dir(path)
    if (layout == "isolated") {
        private <- file.path(path, "evaluator_truth")
        export_benchmark(x, private, include_complete, layout = "legacy", design = design)
        saveRDS(x$parameters, file.path(private, "config.rds"))
        saveRDS(x, file.path(private, "simulation.rds"))
        views <- list()
        for (m in x$parameters$modalities) {
            v <- .public_view(x, m, design);views[[m]] <- v
            md <- file.path(path, "public_inputs", m);.write_matrix_bundle(v$counts, md)
            utils::write.csv(data.frame(unit = v$ids, x = v$coordinates[, 1L], y = v$coordinates[, 2L], batch = v$batch), file.path(md, "coordinates.csv"), row.names = FALSE)
            utils::write.csv(v$mapping, file.path(private, paste0("public_unit_map_", m, ".csv")), row.names = FALSE)
        }
        if (design == "V1") for (nm in names(x$correspondence)) {
            mods <- strsplit(nm, "__", fixed = TRUE)[[1L]];co <- methods::as(x$correspondence[[nm]], "TsparseMatrix");ok <- co@x >= 1 - 1e-12
            pairs <- data.frame(source = views[[mods[1L]]]$ids[co@i[ok] + 1L], target = views[[mods[2L]]]$ids[co@j[ok] + 1L])
            pd <- .safe_dir(file.path(path, "public_inputs", "visible_pairs"));utils::write.csv(pairs, file.path(pd, paste0(nm, ".csv")), row.names = FALSE)
        }
        writeLines(c("HeteroSpaSim isolated export", paste("Design:", design), "Algorithm inputs: public_inputs only. Evaluation/config/full object: evaluator_truth only.", "V1 exposes declared exact pairs; V2 hides matching and stores ID/order maps privately.", "Matrices: units x features. Cross-resolution truth remains soft Jaccard, not obligatory one-to-one."), file.path(path, "README.txt"))
        return(invisible(path))
    }
    for (m in x$parameters$modalities) {
        md <- file.path(path, m)
        .write_matrix_bundle(x$assays[[m]], md)
        xy <- x$spatial$observed_coordinates[[m]]
        txy_obs <- x$spatial$observed_true_coordinates[[m]]
        utils::write.csv(data.frame(unit = rownames(xy), x = xy[, 1L], y = xy[, 2L]),
                         file.path(md, "coordinates.csv"), row.names = FALSE)
        utils::write.csv(data.frame(unit = rownames(txy_obs), x = txy_obs[, 1L], y = txy_obs[, 2L]),
                         file.path(md, "coordinates_true.csv"), row.names = FALSE)
        if (include_complete && length(x$truth$complete_assays) && !is.null(x$truth$complete_assays[[m]])) {
            cd <- file.path(md, "complete_truth")
            .write_matrix_bundle(x$truth$complete_assays[[m]], cd)
            cxy <- x$spatial$complete_coordinates[[m]]
            txy <- x$spatial$true_measurement_coordinates[[m]]
            utils::write.csv(data.frame(unit = rownames(cxy), x = cxy[, 1L], y = cxy[, 2L]),
                             file.path(cd, "coordinates_observed_complete.csv"), row.names = FALSE)
            utils::write.csv(data.frame(unit = rownames(txy), x = txy[, 1L], y = txy[, 2L]),
                             file.path(cd, "coordinates_true.csv"), row.names = FALSE)
        }
    }

    td <- .safe_dir(file.path(path, "truth"))
    utils::write.csv(x$metadata$cells, file.path(td, "latent_cells.csv"), row.names = FALSE)
    utils::write.csv(x$metadata$domain_celltype_composition,
                     file.path(td, "domain_celltype_composition.csv"), row.names = FALSE)
    for (m in x$parameters$modalities) {
        Matrix::writeMM(x$measurement$operators[[m]], file.path(td, paste0("measurement_operator_", m, ".mtx")))
        utils::write.table(data.frame(unit = rownames(x$measurement$operators[[m]])),
                           file.path(td, paste0("measurement_operator_", m, "_rows.tsv")),
                           sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
        utils::write.csv(data.frame(unit = x$masks[[m]]$complete_unit_ids,
                                    observed = x$masks[[m]]$observed),
                         file.path(td, paste0("missing_mask_", m, ".csv")), row.names = FALSE)
        utils::write.csv(x$truth$feature[[m]]$shared_loadings,
                         file.path(td, paste0("shared_loadings_", m, ".csv")), row.names = TRUE)
        utils::write.csv(x$truth$feature[[m]]$specific_loadings,
                         file.path(td, paste0("specific_loadings_", m, ".csv")), row.names = TRUE)
        ann <- x$truth$measurement_annotations[[m]]
        utils::write.csv(data.frame(
            unit = rownames(x$measurement$operators[[m]]),
            domain = ann$domain, cell_type = ann$cell_type,
            batch = if (is.null(ann$batch)) NA_integer_ else ann$batch
        ), file.path(td, paste0("measurement_labels_", m, ".csv")), row.names = FALSE)
        obs_ids <- x$masks[[m]]$observed_unit_ids
        utils::write.csv(data.frame(
            unit = obs_ids,
            domain = ann$observed_domain, cell_type = ann$observed_cell_type,
            batch = if (is.null(ann$observed_batch)) NA_integer_ else ann$observed_batch
        ), file.path(td, paste0("observed_measurement_labels_", m, ".csv")), row.names = FALSE)
        utils::write.csv(ann$domain_composition,
                         file.path(td, paste0("domain_composition_", m, ".csv")), row.names = TRUE)
        utils::write.csv(ann$celltype_composition,
                         file.path(td, paste0("celltype_composition_", m, ".csv")), row.names = TRUE)
    }
    for (nm in names(x$correspondence)) {
        C <- x$correspondence[[nm]]
        Matrix::writeMM(C, file.path(td, paste0("correspondence_", nm, ".mtx")))
        utils::write.table(data.frame(unit = rownames(C)), file.path(td, paste0("correspondence_", nm, "_rows.tsv")),
                           sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
        utils::write.table(data.frame(unit = colnames(C)), file.path(td, paste0("correspondence_", nm, "_cols.tsv")),
                           sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    }
    if (length(x$truth$complete_correspondence)) {
        for (nm in names(x$truth$complete_correspondence)) {
            C <- x$truth$complete_correspondence[[nm]]
            Matrix::writeMM(C, file.path(td, paste0("complete_correspondence_", nm, ".mtx")))
            utils::write.table(data.frame(unit = rownames(C)), file.path(td, paste0("complete_correspondence_", nm, "_rows.tsv")),
                               sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
            utils::write.table(data.frame(unit = colnames(C)), file.path(td, paste0("complete_correspondence_", nm, "_cols.tsv")),
                               sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
        }
    }
    utils::write.csv(x$latent$shared, file.path(td, "shared_factors.csv"), row.names = FALSE)
    for (m in names(x$latent$specific)) {
        utils::write.csv(x$latent$specific[[m]], file.path(td, paste0("specific_factors_", m, ".csv")), row.names = FALSE)
    }
    saveRDS(x$latent$shared, file.path(td, "shared_factors.rds"))
    saveRDS(x$latent$specific, file.path(td, "specific_factors.rds"))
    .write_metadata_sidecar(
        x$transformations,
        file.path(td, "spatial_transformations.json"),
        file.path(td, "spatial_transformations.dput")
    )
    .write_metadata_sidecar(
        x$parameters,
        file.path(path, "config.json"),
        file.path(path, "config.dput")
    )
    writeLines(c(
        "HeteroSpaSim benchmark export",
        "Rows of assay matrices are measurement units; columns are features.",
        "Matrix Market files use 1-based sparse matrix coordinates.",
        "The truth directory contains latent-cell labels, measurement operators, observed and complete correspondence, missingness and feature-loadings ground truth."
    ), file.path(path, "README.txt"))
    invisible(path)
}

#' Convert one modality to SpatialExperiment
#'
#' @param x A HeteroSpaSim object.
#' @param modality Modality name.
#' @return A SpatialExperiment object.
#' @export
as_SpatialExperiment <- function(x, modality, include_truth = TRUE) {
    validate_simulation(x, stop_on_error = TRUE)
    if (!requireNamespace("SpatialExperiment", quietly = TRUE) ||
        !requireNamespace("SummarizedExperiment", quietly = TRUE) ||
        !requireNamespace("SingleCellExperiment", quietly = TRUE) ||
        !requireNamespace("S4Vectors", quietly = TRUE)) {
        stop("Install SpatialExperiment, SingleCellExperiment, SummarizedExperiment and S4Vectors to use this function.", call. = FALSE)
    }
    modality <- match.arg(modality, x$parameters$modalities)
    counts <- Matrix::t(x$assays[[modality]])
    xy <- x$spatial$observed_coordinates[[modality]]
    ann <- x$truth$measurement_annotations[[modality]]
    cd <- S4Vectors::DataFrame(
        domain = ann$observed_domain,
        cell_type = ann$observed_cell_type,
        batch = if (is.null(ann$observed_batch)) NA_integer_ else ann$observed_batch,
        row.names = rownames(xy)
    )
    if (!include_truth) cd <- cd[, "batch", drop = FALSE]
    spe <- SpatialExperiment::SpatialExperiment(
        assays = list(counts = counts),
        colData = cd,
        spatialCoords = xy
    )
    SingleCellExperiment::reducedDim(spe, "spatial") <- xy
    spe
}

#' Convert all modalities to MultiAssayExperiment
#'
#' @param x A HeteroSpaSim object.
#' @return A MultiAssayExperiment containing SpatialExperiment assays.
#' @export
as_MultiAssayExperiment <- function(x, include_truth = TRUE) {
    if (!requireNamespace("MultiAssayExperiment", quietly = TRUE) ||
        !requireNamespace("S4Vectors", quietly = TRUE)) {
        stop("Install MultiAssayExperiment, SpatialExperiment and S4Vectors.", call. = FALSE)
    }
    exps <- lapply(x$parameters$modalities, function(m) as_SpatialExperiment(x, m, include_truth = include_truth))
    names(exps) <- x$parameters$modalities
    sm <- do.call(rbind, lapply(names(exps), function(m) {
        cn <- colnames(exps[[m]])
        data.frame(assay = m, primary = paste(m, cn, sep = "::"), colname = cn,
                   stringsAsFactors = FALSE)
    }))
    prim <- unique(sm$primary)
    cd <- S4Vectors::DataFrame(row.names = prim)
    MultiAssayExperiment::MultiAssayExperiment(
        experiments = exps,
        colData = cd,
        sampleMap = S4Vectors::DataFrame(sm)
    )
}

#' Export each modality as H5AD
#'
#' @param x A HeteroSpaSim object.
#' @param path Destination directory.
#' @return Invisibly, paths to generated H5AD files.
#' @export
export_h5ad <- function(x, path, include_truth = FALSE, design = c("V1", "V2"),
                        backend = c("auto", "python", "zellkonverter")) {
    validate_simulation(x, stop_on_error = TRUE);design <- match.arg(design);backend <- match.arg(backend)
    if (backend == "auto") backend <- if (requireNamespace("reticulate", quietly = TRUE) && reticulate::py_available(initialize = FALSE)) "python" else "zellkonverter"
    if (backend == "zellkonverter" && !requireNamespace("zellkonverter", quietly = TRUE)) {
        stop("Install zellkonverter to export H5AD files.", call. = FALSE)
    }
    path <- .safe_dir(path)
    private <- .safe_dir(file.path(path, "evaluator_truth"))
    public <- .safe_dir(if (include_truth) file.path(private, "annotated_h5ad") else file.path(path, "public_inputs"))
    out <- character()
    for (m in x$parameters$modalities) {
        view <- .public_view(x, m, design)
        f <- file.path(public, paste0(m, ".h5ad"));if (file.exists(f)) stop("Refusing to overwrite H5AD: ", f, call. = FALSE)
        utils::write.csv(view$mapping, file.path(private, paste0(m, "_unit_map.csv")), row.names = FALSE)
        if (backend == "python") {
            if (!requireNamespace("reticulate", quietly = TRUE)) stop("Install reticulate and Python anndata/scipy/pandas.", call. = FALSE)
            ad <- reticulate::import("anndata", convert = FALSE);pd <- reticulate::import("pandas", convert = FALSE)
            obs <- pd$DataFrame(reticulate::dict(batch = as.integer(view$batch)), index = view$ids)
            if (include_truth) {
                obs$`__setitem__`("domain", reticulate::r_to_py(x$truth$measurement_annotations[[m]]$observed_domain[view$index]))
                obs$`__setitem__`("cell_type", reticulate::r_to_py(x$truth$measurement_annotations[[m]]$observed_cell_type[view$index]))
            }
            var <- pd$DataFrame(reticulate::dict(feature_name = colnames(view$counts)), index = colnames(view$counts))
            a <- ad$AnnData(X = reticulate::r_to_py(methods::as(view$counts, "dgCMatrix")), obs = obs, var = var)
            # Python mapping proxies must be mutated explicitly. Nested R [[<-
            # assignment can silently modify a converted temporary instead.
            a$obsm$`__setitem__`("spatial", reticulate::r_to_py(unname(view$coordinates)))
            a$layers$`__setitem__`("counts", a$X$copy())
            a$uns$`__setitem__`("heterospa_design", design)
            a$uns$`__setitem__`("include_truth", include_truth)
            a$write_h5ad(f)
        } else {
            spe <- as_SpatialExperiment(x, m, include_truth = include_truth)
            spe <- spe[, view$index];colnames(spe) <- view$ids
            # Keep an explicit counts layer while selecting an identical X assay.
            SummarizedExperiment::assay(spe, "X") <- SummarizedExperiment::assay(spe, "counts")
            S4Vectors::metadata(spe) <- list(heterospa_design = design, include_truth = include_truth)
            zellkonverter::writeH5AD(spe, f, X_name = "X")
        }
        out[m] <- f
    }
    invisible(out)
}

#' Export a MuData H5MU file
#'
#' @param x A HeteroSpaSim object.
#' @param file Destination `.h5mu` file.
#' @return Invisibly, the destination path.
#' @export
export_h5mu <- function(x, file, include_truth = FALSE, design = c("V1", "V2"),
                        backend = c("auto", "python", "zellkonverter")) {
    design <- match.arg(design);backend <- match.arg(backend)
    if (file.exists(file)) stop("Refusing to overwrite H5MU.", call. = FALSE)
    .safe_dir(dirname(file))
    if (!requireNamespace("reticulate", quietly = TRUE)) {
        stop("Install reticulate; Python packages anndata and mudata and the selected H5AD backend are also required.", call. = FALSE)
    }
    tmp <- tempfile("heterospa_h5ad_")
    dir.create(tmp)
    on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
    paths <- export_h5ad(x, tmp, include_truth = include_truth, design = design, backend = backend)
    ad <- reticulate::import("anndata", delay_load = FALSE, convert = FALSE)
    mudata <- reticulate::import("mudata", delay_load = FALSE, convert = FALSE)
    mods <- list()
    for (m in names(paths)) mods[[m]] <- ad$read_h5ad(paths[[m]])
    md <- mudata$MuData(mods)
    md$write_h5mu(normalizePath(file, winslash = "/", mustWork = FALSE))
    private <- .safe_dir(file.path(dirname(file), paste0(basename(file), "_evaluator_truth")))
    for (p in list.files(file.path(tmp, "evaluator_truth"), full.names = TRUE)) file.copy(p, private, overwrite = FALSE)
    invisible(file)
}

#' Export major matrices to MATLAB MAT format
#'
#' @param x A HeteroSpaSim object.
#' @param file Destination `.mat` file.
#' @return Invisibly, the destination path.
#' @export
export_mat <- function(x, file, include_truth = FALSE, design = c("V1", "V2")) {
    validate_simulation(x, stop_on_error = TRUE);design <- match.arg(design)
    if (file.exists(file)) stop("Refusing to overwrite MAT.", call. = FALSE)
    .safe_dir(dirname(file))
    if (!requireNamespace("R.matlab", quietly = TRUE)) {
        stop("Install R.matlab to export MAT files.", call. = FALSE)
    }
    args <- list()
    for (m in x$parameters$modalities) {
        if (prod(dim(x$assays[[m]])) > 5e7) {
            stop("MAT export would densify a very large assay. Use export_benchmark() or H5AD/H5MU for large simulations.", call. = FALSE)
        }
        v <- .public_view(x, m, design)
        args[[paste0("X_", m)]] <- as.matrix(v$counts)
        args[[paste0("S_", m)]] <- v$coordinates
        args[[paste0("units_", m)]] <- v$ids
        args[[paste0("features_", m)]] <- colnames(v$counts)
        if (design == "V2") {
            pd <- .safe_dir(file.path(dirname(file), paste0(basename(file), "_evaluator_truth")))
            utils::write.csv(v$mapping, file.path(pd, paste0(m, "_unit_map.csv")), row.names = FALSE)
        }
    }
    if (include_truth) {
    args$latent_coordinates <- x$spatial$latent_coordinates
    args$domain <- x$metadata$cells$domain
    args$cell_type <- x$metadata$cells$cell_type
    args$shared_factors <- x$latent$shared
    }
    do.call(R.matlab::writeMat, c(list(con = file), args))
    invisible(file)
}

.public_view <- function(x, modality, design) {
    X <- x$assays[[modality]];original <- rownames(X);n <- nrow(X)
    cfg <- x$parameters;cfg$rng_mode <- "isolated"
    index <- if (design == "V2") .with_rng_stream(cfg, paste0("export/V2/", modality), sample.int(n)) else seq_len(n)
    ids <- if (design == "V2") sprintf(paste0(modality, "_public_%06d"), seq_len(n)) else original
    if (design == "V1" && x$parameters$resolution[[modality]] <= 1.05) {
        codes <- .with_rng_stream(cfg, "export/V1/barcodes", sample.int(cfg$n_cells))
        A <- methods::as(x$measurement$operators[[modality]][x$masks[[modality]]$observed, , drop = FALSE], "TsparseMatrix")
        cell <- integer(n);cell[A@i + 1L] <- A@j + 1L;ids <- sprintf("Spot_%06d", codes[cell])
    }
    counts <- X[index, , drop = FALSE];rownames(counts) <- ids
    xy <- x$spatial$observed_coordinates[[modality]][index, , drop = FALSE];rownames(xy) <- ids
    list(counts = counts, coordinates = xy, ids = ids, index = index, batch = x$truth$measurement_annotations[[modality]]$observed_batch[index], mapping = data.frame(public_unit_id = ids, original_unit_id = original[index], original_row = index))
}
