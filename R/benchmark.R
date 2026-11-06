#' Built-in benchmark definitions
#'
#' @return A named list describing the six standard benchmark axes.
#' @export
benchmark_presets <- function() {
    list(
        B1_shared_strength = list(parameter = "shared_strength", values = c(0.2, 0.4, 0.6, 0.8, 1.0),
                                  description = "Cross-modal shared biological signal"),
        B2_partial_pairing = list(parameter = "pairing", values = c(1.0, 0.75, 0.50, 0.25, 0.0),
                                  description = "Exact pairing fraction"),
        B3_mosaic_overlap = list(parameter = "overlap", values = c(1.0, 0.75, 0.50, 0.25, 0.0),
                                 description = "Spatial coverage overlap"),
        B4_resolution_mismatch = list(parameter = "resolution_secondary", values = c(1, 2, 4, 8, 16),
                                      description = "Cells per measurement unit in non-reference modalities"),
        B5_misalignment = list(parameter = "misalignment", values = c("none", "mild", "moderate", "severe"),
                               description = "Spatial registration difficulty"),
        B6_technical_noise = list(parameter = "technical_severity", values = c("clean", "low", "moderate", "severe"),
                                  description = "Combined technical-noise severity")
    )
}

.apply_benchmark_value <- function(cfg, parameter, value) {
    if (parameter == "shared_strength") {
        cfg$shared_strength[] <- as.numeric(value)
    } else if (parameter == "pairing") {
        cfg$pairing <- as.numeric(value)
    } else if (parameter == "overlap") {
        cfg$overlap <- as.numeric(value)
    } else if (parameter == "resolution_secondary") {
        cfg$resolution[] <- as.numeric(value)
        cfg$resolution[[cfg$modalities[1L]]] <- 1
    } else if (parameter == "misalignment") {
        cfg$misalignment <- as.character(value)
    } else if (parameter == "technical_severity") {
        v <- as.character(value)
        if (v == "clean") {
            cfg$missing_rate[] <- 0; cfg$batch_count <- 1L; cfg$batch_strength <- 0
            cfg$rna_library_sd <- 0.20; cfg$rna_dispersion <- 10; cfg$atac_accessibility <- 0.10; cfg$protein_dispersion <- 14; cfg$protein_ambient <- 0.01
        } else if (v == "low") {
            cfg$missing_rate[] <- 0.05; cfg$batch_count <- 2L; cfg$batch_strength <- 0.10
            cfg$rna_library_sd <- 0.30; cfg$rna_dispersion <- 8; cfg$atac_accessibility <- 0.08; cfg$protein_dispersion <- 12; cfg$protein_ambient <- 0.03
        } else if (v == "moderate") {
            cfg$missing_rate[] <- 0.20; cfg$batch_count <- 2L; cfg$batch_strength <- 0.25
            cfg$rna_library_sd <- 0.45; cfg$rna_dispersion <- 5; cfg$atac_accessibility <- 0.055; cfg$protein_dispersion <- 8; cfg$protein_ambient <- 0.08
        } else {
            cfg$missing_rate[] <- 0.40; cfg$batch_count <- 3L; cfg$batch_strength <- 0.45
            cfg$rna_library_sd <- 0.65; cfg$rna_dispersion <- 3; cfg$atac_accessibility <- 0.035; cfg$protein_dispersion <- 5; cfg$protein_ambient <- 0.18
        }
    }
    validate_config(cfg)
}

.benchmark_base_config <- function(benchmark, modalities, n_cells, seed) {
    cfg <- heterospa_config(
        n_cells = n_cells, modalities = modalities,
        shared_strength = 0.70, resolution = 1, pairing = 1, capture_fraction = NULL,
        overlap = 1, misalignment = "none", missing_rate = 0,
        batch_count = 1, batch_strength = 0, factor_noise = 0.30, seed = seed
    )
    if (benchmark == "B2_partial_pairing") cfg$capture_fraction <- 0.45
    cfg
}

#' Generate a reproducible benchmark suite
#'
#' @param benchmark Name returned by [benchmark_presets()].
#' @param n_replicates Number of independent seeds per parameter setting.
#' @param n_cells Number of latent tissue cells.
#' @param modalities Modalities to simulate.
#' @param base_config Optional baseline configuration.
#' @param out_dir Optional directory. When supplied, each simulation is exported.
#' @param seed_start First random seed.
#' @return A manifest data frame, invisibly if out_dir is supplied; otherwise a list with manifest and simulations.
#' @export
generate_benchmark <- function(benchmark, n_replicates = 10L, n_cells = 3000L,
                               modalities = c("RNA", "ATAC", "Protein"),
                               base_config = NULL, out_dir = NULL, seed_start = 1L) {
    defs <- benchmark_presets()
    if (!benchmark %in% names(defs)) stop("Unknown benchmark: ", benchmark, call. = FALSE)
    def <- defs[[benchmark]]
    cfg0 <- base_config %||% .benchmark_base_config(benchmark, modalities, n_cells, seed_start)
    if (benchmark %in% c("B1_shared_strength", "B2_partial_pairing", "B3_mosaic_overlap", "B4_resolution_mismatch", "B5_misalignment") &&
        length(cfg0$modalities) < 2L) {
        stop(benchmark, " requires at least two modalities.", call. = FALSE)
    }
    manifest <- list(); sims <- list(); k <- 0L
    if (!is.null(out_dir)) .safe_dir(out_dir)

    for (v in def$values) {
        for (r in seq_len(n_replicates)) {
            k <- k + 1L
            cfg <- .apply_benchmark_value(cfg0, def$parameter, v)
            # Use matched random seeds across parameter levels: replicate r
            # therefore shares the same underlying random stream for every
            # difficulty value, reducing nuisance variation in benchmark curves.
            cfg$seed <- as.integer(seed_start + r - 1L)
            sim <- simulate_spatial_multiomics(config = cfg)
            value_tag <- gsub("[^A-Za-z0-9._-]", "_", as.character(v))
            id <- sprintf("%s_%s_rep%02d", benchmark, value_tag, r)
            if (!is.null(out_dir)) {
                export_benchmark(sim, file.path(out_dir, id))
            } else {
                sims[[id]] <- sim
            }
            manifest[[k]] <- data.frame(
                id = id, benchmark = benchmark, parameter = def$parameter,
                value = as.character(v), replicate = r, seed = cfg$seed,
                stringsAsFactors = FALSE
            )
        }
    }
    man <- do.call(rbind, manifest)
    if (!is.null(out_dir)) {
        utils::write.csv(man, file.path(out_dir, "manifest.csv"), row.names = FALSE)
        return(invisible(man))
    }
    list(manifest = man, simulations = sims)
}
