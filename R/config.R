#' Construct a HeteroSpaSim configuration
#'
#' @param n_cells Number of latent tissue cells.
#' @param n_domains Number of spatial domains.
#' @param n_cell_types Number of cell types.
#' @param celltype_composition Optional K-by-C domain-specific cell-type composition matrix.
#' @param modalities Character vector containing any of RNA, ATAC and Protein.
#' @param n_features Number of features per modality; scalar, unnamed vector or named vector.
#' @param shared_dim Shared latent dimension.
#' @param specific_dim Modality-specific latent dimensions.
#' @param shared_strength Effective shared-signal proportions in [0,1].
#' @param resolution Average latent cells per spatial measurement unit.
#' @param pairing Target exact pairing fraction for cell-resolution measurements.
#' @param capture_fraction Optional fixed fraction of eligible latent cells measured per modality; useful for pairing benchmarks without sample-size confounding.
#' @param overlap Target spatial coverage overlap in [0,1].
#' @param misalignment One of none, mild, moderate, severe or a parameter list.
#' @param missing_rate Missing measurement-unit rates per modality.
#' @param batch_count Number of technical batches.
#' @param batch_strength Batch-effect strength.
#' @param seed Random seed.
#' @return A validated list of simulation parameters.
#' @export
heterospa_config <- function(
    n_cells = 3000L,
    n_domains = 4L,
    n_cell_types = 6L,
    celltype_composition = NULL,
    modalities = c("RNA", "ATAC", "Protein"),
    n_features = c(RNA = 1000L, ATAC = 2000L, Protein = 40L),
    geometry = "irregular",
    domain_pattern = "voronoi",
    spatial_strength = 0.85,
    domain_purity = 0.65,
    shared_dim = 6L,
    specific_dim = 4L,
    shared_strength = 0.70,
    shared_length_scale = 0.22,
    specific_length_scale = 0.16,
    n_spatial_basis = 25L,
    factor_noise = 0.30,
    signal_strength = 0.80,
    loading_density = 0.15,
    resolution = 1,
    pairing = 1,
    capture_fraction = NULL,
    overlap = 1,
    misalignment = "none",
    missing_rate = 0,
    missing_mechanism = "mcar",
    batch_count = 1L,
    batch_strength = 0,
    batch_mode = "spatial",
    rna_mean = 1.5,
    rna_dispersion = 8,
    rna_library_sd = 0.35,
    atac_binary = TRUE,
    atac_accessibility = 0.08,
    atac_depth = 1,
    protein_mean = 2.5,
    protein_dispersion = 12,
    protein_ambient = 0.05,
    chunk_size = 500L,
    keep_cell_truth = FALSE,
    keep_complete_assays = TRUE,
    seed = 1L,
    rng_mode = c("isolated", "legacy")
) {
    ints <- list(n_cells = n_cells, n_domains = n_domains, n_cell_types = n_cell_types, n_features = n_features, shared_dim = shared_dim, specific_dim = specific_dim, n_spatial_basis = n_spatial_basis, chunk_size = chunk_size, batch_count = batch_count, seed = seed)
    for (nm in names(ints)) {
        v <- ints[[nm]]
        if (!is.numeric(v) || !length(v) || any(!is.finite(v)) || any(v != round(v)) || any(abs(v) > .Machine$integer.max)) stop(nm, " must contain finite R integers.", call. = FALSE)
    }
    for (v in list(atac_binary, keep_cell_truth, keep_complete_assays)) if (!is.logical(v) || length(v) != 1L || is.na(v)) stop("Boolean configuration fields require TRUE or FALSE.", call. = FALSE)
    modalities <- unique(as.character(modalities))
    allowed <- c("RNA", "ATAC", "Protein")
    if (!length(modalities) || any(!modalities %in% allowed)) {
        stop("modalities must be drawn from RNA, ATAC and Protein.", call. = FALSE)
    }

    n_features <- .recycle_named(n_features, modalities,
        default = c(RNA = 1000L, ATAC = 2000L, Protein = 40L)[modalities])
    specific_dim <- .recycle_named(specific_dim, modalities, default = 4L)
    shared_strength <- .recycle_named(shared_strength, modalities, default = 0.7)
    resolution <- .recycle_named(resolution, modalities, default = 1)
    missing_rate <- .recycle_named(missing_rate, modalities, default = 0)

    cfg <- list(
        n_cells = as.integer(n_cells),
        n_domains = as.integer(n_domains),
        n_cell_types = as.integer(n_cell_types),
        celltype_composition = celltype_composition,
        modalities = modalities,
        n_features = as.integer(n_features),
        geometry = match.arg(geometry, c("irregular", "grid", "hex")),
        domain_pattern = match.arg(domain_pattern, c("voronoi", "layered", "nested", "gradient")),
        spatial_strength = as.numeric(spatial_strength),
        domain_purity = as.numeric(domain_purity),
        shared_dim = as.integer(shared_dim),
        specific_dim = as.integer(specific_dim),
        shared_strength = as.numeric(shared_strength),
        shared_length_scale = as.numeric(shared_length_scale),
        specific_length_scale = as.numeric(specific_length_scale),
        n_spatial_basis = as.integer(n_spatial_basis),
        factor_noise = as.numeric(factor_noise),
        signal_strength = as.numeric(signal_strength),
        loading_density = as.numeric(loading_density),
        resolution = as.numeric(resolution),
        pairing = as.numeric(pairing),
        capture_fraction = if (is.null(capture_fraction)) NULL else as.numeric(capture_fraction),
        overlap = as.numeric(overlap),
        misalignment = misalignment,
        missing_rate = as.numeric(missing_rate),
        missing_mechanism = match.arg(missing_mechanism, c("mcar", "domain", "spatial")),
        batch_count = as.integer(batch_count),
        batch_strength = as.numeric(batch_strength),
        batch_mode = match.arg(batch_mode, c("random", "spatial")),
        rna_mean = as.numeric(rna_mean),
        rna_dispersion = as.numeric(rna_dispersion),
        rna_library_sd = as.numeric(rna_library_sd),
        atac_binary = isTRUE(atac_binary),
        atac_accessibility = as.numeric(atac_accessibility),
        atac_depth = as.numeric(atac_depth),
        protein_mean = as.numeric(protein_mean),
        protein_dispersion = as.numeric(protein_dispersion),
        protein_ambient = as.numeric(protein_ambient),
        chunk_size = as.integer(chunk_size),
        keep_cell_truth = isTRUE(keep_cell_truth),
        keep_complete_assays = isTRUE(keep_complete_assays),
        seed = as.integer(seed),
        rng_mode = match.arg(rng_mode)
    )
    names(cfg$n_features) <- modalities
    names(cfg$specific_dim) <- modalities
    names(cfg$shared_strength) <- modalities
    names(cfg$resolution) <- modalities
    names(cfg$missing_rate) <- modalities
    validate_config(cfg)
}

validate_config <- function(cfg) {
    stopifnot(is.list(cfg))
    cfg$rng_mode <- cfg$rng_mode %||% "legacy"
    if (length(cfg$rng_mode) != 1L || !cfg$rng_mode %in% c("isolated", "legacy")) stop("Invalid rng_mode.", call. = FALSE)
    scalar <- c("n_cells", "n_domains", "n_cell_types", "spatial_strength", "domain_purity", "shared_dim", "shared_length_scale", "specific_length_scale", "n_spatial_basis", "factor_noise", "signal_strength", "loading_density", "pairing", "overlap", "batch_count", "batch_strength", "rna_mean", "rna_dispersion", "rna_library_sd", "atac_accessibility", "atac_depth", "protein_mean", "protein_dispersion", "protein_ambient", "chunk_size", "seed")
    for (nm in scalar) if (length(cfg[[nm]]) != 1L || !is.numeric(cfg[[nm]]) || !is.finite(cfg[[nm]])) stop(nm, " must be one finite numeric value.", call. = FALSE)
    if (!length(cfg$modalities) || anyNA(cfg$modalities) || anyDuplicated(cfg$modalities) || any(!cfg$modalities %in% c("RNA", "ATAC", "Protein"))) stop("Invalid modalities.", call. = FALSE)
    for (nm in c("n_features", "specific_dim", "shared_strength", "resolution", "missing_rate")) {
        v <- cfg[[nm]]
        if (!is.numeric(v) || !identical(names(v), cfg$modalities) || any(!is.finite(v))) stop(nm, " must have finite values named in modality order.", call. = FALSE)
    }
    for (nm in c("n_cells", "n_domains", "n_cell_types", "shared_dim", "n_spatial_basis", "chunk_size", "batch_count", "seed", "n_features", "specific_dim")) if (any(cfg[[nm]] != round(cfg[[nm]]))) stop(nm, " must contain integers.", call. = FALSE)
    for (nm in c("atac_binary", "keep_cell_truth", "keep_complete_assays")) if (!is.logical(cfg[[nm]]) || length(cfg[[nm]]) != 1L || is.na(cfg[[nm]])) stop(nm, " must be TRUE or FALSE.", call. = FALSE)
    opts <- list(geometry = c("irregular", "grid", "hex"), domain_pattern = c("voronoi", "layered", "nested", "gradient"), batch_mode = c("random", "spatial"), missing_mechanism = c("mcar", "domain", "spatial"))
    for (nm in names(opts)) if (length(cfg[[nm]]) != 1L || !is.character(cfg[[nm]]) || is.na(cfg[[nm]]) || !cfg[[nm]] %in% opts[[nm]]) stop("Invalid ", nm, call. = FALSE)
    if (!is.list(cfg$misalignment) && (length(cfg$misalignment) != 1L || !cfg$misalignment %in% c("none", "mild", "moderate", "severe"))) stop("Invalid misalignment.", call. = FALSE)
    if (cfg$n_cells < 20L) stop("n_cells must be at least 20.", call. = FALSE)
    if (cfg$n_domains < 1L || cfg$n_domains > cfg$n_cells) stop("Invalid n_domains.", call. = FALSE)
    if (cfg$n_cell_types < 1L || cfg$n_cell_types > cfg$n_cells) stop("Invalid n_cell_types.", call. = FALSE)
    if (!is.null(cfg$celltype_composition)) {
        cc <- as.matrix(cfg$celltype_composition)
        if (!all(dim(cc) == c(cfg$n_domains, cfg$n_cell_types))) {
            stop("celltype_composition must be n_domains by n_cell_types.", call. = FALSE)
        }
        if (any(!is.finite(cc)) || any(cc < 0) || any(rowSums(cc) <= 0)) {
            stop("celltype_composition must be finite, nonnegative, with positive row sums.", call. = FALSE)
        }
        cfg$celltype_composition <- cc / rowSums(cc)
    }
    if (cfg$shared_dim < 1L) stop("shared_dim must be positive.", call. = FALSE)
    if (cfg$shared_dim + max(cfg$specific_dim) >= cfg$n_cells) stop("Latent dimensions require shared_dim + specific_dim < n_cells.", call. = FALSE)
    if (any(cfg$specific_dim < 1L)) stop("specific_dim must be positive.", call. = FALSE)
    if (any(cfg$n_features < 2L)) stop("Each modality needs at least two features.", call. = FALSE)
    if (cfg$spatial_strength < 0 || cfg$spatial_strength > 1) stop("spatial_strength must lie in [0,1].", call. = FALSE)
    if (cfg$domain_purity < 0 || cfg$domain_purity > 1) stop("domain_purity must lie in [0,1].", call. = FALSE)
    if (any(cfg$shared_strength < 0 | cfg$shared_strength > 1)) stop("shared_strength must lie in [0,1].", call. = FALSE)
    if (any(cfg$resolution < 1)) stop("resolution must be >= 1.", call. = FALSE)
    if (cfg$pairing < 0 || cfg$pairing > 1) stop("pairing must lie in [0,1].", call. = FALSE)
    if (!is.null(cfg$capture_fraction) && (length(cfg$capture_fraction) != 1L || cfg$capture_fraction <= 0 || cfg$capture_fraction > 1)) {
        stop("capture_fraction must be NULL or a scalar in (0,1].", call. = FALSE)
    }
    if (cfg$overlap < 0 || cfg$overlap > 1) stop("overlap must lie in [0,1].", call. = FALSE)
    if (any(cfg$missing_rate < 0 | cfg$missing_rate >= 1)) stop("missing_rate must lie in [0,1).", call. = FALSE)
    if (cfg$loading_density <= 0 || cfg$loading_density > 1) stop("loading_density must lie in (0,1].", call. = FALSE)
    if (cfg$n_spatial_basis < 3L) stop("n_spatial_basis must be at least 3.", call. = FALSE)
    if (cfg$chunk_size < 10L) stop("chunk_size must be at least 10.", call. = FALSE)
    if (cfg$batch_count < 1L) stop("batch_count must be positive.", call. = FALSE)
    if (!is.finite(cfg$shared_length_scale) || cfg$shared_length_scale <= 0) stop("shared_length_scale must be positive.", call. = FALSE)
    if (!is.finite(cfg$specific_length_scale) || cfg$specific_length_scale <= 0) stop("specific_length_scale must be positive.", call. = FALSE)
    if (!is.finite(cfg$factor_noise) || cfg$factor_noise < 0) stop("factor_noise must be nonnegative.", call. = FALSE)
    if (!is.finite(cfg$signal_strength) || cfg$signal_strength < 0) stop("signal_strength must be nonnegative.", call. = FALSE)
    if (!is.finite(cfg$batch_strength) || cfg$batch_strength < 0) stop("batch_strength must be nonnegative.", call. = FALSE)
    if (!is.finite(cfg$rna_mean) || cfg$rna_mean <= 0) stop("rna_mean must be positive.", call. = FALSE)
    if (!is.finite(cfg$rna_dispersion) || cfg$rna_dispersion <= 0) stop("rna_dispersion must be positive.", call. = FALSE)
    if (!is.finite(cfg$rna_library_sd) || cfg$rna_library_sd < 0) stop("rna_library_sd must be nonnegative.", call. = FALSE)
    if (!is.finite(cfg$atac_accessibility) || cfg$atac_accessibility <= 0 || cfg$atac_accessibility >= 1) stop("atac_accessibility must lie in (0,1).", call. = FALSE)
    if (!is.finite(cfg$atac_depth) || cfg$atac_depth <= 0) stop("atac_depth must be positive.", call. = FALSE)
    if (!is.finite(cfg$protein_mean) || cfg$protein_mean <= 0) stop("protein_mean must be positive.", call. = FALSE)
    if (!is.finite(cfg$protein_dispersion) || cfg$protein_dispersion <= 0) stop("protein_dispersion must be positive.", call. = FALSE)
    if (!is.finite(cfg$protein_ambient) || cfg$protein_ambient < 0) stop("protein_ambient must be nonnegative.", call. = FALSE)
    if (length(cfg$seed) != 1L || is.na(cfg$seed)) stop("seed must be a single non-missing integer.", call. = FALSE)
    if (is.list(cfg$misalignment)) {
        req <- c("rotation_deg", "translation", "scale", "warp", "jitter")
        tmp <- utils::modifyList(list(rotation_deg = 0, translation = c(0, 0), scale = 1, warp = 0, jitter = 0), cfg$misalignment)
        if (length(tmp$translation) != 2L || any(!is.finite(unlist(tmp[req])))) stop("Invalid custom misalignment specification.", call. = FALSE)
        if (tmp$scale <= 0 || tmp$warp < 0 || tmp$jitter < 0) stop("misalignment scale must be positive and warp/jitter nonnegative.", call. = FALSE)
    }
    cfg
}

#' Standard scenario presets
#'
#' @param scenario One of easy, moderate or hard.
#' @param modalities Modalities to include.
#' @param n_cells Number of latent cells.
#' @param seed Random seed.
#' @return A HeteroSpaSim configuration list.
#' @export
scenario_preset <- function(scenario = c("easy", "moderate", "hard"),
                            modalities = c("RNA", "ATAC", "Protein"),
                            n_cells = 3000L, seed = 1L) {
    scenario <- match.arg(scenario)
    if (scenario == "easy") {
        heterospa_config(
            n_cells = n_cells, modalities = modalities,
            shared_strength = 0.85, resolution = 1,
            pairing = 1, overlap = 1, misalignment = "none",
            missing_rate = 0, batch_count = 1, batch_strength = 0,
            seed = seed
        )
    } else if (scenario == "moderate") {
        res <- rep(2, length(modalities)); names(res) <- modalities
        res[modalities[1L]] <- 1
        heterospa_config(
            n_cells = n_cells, modalities = modalities,
            shared_strength = 0.65, resolution = res,
            pairing = 0.65, overlap = 0.70, misalignment = "mild",
            missing_rate = 0.10, batch_count = 2, batch_strength = 0.15,
            seed = seed
        )
    } else {
        res <- rep(6, length(modalities)); names(res) <- modalities
        res[modalities[1L]] <- 1
        heterospa_config(
            n_cells = n_cells, modalities = modalities,
            shared_strength = 0.35, resolution = res,
            pairing = 0.30, overlap = 0.35, misalignment = "severe",
            missing_rate = 0.30, batch_count = 3, batch_strength = 0.35,
            spatial_strength = 0.70, factor_noise = 0.45,
            seed = seed
        )
    }
}
