\name{heterospa_config}
\alias{heterospa_config}
\title{Construct a HeteroSpaSim configuration}
\usage{
heterospa_config (n_cells = 3000L, n_domains = 4L, n_cell_types = 6L, 
    celltype_composition = NULL, modalities = c("RNA", "ATAC", 
        "Protein"), n_features = c(RNA = 1000L, ATAC = 2000L, 
        Protein = 40L), geometry = "irregular", domain_pattern = "voronoi", 
    spatial_strength = 0.85, domain_purity = 0.65, shared_dim = 6L, 
    specific_dim = 4L, shared_strength = 0.7, shared_length_scale = 0.22, 
    specific_length_scale = 0.16, n_spatial_basis = 25L, factor_noise = 0.3, 
    signal_strength = 0.8, loading_density = 0.15, resolution = 1, 
    pairing = 1, capture_fraction = NULL, overlap = 1, misalignment = "none", 
    missing_rate = 0, missing_mechanism = "mcar", batch_count = 1L, 
    batch_strength = 0, batch_mode = "spatial", rna_mean = 1.5, 
    rna_dispersion = 8, rna_library_sd = 0.35, atac_binary = TRUE, 
    atac_accessibility = 0.08, atac_depth = 1, protein_mean = 2.5, 
    protein_dispersion = 12, protein_ambient = 0.05, chunk_size = 500L, 
    keep_cell_truth = FALSE, keep_complete_assays = TRUE, seed = 1L, 
    rng_mode = c("isolated", "legacy"))
}
\description{Creates and validates a complete simulation configuration. Scalar modality parameters are recycled to all selected modalities.}
\value{A validated list.}

\arguments{
\item{n_cells}{Number of latent cells used by a scenario preset.}
\item{n_domains}{Number of spatial domains.}
\item{n_cell_types}{Number of cell types.}
\item{celltype_composition}{Optional K-by-C domain-specific cell-type composition matrix.}
\item{modalities}{Modalities used by a scenario preset.}
\item{n_features}{Number of features per modality; scalar, unnamed vector or named vector.}
\item{geometry}{Tissue coordinates: irregular, grid or hex.}
\item{domain_pattern}{Domain partition: voronoi, layered, nested or gradient.}
\item{spatial_strength}{Domain boundary regularity between zero and one.}
\item{domain_purity}{Domain-specific cell-type mixing strength between zero and one.}
\item{shared_dim}{Shared latent dimension.}
\item{specific_dim}{Modality-specific latent dimensions.}
\item{shared_strength}{Effective shared-signal proportions in 0,1.}
\item{shared_length_scale}{Length scale of shared spatial factor basis on normalized coordinates.}
\item{specific_length_scale}{Length scale of modality-specific spatial factor basis.}
\item{n_spatial_basis}{Positive integer number of spatial basis centers.}
\item{factor_noise}{Nonnegative latent factor noise standard deviation.}
\item{signal_strength}{Nonnegative scalar molecular linear predictor strength.}
\item{loading_density}{Feature loading retention probability between zero and one; empty loading rows are rescued.}
\item{resolution}{Average latent cells per spatial measurement unit.}
\item{pairing}{Target exact pairing fraction for cell-resolution measurements.}
\item{capture_fraction}{Optional fixed fraction of eligible latent cells measured per modality; useful for pairing benchmarks without sample-size confounding.}
\item{overlap}{Target spatial coverage overlap in 0,1.}
\item{misalignment}{One of none, mild, moderate, severe or a parameter list.}
\item{missing_rate}{Missing measurement-unit rates per modality.}
\item{missing_mechanism}{Measurement-unit missingness: mcar, domain or spatial.}
\item{batch_count}{Number of technical batches.}
\item{batch_strength}{Batch-effect strength.}
\item{batch_mode}{Batch assignment: spatial or random.}
\item{rna_mean}{Positive RNA baseline count-scale multiplier.}
\item{rna_dispersion}{Positive RNA negative-binomial size; smaller values imply greater conditional variance.}
\item{rna_library_sd}{Nonnegative RNA cell library log-normal standard deviation.}
\item{atac_binary}{Logical; when true use binary accessibility and binarize reference counts before fitting. atac_depth has no effect on binary ATAC generation.}
\item{atac_accessibility}{ATAC baseline accessibility probability strictly between zero and one.}
\item{atac_depth}{Positive Poisson depth multiplier, active only for nonbinary ATAC.}
\item{protein_mean}{Positive protein baseline count-scale multiplier.}
\item{protein_dispersion}{Positive protein negative-binomial size.}
\item{protein_ambient}{Nonnegative additive Poisson background multiplier; it is not identifiable ambient truth in real data.}
\item{chunk_size}{Positive integer computation block size. Isolated streams make count results invariant to this value.}
\item{keep_cell_truth}{Logical; retain cell-level molecular counts.}
\item{keep_complete_assays}{Logical; retain assays before measurement-unit missingness.}
\item{seed}{Finite integer random seed; an explicit simulate_spatial_multiomics seed overrides the configuration seed.}
\item{rng_mode}{Named isolated random streams by default; legacy uses the sequential v0.1.7 random stream. Saved configurations without this field default to legacy.}
}
