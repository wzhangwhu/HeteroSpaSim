\name{simulate_spatial_multiomics}
\alias{simulate_spatial_multiomics}
\title{Simulate heterogeneous spatial multi-omics data}
\usage{
simulate_spatial_multiomics (config = NULL, scenario = NULL, modalities = c("RNA", 
    "ATAC", "Protein"), n_cells = 3000L, seed = NULL, coordinates = NULL, 
    domain_labels = NULL, cell_type_labels = NULL, ...)
}
\description{Generates latent tissue architecture, shared and modality-specific biological factors, molecular observations, heterogeneous spatial measurement operators and complete ground truth.}
\value{An object of class \code{HeteroSpaSim}.}

\arguments{
\item{config}{Optional configuration from heterospa_config().}
\item{scenario}{Optional easy, moderate or hard preset used when config is NULL.}
\item{modalities}{Modalities used by a scenario preset.}
\item{n_cells}{Number of latent cells used by a scenario preset.}
\item{seed}{Finite integer random seed; an explicit simulate_spatial_multiomics seed overrides the configuration seed.}
\item{coordinates}{Optional n-by-2 tissue template coordinates. They are rescaled to 0,1^2 internally.}
\item{domain_labels}{Optional template spatial-domain labels for the latent cells.}
\item{cell_type_labels}{Optional template cell-type labels for the latent cells.}
\item{...}{Additional configuration or method arguments as described for the relevant function.}
}
\examples{
cfg <- heterospa_config(n_cells = 40,
    modalities = c("RNA", "ATAC", "Protein"),
    n_features = c(RNA = 12, ATAC = 15, Protein = 5),
    shared_dim = 2, specific_dim = 2, seed = 42)
sim <- simulate_spatial_multiomics(cfg)
stopifnot(validate_simulation(sim))
summary(sim)
}
