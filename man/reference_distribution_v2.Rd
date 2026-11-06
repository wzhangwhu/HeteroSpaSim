\name{fit_reference_profile}
\alias{fit_reference_profile}
\alias{compare_reference_statistics}
\title{Distribution profiles and reference statistical diagnostics}
\usage{
fit_reference_profile(reference, atac_binary = TRUE, estimate_shared = TRUE)

compare_reference_statistics(reference, simulations, atac_binary = TRUE)
}
\arguments{
\item{reference}{Named RNA, ATAC and/or Protein list, each a units-by-features nonnegative integer matrix, or list with counts and optional coordinates, domain and cell_type. ADT is accepted as a Protein alias. IDs must be unique. Coordinate row IDs must agree when available.}
\item{atac_binary}{Treat ATAC as binary detected peaks; otherwise retain counts and Poisson generation.}
\item{estimate_shared}{Estimate operational shared strength only with ID-established pairs. Equal row counts alone do not establish pairing.}
\item{simulations}{Named list of simulated HeteroSpaSim objects with the same feature universe and ATAC mode as reference.}
}
\description{Fits a serializable profile without retaining reference expression matrices and computes descriptive statistical distances with explicit applicability and calibration roles.}
\details{Profiles retain empirical unit depths, feature mean, variance, detection and shrunken depth-offset NB2 moment estimates. Low detection features are shrunk to a robust center; dispersion is not a uniquely identified biological parameter. Input matrices are never copied as generation templates. Comparator reports exact empirical one-dimensional Wasserstein distances, depth and detected-feature quantiles and sample CV. Independent spatial diagnostics use deterministic subsets of at most 300 units and 60 features, library-normalized log1p counts, directed six-neighbor Moran I, and between-domain over total signal variance. All-zero features or absent spatial labels yield explicit missing/NA entries. Descriptive distances and three computational seeds do not establish general biological fidelity. Protein ambient is a configured assumption, not identified from reference data.}
\value{fit_reference_profile returns a HeteroSpaReferenceProfile list. compare_reference_statistics returns a long data frame with reference/simulated summaries, distances, effective sample sizes, status and calibration role.}
\examples{
set.seed(12)
ref <- list(RNA = matrix(rpois(40 * 20, 3), 40, 20))
profile <- fit_reference_profile(ref)
sim <- simulate_from_reference(ref, calibration = "distribution_v2",
    reference_profile = profile, seed = 42, n_spatial_basis = 4)
head(compare_reference_statistics(ref, list(new = sim)))
}
