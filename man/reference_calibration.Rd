\name{calibrate_reference}
\alias{calibrate_reference}
\alias{simulate_from_reference}
\title{Reference-calibrated spatial multi-omics simulation}
\usage{
calibrate_reference (reference, estimate_shared = TRUE, atac_binary = TRUE,
    calibration = c("moment", "legacy", "distribution_v2"))

simulate_from_reference (reference, preserve_geometry = TRUE, preserve_labels = TRUE, 
    ..., calibration = c("moment", "legacy", "distribution_v2"),
    reference_profile = NULL, depth_model = c("empirical", "lognormal"),
    distribution_control = list())
}
\description{Estimates operational simulation parameters from real units-by-features matrices and uses them to initialize HeteroSpaSim. Shared strength is an operational calibration quantity rather than uniquely identifiable biological truth.}
\value{A calibration list or a HeteroSpaSim object.}

\arguments{
\item{reference}{Reference list accepted by calibrate_reference().}
\item{estimate_shared}{Logical; estimate an operational shared-signal parameter when paired units are available.}
\item{atac_binary}{Logical; when true use binary accessibility and binarize reference counts before fitting. atac_depth has no effect on binary ATAC generation.}
\item{preserve_geometry}{Use coordinates from the first reference modality when available and dimensionally compatible.}
\item{preserve_labels}{Reuse domain and cell-type labels from the first reference modality when available.}
\item{...}{Additional configuration or method arguments as described for the relevant function.}
\item{calibration}{moment performs two bounded median-count corrections for matched baselines; legacy omits these corrections. Explicit distribution_v2 uses a reusable distribution profile and pre-sampling measurement-depth calibration. Defaults and old streams are unchanged. No mode guarantees biological fidelity.}
\item{reference_profile}{Optional HeteroSpaReferenceProfile from fit_reference_profile, fitted to the intended feature universe; no real counts template is stored.}
\item{depth_model}{Empirical measurement-depth sampling (default for distribution_v2), or positive lognormal with a zero mass. Binary ATAC uses feature detection, not fragment depth.}
\item{distribution_control}{Named list: max_mu (default 1e10) caps cell-feature expectations with a reported residual; protocol_policy is warn or error for heterogeneous measurement overrides. New mode requires isolated RNG.}
}
\details{distribution_v2 scales latent-cell intensities before NB/Poisson sampling and retains the original measurement operator. Protein ambient expectation is subtracted from the total signal budget; a target below the ambient floor produces a recorded mismatch. Binary ATAC retains Bernoulli generation and OR aggregation. Only peak-level detection is calibrated; unit detected-peak variation remains a diagnostic. The reference empirical depths describe observed units. Coarse resolution, coverage, pairing or missingness overrides do not establish protocol equivalence. Profiles infer operational shared strength only from ID-established pairing.}
