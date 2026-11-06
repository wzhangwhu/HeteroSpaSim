\name{HeteroSpaSim-metrics}
\alias{score_domains}
\alias{score_alignment}
\alias{score_imputation}
\alias{score_binary_imputation}
\alias{score_correspondence}
\alias{score_latent_recovery}
\title{Standard HeteroSpaSim benchmark metrics}
\usage{
score_domains (truth, predicted)

score_alignment (truth, estimated)

score_imputation (truth, predicted)

score_binary_imputation (truth, score)

score_correspondence (truth, predicted)

score_latent_recovery (truth, estimated)
}
\description{Reference implementations of common benchmark metrics for domain recovery, spatial alignment, continuous and binary imputation, cross-modal correspondence and latent-factor recovery.}
\value{Named numeric vectors of task-specific metrics.}

\arguments{
\item{truth}{Ground truth for the selected task: label vector, coordinate matrix, expression matrix, binary vector, correspondence matrix or latent matrix. Dimensions and unit order must match the estimate.}
\item{predicted}{Task-specific estimate: domain labels, expression matrix or one target match index per source unit. See the corresponding function usage.}
\item{estimated}{Coordinates or latent-factor representation for the same units in the same row order.}
\item{score}{Predicted probabilities or continuous scores.}
}
