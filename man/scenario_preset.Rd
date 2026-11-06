\name{scenario_preset}
\alias{scenario_preset}
\title{Standard simulation scenario}
\usage{
scenario_preset (scenario = c("easy", "moderate", "hard"), modalities = c("RNA", 
    "ATAC", "Protein"), n_cells = 3000L, seed = 1L)
}
\description{Returns a ready-to-run easy, moderate or hard simulation configuration.}
\value{A HeteroSpaSim configuration list.}

\arguments{
\item{scenario}{Optional easy, moderate or hard preset used when config is NULL.}
\item{modalities}{Modalities used by a scenario preset.}
\item{n_cells}{Number of latent cells used by a scenario preset.}
\item{seed}{Finite integer random seed; an explicit simulate_spatial_multiomics seed overrides the configuration seed.}
}
