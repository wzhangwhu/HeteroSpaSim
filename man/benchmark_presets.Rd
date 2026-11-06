\name{benchmark_presets}
\alias{benchmark_presets}
\alias{generate_benchmark}
\title{Standard HeteroSpaSim benchmark suites}
\usage{
benchmark_presets ()

generate_benchmark (benchmark, n_replicates = 10L, n_cells = 3000L, modalities = c("RNA", 
    "ATAC", "Protein"), base_config = NULL, out_dir = NULL, seed_start = 1L)
}
\description{Defines and generates reproducible benchmark axes for shared signal, partial pairing, mosaic overlap, resolution mismatch, spatial misalignment and technical noise.}
\value{\code{benchmark_presets()} returns definitions. \code{generate_benchmark()} returns a manifest and simulations, or writes benchmark directories.}

\arguments{
\item{benchmark}{Name returned by benchmark_presets().}
\item{n_replicates}{Number of independent seeds per parameter setting.}
\item{n_cells}{Number of latent cells used by a scenario preset.}
\item{modalities}{Modalities used by a scenario preset.}
\item{base_config}{Optional baseline configuration.}
\item{out_dir}{Optional directory. When supplied, each simulation is exported.}
\item{seed_start}{First random seed.}
}
