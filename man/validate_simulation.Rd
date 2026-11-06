\name{validate_simulation}
\alias{validate_simulation}
\title{Validate a HeteroSpaSim object}
\usage{
validate_simulation (x, stop_on_error = FALSE)
}
\description{Checks top-level object structure, modality dimensions, coordinate dimensions and measurement-operator dimensions.}
\value{Invisible TRUE on success or a character vector of problems.}

\arguments{
\item{x}{A simulation object.}
\item{stop_on_error}{Stop if validation fails.}
}
