\name{export_benchmark}
\alias{export_benchmark}
\alias{export_h5ad}
\alias{export_h5mu}
\alias{export_mat}
\alias{as_SpatialExperiment}
\alias{as_MultiAssayExperiment}
\title{Export and interoperability helpers}
\usage{
export_benchmark (x, path, include_complete = TRUE, layout = c("isolated", 
    "legacy"), design = c("V1", "V2"))

export_h5ad (x, path, include_truth = FALSE, design = c("V1", "V2"), 
    backend = c("auto", "python", "zellkonverter"))

export_h5mu (x, file, include_truth = FALSE, design = c("V1", "V2"), 
    backend = c("auto", "python", "zellkonverter"))

export_mat (x, file, include_truth = FALSE, design = c("V1", "V2"))

as_SpatialExperiment (x, modality, include_truth = TRUE)

as_MultiAssayExperiment (x, include_truth = TRUE)
}
\description{Exports language-neutral benchmark directories, H5AD/H5MU/MAT files, or converts simulations to Bioconductor spatial containers. Optional formats require their corresponding suggested packages.}
\value{Paths or converted objects, depending on the function.}

\arguments{
\item{x}{A simulation object.}
\item{path}{Destination directory.}
\item{include_complete}{Write complete pre-missingness assays when available.}
\item{layout}{isolated separates public_inputs from evaluator_truth and refuses a nonempty destination. legacy retains the old complete evaluator bundle and must never be passed directly to algorithms.}
\item{design}{V1 exposes legitimate exact pairing through opaque shared spot IDs; V2 permutes each modality independently and uses modality-specific IDs. Private mappings are evaluator-only.}
\item{include_truth}{Include private ground-truth annotations. False by default for file exporters; True by default for analytical Bioconductor containers. Annotated files must be kept evaluator-only.}
\item{backend}{H5AD writing backend. python needs reticulate and Python anndata/scipy/pandas; zellkonverter uses its managed Python environment. auto chooses Python only when already initialized.}
\item{file}{Destination .mat file.}
\item{modality}{Modality name.}
}

\details{H5AD stores observations in rows and features in columns. Both supported
backends write raw counts to X and an identical counts layer, observed coordinates
to obsm spatial, and preserve public observation and feature order. H5MU retains
the corresponding per-modality fields. Truth is excluded by default; annotated
exports belong to evaluator-only directories. A D1 task is configured through
pairing, overlap and resolution, then exported with design V2; D1 is not an export
design argument. Native readback must be tested in the actual optional runtime.
}
