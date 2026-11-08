# HeteroSpaSim

An R framework for heterogeneous spatial multi-omics simulation with explicit
latent biology, measurement operators and multi-layer ground truth.

**Version: 0.1.9**

[Source code](https://github.com/wzhangwhu/HeteroSpaSim) ·
[Download](https://github.com/wzhangwhu/HeteroSpaSim/releases/tag/v0.1.9.9000) ·
[Automated checks](https://github.com/wzhangwhu/HeteroSpaSim/actions) ·
[Report an issue](https://github.com/wzhangwhu/HeteroSpaSim/issues)

Two routes are supported: **de novo simulation** without a real-data reference,
and **reference-calibrated simulation** from raw reference counts and optional
coordinates/labels. Both generate observations with explicit measurement and
ground-truth records.

## Installation

Requires **R >= 4.3.0** and a compatible version of **Matrix**. Core simulation
does not require Python. Install the tagged version from GitHub:

```r
install.packages("remotes")
remotes::install_github(
  "wzhangwhu/HeteroSpaSim",
  ref = "v0.1.9.9000",
  dependencies = NA,  # Required dependencies, not all optional packages.
  upgrade = "never"
)
```

Alternatively, download `HeteroSpaSim_0.1.9.9000.tar.gz` from the release page and
install it from its local path:

```r
install.packages("HeteroSpaSim_0.1.9.9000.tar.gz", repos = NULL, type = "source")
library(HeteroSpaSim)
packageVersion("HeteroSpaSim")
find.package("HeteroSpaSim")
```

From a repository checkout (the checkout root contains DESCRIPTION):

```sh
R CMD INSTALL .
```

## Route 1: de novo simulation

```r
library(HeteroSpaSim)
cfg <- heterospa_config(
  n_cells=120, n_domains=3, n_cell_types=4,
  modalities=c("RNA", "ATAC", "Protein"),
  n_features=c(RNA=40, ATAC=60, Protein=8),
  shared_dim=2, specific_dim=2, n_spatial_basis=6,
  resolution=1, pairing=1, overlap=1, misalignment="none",
  rng_mode="isolated", keep_cell_truth=TRUE, seed=42)
sim <- simulate_spatial_multiomics(config=cfg)
stopifnot(validate_simulation(sim))
export_benchmark(sim, "example_outputs/de_novo", design="V1")
```

The full runnable [de novo example](inst/examples/01_de_novo.R) also demonstrates
partial coverage, different spatial resolutions and coordinate misalignment.

## Route 2: reference-calibrated simulation

Supply a named list with `RNA`, `ATAC` and/or `Protein`. Each entry contains raw
nonnegative integer `counts` (**units × features**), optional `coordinates` (units
× 2), and optional `domain`/`cell_type` labels in the same row order. Real pairing
requires genuine matching unit IDs, not equal row counts or nearest neighbors.

```r
# Requires your own reference.rds in the schema described above.
reference <- readRDS("reference.rds")
profile <- fit_reference_profile(reference)
sim_ref <- simulate_from_reference(reference,
  calibration="distribution_v2", reference_profile=profile,
  depth_model="empirical", rng_mode="isolated", seed=42)
stopifnot(validate_simulation(sim_ref))
```

Start with [02_reference_quickstart.R](inst/examples/02_reference_quickstart.R): it
runs immediately using an explicitly synthetic input fixture. Then use
[03_reference_real_data.R](inst/examples/03_reference_real_data.R) with your own
RDS or Matrix Market inputs.
Each counts matrix needs unique, nonempty unit and feature IDs. Coordinates must
have exactly the same row IDs and order as counts; labels must follow that order.
Use `Protein` for protein/ADT. Do not pass log-normalized, CLR or TF-IDF matrices
as raw reference counts.

For example, a real-reference RDS can be prepared as follows (replace these
variables with your own data):

```r
reference <- list(
  RNA = list(counts = rna_counts, coordinates = rna_xy),
  ATAC = list(counts = atac_counts, coordinates = atac_xy)
)
saveRDS(reference, "reference.rds")
```

The real-data entry point also accepts modality folders containing `counts.mtx`
(units x features), `units.tsv` and `features.tsv` (one ID per line, no header),
and optional `coordinates.csv` (`unit,x,y`) and `labels.csv`
(`unit,domain,cell_type`, either label optional). See the validation in
[03_reference_real_data.R](inst/examples/03_reference_real_data.R).

### Reference calibration parameters

In the example above, `calibration="distribution_v2"` selects distribution-based
reference calibration. It uses reference measurement-unit depths and feature-level
statistics to generate new observations through the latent-cell model and
measurement operators.

| Parameter | Example value | Purpose |
| --- | --- | --- |
| `calibration` | `"distribution_v2"` | Select distribution-based reference calibration |
| `reference_profile` | `profile` | Reuse the summaries fitted by `fit_reference_profile()` |
| `depth_model` | `"empirical"` | Sample target depths from the reference depth distribution |
| `seed` | `42` | Set the simulation seed |

These options belong to the reference simulation route. Save a fitted profile
with `saveRDS(profile, "reference_profile.rds")` to reuse it in subsequent runs.
For the complete parameter list, see `?simulate_from_reference`.

## Running the tutorials

From the checkout root (the folder containing `DESCRIPTION`), after installing
this exact source, run these commands in a **terminal**, not the R console:

```sh
Rscript --vanilla inst/examples/01_de_novo.R example_outputs/de_novo
Rscript --vanilla inst/examples/02_reference_quickstart.R example_outputs/reference_demo
Rscript --vanilla inst/examples/03_reference_real_data.R reference.rds example_outputs/reference_real
```

In RStudio or an interactive R session, run the installed examples as follows:

```r
ex <- system.file("examples", package = "HeteroSpaSim")

# De novo: complete paired and heterogeneous three-modality examples.
source(file.path(ex, "01_de_novo.R"))
run_de_novo_example("example_outputs/de_novo",
  helpers = file.path(ex, "example_helpers.R"))

# Reference tutorial: uses a clearly labeled SYNTHETIC reference fixture.
source(file.path(ex, "02_reference_quickstart.R"))
run_reference_quickstart("example_outputs/reference_demo",
  helpers = file.path(ex, "example_helpers.R"))

# For your own real reference, run this separately after preparing reference.rds:
# source(file.path(ex, "03_reference_real_data.R"))
# run_reference_real_data("reference.rds", "example_outputs/reference_real",
#   helpers = file.path(ex, "example_helpers.R"))
```

Sourcing an example alone does not execute it. Every example refuses nonempty
output directories; use a new path for each run. Tutorials are small, fixed-seed
examples with finite time limits. Save the full object,
profile, session information and warnings for evaluation/reproducibility, but pass
only `public_inputs/` to integration algorithms. Use `find.package()` to locate
the installed package.

## Ground truth and evaluation

The object stores latent factors, domain/cell labels, measurement operators,
observed and complete counts, transformations, missingness masks and support
correspondence. Matrix rows are observed units, which need not be individual cells.
For RNA/Protein/count ATAC, complete counts are `A %*% Xcell`; binary ATAC uses OR
aggregation. Cross-modal correspondence is support-set Jaccard.

`V1` export exposes declared exact pairs; `V2` hides cross-modal ID/order mappings.
Coverage and spatial resolution are configured separately through the simulation
parameters. Binary ATAC row totals represent detected peaks.

## Export formats and optional dependencies

| Format or function | Dependencies |
| --- | --- |
| Matrix Market matrices and CSV coordinates | Core package dependencies |
| JSON metadata | jsonlite; R-readable metadata is written if jsonlite is unavailable |
| H5AD, Python backend | reticulate with Python anndata, scipy and pandas |
| H5AD, zellkonverter backend | zellkonverter and SpatialExperiment |
| H5MU | reticulate with Python anndata and mudata, plus the selected H5AD backend |
| SpatialExperiment / MultiAssayExperiment containers | The corresponding Bioconductor packages |
| MAT files | R.matlab |

To view detailed options for an export function, run its help command in the
**R console**. For example:

```r
help("export_benchmark")
help("export_h5ad")
help("export_h5mu")
help("export_mat")
```

## Contact and license

Software attribution: **HeteroSpaSim Development Team**. MIT license, see
[LICENSE.md](LICENSE.md).

Contact: **wzhang_math@whu.edu.cn**. For bug reports, use
[GitHub Issues](https://github.com/wzhangwhu/HeteroSpaSim/issues) and include the
package version, `sessionInfo()` and a small reproducible example.
