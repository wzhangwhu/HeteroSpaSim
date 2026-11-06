# Explicit optional demo: do not install Python or its modules automatically.
run_optional_h5ad <- function(out_dir="example_outputs/h5ad",
                             helpers="inst/examples/example_helpers.R") {
  if(!requireNamespace("reticulate",quietly=TRUE) || !reticulate::py_module_available("anndata"))
    stop("Optional export needs reticulate and a selected Python with anndata; MTX examples need neither.")
  source(helpers,local=TRUE); out_dir <- example_setup(out_dir)
  on.exit(setTimeLimit(cpu=Inf,elapsed=Inf,transient=FALSE),add=TRUE)
  sim <- simulate_spatial_multiomics(heterospa_config(n_cells=40,n_domains=3,
    n_features=c(RNA=12,ATAC=15,Protein=5),shared_dim=2,specific_dim=2,n_spatial_basis=4,seed=421))
  paths <- export_h5ad(sim,file.path(out_dir,"h5ad"),include_truth=FALSE,design="V2",backend="python")
  ad <- reticulate::import("anndata",convert=FALSE)
  for(p in paths) {
    a <- ad$read_h5ad(p)
    stopifnot(reticulate::py_to_r(a$layers$`__contains__`("counts")),
              reticulate::py_to_r(a$obsm$`__contains__`("spatial")))
  }
  writeLines("PASS: H5AD count/spatial fields read back",file.path(out_dir,"COMPLETE.txt"))
  invisible(paths)
}
if(sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  run_optional_h5ad(if(length(args)) args[1] else "example_outputs/h5ad")
}
