# A small runnable tutorial; never starts an integration algorithm.
run_de_novo_example <- function(out_dir="example_outputs/de_novo",
                               helpers="inst/examples/example_helpers.R") {
  source(helpers, local=TRUE)
  out_dir <- example_setup(out_dir)
  on.exit(setTimeLimit(cpu=Inf,elapsed=Inf,transient=FALSE), add=TRUE)
  with_example_warnings({
    cfg <- heterospa_config(n_cells=120, n_domains=3, n_cell_types=4,
      modalities=c("RNA","ATAC","Protein"), n_features=c(RNA=40,ATAC=60,Protein=8),
      shared_dim=2, specific_dim=2, n_spatial_basis=6,
      resolution=1, pairing=1, overlap=1, misalignment="none",
      rng_mode="isolated", keep_cell_truth=TRUE, seed=42)
    paired <- simulate_spatial_multiomics(config=cfg)
    replay <- simulate_spatial_multiomics(config=cfg)
    stopifnot(identical(paired$assays,replay$assays), identical(paired$latent,replay$latent),
              identical(paired$measurement,replay$measurement))
    save_example(paired,file.path(out_dir,"paired"),"V1")
    cfg_het <- cfg
    cfg_het$resolution <- c(RNA=1,ATAC=3,Protein=2)
    cfg_het$pairing <- .5; cfg_het$capture_fraction <- .7
    cfg_het$overlap <- .6; cfg_het$misalignment <- "mild"
    heterogeneous <- simulate_spatial_multiomics(config=cfg_het)
    save_example(heterogeneous,file.path(out_dir,"heterogeneous"),"V2")
    writeLines(c("PASS: deterministic scientific replay and A/count/Jaccard GT checks",
                 "Small tutorials only. No model training or paper Benchmark reproduced."),
               file.path(out_dir,"COMPLETE.txt"))
    invisible(list(paired=paired, heterogeneous=heterogeneous))
  },file.path(out_dir,"warnings.txt"))
}
if(sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  run_de_novo_example(if(length(args)) args[1] else "example_outputs/de_novo")
}
