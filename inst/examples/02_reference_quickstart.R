# Synthetic reference fixture for learning the interface, not real biological data.
make_reference_demo_fixture <- function() {
  set.seed(190019)
  xy <- as.matrix(expand.grid(x=seq(0,1,length.out=10),y=seq(0,1,length.out=8)))
  ids <- sprintf("paired_unit_%03d",seq_len(nrow(xy))); rownames(xy) <- ids
  domain <- ifelse(xy[,1]<.5,"left","right")
  depth <- exp(rnorm(nrow(xy),0,.6))
  make <- function(m,p) {
    mu <- outer(depth,seq(.5,3,length.out=p))
    counts <- if(m=="ATAC") matrix(rbinom(length(mu),1,plogis(log(mu)-3)),nrow(mu)) else
      matrix(rnbinom(length(mu),mu=as.numeric(mu),size=8),nrow(mu))
    dimnames(counts) <- list(ids,paste0(m,"_",seq_len(p)))
    list(counts=Matrix::Matrix(counts,sparse=TRUE),coordinates=xy,domain=domain)
  }
  list(RNA=make("RNA",40),ATAC=make("ATAC",50),Protein=make("Protein",8))
}

run_reference_quickstart <- function(out_dir="example_outputs/reference_demo",
                                    helpers="inst/examples/example_helpers.R") {
  source(helpers,local=TRUE)
  out_dir <- example_setup(out_dir)
  on.exit(setTimeLimit(cpu=Inf,elapsed=Inf,transient=FALSE),add=TRUE)
  with_example_warnings({
    reference <- make_reference_demo_fixture()
    example_validate_reference(reference)
    profile <- fit_reference_profile(reference)
    stopifnot(!profile$provenance$full_counts_stored)
    shared_args <- list(reference=reference, seed=42, shared_dim=2,specific_dim=2,
                        n_spatial_basis=6,rng_mode="isolated",keep_cell_truth=TRUE)
    moment <- do.call(simulate_from_reference,c(shared_args,list(calibration="moment")))
    calibrated <- do.call(simulate_from_reference,c(shared_args,list(calibration="distribution_v2",
                      reference_profile=profile,depth_model="empirical")))
    replay <- do.call(simulate_from_reference,c(shared_args,list(calibration="distribution_v2",
                      reference_profile=profile,depth_model="empirical")))
    stopifnot(identical(calibrated$assays,replay$assays))
    save_example(moment,file.path(out_dir,"moment"),"V1")
    save_example(calibrated,file.path(out_dir,"distribution_v2"),"V1")
    saveRDS(reference,file.path(out_dir,"SYNTHETIC_reference.rds"))
    saveRDS(profile,file.path(out_dir,"reference_profile.rds"))
    metrics <- compare_reference_statistics(reference,list(moment=moment,distribution_v2=calibrated))
    utils::write.csv(metrics,file.path(out_dir,"statistical_comparison.csv"),row.names=FALSE,na="NA")
    writeLines(c("PASS: three-modality synthetic input, saved profile, GT and replay checks",
                 "This reference is SYNTHETIC. It does not demonstrate real-data fidelity."),
               file.path(out_dir,"COMPLETE.txt"))
    invisible(list(reference=reference,profile=profile,simulation=calibrated))
  },file.path(out_dir,"warnings.txt"))
}
if(sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  run_reference_quickstart(if(length(args)) args[1] else "example_outputs/reference_demo")
}
