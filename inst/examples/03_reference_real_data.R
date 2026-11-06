read_reference_folder <- function(path) {
  mods <- c("RNA","ATAC","Protein")
  mods <- mods[dir.exists(file.path(path,mods))]
  if(!length(mods)) stop("Folder needs RNA, ATAC and/or Protein subdirectories.")
  out <- lapply(mods,function(m) {
    d <- file.path(path,m)
    x <- Matrix::readMM(file.path(d,"counts.mtx"))
    # Matrix Market may store binary counts as an implicit-one pattern matrix.
    # Materialize numeric ones without changing dimensions, sparsity or values.
    if(inherits(x,"nMatrix") || inherits(x,"lMatrix")) x <- methods::as(x,"dMatrix")
    units <- readLines(file.path(d,"units.tsv"),warn=FALSE)
    feats <- readLines(file.path(d,"features.tsv"),warn=FALSE)
    if(length(units)!=nrow(x)||length(feats)!=ncol(x)) stop(m,": counts must be units x features")
    dimnames(x) <- list(units,feats); z <- list(counts=x)
    coord_file <- file.path(d,"coordinates.csv")
    if(file.exists(coord_file)) {
      xy <- utils::read.csv(coord_file,check.names=FALSE,colClasses=c(unit="character"))
      if(!all(c("unit","x","y")%in%names(xy))||!identical(xy$unit,units)) stop(m,": coordinate IDs/order differ")
      z$coordinates <- as.matrix(xy[,c("x","y")]); rownames(z$coordinates) <- units
    }
    label_file <- file.path(d,"labels.csv")
    if(file.exists(label_file)) {
      ann <- utils::read.csv(label_file,check.names=FALSE,colClasses=c(unit="character"))
      if(!identical(ann$unit,units)) stop(m,": label IDs/order differ")
      for(k in c("domain","cell_type")) if(k%in%names(ann)) z[[k]] <- ann[[k]]
    }
    z
  }); names(out) <- mods; out
}

run_reference_real_data <- function(input, out_dir="example_outputs/reference_real",
                                   helpers="inst/examples/example_helpers.R") {
  source(helpers,local=TRUE)
  if(missing(input)||!file.exists(input)) stop("Provide a real reference RDS file or input folder; see input schema.")
  reference <- if(dir.exists(input)) read_reference_folder(input) else readRDS(input)
  example_validate_reference(reference)
  too_large <- any(vapply(reference,function(z) nrow(z$counts)>5000 || ncol(z$counts)>10000,logical(1)))
  if(too_large && Sys.getenv("HETEROSPASIM_ALLOW_LARGE_REFERENCE")!="true")
    stop("Large input: review resource needs, then explicitly opt in; no automatic subsetting is performed.")
  out_dir <- example_setup(out_dir)
  on.exit(setTimeLimit(cpu=Inf,elapsed=Inf,transient=FALSE),add=TRUE)
  with_example_warnings({
    files <- if(dir.exists(input)) list.files(input,recursive=TRUE,full.names=TRUE) else input
    files <- files[!file.info(files)$isdir]
    utils::write.csv(data.frame(path=normalizePath(files,winslash="/"),md5=unname(tools::md5sum(files))),
                     file.path(out_dir,"private_input_provenance.csv"),row.names=FALSE)
    profile <- fit_reference_profile(reference)
    saveRDS(profile,file.path(out_dir,"reference_profile.rds"))
    sim <- simulate_from_reference(reference,calibration="distribution_v2",reference_profile=profile,
      preserve_geometry=!is.null(reference[[1]]$coordinates),preserve_labels=TRUE,
      depth_model="empirical",seed=42,rng_mode="isolated",keep_cell_truth=TRUE)
    save_example(sim,file.path(out_dir,"simulation"),"V2")
    metrics <- compare_reference_statistics(reference,list(distribution_v2=sim))
    utils::write.csv(metrics,file.path(out_dir,"statistical_comparison.csv"),row.names=FALSE,na="NA")
    writeLines(c("PASS: user reference validation, simulation, aggregation and Jaccard checks",
                 "V2 hides public pairing IDs; it does not itself impose partial support coverage."),
               file.path(out_dir,"COMPLETE.txt"))
    invisible(sim)
  },file.path(out_dir,"warnings.txt"))
}
if(sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  if(length(args)<1L) stop("Usage: Rscript inst/examples/03_reference_real_data.R INPUT [NEW_OUTPUT]")
  run_reference_real_data(args[1],if(length(args)>1L) args[2] else "example_outputs/reference_real")
}
