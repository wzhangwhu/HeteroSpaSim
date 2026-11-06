example_setup <- function(out_dir) {
  library(HeteroSpaSim)
  if (as.character(utils::packageVersion("HeteroSpaSim")) != "0.1.9.9000")
    stop("Install/load the prepared 0.1.9.9000 library first; inspect find.package().")
  if (dir.exists(out_dir) && length(list.files(out_dir, all.files=TRUE, no..=TRUE)))
    stop("Output directory is not empty; use a new path.")
  dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
  setTimeLimit(cpu=Inf, elapsed=180, transient=FALSE)
  writeLines(c(capture.output(sessionInfo()), paste("Loaded package:", find.package("HeteroSpaSim"))),
             file.path(out_dir,"sessionInfo.txt"))
  invisible(normalizePath(out_dir, winslash="/", mustWork=TRUE))
}

with_example_warnings <- function(expr, path) {
  seen <- character()
  on.exit(writeLines(if(length(seen)) seen else "No warnings", path), add=TRUE)
  withCallingHandlers(expr, warning=function(w) {
    seen <<- c(seen, conditionMessage(w)) # Keep the original warning visible as well.
  })
}

example_gt_checks <- function(sim) {
  stopifnot(isTRUE(HeteroSpaSim::validate_simulation(sim)))
  rows <- lapply(names(sim$assays), function(m) {
    a <- sim$measurement$operators[[m]]
    agg <- a %*% sim$truth$cell_assays[[m]]
    if(m=="ATAC" && sim$parameters$atac_binary) agg <- 1*(agg>0)
    stopifnot(isTRUE(all.equal(as.matrix(agg), as.matrix(sim$truth$complete_assays[[m]]), check.attributes=FALSE)))
    mask <- sim$masks[[m]]$observed
    stopifnot(isTRUE(all.equal(as.matrix(sim$truth$complete_assays[[m]][mask,,drop=FALSE]),
                               as.matrix(sim$assays[[m]]), check.attributes=FALSE)))
    data.frame(modality=m, latent_cells=ncol(a), complete_units=nrow(a),
               observed_units=nrow(sim$assays[[m]]), features=ncol(sim$assays[[m]]),
               cells_per_unit_mean=mean(Matrix::rowSums(a)), aggregation_matches=TRUE)
  })
  for(nm in names(sim$correspondence)) {
    mm <- strsplit(nm,"__",fixed=TRUE)[[1]]
    a <- sim$measurement$operators[[mm[1]]][sim$masks[[mm[1]]]$observed,,drop=FALSE]
    b <- sim$measurement$operators[[mm[2]]][sim$masks[[mm[2]]]$observed,,drop=FALSE]
    overlap <- as.matrix(Matrix::tcrossprod(a,b))
    union <- outer(as.numeric(Matrix::rowSums(a)),as.numeric(Matrix::rowSums(b)),"+")-overlap
    gt <- overlap/pmax(union,1)
    stopifnot(isTRUE(all.equal(gt,as.matrix(sim$correspondence[[nm]]),check.attributes=FALSE)))
  }
  do.call(rbind,rows)
}

save_example <- function(sim, path, design) {
  dir.create(path, recursive=TRUE, showWarnings=FALSE)
  checks <- example_gt_checks(sim)
  utils::write.csv(checks,file.path(path,"GT_checks.csv"),row.names=FALSE)
  private <- file.path(path,"private"); dir.create(private)
  saveRDS(sim,file.path(private,"simulation.rds"))
  saveRDS(sim$parameters,file.path(private,"config.rds"))
  HeteroSpaSim::export_benchmark(sim,file.path(path,"benchmark"),design=design)
  stopifnot(dir.exists(file.path(path,"benchmark","public_inputs")),
            dir.exists(file.path(path,"benchmark","evaluator_truth")))
  invisible(checks)
}

example_validate_reference <- function(reference) {
  if(!is.list(reference) || is.null(names(reference)) || anyDuplicated(names(reference)) ||
     !length(reference) || !all(names(reference)%in%c("RNA","ATAC","Protein")))
    stop("Use a uniquely named RNA/ATAC/Protein list.")
  for(m in names(reference)) {
    z <- reference[[m]]; x <- z$counts
    if(is.null(x) || length(dim(x))!=2L || any(dim(x)<2L)) stop(m,": invalid counts dimensions")
    for(ids in list(rownames(x),colnames(x)))
      if(is.null(ids)||anyNA(ids)||any(ids=="")||anyDuplicated(ids)) stop(m,": unique row/feature IDs required")
    v <- if(inherits(x,"sparseMatrix")) x@x else as.numeric(x)
    if(any(!is.finite(v))||any(v<0)||any(v!=round(v))) stop(m,": raw nonnegative integer counts required")
    if(!is.null(z$coordinates)) {
      xy <- as.matrix(z$coordinates)
      if(!is.numeric(xy)||nrow(xy)!=nrow(x)||ncol(xy)!=2||any(!is.finite(xy))||
         !identical(rownames(xy),rownames(x))) stop(m,": coordinates must exactly match counts IDs/order")
    }
    for(k in c("domain","cell_type"))
      if(!is.null(z[[k]]) && (length(z[[k]])!=nrow(x)||anyNA(z[[k]]))) stop(m,": invalid ",k," labels")
  }
  invisible(TRUE)
}
