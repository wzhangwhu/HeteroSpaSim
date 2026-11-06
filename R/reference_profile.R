.distribution_summary <- function(x) {
    x <- as.numeric(x); x <- x[is.finite(x)]
    if (!length(x)) return(c(mean=NA,median=NA,P10=NA,P25=NA,P75=NA,P90=NA,CV=NA))
    c(mean=mean(x),median=stats::median(x),setNames(as.numeric(stats::quantile(x,c(.1,.25,.75,.9),names=FALSE)),c("P10","P25","P75","P90")),
      CV=if (length(x)>1 && mean(x)>0) stats::sd(x)/mean(x) else NA_real_)
}

.profile_feature_stats <- function(x, count_model=TRUE) {
    n <- nrow(x); mu <- as.numeric(Matrix::colMeans(x)); sx <- as.numeric(Matrix::colSums(x))
    sx2 <- as.numeric(Matrix::colSums(x^2)); va <- pmax(0,(sx2-n*mu^2)/max(n-1,1))
    detected <- as.numeric(Matrix::colSums(x>0)); depth <- as.numeric(Matrix::rowSums(x))
    # NB2 moment residual conditional on library size; no dense units x features copy.
    q <- if (sum(depth)>0) sx/sum(depth) else rep(0,ncol(x))
    dx <- as.numeric(Matrix::crossprod(Matrix::Matrix(depth,ncol=1),x))
    denominator <- q^2 * sum(depth^2)
    alpha <- (sx2-2*q*dx+denominator-sx)/pmax(denominator,1e-12)
    usable <- is.finite(alpha)&alpha>0&detected>=3&mu>0
    center <- if (any(usable)) stats::median(.clip(alpha[usable],1e-6,10)) else .1
    raw <- .clip(ifelse(is.finite(alpha),pmax(alpha,0),center),1e-6,10)
    weight <- detected/(detected+20)
    theta <- 1/(weight*raw+(1-weight)*center)
    theta <- .clip(theta,.1,1e6)
    list(mean=mu,variance=va,detected_fraction=detected/n,zero_fraction=1-detected/n,
         dispersion=if(count_model) theta else rep(NA_real_,ncol(x)),
         dispersion_method=if(count_model) "depth-offset NB2 residual moments; detection-count shrinkage; theta [0.1,1e6]" else "NOT_APPLICABLE_BERNOULLI",
         dispersion_unstable_features=sum(!usable),all_zero_features=sum(sx==0))
}

#' Fit a reusable distribution-level reference calibration profile
#'
#' @param reference Named RNA/ATAC/Protein list of units-by-features integer counts,
#'   optionally with coordinates, domain and cell_type. No matrix is stored in the profile.
#' @param atac_binary Binarize integer ATAC input before fitting detection statistics.
#' @param estimate_shared Estimate operational shared strength only for ID-established pairs.
#' @return A serializable HeteroSpaReferenceProfile. Feature summaries and empirical
#'   measurement-depth samples are retained, not real expression templates.
#' @export
fit_reference_profile <- function(reference, atac_binary=TRUE, estimate_shared=TRUE) {
    # Operational fitting may subsample large references. Isolate that RNG as well.
    had_seed<-exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
    old_seed<-if(had_seed) get(".Random.seed",envir=.GlobalEnv) else NULL
    on.exit({if(had_seed) assign(".Random.seed",old_seed,envir=.GlobalEnv) else if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)) rm(".Random.seed",envir=.GlobalEnv)},add=TRUE)
    set.seed(190019L)
    # Validate coordinate IDs before the old preparation function drops no information.
    for (m in names(reference)) {
        v <- reference[[m]];x <- .reference_counts(v);xy <- .reference_coordinates(v)
        if (!is.null(xy) && !is.null(rownames(x)) && !is.null(rownames(xy)) && !identical(rownames(x),rownames(xy)))
            stop(m," coordinate IDs/order must match counts rows.",call.=FALSE)
        values <- if (inherits(x,"sparseMatrix")) x@x else as.numeric(x)
        if (any(is.finite(values)&values>2^53)) stop("Counts exceed exact integer precision (2^53).",call.=FALSE)
    }
    ref <- .prepare_reference(reference,atac_binary=atac_binary);mods <- names(ref)
    paired <- list(); reliable <- TRUE
    if(length(mods)>1) for(i in seq_len(length(mods)-1)) for(j in (i+1):length(mods)) {
        a <- rownames(.reference_counts(ref[[mods[i]]]));b <- rownames(.reference_counts(ref[[mods[j]]]))
        n <- if (!is.null(a)&&!is.null(b)) length(intersect(a,b)) else 0L
        established <- n>=3;reliable <- reliable&&established
        paired[[paste(mods[c(i,j)],collapse="__")]] <- list(ID_established=established,common_IDs=n,
            note="Equal row counts alone are not evidence of pairing; simulated pairing remains an explicit design parameter.")
    }
    fit <- calibrate_reference(ref,estimate_shared=isTRUE(estimate_shared)&&reliable,atac_binary=atac_binary)
    per <- lapply(mods,function(m) {
        x <- .reference_counts(ref[[m]]);dep <- as.numeric(Matrix::rowSums(x));nz <- as.numeric(Matrix::rowSums(x>0))
        binary <- m=="ATAC"&&isTRUE(atac_binary);fs <- .profile_feature_stats(x,!binary)
        list(modality=m,n_units=nrow(x),n_features=ncol(x),feature_ids=colnames(x),
            depth=dep,depth_summary=.distribution_summary(dep),log1p_depth_summary=.distribution_summary(log1p(dep)),
            log_depth_model=list(mean=mean(log(pmax(dep[dep>0],1))),sd=if(sum(dep>0)>1) stats::sd(log(dep[dep>0])) else 0,zero_fraction=mean(dep==0)),
            detected_features=nz,detected_summary=.distribution_summary(nz),feature=fs,
            depth_semantics=if(binary) "detected_peaks_NOT_fragment_sequencing_depth" else "total_counts_per_reference_observed_unit",
            zero_rows=sum(dep==0),global_zero_fraction=1-sum(x>0)/(nrow(x)*ncol(x)),
            coordinates_available=!is.null(.reference_coordinates(ref[[m]])),
            coordinate_order_evidence=if(!is.null(rownames(.reference_coordinates(ref[[m]])))) "row_IDs_checked" else "row_order_assumed_when_coordinates_present",
            ambient=if(m=="Protein") "Not identified; configured protein_ambient retained as explicit assumption" else NULL)
    });names(per)<-mods
    result <- list(protocol="reference-distribution-v2.1",profile_version=1L,modalities=mods,atac_binary=isTRUE(atac_binary),
        per_modality=per,operational_calibration=fit,pairing_diagnostics=paired,
        calibration_targets=c("measurement_depth","feature_abundance","ATAC_peak_detection","depth_offset_dispersion"),
        independent_diagnostics=c("detected_feature_distribution","spatial_autocorrelation","within_between_domain_signal","GT_consistency"),
        provenance=list(created_at=as.character(Sys.time()),input_orientation="units_by_features",full_counts_stored=FALSE,
            pairing_inferred_from_equal_rows=FALSE,shared_strength_estimated=isTRUE(estimate_shared)&&reliable))
    class(result)<-c("HeteroSpaReferenceProfile","list");result
}
