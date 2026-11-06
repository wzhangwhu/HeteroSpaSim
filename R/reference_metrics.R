# Exact one-dimensional empirical Wasserstein distance, including unequal sizes.
.reference_w1 <- function(a,b) {
    a<-sort(as.numeric(a));b<-sort(as.numeric(b));a<-a[is.finite(a)];b<-b[is.finite(b)]
    if(!length(a)||!length(b)) return(NA_real_)
    breaks<-sort(unique(c(0,seq_along(a)/length(a),seq_along(b)/length(b))))
    mid<-(utils::head(breaks,-1)+utils::tail(breaks,-1))/2
    sum(diff(breaks)*abs(a[pmax(1,ceiling(mid*length(a)))]-b[pmax(1,ceiling(mid*length(b)))]))
}

.reference_spatial_stats <- function(x,xy,labels,max_units=300L,max_features=60L) {
    # Deterministic equally spaced row/feature subsets: no RNG and no full dense matrix.
    ii<-unique(round(seq(1,nrow(x),length.out=min(nrow(x),max_units))))
    jj<-unique(round(seq(1,ncol(x),length.out=min(ncol(x),max_features))))
    dep<-as.numeric(Matrix::rowSums(x));v<-log1p(as.matrix(x[ii,jj,drop=FALSE])/pmax(dep[ii],1)*1e4)
    moran<-rep(NA_real_,ncol(v));domain<-rep(NA_real_,ncol(v))
    if(!is.null(xy)&&length(ii)>3) {
        ds<-as.matrix(stats::dist(xy[ii,,drop=FALSE]));diag(ds)<-Inf;k<-min(6L,length(ii)-1L)
        W<-matrix(0,length(ii),length(ii));for(i in seq_along(ii)) W[i,order(ds[i,])[seq_len(k)]]<-1
        z<-sweep(v,2,colMeans(v));den<-colSums(z^2)
        moran<-length(ii)/sum(W)*colSums(z*(W%*%z))/den;moran[den<=0]<-NA_real_
    }
    if(!is.null(labels)&&length(labels)==nrow(x)) {
        lab<-factor(labels[ii]);grand<-colMeans(v);total<-colSums(sweep(v,2,grand)^2)
        within<-numeric(ncol(v));for(lev in levels(lab)) {g<-v[lab==lev,,drop=FALSE];within<-within+colSums(sweep(g,2,colMeans(g))^2)}
        domain<-(total-within)/total;domain[total<=0]<-NA_real_
    }
    list(Moran_I=moran,domain_between_total_fraction=domain,units=length(ii),features=length(jj))
}

#' Compare reference and simulated measurement statistics
#'
#' Descriptive distances are not biological fidelity or formal hypothesis tests.
#' Uses the same feature universe and binary/count ATAC interpretation for each input.
#' @param reference Named reference list accepted by fit_reference_profile.
#' @param simulations Named list of HeteroSpaSim simulation objects.
#' @param atac_binary Whether ATAC counts represent detection indicators.
#' @return Long data frame with source value, reference value, distance, applicability,
#' calibration role and effective sample sizes. Spatial diagnostics use deterministic
#' subsets of at most 300 units and 60 features and library-normalized log1p counts.
#' @export
compare_reference_statistics <- function(reference,simulations,atac_binary=TRUE) {
    ref<-.prepare_reference(reference,atac_binary=atac_binary);rows<-list();z<-0L
    for(m in names(ref)) {
        rx<-.reference_counts(ref[[m]]);rxy<-.reference_coordinates(ref[[m]]);rlab<-.reference_labels(ref[[m]])$domain
        stat<-function(x,xy,lab) {
            f<-.profile_feature_stats(x,!(m=="ATAC"&&atac_binary));d<-as.numeric(Matrix::rowSums(x));nz<-as.numeric(Matrix::rowSums(x>0));s<-.reference_spatial_stats(x,xy,lab)
            ok<-f$mean>0&f$variance>0
            mv<-if(sum(ok)>2) suppressWarnings(stats::cor(log1p(f$mean[ok]),log1p(f$variance[ok]))) else NA_real_
            list(depth=d,log1p_depth=log1p(d),detected_features=nz,feature_mean=f$mean,feature_detection=f$detected_fraction,
                 feature_zero=f$zero_fraction,log1p_feature_variance=log1p(f$variance),log1p_dispersion=log1p(f$dispersion),
                 mean_variance_log_correlation=mv,Moran_I=s$Moran_I,domain_between_total_fraction=s$domain_between_total_fraction)
        }
        rs<-stat(rx,rxy,rlab)
        for(nm in names(simulations)) {
            sim<-simulations[[nm]];x<-sim$assays[[m]]
            if(is.null(x)) stop("Simulation is missing reference modality: ",m,call.=FALSE)
            if(ncol(x)!=ncol(rx)) stop("Comparison requires the same feature universe (no implicit feature resampling).",call.=FALSE)
            if(m=="ATAC"&&!identical(isTRUE(sim$parameters$atac_binary),isTRUE(atac_binary))) stop("ATAC mode differs between reference and simulation.",call.=FALSE)
            ss<-stat(x,sim$spatial$observed_coordinates[[m]],sim$truth$measurement_annotations[[m]]$observed_domain)
            for(metric in names(rs)) {
                r<-rs[[metric]];v<-ss[[metric]];role<-if(metric%in%c("depth","log1p_depth","feature_mean","feature_detection","feature_zero","log1p_dispersion")) "calibration_related" else "independent_diagnostic"
                if(metric%in%c("feature_detection","feature_zero")&&m!="ATAC") role<-"independent_diagnostic"
                if(metric=="log1p_dispersion"&&m=="ATAC") role<-"diagnostic_not_Poisson_calibration_target"
                if(metric%in%c("depth","log1p_depth")&&m=="ATAC"&&atac_binary) role<-"independent_detected_peaks_distribution"
                z<-z+1L;distance<-.reference_w1(r,v)
                rows[[z]]<-data.frame(simulation=nm,modality=m,metric=metric,calibration_role=role,statistic="empirical_W1",reference_value=NA_real_,simulation_value=NA_real_,distance=distance,n_reference=sum(is.finite(r)),n_simulation=sum(is.finite(v)),status=if(is.finite(distance)) "DEFINED" else "NOT_APPLICABLE_OR_MISSING",stringsAsFactors=FALSE)
                if(metric%in%c("feature_zero","feature_detection","mean_variance_log_correlation")) {
                    z<-z+1L;ar<-if(any(is.finite(r))) mean(r[is.finite(r)]) else NA_real_;av<-if(any(is.finite(v))) mean(v[is.finite(v)]) else NA_real_
                    rows[[z]]<-data.frame(simulation=nm,modality=m,metric=metric,calibration_role=role,statistic="mean",reference_value=ar,simulation_value=av,distance=abs(ar-av),n_reference=sum(is.finite(r)),n_simulation=sum(is.finite(v)),status=if(is.finite(ar)&&is.finite(av)) "DEFINED" else "NOT_APPLICABLE_OR_MISSING",stringsAsFactors=FALSE)
                }
                if(metric%in%c("depth","detected_features")) {
                    a<-.distribution_summary(r);b<-.distribution_summary(v)
                    for(t in names(a)) {z<-z+1L;rows[[z]]<-data.frame(simulation=nm,modality=m,metric=metric,calibration_role=role,statistic=t,reference_value=a[[t]],simulation_value=b[[t]],distance=abs(a[[t]]-b[[t]]),n_reference=length(r),n_simulation=length(v),status=if(is.finite(a[[t]])&&is.finite(b[[t]])) "DEFINED" else "NOT_APPLICABLE_OR_MISSING",stringsAsFactors=FALSE)}
                }
            }
        }
    }
    rownames_out<-do.call(rbind,rows);rownames(rownames_out)<-NULL;rownames_out
}
