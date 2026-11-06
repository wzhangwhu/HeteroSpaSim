.reference_depth_draw <- function(s,n,cfg,m,model) {
    .with_rng_stream(cfg,paste0("reference_distribution_v2/depth/",m),{
        if(model=="empirical") s$depth[sample.int(length(s$depth),n,replace=TRUE)] else {
            z <- s$log_depth_model;mu<-z$mean;if(!is.finite(mu)) mu<-0
            ans<-stats::rlnorm(n,mu,z$sd);ans[stats::runif(n)<z$zero_fraction]<-0;ans
        }
    })
}

.reference_eta_block <- function(idx,H0,Hs,W0,Ws,rho,cfg,batch,be) {
    u0<-.standardize_cols(H0%*%t(W0[idx,,drop=FALSE]));us<-.standardize_cols(Hs%*%t(Ws[idx,,drop=FALSE]))
    eta<-cfg$signal_strength*(sqrt(rho)*u0+sqrt(1-rho)*us)
    if(!is.null(be)) eta<-eta+be[batch,idx,drop=FALSE]
    .clip(eta,-6,6)
}

.reference_binary_baseline <- function(eta,A,target) {
    r<-mean(as.numeric(Matrix::rowSums(A)))
    b<-stats::qlogis(.clip(-expm1(log1p(-.clip(target,1e-12,1-1e-12))/r),1e-12,1-1e-12))
    # Bounded feature-wise Newton solve of mean_j[1-product_i(1-p_ig)].
    for(iter in seq_len(20L)) {
        logits<-sweep(eta,2,b,"+");prob<-stats::plogis(logits)
        log0<--log1p(exp(.clip(logits,-40,40)));logspot<-as.matrix(A%*%log0)
        expected<--expm1(logspot);grad<-exp(logspot)*as.matrix(A%*%prob)
        step<-(colMeans(expected)-target)/pmax(colMeans(grad),1e-10)
        b<-.clip(b-.clip(step,-2,2),-35,35)
    }
    b[target==0]<--Inf;b[target==1]<-Inf;b
}

.simulate_reference_assay <- function(modality,cfg,H0,Hs,batch,A,keep_cell) {
    if(is.null(A)) stop("distribution_v2 requires an actual measurement operator.",call.=FALSE)
    if(any(A@x!=1)||any(Matrix::colSums(A)>1)) stop("distribution_v2 currently requires disjoint binary support rows.",call.=FALSE)
    s<-cfg$reference_distribution_profile$per_modality[[modality]];p<-cfg$n_features[[modality]];n<-nrow(H0);nr<-nrow(A)
    W0<-.simulate_loadings(p,ncol(H0),cfg$loading_density,cfg=cfg,key=paste0("loadings/",modality,"/shared"))
    Ws<-.simulate_loadings(p,ncol(Hs),cfg$loading_density,cfg=cfg,key=paste0("loadings/",modality,"/specific"))
    rho<-cfg$shared_strength[[modality]]
    source_idx<-if(p==s$n_features) seq_len(p) else .with_rng_stream(cfg,paste0("reference_distribution_v2/features/",modality),sample.int(s$n_features,p,replace=TRUE))
    target_mean<-s$feature$mean[source_idx];theta<-s$feature$dispersion[source_idx]
    baseline<-log(target_mean);be<-NULL
    if(cfg$batch_count>1&&cfg$batch_strength>0) be<-.with_rng_stream(cfg,paste0("batch_effects/",modality),matrix(stats::rnorm(cfg$batch_count*p,0,cfg$batch_strength),cfg$batch_count,p))
    mapping<-integer(n);aT<-methods::as(A,"TsparseMatrix");mapping[aT@j+1L]<-aT@i+1L
    captured<-mapping>0;starts<-seq.int(1L,p,by=cfg$chunk_size);denom<-numeric(nr)
    binary<-modality=="ATAC"&&cfg$atac_binary
    # Pass1 calibrates feature intercepts and accumulates expected spot intensity.
    for(start in starts) {
        idx<-start:min(p,start+cfg$chunk_size-1L)
        eta<-.reference_eta_block(idx,H0,Hs,W0,Ws,rho,cfg,batch,be)
        if(binary) baseline[idx]<-.reference_binary_baseline(eta,A,s$feature$detected_fraction[source_idx[idx]]) else {
            # Analytic centering preserves feature heterogeneity, not exact final column sums.
            baseline[idx]<-log(target_mean[idx])-log(colMeans(exp(eta[captured,,drop=FALSE])))
            h<-exp(sweep(eta,2,baseline[idx],"+"));denom<-denom+as.numeric(A%*%rowSums(h))
        }
    }
    target<-if(binary) rep(NA_real_,nr) else .reference_depth_draw(s,nr,cfg,modality,cfg$reference_depth_model)
    ambient<-if(modality=="Protein") as.numeric(Matrix::rowSums(A))*p*cfg$protein_ambient else rep(0,nr)
    scale<-rep(1,n);ambient_floor<-if(binary) rep(FALSE,nr) else target<ambient
    if(!binary) {
        if(any(ambient_floor)) warning("Protein target below configured ambient floor; signal expectation is zero and mismatch is recorded.",call.=FALSE)
        if(any(denom==0&target>ambient)) warning("Zero feature intensity cannot meet positive target depth; residual recorded.",call.=FALSE)
        spot_scale<-pmax(target-ambient,0)/pmax(denom,.Machine$double.xmin)
        spot_scale[denom==0]<-0;scale[captured]<-spot_scale[mapping[captured]]
    }
    blocks<-vector("list",length(starts));cells<-if(keep_cell) vector("list",length(starts)) else NULL
    achieved<-numeric(nr);clipped<-0;max_mu<-cfg$reference_distribution_control$max_mu
    for(b in seq_along(starts)) {
        idx<-starts[b]:min(p,starts[b]+cfg$chunk_size-1L)
        eta<-.reference_eta_block(idx,H0,Hs,W0,Ws,rho,cfg,batch,be)
        if(binary) {
            logits<-sweep(eta,2,baseline[idx],"+");prob<-stats::plogis(logits)
            expected<--expm1(as.matrix(A%*%(-log1p(exp(.clip(logits,-40,40))))))
            expected[,is.infinite(baseline[idx])&baseline[idx]<0]<-0
            expected[,is.infinite(baseline[idx])&baseline[idx]>0]<-1
            achieved<-achieved+rowSums(expected)
            x<-vapply(seq_along(idx),function(j) .with_rng_stream(cfg,paste0("reference_distribution_v2/counts/",modality,"/",idx[j]),stats::rbinom(n,1,prob[,j])),numeric(n))
        } else {
            mu<-exp(sweep(eta,2,baseline[idx],"+"))*scale
            clipped<-clipped+sum(mu>max_mu);mu<-pmin(mu,max_mu)
            if(any(!is.finite(mu))) stop("Nonfinite calibrated expectation; check extreme reference depths.",call.=FALSE)
            achieved<-achieved+as.numeric(A%*%rowSums(mu))
            x<-vapply(seq_along(idx),function(j) .with_rng_stream(cfg,paste0("reference_distribution_v2/counts/",modality,"/",idx[j]),{
                if(modality=="ATAC") stats::rpois(n,mu[,j]) else stats::rnbinom(n,mu=mu[,j],size=theta[idx[j]])
            }),numeric(n))
            if(modality=="Protein") x<-x+vapply(idx,function(g) .with_rng_stream(cfg,paste0("reference_distribution_v2/ambient/",modality,"/",g),stats::rpois(n,cfg$protein_ambient)),numeric(n))
        }
        xb<-Matrix::Matrix(x,sparse=TRUE);agg<-A%*%xb
        if(binary) agg<-Matrix::Matrix(1*(agg>0),sparse=TRUE)
        blocks[[b]]<-agg;if(keep_cell) cells[[b]]<-xb
    }
    if(clipped>0) warning("distribution_v2 cell expectations reached max_mu; expected depth residual recorded.",call.=FALSE)
    feat<-sprintf("%s_%06d",.feature_prefix(modality),seq_len(p));measured<-.sparse_cbind(blocks)
    dimnames(measured)<-list(rownames(A),feat);cell_counts<-if(keep_cell) .sparse_cbind(cells) else NULL
    if(keep_cell) dimnames(cell_counts)<-list(sprintf("Cell_%06d",seq_len(n)),feat)
    rownames(W0)<-rownames(Ws)<-feat
    list(counts=measured,cell_counts=cell_counts,shared_loadings=W0,specific_loadings=Ws,baseline=baseline,dispersion=theta,
        feature_names=feat,shared_strength=rho,library_factors=scale,batch_effects=be,
        distribution_calibration=list(protocol="reference-distribution-v2.1",depth_model=cfg$reference_depth_model,
            target_expected_depth=target,achieved_expected_depth=if(binary) achieved else achieved+ambient,
            ambient_expected=ambient,ambient_floor_exceeds_target=sum(ambient_floor),clipped_cell_feature_expectations=clipped,
            reference_feature_index=source_idx,uncaptured_cell_count=sum(!captured),uncaptured_rule="scale=1, uncalibrated cell molecules; never enter A aggregation",
            binary_unit_offsets=FALSE,binary_remaining_error="Only feature mean detection calibrated in expectation; detected-peaks distribution is independent diagnostic"))
}

.simulate_distribution_reference <- function(reference,preserve_geometry,preserve_labels,args,profile,depth_model,control) {
    binary<-args$atac_binary%||%TRUE
    if(is.null(profile)) profile<-fit_reference_profile(reference,atac_binary=binary)
    if(!inherits(profile,"HeteroSpaReferenceProfile")) stop("reference_profile must come from fit_reference_profile().",call.=FALSE)
    for(m in names(reference)) {
        x<-.reference_counts(reference[[m]]);xy<-.reference_coordinates(reference[[m]])
        if(!is.null(xy)&&!is.null(rownames(x))&&!is.null(rownames(xy))&&!identical(rownames(x),rownames(xy))) stop(m," coordinate IDs/order must match counts rows.",call.=FALSE)
        values<-if(inherits(x,"sparseMatrix")) x@x else as.numeric(x)
        if(any(is.finite(values)&values>2^53)) stop("Counts exceed exact integer precision (2^53).",call.=FALSE)
    }
    reference<-.prepare_reference(reference,atac_binary=binary);fit<-profile$operational_calibration
    if(!identical(profile$modalities,names(reference))||!identical(profile$atac_binary,isTRUE(binary))) stop("Reference profile modality/mode mismatch.",call.=FALSE)
    for(m in names(reference)) {
        x<-.reference_counts(reference[[m]]);s<-profile$per_modality[[m]]
        if(ncol(x)!=s$n_features||(!is.null(s$feature_ids)&&!identical(colnames(x),s$feature_ids))) stop("Reference profile feature universe/order mismatch.",call.=FALSE)
    }
    base<-list(n_cells=fit$suggested_n_cells,modalities=fit$modalities,n_features=setNames(vapply(profile$per_modality,function(s)s$n_features,integer(1)),fit$modalities),
        shared_strength=fit$shared_strength,shared_length_scale=fit$shared_length_scale,specific_length_scale=max(.03,.75*fit$shared_length_scale))
    if(!is.null(fit$domain_celltype_composition)) {base$n_domains<-fit$n_domains;base$n_cell_types<-fit$n_cell_types;base$celltype_composition<-fit$domain_celltype_composition}
    if(("n_domains"%in%names(args)||"n_cell_types"%in%names(args))&&!"celltype_composition"%in%names(args)) base$celltype_composition<-NULL
    xy<-.reference_coordinates(reference[[fit$modalities[1]]]);lab<-if(preserve_labels) .reference_labels(reference[[fit$modalities[1]]]) else NULL
    if(preserve_geometry&&!is.null(xy)&&!"n_cells"%in%names(args)) base$n_cells<-nrow(xy)
    cfg<-do.call(heterospa_config,utils::modifyList(base,args))
    if(cfg$rng_mode!="isolated") stop("distribution_v2 requires rng_mode='isolated' for named stream reproducibility.",call.=FALSE)
    if(!is.list(control)||length(setdiff(names(control),c("max_mu","protocol_policy")))) stop("Unknown distribution_control field.",call.=FALSE)
    control<-utils::modifyList(list(max_mu=1e10,protocol_policy="warn"),control)
    if(length(control$max_mu)!=1||!is.finite(control$max_mu)||control$max_mu<=0) stop("max_mu must be positive finite.",call.=FALSE)
    control$protocol_policy<-match.arg(control$protocol_policy,c("warn","error"))
    feature_override<-any(cfg$n_features!=vapply(profile$per_modality,function(s)s$n_features,integer(1)))
    hetero<-feature_override||any(cfg$resolution>1.05)||cfg$pairing<1||cfg$overlap<1||any(cfg$missing_rate>0)||(!is.null(cfg$capture_fraction)&&cfg$capture_fraction<1)
    if(hetero) {
        message<-"distribution_v2 heterogeneous measurement/feature override: reference observed-unit depths applied to synthetic complete units before missingness; protocol equivalence is not established."
        if(control$protocol_policy=="error") stop(message,call.=FALSE) else warning(message,call.=FALSE)
    }
    if(preserve_geometry&&!is.null(xy)&&nrow(xy)!=cfg$n_cells) {warning("Reference coordinates not used: n_cells mismatch.",call.=FALSE);xy<-NULL}
    if(!preserve_geometry) xy<-NULL
    if(!is.null(lab)) for(nm in names(lab)) if(length(lab[[nm]])!=cfg$n_cells) lab[[nm]]<-NULL
    cfg$reference_distribution_profile<-profile;cfg$reference_depth_model<-depth_model;cfg$reference_distribution_control<-control
    sim<-simulate_spatial_multiomics(cfg,coordinates=xy,domain_labels=lab$domain,cell_type_labels=lab$cell_type)
    sim$reference_calibration<-list(mode="distribution_v2",profile=profile,protocol_equivalence=if(hetero) "NOT_ESTABLISHED_HETEROGENEOUS_OVERRIDE" else "baseline_reference_observed_units",
        depth_model=depth_model,ambient_assumption="Configured protein_ambient, not independently identified",counts_postprocessed=FALSE,
        unit_assumption="At baseline each reference measurement unit is modeled as one synthetic latent cell; physical cell-support size is not inferred from counts")
    sim$provenance$reference_protocol<-"reference-distribution-v2.1";sim
}
