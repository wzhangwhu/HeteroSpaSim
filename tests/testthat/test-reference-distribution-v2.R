v2_fixture <- function() {
    set.seed(990)
    mods<-c("RNA","ATAC","Protein");p<-c(24L,28L,8L)
    ans<-lapply(seq_along(mods),function(i) {
        x<-matrix(rpois(60*p[i],if(i==2) .2 else 2),60,p[i]);dimnames(x)<-list(paste0("u",1:60),paste0(mods[i],1:p[i]))
        xy<-cbind(seq(0,1,length.out=60),rep(0:5,each=10)/5);rownames(xy)<-rownames(x)
        list(counts=Matrix::Matrix(x,sparse=TRUE),coordinates=xy,domain=rep(1:3,each=20))
    });names(ans)<-mods;ans
}
v2_run <- function(ref=v2_fixture(),...) simulate_from_reference(ref,calibration="distribution_v2",shared_dim=2,specific_dim=2,n_spatial_basis=4,keep_cell_truth=TRUE,seed=42,...)

test_that("distribution profiles validate counts, IDs, coordinates and pairing", {
    r<-v2_fixture();p<-fit_reference_profile(r)
    expect_s3_class(p,"HeteroSpaReferenceProfile");expect_false(p$provenance$full_counts_stored)
    expect_false(any(vapply(p$per_modality,function(x) "counts"%in%names(x),logical(1))))
    for(mod in names(r)) expect_s3_class(fit_reference_profile(r[mod]),"HeteroSpaReferenceProfile")
    z<-r;z$RNA$counts[1,1]<-.5;expect_error(fit_reference_profile(z),"integer")
    z<-r;rownames(z$RNA$counts)[1]<-rownames(z$RNA$counts)[2];z$RNA$coordinates<-NULL;expect_error(fit_reference_profile(z),"unique")
    z<-r;z$RNA$coordinates<-z$RNA$coordinates[60:1,];expect_error(fit_reference_profile(z),"order")
    z<-r;z$ATAC$counts<-z$ATAC$counts[1:40,];z$ATAC$coordinates<-z$ATAC$coordinates[1:40,];z$ATAC$domain<-z$ATAC$domain[1:40]
    rownames(z$ATAC$counts)<-rownames(z$ATAC$coordinates)<-paste0("other",1:40)
    q<-fit_reference_profile(z);expect_false(q$pairing_diagnostics$RNA__ATAC$ID_established)
    expect_false(q$provenance$shared_strength_estimated)
    z<-r;z$RNA$counts[1,]<-0;z$RNA$counts[,1]<-0;q<-fit_reference_profile(z)
    expect_gte(q$per_modality$RNA$zero_rows,1);expect_gte(q$per_modality$RNA$feature$all_zero_features,1)
    expect_true(all(is.finite(q$per_modality$RNA$feature$dispersion)))
    z<-r;z$RNA$counts[1,1]<-2^54;expect_error(fit_reference_profile(z),"precision")
    z<-r;z$RNA$counts<-as.matrix(z$RNA$counts);expect_s3_class(fit_reference_profile(z),"HeteroSpaReferenceProfile")
})

test_that("new mode is reproducible, restores RNG, preserves biology and supports chunks", {
    r<-v2_fixture();set.seed(812);before<-.Random.seed
    p<-fit_reference_profile(r);expect_identical(.Random.seed,before)
    a<-v2_run(r,reference_profile=p,chunk_size=10);expect_identical(.Random.seed,before)
    b<-v2_run(r,reference_profile=p,chunk_size=17)
    expect_identical(a$assays,b$assays);expect_identical(a$latent,b$latent)
    oldcfg<-a$parameters;oldcfg$reference_distribution_profile<-NULL;oldcfg$reference_depth_model<-NULL;oldcfg$reference_distribution_control<-NULL
    old<-simulate_spatial_multiomics(oldcfg,coordinates=r$RNA$coordinates,domain_labels=r$RNA$domain)
    expect_identical(a$latent,old$latent);expect_identical(a$measurement,old$measurement)
    for(m in names(r)) expect_identical(a$truth$feature[[m]]$shared_loadings,old$truth$feature[[m]]$shared_loadings)
    rm(".Random.seed",envir=.GlobalEnv);invisible(v2_run(r));expect_false(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))
    expect_error(v2_run(r,rng_mode="legacy"),"isolated")
    expect_error(v2_run(r,distribution_control=list(made_up=1)),"Unknown")
    expect_warning(v2_run(r,n_features=c(RNA=30,ATAC=28,Protein=8)),"protocol equivalence")
    z<-r;colnames(z$RNA$counts)<-rev(colnames(z$RNA$counts));expect_error(v2_run(z,reference_profile=p),"feature universe")
    a<-v2_run(r);bad<-a$measurement$operators$RNA;bad[2,1]<-1
    expect_error(HeteroSpaSim:::.simulate_reference_assay("RNA",a$parameters,a$latent$shared,a$latent$specific$RNA,a$metadata$batch,bad,TRUE),"disjoint")
})

test_that("GT and cell aggregation remain valid for heterogeneous measurement overrides", {
    r<-v2_fixture()
    conditions<-list(cell=list(),coarse=list(resolution=c(RNA=1,ATAC=3,Protein=2)),
        partial=list(pairing=.5),mosaic=list(overlap=.5),shift=list(misalignment="severe"),
        missing=list(missing_rate=.2),batch=list(batch_count=2,batch_strength=.3))
    for(cc in conditions) {
        x<-suppressWarnings(do.call(v2_run,c(list(ref=r),cc)));expect_true(validate_simulation(x))
        for(m in names(r)) {
            agg<-x$measurement$operators[[m]]%*%x$truth$cell_assays[[m]];if(m=="ATAC") agg<-1*(agg>0)
            expect_equal(as.matrix(agg),as.matrix(x$truth$complete_assays[[m]]),ignore_attr=TRUE)
            # validate_simulation also checks observed/complete correspondence and masks.
            z<-x$truth$feature[[m]]$distribution_calibration
            if(m!="ATAC") expect_equal(z$achieved_expected_depth,z$target_expected_depth,tolerance=1e-9)
        }
    }
    expect_warning(v2_run(r,resolution=2),"protocol equivalence")
    expect_error(v2_run(r,resolution=2,distribution_control=list(protocol_policy="error")),"protocol equivalence")
})

test_that("zero depth, ambient floor, count ATAC and numerical caps are explicit", {
    r<-v2_fixture();r$RNA$counts[,]<-0
    x<-v2_run(r);expect_equal(sum(x$assays$RNA),0)
    r$Protein$counts[,]<-0
    expect_warning(x<-v2_run(r),"ambient floor")
    expect_gt(x$truth$feature$Protein$distribution_calibration$ambient_floor_exceeds_target,0)
    x<-suppressWarnings(v2_run(v2_fixture(),atac_binary=FALSE));expect_true(validate_simulation(x))
    z<-x$truth$feature$ATAC$distribution_calibration;expect_equal(z$achieved_expected_depth,z$target_expected_depth,tolerance=1e-9)
    warnings_seen<-character()
    x<-withCallingHandlers(v2_run(v2_fixture(),distribution_control=list(max_mu=.001)),warning=function(w) {warnings_seen<<-c(warnings_seen,conditionMessage(w));invokeRestart("muffleWarning")})
    expect_true(any(grepl("max_mu",warnings_seen)))
    expect_gt(x$truth$feature$RNA$distribution_calibration$clipped_cell_feature_expectations,0)
    expect_true(all(is.finite(x$assays$RNA@x)))
})

test_that("statistical comparator distinguishes distances, NA and calibration roles", {
    r<-v2_fixture();x<-v2_run(r)
    d<-compare_reference_statistics(r,list(new=x));expect_true(all(c("calibration_role","n_reference","status")%in%names(d)))
    expect_true(any(d$metric=="Moran_I"&d$status=="DEFINED"))
    expect_true(any(d$metric=="domain_between_total_fraction"&d$status=="DEFINED"))
    expect_true(all(d$status[d$metric=="log1p_dispersion"&d$modality=="ATAC"]=="NOT_APPLICABLE_OR_MISSING"))
    w<-HeteroSpaSim:::.reference_w1;expect_equal(w(c(0,2),1),1);expect_equal(w(0:2,0:2),0);expect_true(is.na(w(NA,1)))
})
