#!/usr/bin/env Rscript
# MOIS: actual Bayesian HMSC probit fourth-corner of PHOTO colour × environment
# against source-authenticated GBIF recording in WCVP broad-native cells.
# Do not interpret GBIF unrecorded as field-verified ecological absence.
# Do not interpret species-modal photographic hue as field-observed petal colour.
# No genus effect, no genus weighting. The site random effect/latent factor is
# present at TRAINING sites only. External test is fixed-only prediction.
suppressPackageStartupMessages({
  library(Hmsc)
  library(jsonlite)
  library(coda)
})
options(warn=1)
args<-commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==3)
folder<-args[[1]]
out<-args[[2]]
panel_name<-args[[3]]
if (dir.exists(out)) stop("Never overwrite original JSDM posterior results")
dir.create(out,recursive=TRUE)
set.seed(20261117)
X<-read.csv(file.path(folder,"X_true_GBIF_100km_environment.csv"),
            check.names=FALSE)
Tr<-read.csv(file.path(folder,"TR_photo9_species_traits.csv"),
             check.names=FALSE,colClasses="character")
Y<-as.matrix(read.csv(file.path(folder,"Y_apparent_recorded_not_absence.csv"),
                      check.names=FALSE,na.strings="NA"))
native<-as.matrix(read.csv(file.path(folder,"Y_species_native_TDWG3_eligibility.csv"),
                           check.names=FALSE))
storage.mode(Y)<-"numeric"
storage.mode(native)<-"numeric"
stopifnot(nrow(X)==300, nrow(Y)==nrow(X), ncol(Y)==nrow(Tr),
          all(colnames(Y)==Tr$key), all(is.na(Y)==(native==0)),
          all(Y[!is.na(Y)] %in% c(0,1)),
          "GBIF_all_selected_species_record_count" %in% names(X),
          !any(grepl("genus",names(X),ignore.case=TRUE)))
# One prespecified geographical split. Each test 500km block is entirely held out.
blocks<-sort(unique(as.character(X$block_500km)))
hold_blocks<-blocks[seq_along(blocks) %% 5L == 0L]
test<-which(as.character(X$block_500km) %in% hold_blocks)
train<-setdiff(seq_len(nrow(Y)),test)
stopifnot(length(hold_blocks)>=6L,length(train)>180L,length(test)>25L,
          !length(intersect(as.character(X$block_500km[train]),
                            as.character(X$block_500km[test]))))
# Training-only species support. Initial roster was sampled from GBIF
# recorded occurrence opportunity across all sites; that earlier
# screening is not fully held out. Do not call validation fully external.
eligible<-colSums(!is.na(Y[train,,drop=FALSE]))
positives<-colSums(Y[train,,drop=FALSE]==1,na.rm=TRUE)
negatives<-colSums(Y[train,,drop=FALSE]==0,na.rm=TRUE)
can<-which(positives>=5L & negatives>=5L & eligible>=25L)
# Select by frozen alphabetic source key, never white status or fitted effects.
can<-can[seq_len(min(length(can),80L))]
stopifnot(length(can)>=40L)
Tr<-Tr[can,,drop=FALSE]
Y<-Y[,can,drop=FALSE]
native<-native[,can,drop=FALSE]
trait<-data.frame(white_proxy=as.numeric(Tr$white_proxy))
rownames(trait)<-Tr$key
colnames(Y)<-Tr$key
stopifnot(sum(trait$white_proxy==1)>8L,
          sum(trait$white_proxy==0)>8L,
          colSums(!is.na(Y[train,,drop=FALSE])) |> min()>=25L)
raw<-data.frame(
    current_temperature_C=as.numeric(X$present_temperature_C),
    log_current_precip=log1p(as.numeric(X$present_precip_mm)),
    NPP_2001_2023=as.numeric(X$NPP_2001_2023_mean_kgC_m2_yr),
    log_woody_AGB=log1p(as.numeric(X$ESA_AGB_woody_Mg_ha_2020)),
    NDVI_2001_2021=as.numeric(X$MODIS_NDVI_2001_2021_mean),
    LGM_warming=as.numeric(X$threeGCM_delta_temperature_C_median),
    LGM_precip_change=as.numeric(X$threeGCM_delta_log_precip_median),
    log_GBIF_record_count=log1p(as.numeric(
        X$GBIF_all_selected_species_record_count))
)
stopifnot(all(is.finite(as.matrix(raw))))
# Scale only using TRAINING site distributions, never test.
smean<-vapply(raw[train,,drop=FALSE],mean,numeric(1))
ssd<-vapply(raw[train,,drop=FALSE],sd,numeric(1))
stopifnot(all(is.finite(ssd)),all(ssd>1e-9))
scaled<-as.data.frame(sweep(sweep(as.matrix(raw),2,smean,"-"),2,ssd,"/"))
stopifnot(all(is.finite(as.matrix(scaled))))
xtrain<-scaled[train,,drop=FALSE]
xtest<-scaled[test,,drop=FALSE]
ytrain<-Y[train,,drop=FALSE]
ytest<-Y[test,,drop=FALSE]
stopifnot(sum(!is.na(ytest))>350,sum(ytest==1,na.rm=TRUE)>100L,
          sum(ytest==0,na.rm=TRUE)>100L)

variants<-list(
  "current_climate_only"=c("current_temperature_C",
        "log_current_precip","log_GBIF_record_count"),
  "plus_NPP_AGB_NDVI"=c("current_temperature_C",
        "log_current_precip","log_GBIF_record_count",
        "NPP_2001_2023","log_woody_AGB","NDVI_2001_2021"),
  "plus_LGM_history"=names(raw)
)
# Probit pooled sample prediction uses posterior Beta fixed component ONLY.
# New sites do not get latent site effects inferred from their outcome Y.
# HMSC XScale=FALSE is essential; data already standardized from train.
fixed_only_prediction<-function(fit,Xnew) {
  A<-model.matrix(fit$XFormula,data=Xnew)
  stopifnot(nrow(A)==nrow(Xnew), ncol(A)==fit$nc,
            all(colnames(A)==fit$covNames))
  cols<-ncol(Y)
  total<-matrix(0,nrow(A),cols)
  seen<-0L
  for (chain in fit$postList) {
    for (draw in chain) {
      beta<-draw$Beta
      stopifnot(all(dim(beta)==c(ncol(A),cols)))
      p<-pnorm(A %*% beta)
      total<-total+p
      seen<-seen+1L
    }
  }
  stopifnot(seen>100L)
  total/seen
}
metric<-function(prob,truth,who) {
  ok<-is.finite(truth)
  p<-pmin(pmax(prob[ok],1e-6),1-1e-6)
  y<-truth[ok]
  stopifnot(length(p)>350L,all(y %in% c(0,1)))
  auc<-unname(as.numeric(pROC::auc(pROC::roc(
       y,p,direction="<",quiet=TRUE))))
  data.frame(variant=who,n_test_species_cells=length(y),
    n_recorded=sum(y),n_unrecorded=sum(1-y),
    log_loss=mean(-(y*log(p)+(1-y)*log(1-p))),
    brier=mean((y-p)^2),
    roc_auc=auc)
}
pres<-colSums(ytrain==1,na.rm=TRUE)
obs<-colSums(!is.na(ytrain))
prevalence<-(pres+0.5)/(obs+1)
base<-matrix(rep(prevalence,each=length(test)),length(test),length(prevalence))
scores<-metric(base,ytest,"train_species_prevalence_baseline")
gc()
contrasts<-list()
posterior_summary<-list()
for (variant in names(variants)) {
  print(paste("MOIS_HMSC_REGION_FIT_BEGIN",panel_name,variant,
              "train",length(train),"test",length(test),
              "species",ncol(ytrain)),quote=FALSE)
  cols<-variants[[variant]]
  formula<-reformulate(cols)
  training<-data.frame(site=factor(sprintf("site_%03d",train)))
  siteR<-HmscRandomLevel(units=training$site)
  siteR$nfMax<-2L
  m<-Hmsc(Y=ytrain,XData=xtrain,
          XFormula=formula,XScale=FALSE,
          TrData=trait,TrFormula=~white_proxy,TrScale=FALSE,
          studyDesign=training,ranLevels=list(site=siteR),
          distr="probit")
  stopifnot(m$nc==length(cols)+1L,m$nt==2L)
  fit<-sampleMcmc(m,samples=180,transient=250,thin=2,
                  nChains=2,nParallel=1,verbose=0,
                  engine="R")
  stopifnot(length(fit$postList)==2,
            all(vapply(fit$postList,length,integer(1))==180L))
  p<-fixed_only_prediction(fit,xtest)
  scores<-rbind(scores,metric(p,ytest,variant))
  write.csv(data.frame(site_row=rep(test,times=ncol(ytest)),
       species_key=rep(colnames(ytest),each=length(test)),
       recorded=as.vector(ytest),p_fixed_only=as.vector(p),
       variant=variant),
       file.path(out,paste0(variant,"_spatial_holdout_predictions.csv")),
       row.names=FALSE,na="")
  post<-getPostEstimate(fit,"Gamma",q=c(.025,.975))
  # Gamma rows = environmental effects, columns = trait intercept and white.
  # Verify dimension/name mapping rather than hard-code an index.
  stopifnot(all(dim(post$mean)==c(m$nc,m$nt)),
            all(dim(post$support)==c(m$nc,m$nt)))
  wcol<-match("white_proxy",fit$trNames)
  if (is.na(wcol)) wcol<-match("white_proxy",colnames(post$mean))
  stopifnot(length(wcol)==1,!is.na(wcol))
  result<-list()
  for (v in cols) {
    ix<-match(v,fit$covNames)
    stopifnot(length(ix)==1,!is.na(ix))
    bychain<-lapply(fit$postList,function(ch)
       vapply(ch,function(draw)as.numeric(draw$Gamma[ix,wcol]),numeric(1)))
    postchain<-coda::mcmc.list(lapply(bychain,coda::mcmc))
    rh<-tryCatch(as.numeric(coda::gelman.diag(postchain,
          autoburnin=FALSE,multivariate=FALSE)$psrf[1,1]),
         error=function(e)NaN)
    ess<-tryCatch(as.numeric(coda::effectiveSize(postchain)[1]),
         error=function(e)NaN)
    value<-as.numeric(post$mean[ix,wcol])
    lo<-as.numeric(post$quantiles[ix,wcol,1])
    hi<-as.numeric(post$quantiles[ix,wcol,2])
    pos<-as.numeric(post$support[ix,wcol])
    # Do NOT call posterior sign/probability a causal selection result.
    result[[v]]<-list(gamma_trait_environment_probit_mean=value,
           CrI_95=c(lo,hi),posterior_support_positive=pos,
           Rhat=rh,ESS=ess,
           convergence_screen_1point1_and_ESS100=
               is.finite(rh)&&rh<=1.1&&ess>=100)
  }
  contrasts[[variant]]<-result
  posterior_summary[[variant]]<-list(
     n_chains=2L,posterior_draws_per_chain=180L,
     transient=250L,thin=2L,
     Hmsc_version=as.character(packageVersion("Hmsc")),
     all_reported_chains_converged=
        all(vapply(result,function(z)
            z$convergence_screen_1point1_and_ESS100,logical(1))))
  print(paste("MOIS_HMSC_FOURTHCORNER_MODEL_COMPLETED",
      panel_name,variant,toJSON(result,
        auto_unbox=TRUE,digits=8)),quote=FALSE)
  rm(fit,m);gc()
}
write.csv(scores,file.path(out,"500km_unseen_block_model_comparison.csv"),
          row.names=FALSE)
writeLines(toJSON(contrasts,pretty=TRUE,auto_unbox=TRUE,digits=8),
           file.path(out,"posterior_flower_white_fourthcorner_gamma.json"))
receipt<-list(
    schema="MOIS_true_GBIF_WCVP_flower_white_Hmsc_JSDM_spatial_holdout_v1",
    region=panel_name,
    source_panel=basename(folder),
    Hmsc_version=as.character(packageVersion("Hmsc")),
    n_train_cells=length(train),
    n_test_cells=length(test),
    n_train_500km_blocks=length(unique(X$block_500km[train])),
    n_test_500km_blocks=length(hold_blocks),
    n_included_species=ncol(Y),
    n_photo_white_species=sum(trait$white_proxy==1),
    n_photo_nonwhite_species=sum(trait$white_proxy==0),
    white_is_species_photo9_proxy_not_verified_petal=TRUE,
    GBIF_unrecorded_is_not_true_absence=TRUE,
    WCVP_native_is_broad_TDWG3_not_individual_wild=TRUE,
    geographic_test_blocks_not_train=TRUE,
    site_latent_effect_not_estimated_using_test_Y=TRUE,
    test_prediction_fixed_species_beta_only=TRUE,
    GBIF_effort_proxy_not_survey_effort=TRUE,
    panel_species_screened_using_GBIF_full_matrix_before_split=TRUE,
    model_interactions_are_associational_not_adaptation=TRUE,
    true_LGM_refugia_not_identified=TRUE,
    photo_hue_selected_as_outcome=FALSE,
    latent_site_random_level_used_in_train=TRUE,
    X_standardization_training_only=TRUE,
    sample_diagnostics=posterior_summary,
    heldout_prediction_metrics=scores,
    fourthcorner_white_photo_gamma=contrasts)
writeLines(toJSON(receipt,pretty=TRUE,auto_unbox=TRUE,digits=8),
           file.path(out,"receipt.json"))
cat("MOIS_HMSC_TRUE_JSDM_FLOWER_COLOUR_FOURTHCORNER_SPATIAL_HOLDOUT_SUCCESS",
    toJSON(list(region=panel_name,n_species=ncol(Y),
       sites_train=length(train),sites_test=length(test),
       performance=scores,trait_effects=contrasts[["plus_LGM_history"]]),
       auto_unbox=TRUE,digits=8),"\n")
