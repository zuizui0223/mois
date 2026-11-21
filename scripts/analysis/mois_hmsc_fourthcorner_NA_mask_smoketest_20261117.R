#!/usr/bin/env Rscript
# MOIS TRUE GBIF native-range recorded-assemblage HMSC source / NA smoke.
# A strict technical initialization and 5+5 MCMC diagnostic only.
# NOT valid posterior inference or a reported flower-colour association.
# No nonnative species/cell combo is relabeled "true ecological absence".
suppressPackageStartupMessages({
    library(Hmsc)
    library(jsonlite)
})
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
folder<-args[[1]]
out<-args[[2]]
if(dir.exists(out))stop("Refuse overwrite previous source signed Hmsc test")
dir.create(out,recursive=TRUE)
X<-read.csv(file.path(folder,"X_true_GBIF_100km_environment.csv"),
            check.names=FALSE)
Tr<-read.csv(file.path(folder,"TR_photo9_species_traits.csv"),
             check.names=FALSE,colClasses="character")
Y<-as.matrix(read.csv(file.path(folder,"Y_apparent_recorded_not_absence.csv"),
                       check.names=FALSE,na.strings="NA"))
Native<-as.matrix(read.csv(file.path(folder,"Y_species_native_TDWG3_eligibility.csv"),
                            check.names=FALSE))
storage.mode(Y)<-"numeric"
storage.mode(Native)<-"numeric"
stopifnot(nrow(X)==nrow(Y),ncol(Y)==nrow(Tr),
          all(colnames(Y)==Tr$key),
          all(is.na(Y)==(Native==0)),
          all(is.finite(Y[!is.na(Y)])),
          "GBIF_all_selected_species_record_count" %in% colnames(X),
          !any(grepl("genus",names(X),ignore.case=TRUE)))
# Select species only by frozen input order (original stable hash and support).
keep<-which(colSums(Y==1,na.rm=TRUE)>=5 &
            colSums(Y==0,na.rm=TRUE)>=5 &
            colSums(!is.na(Y))>=25)
keep<-keep[seq_len(min(80L,length(keep)))]
stopifnot(length(keep)>=40)
Y<-Y[,keep,drop=FALSE]
Tr<-Tr[keep,,drop=FALSE]
trait<-data.frame(white_proxy=as.numeric(Tr$white_proxy))
rownames(trait)<-Tr$key
stopifnot(sum(trait$white_proxy==1)>5,
          sum(trait$white_proxy==0)>5)
raw<-data.frame(
    current_temp=as.numeric(X$present_temperature_C),
    log_precip=log1p(as.numeric(X$present_precip_mm)),
    NPP=as.numeric(X$NPP_2001_2023_mean_kgC_m2_yr),
    LGM=as.numeric(X$threeGCM_delta_temperature_C_median),
    log_AGB=log1p(as.numeric(X$ESA_AGB_woody_Mg_ha_2020)),
    NDVI=as.numeric(X$MODIS_NDVI_2001_2021_mean),
    log_GBIF_record_effort=log1p(as.numeric(
        X$GBIF_all_selected_species_record_count)))
stopifnot(all(is.finite(as.matrix(raw))))
XData<-as.data.frame(scale(raw))
stopifnot(all(is.finite(as.matrix(XData))))
studyDesign<-data.frame(site=factor(sprintf("site_%03d",seq_len(nrow(Y)))))
siteR<-HmscRandomLevel(units=studyDesign$site)
siteR$nfMax<-2L
model<-Hmsc(
    Y=Y, XData=XData,
    XFormula=~current_temp + log_precip + NPP + LGM +
              log_AGB + NDVI + log_GBIF_record_effort,
    TrData=trait,TrFormula=~white_proxy,
    studyDesign=studyDesign,ranLevels=list(site=siteR),
    distr="probit")
# Deliberately insufficient MCMC to prevent unwarranted p/credible intervals.
# This only checks that the correct source matrix fits through the sampler.
fit<-sampleMcmc(model,samples=5,transient=5,thin=1,
                nChains=1,nParallel=1,verbose=0,engine="R")
stopifnot(inherits(fit,"Hmsc"),length(fit$postList)==1L,
          length(fit$postList[[1]])==5L)
s <- list(
    schema="mois_hmsc_true_GBIF_native_range_photo9_technical_smoke_v1",
    source_panel=basename(folder),
    Hmsc_version=as.character(packageVersion("Hmsc")),
    n_sites=nrow(Y),n_species=ncol(Y),
    n_white_species=sum(trait$white_proxy==1),
    n_native_eligible_site_species=sum(!is.na(Y)),
    n_native_GBIF_recorded_site_species=sum(Y==1,na.rm=TRUE),
    n_unknown_not_recorded_but_native_eligible=sum(Y==0,na.rm=TRUE),
    n_outside_native_range_missing=sum(is.na(Y)),
    uses_WCVP_native_range_NA_mask=TRUE,
    zeros_are_GBIF_unrecorded_not_true_absence=TRUE,
    trait_is_photo9_proxy_not_verified_petal=TRUE,
    current_climate_NPP_LGM_AGB_NDVI_all_in_model=TRUE,
    genus_as_model_covariate=FALSE,
    JSDM_has_site_random_level=TRUE,
    mcmc_only_5_samples_after_5_transient=TRUE,
    model_convergence_and_posterior_inference_NOT_assessed=TRUE,
    true_ecological_occupancy_or_refugium_NOT_established=TRUE)
writeLines(toJSON(s,pretty=TRUE,auto_unbox=TRUE),file.path(out,"receipt.json"))
cat("MOIS_TRUE_GBIF_HMSC_FOURTHCORNER_NA_SMOKE_SUCCESS",
    toJSON(s[c("n_sites","n_species","n_native_eligible_site_species",
                "n_outside_native_range_missing","Hmsc_version")],
           auto_unbox=TRUE),"\n")
