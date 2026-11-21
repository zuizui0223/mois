#!/usr/bin/env Rscript
# MOIS: actual gllvm fourth-corner JSDM pilot for *GBIF recording events*
# (NOT true occurrence/absence). WCVP native-range eligibility is an NA mask.
#
# Response: species×site GBIF-recorded 1 vs native-eligible but not recorded 0.
# Excluded: outside WCVP TDWG3 native range = NA (no invented 0).
# Trait: observed photograph9 white-proxy, not independently verified petals.
# Predictor: contemporary temp / precipitation / NDVI / woody AGB,
#            MOD17 NPP, LGM climate change, GBIF record effort proxy.
# 500km site-block heldout. Predictions at new sites use level=0 and do not
# use held-out site's Y to infer a latent factor. NO genus fixed effect.
# Test is about recording-pattern prediction; no inference of local adaptation.
suppressPackageStartupMessages({
  library(gllvm)
  library(jsonlite)
})
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
panel <- args[[1]]
output <- args[[2]]
if (dir.exists(output)) stop("refusing to replace prior gllvm results")
dir.create(output, recursive=TRUE)
v <- as.character(packageVersion("gllvm"))
stopifnot(utils::compareVersion(v,"2.0.0")>=0)
X <- read.csv(file.path(panel,"X_true_GBIF_100km_environment.csv"),
              check.names=FALSE)
TR <- read.csv(file.path(panel,"TR_photo9_species_traits.csv"),
               colClasses="character",check.names=FALSE)
Y <- as.matrix(read.csv(file.path(panel,"Y_apparent_recorded_not_absence.csv"),
                        check.names=FALSE,na.strings="NA"))
A <- as.matrix(read.csv(file.path(panel,"Y_species_native_TDWG3_eligibility.csv"),
                        check.names=FALSE))
storage.mode(Y) <- "numeric"
storage.mode(A) <- "numeric"
stopifnot(nrow(Y)==nrow(X),ncol(Y)==nrow(TR),all(colnames(Y)==TR$key),
          all(colnames(A)==TR$key),
          all(is.na(Y)==(A==0)),all(Y[!is.na(Y)] %in% c(0,1)),
          all(TR$source_hue9 %in% c("white","yellow","orange","red",
                                    "pink","magenta","purple","blue","bronze")))
TR$white_proxy <- as.integer(TR$white_proxy)
stopifnot(length(unique(TR$white_proxy))==2,
          ncol(Y)>=100,nrow(Y)==300,
          "GBIF_all_selected_species_record_count" %in% colnames(X))
# Geographic split fixed by stable order of 500-km blocks, not flower colour.
blocks <- sort(unique(as.character(X$block_500km)))
held <- blocks[seq_along(blocks) %% 5L == 0L]
te <- which(X$block_500km %in% held)
tr <- setdiff(seq_len(nrow(X)),te)
stopifnot(length(te)>25,length(tr)>100,length(unique(X$block_500km[tr]))>25,
          !any(X$block_500km[tr] %in% X$block_500km[te]))
keep <- colSums(Y[tr,,drop=FALSE]==1,na.rm=TRUE)>=4 &
        colSums(Y[tr,,drop=FALSE]==0,na.rm=TRUE)>=4 &
        colSums(!is.na(Y[tr,,drop=FALSE]))>=30
# Model stability gate, not independent colour/response cherry-picking:
# retain up to 110 already trait-blind SHA-ranked accepted species,
# without looking at their white_proxy.
k <- which(keep)
if(length(k)<60)stop("insufficient survey-eligible recording data for JSDM")
k <- k[seq_len(min(length(k),110L))]
Y <- Y[,k,drop=FALSE]
TR <- TR[k,,drop=FALSE]
stopifnot(ncol(Y)==nrow(TR),sum(TR$white_proxy)>5,
          sum(TR$white_proxy==0)>5)
R <- data.frame(
  current_temp = as.numeric(X$present_temperature_C),
  log_precip = log1p(as.numeric(X$present_precip_mm)),
  NDVI = as.numeric(X$MODIS_NDVI_2001_2021_mean),
  log_AGB = log1p(as.numeric(X$ESA_AGB_woody_Mg_ha_2020)),
  NPP = as.numeric(X$NPP_2001_2023_mean_kgC_m2_yr),
  LGM = as.numeric(X$threeGCM_delta_temperature_C_median),
  log_effort = log1p(as.numeric(X$GBIF_all_selected_species_record_count)))
stopifnot(all(is.finite(as.matrix(R))))
mu <- vapply(R[tr,,drop=FALSE],mean,numeric(1))
sig <- vapply(R[tr,,drop=FALSE],sd,numeric(1))
stopifnot(all(sig>1e-9))
Xscaled <- data.frame(scale(R,center=mu,scale=sig))
stopifnot(all(is.finite(as.matrix(Xscaled))))
trainY <- Y[tr,,drop=FALSE]
testY <- Y[te,,drop=FALSE]
trainX <- Xscaled[tr,,drop=FALSE]
testX <- Xscaled[te,,drop=FALSE]
traits <- data.frame(white_proxy=TR$white_proxy)
rownames(traits) <- TR$key

# Explicit formulas as in the CRAN gllvm fourth-corner vignette. Species
# intercepts are internally fitted. 1 latent factor is JSDM structure.
form_current <- as.formula(
  "y ~ current_temp + log_precip + NDVI + log_AGB + log_effort +
    (current_temp + log_precip + NDVI + log_AGB):white_proxy")
form_full <- as.formula(
  "y ~ current_temp + log_precip + NDVI + log_AGB + log_effort + NPP + LGM +
    (current_temp + log_precip + NDVI + log_AGB + NPP + LGM):white_proxy")

fit_one <- function(formula,tag,seed) {
  m <- gllvm(y=trainY,X=trainX,TR=traits,
             family="binomial",num.lv=1,method="VA",
             formula=formula,seed=seed,row.eff=FALSE,
             control.start=list(n.init=1))
  p <- predict(m,newX=testX,newTR=traits,type="response",level=0)
  p <- as.matrix(p)
  stopifnot(identical(dim(p),dim(testY)),
            all(is.finite(p)),all(p>=0),all(p<=1))
  use <- is.finite(testY)
  obs <- testY[use]
  pred <- p[use]
  stopifnot(length(obs)>500,
            sum(obs==1)>50,sum(obs==0)>50)
  pred <- pmax(1e-7,pmin(1-1e-7,pred))
  log_loss <- -mean(obs*log(pred)+(1-obs)*log(1-pred))
  brier <- mean((pred-obs)^2)
  pr <- sum(obs==1)/length(obs)
  rank_auc <- sum(rank(pred)[obs==1])
  n1 <- sum(obs==1)
  n0 <- sum(obs==0)
  auc <- (rank_auc-n1*(n1+1)/2)/(n1*n0)
  fourth <- m$fourth.corner
  fc <- if(!is.null(fourth)) as.matrix(fourth) else matrix(numeric(0),0,0)
  stat <- list(
    model=tag,one_latent_factor=TRUE,
    n_train_sites=nrow(trainY),n_test_sites=nrow(testY),
    n_species=ncol(trainY),
    n_eligible_test_species_sites=length(obs),
    n_recorded_test_species_sites=n1,
    test_fraction_recorded=pr,
    heldout_500km_blocks=as.character(held),
    heldout_500km_block_count=length(held),
    test_logloss=as.numeric(log_loss),
    test_brier=as.numeric(brier),
    test_AUC=as.numeric(auc),
    fourthcorner_coefficients=unclass(as.data.frame(fc)),
    fourthcorner_rownames=rownames(fc),
    fourthcorner_colnames=colnames(fc),
    predictor_scaling_sd=as.list(sig),
    target_is_GBIF_recorded_NOT_true_occupancy=TRUE,
    unrecorded_does_not_mean_absent=TRUE,
    native_is_WCVP_TDWG3_not_individual=TRUE,
    colour_is_photo9_trait_proxy=TRUE
  )
  saveRDS(m,file.path(output,paste0(tag,"_gllvm.rds")))
  cat("MOIS_GLLVM_TRUE_RECORDING_FOURTHCORNER_MODEL",
      jsonlite::toJSON(stat[c("model","n_train_sites","n_test_sites",
         "n_species","test_logloss","test_brier","test_AUC",
         "n_eligible_test_species_sites")],auto_unbox=TRUE),"\n")
  stat
}
current <- fit_one(form_current,"current_environment",20261117)
full <- fit_one(form_full,"plus_NPP_LGM_fourthcorner",20261118)
receipt <- list(
 source="signed_GBIF_WCVP_photo9_native_100km_20261117",
 source_panel=basename(panel),
 gllvm_version=v,
 n_sites=nrow(Y),
 n_species=ncol(Y),
 n_white_proxy_species=sum(TR$white_proxy==1),
 n_nonwhite_proxy_species=sum(TR$white_proxy==0),
 n_total_recorded_species_sites=sum(Y==1,na.rm=TRUE),
 response_is_opportunistic_GBIF_recording=TRUE,
 model_is_not_true_SD_or_refugia=TRUE,
 genus_effect_in_model=FALSE,
 no_photographic_petal_truth=TRUE,
 prediction_heldout_block_and_level0=TRUE,
 models=list(current,full),
 heldout_delta_logloss_full_minus_current=
     full$test_logloss-current$test_logloss,
 heldout_delta_brier_full_minus_current=
     full$test_brier-current$test_brier,
 heldout_delta_AUC_full_minus_current=
     full$test_AUC-current$test_AUC,
 result_is_descriptive_prediction_not_occupancy_causality=TRUE)
writeLines(toJSON(receipt,pretty=TRUE,auto_unbox=TRUE,digits=9),
           con=file.path(output,"receipt.json"))
cat("MOIS_GLLVM_RECORDING_FOURTHCORNER_PILOT_SUCCESS",
    toJSON(receipt[c("gllvm_version","source_panel",
       "n_species","heldout_delta_logloss_full_minus_current",
       "heldout_delta_AUC_full_minus_current")],auto_unbox=TRUE),"\n")
