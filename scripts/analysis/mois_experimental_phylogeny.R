# Comparative sensitivity using full first-stage sampling covariance.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==3)
root <- args[[1]]; treefile <- args[[2]]; out <- args[[3]]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
library(ape); library(metafor)
p <- read.csv(file.path(root,'adjusted_species_profiles.csv'),check.names=FALSE)
V <- as.matrix(read.csv(file.path(root,'adjusted_species_full_sampling_covariance.csv'),row.names=1,check.names=FALSE))
asymmetry <- max(abs(V-t(V)))
stopifnot(identical(rownames(V),colnames(V)),asymmetry<=1e-10*max(1,max(abs(V))))
# Remove only verified floating-point asymmetry; no eigenvalue/variance repair.
V <- (V+t(V))/2
write.csv(data.frame(max_original_asymmetry=asymmetry),file.path(out,'sampling_covariance_roundoff_audit.csv'),row.names=FALSE)
e <- new.env(); load(treefile,envir=e); stopifnot(exists('diel.tree.out',envir=e))
tree <- e$diel.tree.out
stopifnot(!anyDuplicated(p$species),!anyDuplicated(p$phylo),all(p$phylo %in% tree$tip.label))
C <- vcv(tree,cor=TRUE)
rows <- list(); diagnostics <- list(); failures <- list()
for (definition in c('night_share','raw_control_matched','net_night_relative_outcome')) {
  d <- p[!is.na(p[[definition]]),]
  d$ratio_centered <- d[[definition]]-.5
  d$species <- as.character(d$species);d$phylo <- as.character(d$phylo)
  sampling <- V[d$species,d$species,drop=FALSE]
  phy <- C[d$phylo,d$phylo,drop=FALSE]
  stopifnot(min(eigen(sampling,symmetric=TRUE,only.values=TRUE)$values)>0)
  for (kind in c('independent_species','phylogenetic_plus_species')) {
    label <- paste(definition,kind,sep='__')
    result <- tryCatch({
      if(kind=='independent_species') {
        rma.mv(yi=adjusted_log_alan,V=sampling,mods=~white*ratio_centered,
               random=~1|species,data=d,method='REML',test='t',dfs='residual')
      } else {
        rma.mv(yi=adjusted_log_alan,V=sampling,mods=~white*ratio_centered,
               random=list(~1|species,~1|phylo),R=list(phylo=phy),
               data=d,method='REML',test='t',dfs='residual')
      }
    },error=function(error) error)
    if(inherits(result,'error')) {
      failures[[label]] <- data.frame(analysis=label,error=conditionMessage(result))
      next
    }
    rows[[label]] <- data.frame(analysis=label,term=rownames(result$b),estimate=as.numeric(result$b),
                               se=result$se,ci95_low=result$ci.lb,ci95_high=result$ci.ub,
                               n_species=nrow(d),df=result$ddf)
    diagnostics[[label]] <- data.frame(analysis=label,n_species=nrow(d),logLik=as.numeric(logLik(result)),
                                      species_variance=result$sigma2[1],phylo_variance=if(length(result$sigma2)>1)result$sigma2[2] else NA)
    saveRDS(result,file.path(out,paste0(label,'.rds')))
  }
}
write.csv(if(length(failures))do.call(rbind,failures) else data.frame(analysis=character(),error=character()),file.path(out,'unsupported_phylogenetic_models.csv'),row.names=FALSE)
stopifnot(length(rows)>0)
write.csv(do.call(rbind,rows),file.path(out,'phylogenetic_sensitivity_coefficients.csv'),row.names=FALSE)
write.csv(do.call(rbind,diagnostics),file.path(out,'phylogenetic_model_diagnostics.csv'),row.names=FALSE)
write.csv(if(length(failures))do.call(rbind,failures) else data.frame(analysis=character(),error=character()),file.path(out,'unsupported_phylogenetic_models.csv'),row.names=FALSE)
write.csv(data.frame(species=p$species,phylo_tip=p$phylo),file.path(out,'verified_species_tip_mapping.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'R_sessionInfo.txt'))
cat('phylogenetic sensitivity successful models:',length(rows),'failures:',length(failures),'\n')
