# Run from repository root after verifying the v1 download receipt.
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1]] else 'signed/diel_figshare_v1'
e <- new.env()
load(file.path(root, 'DvN_analysis_dataframe final.RData'), envir=e)
stopifnot(exists('diel.all.diffs', envir=e), is.data.frame(e$diel.all.diffs))
write.csv(e$diel.all.diffs, file.path(root, 'published_analysis_export.csv'), row.names=FALSE)
