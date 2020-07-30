#find and remove features that only showed up in one cohort (for multi cohort phenotypes only) after meta-analysis

suppressMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
args <- c('2020_T2D_Data/Karlsson_metaanlysis.rds', '2020_T2D_Data/Karlsson_metaanalysis_filtered.rds')

d=readRDS(args[[1]])

#get number of datasets
dsetnum=max(unlist(unname(lapply(d,function(x) nrow(x[[1]]$data)))))

if(is.finite(dsetnum)==FALSE){
  dsetnum=1
}

if(dsetnum>1){
#remove rows with not enough datasets
lengths=sapply(d,function(x) nrow(x[[1]]$data))
lengths=unname(unlist(lapply(lengths, (function(x) {if (is.null(x) | length(x) == 0) {0} else { x }}))))
toremove=colnames(d[,lengths<=1])
d=d[,lengths>1]

#save output file
saveRDS(d,args[[2]])
#saveRDS(toremove,args[[3]])
}
if(dsetnum==1){
  saveRDS(d,args[[2]])
}