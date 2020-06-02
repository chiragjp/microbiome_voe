library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)
library(magrittr)
library(broom)
library(stringr) # regex pattern matching
library(rlang) # rlang::duplicate()
library(cowplot)
library(rje) # powerSet()
library(meta)
library(metafor)

get_good_metadfs <- function(meta_df) {
  toremove=list()
  count=0
  for(x in 1:length(meta_df)){
    if(class(meta_df[[x]][[1]])[[1]]!='metagen'){
      toremove[as.character(count)]=x
      count=count+1
    }
  }
  if(length(toremove)>0){
    meta_df=meta_df[-unname(unlist(toremove))]
  }
  return(meta_df)
}

get_summary_stats <- function(input_meta_df, cohort_type, mapping) {
  if (cohort_type == "multi") {
    meta_df <- get_good_metadfs(input_meta_df) # for meta-analyzed outputs, filter for good outputs first
    return(
      tibble(
        feature = colnames(meta_df),
        feature_name = map_chr(feature, ~mapping[[as.integer(str_split(., "_")[[1]][[2]]), 1]]),
        estimate = map_dbl(meta_df, ~.[[1]]$TE.random),
        p.val = map_dbl(meta_df, ~.[[1]]$pval.random),
        bonferroni.p.val = p.adjust(p.val, method = "bonferroni"),
        fdr.p.val = p.adjust(p.val, method = "fdr"),
        by.p.val = p.adjust(p.val, method = "BY"),
        CI_95_lower = map_dbl(meta_df, ~.[[1]]$lower.random),
        CI_95_upper = map_dbl(meta_df, ~.[[1]]$upper.random)
      )
    )
  } else if (cohort_type == "single") {
    meta_df <- input_meta_df # skip get_good_metadfs step
    rows=list()
    count=0
    for(x in meta_df){
      count=count+1
      rows[[count]]=nrow(x[[1]])
    }
    rows=unname(unlist(rows))
    clean_meta_df=meta_df[rows==1]
    return(
      tibble(
        feature = colnames(clean_meta_df),
        feature_name = map_chr(feature, ~mapping[[as.integer(str_split(., "_")[[1]][[2]]), 1]]),
        estimate = map_dbl(clean_meta_df, ~.[[1]]$estimate),
        p.val = map_dbl(clean_meta_df, ~.[[1]]$p.value),
        bonferroni.p.val = p.adjust(p.val, method = "bonferroni"),
        fdr.p.val = p.adjust(p.val, method = "fdr"),
        by.p.val = p.adjust(p.val, method = "BY"),
        CI_95_lower = map_dbl(clean_meta_df, ~(.[[1]]$estimate - 2 * .[[1]]$std.error)),
        CI_95_upper = map_dbl(clean_meta_df, ~(.[[1]]$estimate + 2 * .[[1]]$std.error))
      )
    )
  }
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  meta_outputs <- as_tibble(readRDS(args[[1]]))
  mapping <- as_tibble(readRDS(args[[2]]))
  cohort_type <- args[[3]]
  outputname <- args[[4]]

  output <- get_summary_stats(meta_outputs, cohort_type, mapping) 
  saveRDS(output, 
          outputname)
}

main()







