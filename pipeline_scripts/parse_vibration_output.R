library(tidyverse)
library(stringr)
library(cowplot)
library(ggplot2)
library(tidyr)
library(magrittr)
library(broom)
library(meta)
library(stringr) 
library(rlang) 
library(lme4)
library(lmerTest)
library(broom.mixed)

args <- commandArgs(trailingOnly = TRUE)

filter_unnest_feature_vib <- function(vib_df, feature_num) {
  vib_df = vib_df %>% filter(feature == feature_num)
  vib_df = vib_df %>% slice(which(map_lgl(vib_df$feature_fit, ~class(.)[[1]] == "tbl_df")))
  vib_df= vib_df %>% filter(unlist(map(feature_fit,nrow))==1)
  temp=vib_df$feature_fit
  temp=bind_rows(temp)
  vib_df = vib_df %>% select(-feature_fit) %>% bind_cols(temp)
  return(vib_df)
}

get_adjuster_expanded_vibrations <- function(voe_df,feature) {
  if(nrow(voe_df)>50000){
    voe_df=sample_n(voe_df,50000)
  }
  copy_voe_df <- duplicate(voe_df, shallow = FALSE)
  adjusters <- get_all_adjusters(voe_df,feature)
  for (variable in adjusters) {
    copy_voe_df %<>% mutate(newcol = map_int(copy_voe_df$vars, ~(variable %in% .)))
    colnames(copy_voe_df)[length(colnames(copy_voe_df))] <- variable
  }
  copy_voe_df=filter_unnest_feature_vib(copy_voe_df,feature)
  return(copy_voe_df)
}

get_all_adjusters <- function(vib_df,feature_num) {
  adjusters=vib_df %>% filter(feature==feature_num) %>% select(vars) %>% unlist %>% unique
  return(adjusters)
}

voe_summary <- function(df){
  estimates=quantile(voe_adjust$estimate,c(.01,.5,.99))
  pvals=quantile(-log10(voe_adjust$p.value),c(.01,.5,.99))
  estimate_diff=unname(estimates[3]-estimates[1])
  pvals_diff=unname(pvals[3]-pvals[1])
  num=nrow(df)
  output=unname(c(estimates,estimate_diff,pvals,pvals_diff,num))
  names(output)=c('estimate_1%','estimate_50%','estimate_99%','estimate_range','p.value_1%','p.value_50%','p.value_99%','p.value_range','number_of_models')
  return(output)
}

metaphlan_mapping=readRDS('metaphlan_mapping') %>% mutate(feature_num=paste('feature',num,sep='_'),datatype='species')
pathway_mapping=readRDS('pathways_mapping') %>% mutate(feature_num=paste('feature',num,sep='_'),datatype='pathways')
genefamily_map=readRDS('genefamilies_stool_mapping') %>% mutate(feature_num=paste('feature',num,sep='_'),datatype='genefamilies')

mapping=rbind(metaphlan_mapping,pathway_mapping,genefamily_map)

#files='101589_T2D_genefamilies_vibration_output.rds'
files=args[[1]]

summary_output=list()
estimate_models=list()
pvalue_models=list()
voe_list_for_reg=list()
for(f in files){
  print(f)
  voe=readRDS(f)
  voe=voe[unlist(lapply(voe$vars,function(x) !("dataset_name" %in% x))),]
  voe=voe[unlist(lapply(voe$vars,function(x) !("1" %in% x))),]
  voe=drop_na(voe)
  voe=voe[unlist(lapply(voe$feature_fit,function(x) length(x)!=2)),]
  features=unique(voe$feature)[!is.na(unique(voe$feature))]
  for(feat in features){
    print(feat)
    #generate summary statistics
    voe2= voe %>% filter(feature==feat)
    voe_adjust=get_adjuster_expanded_vibrations(voe2,feat) %>% drop_na()
    if(nrow(voe_adjust)==0){
      next
    }
    temp=strsplit(f,'_')[[1]][2]
    if(grepl('metaphlan',temp)==TRUE){
      dtype='species'
    }
    if(grepl('pathway',temp)==TRUE){
      dtype='pathways'
    }
    voe_adjust$phenotype=gsub('pathways-','',gsub('-meta-stats-to','',gsub('metaphlan-','',temp)))
    voe_list_for_reg[[feat]]=voe_adjust
    foo=strsplit(f,'_')[[1]][4][[1]][1]
    if(exists("dtype")==FALSE){
      dtype='genefamilies'
    }
    summary_name=mapping %>% filter(datatype==dtype, feature_num==feat) %>% select(feature_names) %>% as.character()
    sumout=voe_summary(voe_adjust)
    sumout[['feature']]=feat
    sumout[['datatype']]=dtype
    sumout[['id']]=paste(paste(strsplit(f,'_')[[1]][2:4],collapse='_'),'_',feat,sep='')
    sumout[['fraction_positive']]=(voe_adjust %>% filter(estimate>0) %>% nrow())/nrow(voe_adjust)
    summary_output[[summary_name]]=sumout
  }
}

summary_output=as.data.frame(do.call('rbind',summary_output),stringsAsFactors = FALSE)
summary_output=summary_output %>% rownames_to_column() %>%arrange(desc(estimate_range))
summary_output$estimate_range <- as.numeric(summary_output$estimate_range)
summary_output$p.value_range <-  as.numeric(summary_output$p.value_range)

write.csv(summary_output,paste('summary_output_vibrations',args[[1]],'.csv',sep=''))

#use this in the next step, where we regress confounder on effect size
voe_list_for_reg=bind_rows(voe_list_for_reg)
voe_list_for_reg[is.na(voe_list_for_reg)] <- 0
voe_list_for_reg= voe_list_for_reg %>% select(-vars,-full_fits)

write.csv(voe_list_for_reg,paste('voe_output_disease_regression',args[[1]],'.csv',sep=''))
