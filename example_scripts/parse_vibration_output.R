suppressMessages(library(tidyverse))
suppressMessages(library(stringr))
suppressMessages(library(cowplot))
suppressMessages(library(ggplot2))
suppressMessages(library(tidyr))
suppressMessages(library(magrittr))
suppressMessages(library(broom))
suppressMessages(library(meta))
suppressMessages(library(stringr) )
suppressMessages(library(rlang) )
suppressMessages(library(lme4))
suppressMessages(library(lmerTest))
suppressMessages(library(broom.mixed))



theme_set(theme_cowplot())

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
  estimates=quantile(df$estimate,c(.01,.5,.99))
  pvals=quantile(-log10(df$p.value),c(.01,.5,.99))
  estimate_diff=unname(estimates[3]-estimates[1])
  pvals_diff=unname(pvals[3]-pvals[1])
  num=nrow(df)
  output=unname(c(estimates,estimate_diff,pvals,pvals_diff,num))
  names(output)=c('estimate_1%','estimate_50%','estimate_99%','estimate_range','p.value_1%','p.value_50%','p.value_99%','p.value_range','number_of_models')
  return(output)
}

find_confounders <- function(voe_list_for_reg,ptype,cohorttype){
  voe_adjust_for_reg_ptype <- voe_list_for_reg %>% select(-dataset_name) %>%select_if(~ length(unique(.)) > 1)
  voe_adjust_for_reg_ptype[is.na(voe_adjust_for_reg_ptype)]=0
  voe_adjust_for_reg_ptype$estimate=abs(voe_adjust_for_reg_ptype$estimate)
  if(cohorttype=='multi'){
    fit_estimate=lmer(data=voe_adjust_for_reg_ptype,as.formula(estimate ~ . +(1|feature)-feature - estimate  - std.error - p.value  - statistic),control = lmerControl(optimizer = "bobyqa"))
  }
  else{
    fit_estimate=lmer(data=voe_adjust_for_reg_ptype,as.formula(estimate ~ . +(1|feature) -feature - estimate - std.error - p.value - statistic),control = lmerControl(optimizer = "bobyqa"))
  }

  fit_estimate_forplot=tidy(fit_estimate) %>% mutate(sdmin=(estimate-std.error),sdmax=(estimate+std.error))
  outname=paste('example_pipeline_output/disease_specific_confounders_estimate_',ptype,".csv",sep='')
  write.csv(fit_estimate_forplot,outname)
}

plot_voe <- function(feature_vib_df, feature_num,title,fdr,y,x) {
  if(nrow(feature_vib_df)==0){
    return('')
  }
  plot <- ggplot(data = feature_vib_df, aes(x = estimate, y = -log10(p.value))) +geom_point(alpha=.5,aes(color = dataset_name)) +geom_hline(yintercept = -log10(fdr),linetype='dashed') +theme(legend.position="bottom") +theme(legend.title = element_blank()) +theme(plot.title = element_text(size=8))+ xlab('Estimate') +xlim(-3,3)+ theme(legend.text=element_text(size=6))+ggtitle(title)+ theme(plot.title = element_text(size=12)) + geom_hline(yintercept = -log10(0.05)) + geom_point(aes(x=x, y=-log10(y)),shape=23, fill="blue", color="darkred", size=3)
  return(plot)
}


process_voe <- function(f,meta_analysis_output){
  summary_output=list()
  voe_list_for_reg=list()
    print(f)
    voe=readRDS(f)
    meta_analysis_output=readRDS(meta_analysis_output)
    voe=drop_na(voe)
    voe=voe[unlist(lapply(voe$feature_fit,function(x) length(x)!=2)),]
    features=unique(voe$feature)[!is.na(unique(voe$feature))]
    for(feat in features){
      print(feat)
      #generate summary statistics
      voe2= voe %>% filter(feature==feat)
      voe_adjust=get_adjuster_expanded_vibrations(voe2,feat) %>% drop_na()
      voe_list_for_reg[[feat]]=voe_adjust
      plot_title = mapping %>% filter(feature_num==feat) %>% select(feature_names) %>% as.character()
      plot_title=paste(strsplit(plot_title,'\\|')[[1]][length(strsplit(plot_title,'\\|')[[1]])],'vibration output by dataset')
      max_p_value_still_sig = meta_analysis_output %>% filter(by.p.val<=.05) %>% select(p.val)%>% max
      ma_pval= meta_analysis_output %>% filter(feature==feat) %>% select(p.val) %>% unname %>% unlist %>% as.numeric
      ma_estimate =  meta_analysis_output %>% filter(feature==feat) %>% select(estimate) %>% unname %>% unlist %>% as.numeric
      outplot=plot_voe(voe_adjust,feat,plot_title,max_p_value_still_sig,ma_pval,ma_estimate)
      ggsave(paste('example_pipeline_output/vibration_plot_for_',plot_title,'_',feat,'.png',sep=''))
      if(nrow(voe_adjust)==0){
        next
      }
      summary_name=mapping %>% filter(feature_num==feat) %>% select(feature_names) %>% as.character()
      sumout=voe_summary(voe_adjust)
      sumout[['feature']]=feat
      sumout[['id']]=paste(paste(strsplit(f,'_')[[1]][2:4],collapse='_'),'_',feat,sep='')
      sumout[['fraction_positive']]=(voe_adjust %>% filter(estimate>0) %>% nrow())/nrow(voe_adjust)
      summary_output[[summary_name]]=sumout
    }
 
    summary_output=as.data.frame(do.call('rbind',summary_output),stringsAsFactors = FALSE)
    summary_output=summary_output %>% rownames_to_column() %>%arrange(desc(estimate_range))
    summary_output$estimate_range <- as.numeric(summary_output$estimate_range)
    summary_output$p.value_range <-  as.numeric(summary_output$p.value_range)
    write.csv(summary_output,paste('example_pipeline_output/summary_output_vibrations_',gsub('rds','csv',args[[1]]),sep=''))

    voe_list_for_reg=bind_rows(voe_list_for_reg)
    voe_list_for_reg[is.na(voe_list_for_reg)] <- 0
    voe_list_for_reg= voe_list_for_reg %>% select(-vars,-full_fits)
    find_confounders(voe_list_for_reg,args[[3]],args[[4]])
    write.csv(summary_output,paste('example_pipeline_output/raw_voe_adjuster_data_for_confounder_regression_',gsub('rds','csv',args[[1]]),sep=''))

    return('Finished parsing VoE output. You can find the results in the `example_pipeline_output` folder.')

}


file=args[[1]]
meta_analysis_output=args[[2]]
mapping=readRDS(args[[3]]) %>% mutate(feature_num=paste('feature',num,sep='_'),datatype='species')

process_voe(file,meta_analysis_output)




