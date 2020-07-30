
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(purrr))
suppressMessages(library(tidyr))
suppressMessages(library(magrittr))
suppressMessages(library(broom))
suppressMessages(library(stringr)) # regex pattern matching)
suppressMessages(library(rlang)) # rlang::duplicate())
suppressMessages(library(tidyverse))
suppressMessages(library(rje)) # for function powerSet())


# Workflow:
# 1) get features want to vibrate over, stored as a vector in an RDS (from and lit review)
# 2) get maximal model dataset (modeldfs from prune_metadata.R)
# 3) for each desired feature and for each dataset/cohort, perform vibrations
# 4) Aggregate outputs & visualize (in analyze_vibrations.R)

# Finds smallest non-zero value in input df
get_fudge_factor <- function(metaphlan_df) {
  datatype_df <- metaphlan_df %>% select_if(is.numeric)
  return(min(map_dbl(map(datatype_df, ~.[. > 0]), ~min(.))))
}
# print(get_fudge_factor(t_metaphlan_df))


# Vibrates over a single maximal model for a single feature
vibrate <- function(model_df, feature.num, dataset.name, metaphlan_df,pheno) {
  fudge_factor <- get_fudge_factor(metaphlan_df) 
  new_modeldf <- data.frame(model_df) %>% mutate(dataset_name = dataset.name) 
  feature_num <- as.integer(str_split(feature.num, "_")[[1]][[2]])
  regression_df <- left_join(new_modeldf %>% mutate_if(is.factor, as.character), metaphlan_df %>% mutate(dataset_name=str_replace(dataset_name, "gene_families_", "")) %>%select(c(1,2, (feature_num + 2))), by = c("dataset_name", "sampleID")) %>% select(-dataset_name) %>% mutate_if(is.character, as.factor)
  ################averaging logic
  subnum=length(unique(regression_df$subjectID))
  samnum=length(unique(regression_df$sampleID))
  if(subnum!=samnum){
      #get case subjects and remove ones that are pre disease (this can be changed if we want to compare health controls to pre disease individuals)
      cases=regression_df %>% filter(study_condition==1) %>% select(subjectID) %>% mutate_if(is.factor, as.character) %>% unique() %>% unlist %>% unname()
      regression_df_cases = regression_df %>% filter(subjectID %in% cases) %>% filter(study_condition != 0) %>% select(-c(as.character(feature_num)))
      case_samples=regression_df_cases %>% select(sampleID) %>% mutate_if(is.factor,as.character) %>% unlist %>% unname
      #average 
      regression_df_cases_numeric_averaged = regression_df_cases %>%group_by(subjectID) %>% select_if(is.numeric) %>% summarise_each(mean)
      regression_df_cases_character_averaged = regression_df_cases %>% select_if(negate(is.numeric)) %>% select(-sampleID) %>% distinct() %>% mutate_if(is.factor, as.character)
      subjects=unique(regression_df_cases_character_averaged$subjectID)
      toremove=list()
      for(s in subjects){
        regression_df_cases_character_averaged_sub = regression_df_cases_character_averaged %>% filter(subjectID==s)
        temp=apply(regression_df_cases_character_averaged_sub, 2, function(x) length(unique(x)))
        temp=names(temp[temp>1])
        if(length(temp)>0){
          toremove[s]=temp
        }            
      }
      toremove=unique(unname(unlist(toremove)))
      if(length(toremove)>0){
        regression_df_cases_character_averaged=regression_df_cases_character_averaged %>% select(-c(toremove)) %>% distinct()
      }
      cases_df = full_join(regression_df_cases_numeric_averaged,regression_df_cases_character_averaged,on=subjectID) %>% mutate_if(is.character,as.factor)
      
      #get control subjects
      controls=regression_df %>% filter(study_condition==0) %>% select(subjectID) %>% mutate_if(is.factor, as.character) %>% unique() %>% unlist %>% unname()
      regression_df_controls = regression_df %>% filter(subjectID %in% controls) %>% filter(study_condition != 1) %>% select(-c(as.character(feature_num)))
      control_samples=regression_df_controls %>% select(sampleID) %>% mutate_if(is.factor,as.character) %>% unlist %>% unname
      regression_df_controls_numeric_averaged = regression_df_controls %>%group_by(subjectID) %>% select_if(is.numeric) %>% summarise_each(mean)
      
      regression_df_controls_character_averaged = regression_df_controls %>% select_if(negate(is.numeric)) %>% select(-sampleID) %>% distinct() %>% mutate_if(is.factor, as.character)
      subjects=unique(regression_df_controls_character_averaged$subjectID)
      toremove=list()
      for(s in subjects){
        regression_df_controls_character_averaged_sub = regression_df_controls_character_averaged %>% filter(subjectID==s)
        temp=apply(regression_df_controls_character_averaged_sub, 2, function(x) length(unique(x)))
        temp=names(temp[temp>1])
        if(length(temp)>0){
          toremove[s]=temp
        }
      }
      toremove=unique(unname(unlist(toremove)))
      if(length(toremove)>0){
        regression_df_controls_character_averaged=regression_df_controls_character_averaged %>% select(-c(toremove)) %>% distinct()     
      }
      controls_df = full_join(regression_df_controls_numeric_averaged,regression_df_controls_character_averaged,on=subjectID) %>% mutate_if(is.character,as.factor)
      
      #average across microbial feature
      regression_df_mf = regression_df %>% filter(sampleID %in% case_samples | sampleID %in% control_samples) %>% select(subjectID,as.character(feature_num)) %>% group_by(subjectID) %>% summarise_all(funs(mean=mean,sd=sd))            
      #write.csv(regression_df_mf,paste('averaging_output/',phenotype,'_',datatype,'_',new_feature_name,'.csv',sep=''))
      regression_df_mf=regression_df_mf %>% select(subjectID,mean)
      colnames(regression_df_mf)[2]=feature_num
      
      #merge
      regression_df = bind_rows(cases_df,controls_df) %>% mutate_if(is.character,as.factor)
      regression_df=full_join(regression_df,regression_df_mf,on='subjectID')
  }
  else{
    regression_df %<>% select(-sampleID)
  }
  
  model_with_feature_col = regression_df %>% select(-subjectID)
  
    varset=powerSet(colnames(regression_df %>% select(-subjectID, -study_condition, -c(as.character(feature_num)))))
    if(length(varset)>50000){
      varset=sample(varset,50000)
    }
    return(tibble(
      dataset_name = dataset.name,
      feature = feature.num,
      vars = varset,
      full_fits = map(vars, function(x) tryCatch(tidy(stats::lm(as.formula(str_c("I(log10(`", str_split(feature.num, "_")[[1]][[2]], "` + ", toString(fudge_factor), ")) ~ ",paste(x,'study_condition',collapse='+',sep='+'))), 
                                                                data = model_with_feature_col
      ) 
      ),
      warning = function(w) w,
      error = function(e) e
      )
      ),
      feature_fit = map(full_fits, function(x) tryCatch(filter(x, term == "study_condition"),
                                                        warning = function(w) w, 
                                                        error = function(e) e)
      )
    )
    )
}

compute_vibrations <- function(disease_all_modeldfs, all_feature_nums, metaphlan_df, pheno) {
  print('Computing vibrations for:')
  print(unlist(unname(all_feature_nums)))
  reduce(  # map over all feature's want to vibrate for
    map(all_feature_nums, function(x)
      reduce( # get single feature's vibration across all cohorts
        map(c(1:nrow(disease_all_modeldfs)), function(y) vibrate(disease_all_modeldfs[[y,2]], 
                                                                 x, 
                                                                 disease_all_modeldfs[[y,1]],
                                                                 metaphlan_df,pheno)),
        rbind #combine all cohorts for one feature
      )
    ),
    rbind, #combine all features
    .init = NA_real_ # .init is supplied as the first value to start the accumulation in reduce, as o.w. reduce with throw error for empty starting value
  )
}

main <- function() {
  
  args <- commandArgs(trailingOnly = TRUE)
  model_dfs <- as_tibble(readRDS(args[[1]]))
  features <- readRDS(args[[2]])
  t_data_type <- as_tibble(readRDS(args[[3]]))
  outputname <- args[[4]]
  pheno <- args[[5]]
  output <- compute_vibrations(model_dfs, features, t_data_type,pheno)
  saveRDS(output, outputname)
}

main()