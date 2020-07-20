suppressMessages(library(dplyr))
suppressMessages(library(purrr))
suppressMessages(library(tidyr))
suppressMessages(library(magrittr))
suppressMessages(library(broom))
suppressMessages(library(stringr))



# Workflow
# For each microbial feature (col in data type's tibble): 
#   For each phenotype's modeldfs:
#     For each dataset:
#     1) add dataset_name to all tibbles
#     2) left_join matching ID's abundances using dataset_name and IDs
#     3) If proportion of non-zero values for abundances is at least 1,
#          mutate new col with tidied glm fit summary on model df that now has 1 feature tacked on
#     4) rename column as desired

# Performs above workflow on df outputted by cmd_metadata_metaphlan_v2.Rmd, and transposed, separated metaphlan df outputted
# by same Rmd file.


# For a single feature and single modeldf outputted by prune_metadata.R, computes the proportion of non-zero values per cohort.
compute_microbial_feature_proportions_helper <- function(model_ready_df, metaphlan_df, feature_num) {
  return(
    length(
      which((left_join(model_ready_df %>% mutate_if(is.factor, as.character),
                       metaphlan_df %>% mutate(dataset_name=str_replace(dataset_name, "gene_families_", "")) %>% select(c(1,2, (feature_num + 2))), 
                       by = c("dataset_name", "sampleID")))[[length(model_ready_df) + 1]] > 0)
    ) / nrow(model_ready_df)
  )
}

compute_microbial_feature_proportions <- function(df, metaphlan_df){
  # adding dataset name to nested tibble
  df$model_selected_df <- map(c(seq(from = 1, to = nrow(df), by = 1)), function(x) as.data.frame(df[[x,2]]) %>% mutate(dataset_name = df[[x,1]]))
  for (j in seq_along(metaphlan_df %>% select(-1, -2))) {
    # adding 1 feature, fitting, and tidying summary
    new_feature_name <- colnames(metaphlan_df %>% select(-1, -2) %>% select(j))
    df %<>% mutate(
      new_feature = map_dbl(df$model_selected_df, function(x)
        compute_microbial_feature_proportions_helper(x, metaphlan_df, j)
      )
    )
    colnames(df)[length(colnames(df))] <- paste0("feature_", toString(new_feature_name)) # last step to rename feature to desired variable input name
  }
  return(df)
}

# Takes model_ready_df (output of prune_metadata.R) and a transformed data_type file (outputted from pre_data_for_pipieline.R)
# and computes associations one-by-one and outputs in a tibble.
compute_lm_associations <- function(df, metaphlan_df,phenotype){
  # find fudge factor for logging
  fudge_factor <- min((metaphlan_df %>% select(-1, -2) %>% unlist())[which((metaphlan_df %>% select(-1, -2) %>% unlist()) > 0)]) 
  # clean dataset name folder if working with gene families
  #iterate through each cohort for each feature
  
  #This line just cleans the DF from dataset names that aren't in the metadata but are in the CAGS (we shouldn't have that issue)
  #df=df %>% filter(dataset_name %in% unique(gsub('gene_families_','',metaphlan_df$dataset_name)))
  df$model_selected_df <- map(c(seq(from = 1, to = nrow(df), by = 1)), function(x) data.frame(df[[x,2]]) %>% mutate(dataset_name = df[[x,1]]))
  for (j in seq_along(metaphlan_df %>% select(-1, -2))) {
    # adding 1 feature, fitting, and tidying summary
    new_feature_name <- colnames(metaphlan_df %>% select(-1, -2) %>% select(j))
    print(str_c("feature no.: ", new_feature_name)) # print feature no. currently working on
    
    df %<>% mutate(new_feature = map(df$model_selected_df, function(x)
      if (compute_microbial_feature_proportions_helper(x, metaphlan_df, j) > 0) { # prorportion must be non-zero or else no response to model
        
        regression_df=left_join(as.data.frame(x) %>% mutate_if(is.factor, as.character), metaphlan_df %>% mutate(dataset_name=str_replace(dataset_name, "gene_families_", "")) %>% select(c(1,2, (j + 2))),by = c("dataset_name", "sampleID")) %>% select(-dataset_name) %>% mutate_if(is.character, as.factor) %>% drop_na()
        
        ################averaging logic
        subnum=length(unique(regression_df$subjectID))
        samnum=length(unique(regression_df$sampleID))
        #check if the number of cases vs subjects differs and if so, average across samples, dropping columns that have multiple values for the same subject or, as a result of the averaging, contain only 1 value.
        if(subnum!=samnum){
          #get case subjects and remove ones that are pre disease (this can be changed if we want to compare health controls to pre disease individuals)
          cases=regression_df %>% filter(study_condition==1) %>% select(subjectID) %>% mutate_if(is.factor, as.character) %>% unique() %>% unlist %>% unname()
          regression_df_cases = regression_df %>% filter(subjectID %in% cases) %>% filter(study_condition != 0) %>% select(-new_feature_name)
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
            toremove=unique(unname(unlist(toremove)))
            if(length(toremove)>0){
              regression_df_cases_character_averaged=regression_df_cases_character_averaged %>% select(-c(toremove)) %>% distinct()
            }
            cases_df = full_join(regression_df_cases_numeric_averaged,regression_df_cases_character_averaged,on=subjectID) %>% mutate_if(is.character,as.factor)
            
            #get control subjects
            controls=regression_df %>% filter(study_condition==0) %>% select(subjectID) %>% mutate_if(is.factor, as.character) %>% unique() %>% unlist %>% unname()
            regression_df_controls = regression_df %>% filter(subjectID %in% controls) %>% filter(study_condition != 1) %>% select(-new_feature_name) 
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
            regression_df_mf = regression_df %>% filter(sampleID %in% case_samples | sampleID %in% control_samples) %>% select(subjectID,new_feature_name) %>%group_by(subjectID) %>% summarise_all(funs(mean=mean,sd=sd))            
            #write.csv(regression_df_mf,paste('averaging_output/',phenotype,'_',datatype,'_',new_feature_name,'.csv',sep=''))
            regression_df_mf=regression_df_mf %>% select(subjectID,mean)
            colnames(regression_df_mf)[2]=new_feature_name
            
            #merge
            regression_df = bind_rows(cases_df,controls_df) %>% mutate_if(is.character,as.factor)
            regression_df=full_join(regression_df,regression_df_mf,on='subjectID')
          }
        }
        else{
          regression_df %<>% select(-sampleID)
        }
        
        regression_df %<>% select(-subjectID)
        #remove HLA variable from initial maximal models
        tryCatch({
          regression_df %<>% select(-HLA)
          
        },error=function(e)
        {
          regression_df=regression_df
        })
        #run regression
        regression_df=regression_df %>% drop_na %>% select_if(~ length(unique(.)) > 1)
        
        # Added to change column name study_condition of example data to "disease" so it matches iyr data 
        regression_df <- regression_df %>% 
          rename(
            disease = study_condition
          )
        return(tryCatch(tidy(stats::lm(as.formula(str_c("I(log10(`", new_feature_name, "` + ", toString(fudge_factor), ")) ~ disease")), # fit lm and get tidy outputs, excluding ID vars (sampleID and subjectID)
                                       data = regression_df))
                        
                        ,
                        # NOTE: tidy() ends up removing and singularity fits (e.g. NA fits from a feature that, for the dataset's samples all have 0 relative abundance)
                        warning = function(w) w, # if warning or error, just return them instead of output
                        error = function(e) e
        )
        )
      } else {
        return(simpleError(str_c("No non-zero values for ", new_feature_name, " to be considered in ", x$dataset_name[[1]], "."))
        )
      }
    )
    )
    colnames(df)[length(colnames(df))] <- paste0("feature_", toString(new_feature_name)) # last step to rename feature to desired variable input name
  }
  return(df)
}

main <- function() {
  #2020 T2D Data
  #args <- c('2020_T2D_Data/metadata.rds', '2020_T2D_Data/t_metaphlan_abundance_cmd_example_CRC.rds', 'CRC_example_associations.rds','CRC')
  
  # Example Data 
  args <- c('Example_Data/CRC_example_modeldfs.rds', 'Example_Data/t_metaphlan_abundance_cmd_example_CRC.rds', 'CRC_example_associations.rds','CRC')
  
  cmd_file <- as_tibble(readRDS(args[[1]]))
  datatype_file <- as_tibble(readRDS(args[[2]]))
  outputname <- args[[3]]
  phenotype=args[[4]]
  output <- compute_lm_associations(cmd_file, datatype_file,phenotype)
  saveRDS(output, outputname)
}

main()
