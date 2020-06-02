library(curatedMetagenomicData) 
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)
library(magrittr)
library(tidyverse
        
        )
library(broom)
library(stringr) 

# # Workflow:
# 1) Get names of datasets for a phenotype of interest, typically represented by study_condition and "control"
# 2) Create df of just those datasets, filtered for study_condition and controls only
#     Optional: Do brief exploration 
# 3) Nest data for each dataset to:
#     select vars,
#     prune df,
#     prune disease
#     prune control
#     model df
#     model selected df
#     tidied fits
#   Then unnest
#   Output both tibble of datasets, as well as tidied fits only

# # # Helper functions

# Selects variables with proportion of non-NAs >= threshold, returns selected df
# Args: df: typically a separated out and filtered (using filter_disease)  
# threshold: number between (0-1) indicating proportion
select_conserved <- function(df, threshold) {
  row_count <- nrow(df)
  var_counts <- df %>% 
    select(everything()) %>% 
    summarise_all(funs(list(sum(!is.na(.)))))
  selected_cons_vars_names <- as.data.frame(unlist(var_counts)) %>% rownames_to_column %>% filter(unlist(var_counts) >= (threshold * row_count)) %>% select(rowname) %>% unlist %>% unname
  df_selected <- df %>% select(selected_cons_vars_names)
  return(df_selected)
}

remove_vars <- function(df, vars) {
  return(df %>% select(-which(names(df) %in% vars)))
}


# For given dataset and study condition, filters only samples that are of the desired study_condition or a related intra-dataset control
# args: df, CMD type dataframe
# study_cond, string with study_condition name e.g. "CRC"
get_study_condition_datasets <- function(df, study_cond) {
  matches <- df %>% filter(study_condition == study_cond)
  dset_names <- unique(matches$dataset_name)
  all_study_cond_dsets <- df %>% filter(dataset_name %in% dset_names)
  return(all_study_cond_dsets %>% filter(study_condition == study_cond | study_condition == "control"))
}
# # test case:
# identical(get_study_condition_datasets(clean_metadata_metaphlan, "CRC"), all_crc_control)

get_study_condition_datasets_non_disease <- function(df, phenotype_column) {
  matches <- df %>% drop_na(phenotype_column)
  dset_names <- unique(matches$dataset_name)
  all_study_cond_dsets <- df %>% filter(dataset_name %in% dset_names)
  return(all_study_cond_dsets %>% filter(study_condition == "control"))
}

# Takes input df, converts all chr columns to factors, makes study_condition of interest 0/1 (1 if experimental 0 if control), and scales all numeric columns
# then moves study_condition to the first column
# args: df, CMD like dataframe
# disease, string rep. of study_condition 
convert_factors_and_scale_helper <- function(df, disease) {
  df$study_condition <- as.factor(map_int(df$study_condition, ~ifelse(. == disease, 1L, 0L))) # disease as 0/1 numerical factor, at end will turn back to integer 0/1 (R auto does this, but disease not guaranteed to get 1?)
  # print("encoding disease good")
  factored_df <- df %>% mutate_if(is_character, as.factor) # convert all chr vars to factors
  # print("converting to factors good")
  final_df <- factored_df %>% mutate_if(is.numeric, scale) # scale numeric vars, but not disease
  # print("scaling good")
  final_df$study_condition <- as.integer(as.character(final_df$study_condition)) # convert disease factor back to 0/1 integer
  return(final_df %<>% select(study_condition, everything())) # move disease to the front
}


# For CMD dataframe with 1 or more dataset_names, nests by dataset_name, and for each dataset_name, selects/prunes variables based on user-defined threshold
# args: df, CMD like dataframe (can contain variables that want to remove, just specify which in vars_to_remove)
# generally, this df should be an output of get_study_condition_datasets() which contains data only fora specific study_condition and its related controls
# vars_to_remove, character vector specifying variables to remove from df
# threshold, [0,1] representing proportion of non-na's in column needed to select a variable
# output: nested tibble with dataset_name, data, clean_data (data but with specified vars removed), and selected variable names
nest_compute_helper <- function(df, threshold, vars_to_remove) {
  return(df %>% 
           nest(-dataset_name) %>% 
           mutate(
             clean_data = map(data, ~remove_vars(., vars_to_remove)),
             pruned_vars = map(clean_data, function(x) return(names(select_conserved(x, threshold)))
             )
           )
  )
}

# # test case:
# identical(nest_compute_helper(get_study_condition_datasets(clean_metadata_metaphlan, "CRC"), 0.8), all_crc_control_selected)

# Given a CMD like dataframe, for each dataset nests all data, cleans the data via remove_vars(),
# select vars, prune df, prune disease, prune control, model df, model selected df (further selects out vars that will be problematic when fitting), tidied fits 
# -- maintaining all in a tibble (example of nest, map, unnest workflow)
# Step 1: generate a df with only study_condition of interest and corresponding controls (via get_study_condition_datasets())
# Step 2: nest and generate clean data (remove uninteresting vars) and set of conserved variables based on threshold
# Step 3: compute all remaining data, prunes, and fits listed
# args: df, CMD like dataframe (can contain variables that want to remove, just specify which in vars_to_remove)
# generally, this df should be an output of get_study_condition_datasets() which contains data only fora specific study_condition and its related controls
# study_cond, string rep. of study_condition of interest
# threshold, [0,1] representing proportion of non-na's in column needed to select a variable
# vars_to_remove, character vector specifying variables to remove from df
nest_compute_CMD <- function(df, study_cond, threshold, vars_to_remove) {
    filtered_studycond_df <- get_study_condition_datasets(df, study_cond)
    nested_selected <- nest_compute_helper(filtered_studycond_df, threshold, vars_to_remove)
    return(nested_selected %>% 
             mutate(pruned_data = map(row.names(nested_selected), function(x) return(select(nested_selected[[as.numeric(x), 3]], nested_selected[[as.numeric(x), 4]]) %>% 
                                                                                       drop_na())), # selects vars in data based on pruned_vars, and removes any rows with na
                    model_df = map(pruned_data, function(x) return(convert_factors_and_scale_helper(x, study_cond))), # convert and scale,
                    model_selected_df = map(model_df, function(x) return(x %>%
                                                                           select_if(~!anyNA(.)) %>% # select out any vars that failed conversion and scale (generally bc of only having 1 value repeated)
                                                                           select_if(~ (length(levels(unlist(.))) > 1) | (is.numeric(.) & n_distinct(.) > 1)) # select out any categorical vars that only have 1 level
                    )
                    ),
                    tidy_fits = map(model_selected_df, function(x) tidy(stats::glm(study_condition ~ . - sampleID - subjectID, # fit logistic regression model and get tidy outputs, excluding ID vars (sampleID and subjectID)
                                                                                   data = x,
                                                                                   family = binomial(link = 'logit'),
                                                                                   control = glm.control(maxit=50)) # doubling max iterations to guarantee convergence
                    )
                    
                    )
             )
    )
}

# # # Function for interacting from the command line 
# (simply call Rscript prune_metadata.R [metadatafilename e.g. cmd_all_metadata.rds] [phenotype e.g. "CRC"] [threshold e.g. 0.9] [output filename e.g. "crc_modeldfs.rds"])
main <- function(vars_to_remove) {
  args <- commandArgs(trailingOnly = TRUE)
  cmd_file <- as_tibble(readRDS(args[[1]]))
  phenotype <- args[[2]]
  threshold <- as.double(args[[3]])
  outputname <- args[[4]]
  output <- nest_compute_CMD(cmd_file, phenotype, threshold, vars_to_remove)
  print(str_c("Full metadata pruning process output for ", phenotype, ".", " Additionally add/remove variables from maximal model as necessary."))
  print(output)
  saveRDS(output %>% select(dataset_name, model_selected_df, tidy_fits), 
          outputname)
}
to_remove2 <- c("PMID", "NCBI_accession", "number_bases", "curator", 
                "sequencing_platform", "DNA_extraction_kit", 
                "body_site", "disease", "country", "non_westernized", "location",
                # below are likely culprits of across the board 0/1 fitted probs
                "minimum_read_length", "median_read_length", "antibiotics_current_use")  #including sampleID and subjectID, just excluding in model formula
main(to_remove2)

