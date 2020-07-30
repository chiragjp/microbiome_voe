suppressMessages(library(dplyr))
suppressMessages(library(purrr))
suppressMessages(library(tidyr))
suppressMessages(library(magrittr))
suppressMessages(library(broom))
suppressMessages(library(stringr) )
suppressMessages(library(rlang) )
suppressMessages(library(meta))
suppressMessages(library(metafor))



# # # Meta-analysis and Data combining

# Workflow:
# For both single-cohort and multi-cohort study conditions, 
# for each tibble of tidy fit outputs, filter out microbial feature stats and combine them across cohorts into 1 tibble
# For only multi-cohort study conditions, then map over all microbial features' tibbles and compute meta-analysis (metagen())

# takes in df outputted by compute_microbial_associations, currently of form where
# each row is grouped by dataset_name, and all data is nested, with microbial feature fits starting in the 4th col (first 3 are dataset_name, model_selected_df, and tidy_fits)
# goes through and filters each nested tibble of summary outputs for microbial features (since some fit with singularity and are removed)
# NOTE: features must currently start from col 4 as first 3 are ommitted
# NEW BEHAVIOR: error handling! now for any assocs that are not tibbles (e.g. error handled assocs passed on from comput_microbial_associations.R) will be labeled in this step as NA, and removed in get_combined_feature_fits()
get_filtered_associations <- function(df,ci) {
  copy <- duplicate(df, shallow = FALSE)
  copy[4:ncol(copy)] <- map(copy[4:ncol(copy)], function(x) map(x, function(y) return(tryCatch(filter(y, term == ci), # filtering based on existence of "study_condition" term
                                                                                               warning = function(w) NA_real_, # otherwise return NA which will be removed in the next helper
                                                                                               error = function(e) NA_real_
                                                                                               )
                                                                                             )
                                                                )
                            )
  return(copy)
}

get_combined_feature_fits <- function(df) {
  new_df <- tibble(analysis = "pre-meta-analysis-grouping") # create a new tibble to work with
  n <- 0
  for (i in c(seq(from = 4, to = ncol(df), by = 1))) { # for each feature column
    n <- n + 1
    print(paste0("iteration no: ", toString(n)))
    
    new_colname <- colnames(df %>% select(i)) # get the colum name
    new_df %<>% mutate(new = list(df %>% select(1, i) %>% filter(!is.na(eval(parse(text = new_colname)))) %>% unnest(cols=new_colname))) # unnest the dataset_name and feature columns, and make that a
                                                                                                                                  # a new column in the new tibble
    colnames(new_df)[length(colnames(new_df))] <- new_colname # correctly name the new column 
  }
  return(new_df %>% select(-1)) # remove the placeholder column used to create the tibble in line 1 of this function
}

get_metaanalysis_dfs <- function(df,ci) {
  filtered <- get_filtered_associations(df,ci)
  print("filtering successful")
  browser()
  combined <- get_combined_feature_fits(filtered)
  print("combining successful")
  return(combined)
}

compute_metaanalysis <- function(df) {
  new_df <- tibble(analysis = "meta-analysis") # create new tibble with placeholder column
  df %<>% select_if(~nrow(.[[1]]) > 0) # removing all features that yielded singularities for all cohorts. Failed regressions are caught in metagen() tryCatch() step next.
  n <- 0
  for (i in seq_along(df)) {
    new_colname <- colnames(df %>% select(i))
    n <- n + 1
    print(paste0("iteration no: ", toString(n), " ", new_colname))
    new_df %<>% mutate(new = list(tryCatch(metagen(estimate, 
                                                   std.error,
                                                   data = df[[i]][[1]],
                                                   studlab = dataset_name,
                                                   comb.fixed = FALSE,
                                                   comb.random = TRUE, # random effects model
                                                   method.tau = 'REML', # using REML estimator for tau heterogeneity parameter
                                                   hakn = FALSE, # not using that conserative estimator adjuster
                                                   prediction = TRUE,
                                                   sm = "SMD",
                                                   control=list(maxiter=1000)), # using 10x default max iteractions
                                           warning = function(w) w, # if warning return warning, don't want results at all to keep data and outputs as clean as possible (previously outputted list of 2, results and warnings)
                                           error = function(e) e))  # if results in error, return error
    )
    colnames(new_df)[length(colnames(new_df))] <- new_colname #properly name new column
  }
  return(new_df %>% select(-1)) # remove placeholder column
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  #args <- c("CRC_example_associations.rds", "CRC_example_meta_analysis.rds", "study_condition" )
  args <- c("Karlsson_associations.rds","2020_T2D_Data/Karlsson_metaanlysis.rds", "disease1" )
  association_outputs <- as_tibble(readRDS(args[[1]]))
  association_outputs <- association_outputs %>%
    filter(dataset == "Karlsson")
  #check if phenotype has multiple cohorts, if not don't meta-analyze and just filter out failed regressions
  run_metaanalysis = "FALSE"
    if(nrow(association_outputs)>1){
    run_metaanalysis = "TRUE"
  }
  outputname <- args[[2]]
  column_of_interest <-args[[3]]
  if (run_metaanalysis == "TRUE") {
    output <- compute_metaanalysis(get_metaanalysis_dfs(association_outputs,column_of_interest)) 
    saveRDS(output, 
            outputname)
  } else if (run_metaanalysis == "FALSE") {
    output <- get_metaanalysis_dfs(association_outputs,column_of_interest)
    saveRDS(output, 
            outputname)
  }
}

main()
