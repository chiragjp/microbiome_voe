library(dplyr)
library(purrr)
library(tidyr)
library(magrittr)
library(broom)
library(stringr) # regex pattern matching
library(rlang) # rlang::duplicate()

args = commandArgs(trailingOnly=TRUE)

# purpose of this script is to prep an abundance file (e.g. metaphlan abundances for regression)
# it transposes the data and build a mapping between row number and microbial feature name (as tibbles tend to throw out rownames)

df <- as_tibble(readRDS(args[[1]])) #throwing out row names for now, mapping will come back
df = df %>% select(-c(feature_names))
metaphlan_mapping <- as_tibble(readRDS(args[1])) %>% select(feature_names)
metaphlan_mapping %<>% mutate(num = row_number())


# Takes in CMD metaphlan data file, df, gets dataset_name and sampleID from column names
# and outputs a 2 x ncol(df) tibble with dataset name and sample ID as rows, ready for cbind remaining data
# args: df, CMD metaphlan data file loaded in as a dataframe
generate_name_sampleID_df <- function(df) {
  return(
    tibble(
      dataset_name = map_chr(names(df), ~ unlist(str_split(string = ., pattern = "\\."))[[1]]),
      sampleID = map_chr(names(df), ~ unlist(str_split(string = ., pattern = "\\:"))[[2]])
    )
  )
}

# # Takes metaphlan dataframe, and a list of rows (microbial features in metaphlan)
# # to eventually merge with metadata, generates base name + sampleID df, cbinds with the name_sampleID_df, then merges with the metadata
# # args: df, metaphlan dataframe
# # rows, vector of row numbers to take and perform this transpose + merge operation
# # outputs: ncol(df) x (2 + length(rows)) df with each row being a sample beginning with dataset_name and sampleID
t_metaphlan <- function(df, rows) {
  base_df <- generate_name_sampleID_df(df)
  base_df <- cbind(base_df, t(slice(df, rows)))
  return(as_tibble(base_df))
}

t_df <- t_metaphlan(df, c(seq(from = 1, to = nrow(metaphlan), by = 1)))

saveRDS(t_df, args[2]) 
saveRDS(metaphlan_mapping, args[3])

