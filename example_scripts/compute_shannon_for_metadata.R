suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(cowplot))
suppressMessages(library(vegan)) # for diversity() function)
suppressMessages(library(ggplot2))
suppressMessages(library(stringr))

#ARG 1 = metadata, ARG 2 = output name

args = commandArgs(trailingOnly=TRUE)

metadata <- as_tibble(readRDS(args[[1]]))

df <- as_tibble(readRDS((args[[2]])), rownames = "feature_names") %>% select(feature_names, everything()) #keeping row names, making it a new col placed at the front

##### Computing diversity and genus count

# NOTE: data must be read in with feature_names (row names) kept in order for comptue_metaphlan_diversity() to work

# Computes shannon diversity using vegan::diversity() on  data in the form of dataframe or tibble 
# with sampleIDs as columns and species as rows.
# args: df, dataframe/tibble with sampleIDs as columns and microbial features as rows in 
# Returns: a ncol(df) x 4 tibble containing dataset name, sampleID, and shannon diversity inde, and genus count for each sample (col) from df
compute_diversity <- function(df) {
  no_feature_df <- df %>% select(-feature_names)
  genera_indices <- str_which(df$feature_names, "^k__Bacteria.+g__[^|]+$") # regex guarantees that we are selecting only Bacterial genera that end with genus
   if(length(genera_indices)==0){
    genera_indices=seq(nrow(df))
   }
  return(
    tibble(
      dataset_name = map_chr(names(no_feature_df), ~ unlist(str_split(string = ., pattern = "\\."))[[1]]), # getting dataset_name
      sampleID = map_chr(names(no_feature_df), ~ unlist(str_split(string = ., pattern = "\\:"))[[2]]), # getting sampleID
      shannon_diversity_idx = map_dbl(no_feature_df, ~ diversity(.)), # shannon diversity measure
      genus_count = map_int(no_feature_df, ~sum(.[genera_indices] > 0)) # genus richness by count
    )
  )
}


df_diversity <- compute_diversity(df)
merged_metadata_metaphlan_diversity <- left_join(metadata, df_diversity, by = c("dataset_name", "sampleID"))   # 7152 x 107, expected dimensions
write_rds(merged_metadata_metaphlan_diversity, args[[3]]) 

