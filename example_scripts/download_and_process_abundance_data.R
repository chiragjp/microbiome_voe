suppressMessages(library(curatedMetagenomicData))
suppressMessages(library(tidyverse))

args = commandArgs(trailingOnly=TRUE)

######first, download metadata and raw data

crc_datasets = combined_metadata %>% filter(study_condition=='CRC') %>% select(dataset_name) %>% distinct %>% unname %>% unlist
combined_metadata_subset = combined_metadata %>% filter(dataset_name %in% crc_datasets)

cmd_data=curatedMetagenomicData('*metaphlan*',dryrun=T)
cmd_data=unlist(unname(sapply(crc_datasets, function(x) cmd_data[grepl(x, cmd_data)])))
cmd_data=curatedMetagenomicData(cmd_data,dryrun=F)
cmd_data=mergeData(cmd_data)
cmd_data=exprs(cmd_data)

#the subset data loaded here contains F. nucleatum abundance, which we know will be positively associated with CRC
cmd_data_subset_for_example = cmd_data[readRDS('features_to_sample.rds'),]

######second, compute shannon diversity 

#some of these packages need to be loaded after CMD download is completed due to function name conflicts
suppressMessageslibrary(dplyr))
suppressMessageslibrary(purrr))
suppressMessageslibrary(tidyr))
suppressMessageslibrary(magrittr))
suppressMessageslibrary(broom))
suppressMessageslibrary(stringr)) # regex pattern matching
suppressMessageslibrary(rlang)) # rlang::duplicate()
suppressMessageslibrary(vegan))

compute_diversity <- function(df) {
  no_feature_df <- df %>% select(-feature_names)
  genera_indices <- str_which(df$feature_names, "^k__Bacteria.+g__[^|]+$") # regex guarantees that we are selecting only Bacterial genera that end with genus
   if(length(genera_indices)==0){
    genera_indices=seq(nrow(df))
   }
  return(
    tibble(
      dataset_name = map_chr(names(no_feature_df), ~ unlist(str_split(string = ., pattern = "\\."))[[1]]), # getting dataset_name
      sampleID = map_chr(names(no_feature_df), ~ unlist(gsub('\\.','-',str_split(string = ., pattern = "\\.stool\\.")[[1]][2]))), # getting sampleID
      shannon_diversity_idx = map_dbl(no_feature_df, ~ diversity(.)), # shannon diversity measure
      genus_count = map_int(no_feature_df, ~sum(.[genera_indices] > 0)) # genus richness by count
    )
  )
}

cmd_data_subset_for_example_diversity <- compute_diversity(data.frame(cmd_data_subset_for_example) %>% rownames_to_column('feature_names'))
merged_metadata_metaphlan_diversity <- left_join(combined_metadata_subset, cmd_data_subset_for_example_diversity, by = c("dataset_name", "sampleID"))   
saveRDS(merged_metadata_metaphlan_diversity, args[[1]]) 

# purpose of the next portion of this script is to prep an abundance file (e.g. metaphlan abundances for regression)
# it transposes the data and build a mapping between row number and microbial feature name (as tibbles tend to throw out rownames)

mapping <- data.frame(cmd_data_subset_for_example)  %>% rownames_to_column('feature_names') %>% select(feature_names)
mapping  = data.frame(mapping) %>% mutate(num = row_number())
cmd_data_subset_for_example = as_tibble(cmd_data_subset_for_example)

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

# # Takes abundance data dataframe, and a list of rows (microbial features in metaphlan)
# # to eventually merge with metadata, generates base name + sampleID df, cbinds with the name_sampleID_df, then merges with the metadata
# # args: df, metaphlan dataframe
# # rows, vector of row numbers to take and perform this transpose + merge operation
# # outputs: ncol(df) x (2 + length(rows)) df with each row being a sample beginning with dataset_name and sampleID
t_metaphlan <- function(df, rows) {
  base_df <- generate_name_sampleID_df(df)
  base_df <- cbind(base_df, t(slice(df, rows)))
  return(as_tibble(base_df))
}

t_cmd_data_subset_for_example <- t_metaphlan(cmd_data_subset_for_example, c(seq(from = 1, to = nrow(cmd_data_subset_for_example), by = 1)))

saveRDS(t_cmd_data_subset_for_example, args[[2]]) 
saveRDS(mapping, args[[3]])
