suppressMessages(library(tidyverse))
suppressMessages(library(readxl))

# The combined Metadata Excel sheet should:
#   1. Have a column called "dataset_name" with a reference to the researcher or name of the study associated
#   2. Have a column named "disease" (where 0 indicates control or lack of T2D, 1 indicates prediabetes, and 2 indicates diabetes)


main <- function() {
  # Replace the arguments with the names of your files (save the files in the same directory)
  args <- c('combinedMetadata.xlsx', 't_metaphlan_abundance_cmd_example_CRC.rds', 'CRC_example_associations.rds','CRC')
  
  # * PREPARE METADATA * 
  metadata_df <- read_excel(args[[1]], sheet = 2)
  # Build format of the summary tibble that will be used as input for 'compute associations.R'
  outputDF <- metadata_df %>% nest(-dataset_name) 
  saveRDS(outputDF, file = "metadata.rds")
  
  # * PREPARE CAGS * 
  # cags_df <- read_csv(args[[2]])
  # cagsDF <- readRDS('t_metaphlan_abundance_cmd_example_CRC.rds')
  # cagsDF$dataset_name <- gsub('FengQ_2015', 'Qin', cagsDF$dataset_name)
  # cagsDF$dataset_name <- gsub('HanniganGD_2017', 'Qin', cagsDF$dataset_name)
  # cagsDF$dataset_name <- gsub('ThomasAM_2018a', 'HMP', cagsDF$dataset_name)
  # cagsDF$dataset_name <- gsub('ThomasAM_2018b', 'HMP', cagsDF$dataset_name)
  # cagsDF$dataset_name <- gsub('VogtmannE_2016', 'HMP', cagsDF$dataset_name)
  # cagsDF$dataset_name <- gsub('YuJ_2015', 'Karlsson', cagsDF$dataset_name)
  # cagsDF$dataset_name <- gsub('ZellerG_2014', 'Karlsson', cagsDF$dataset_name)
  # browser()
  
  
}

main()


