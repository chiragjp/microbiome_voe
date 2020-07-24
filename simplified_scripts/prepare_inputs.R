suppressMessages(library(tidyverse))
suppressMessages(library(readxl))

library(BiocFileCache)
library(curatedMetagenomicData)
library(readxl)
library(skimr)

# The combined Metadata Excel sheet should:
#   1. Have a column called "dataset_name" with a reference to the researcher or name of the study associated
#   2. Have a column named "disease" (where 0 indicates control or lack of T2D, and 1 indicates prediabetes and diabetes

main <- function() {
   
  # This is if you have your metadata saved as an excel file.
  # Replace the arguments with the names of your files (save the files in the same directory)
  args <- c('2020_T2D_Data/final_metadata.rds','Example_Data/t_metaphlan_abundance_cmd_example_CRC.rds', '1_prepared_metadata.rds')
  # * PREPARE METADATA *
  metadata <- readRDS(args[[1]])
  metadata <- metadata %>% nest(-dataset) 
  saveRDS(metadata, file = args[[3]])

  
}

main()


