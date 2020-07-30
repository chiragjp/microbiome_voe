suppressMessages(library(tidyverse))
suppressMessages(library(readxl))
library(readxl)

# The combined Metadata Excel sheet should:
#   1. Have a column called "dataset_name" with a reference to the researcher or name of the study associated
#   2. Have a column named "disease" (where 0 indicates control or lack of T2D, and 1 indicates prediabetes and diabetes

main <- function() {
   
  # * PREPARE METADATA *
  # This is if you have your metadata saved as an excel file.
  # Replace the arguments with the names of your files (save the files in the same directory)
  args <- c('2020_T2D_Data/final_metadata.rds','2020_T2D_Data/raw_abundance_data.rds', '1_prepared_metadata.rds', '2020_T2D_Data/Karlsson_Mapping.rds')
  #args <- commandArgs(trailingOnly = TRUE)
  # * PREPARE METADATA *
  metadata <- readRDS(args[[1]]) %>%
    rename(runID = sampleID) %>%
    rename(sampleID = PatientID) %>%
    rename(study_condition = disease) %>%
    mutate(subjectID = sampleID)
  metadata <- metadata %>% nest(-dataset)  
  saveRDS(metadata, file = args[[3]])
  
  abundance_data <- readRDS(args[[2]])
  geneNames <- abundance_data[[1]]
  newNames <-  seq(1,length(geneNames),by=1)
  abundance_data[[1]] = newNames
  
  
  abundance_data <- as.data.frame(t(abundance_data))
  names(abundance_data) <- as.character(unlist(abundance_data[1,]))
  abundance_data <- abundance_data[-1,]
  my_row_names <- rownames(abundance_data)
  abundance_data <- as_tibble(abundance_data) %>%
    mutate(sampleID = my_row_names)
  #abundance_data <- read_tsv(args[[2]])
  # here is where it transposes
  #abundance_data <- abundance_data %>%
  #  rownames_to_column() %>%
  #  pivot_longer(-GeneName, 'variable', 'value') %>%
  #  pivot_wider(variable, GeneName)

  browser()
  mapping <- cbind(geneNames, newNames)
  mapping <- as_tibble(mapping)
  mapping <- mapping %>%
    rename(feature_names = geneNames) %>%
    rename(num = newNames)
  
  saveRDS(mapping, args[[4]])
  
  abundance_data <- abundance_data %>% 
    mutate(dataset_name = rep("Karlsson", nrow(abundance_data)))  %>%
    select(dataset_name, sampleID, everything()) 
  
  saveRDS(abundance_data, '2020_T2D_Data/abundance_data.rds')
  print('finished')
  
  
}

main()


