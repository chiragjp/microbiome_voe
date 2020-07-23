suppressMessages(library(tidyverse))
suppressMessages(library(readxl))

library(BiocFileCache)
library(curatedMetagenomicData)
library(readxl)
library(skimr)

# The combined Metadata Excel sheet should:
#   1. Have a column called "dataset_name" with a reference to the researcher or name of the study associated
#   2. Have a column named "disease" (where 0 indicates control or lack of T2D, 1 indicates prediabetes, and 2 indicates diabetes)


main <- function() {
  
  
  suppressMessages(library(httr))
  # Sample Data
  #sample_data <- read_tsv('https://storage.googleapis.com/gbsc-gcp-project-ipop_public/HMP/clinical_tests/clinical_tests.txt')
  #sample_data$VisitID <- as.character(sample_data$VisitID)
  #sample_data <- sample_data %>% separate(VisitID, into = c('PatientID', "VisitNum"), sep = '-')
  
  # Patient Data  
  #GET('https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-019-1236-x/MediaObjects/41586_2019_1236_MOESM3_ESM.xlsx', write_disk(tf <- tempfile(fileext = ".xlsx")))
  #patient_data <- read_excel(tf, 1L)
  #patient_data <- patient_data %>% rename('PatientID' = 'SubjectID')
  
  #Combining Patient with Sample Data 
  #combined_HMP_data <- full_join(sample_data, patient_data, by = 'PatientID')
  #saveRDS(combined_HMP_data, 'combined_HMP.rds') # saving for debugging purposes 
  
  # Convert Columns to correct values 
  combined_HMP_data <- readRDS('combined_HMP.rds')
  combined_HMP_data <- combined_HMP_data %>% rename(HDL_old = HDL) %>% 
                                              mutate(HDL = HDL_old * 0.0555) %>%
                                              mutate(dataset = rep('HMP', nrow(combined_HMP_data))) %>%
                                              mutate(country = rep('USA', nrow(combined_HMP_data))) %>%
                                              rename(sampleID = VisitNum) %>%
                                              rename(disease = Class)%>%
                                              rename(age = Adj.age)
  
  combined_HMP_data$disease[combined_HMP_data$disease == "Diabetic"] <- 1
  combined_HMP_data$disease[combined_HMP_data$disease == "Prediabetic"] <- 1
  combined_HMP_data$disease[combined_HMP_data$disease == "Crossover"] <- 1
  combined_HMP_data$disease[combined_HMP_data$disease == "Control"] <- 0
  
  final_metadata <- combined_HMP_data %>% select(dataset, PatientID, sampleID, disease, age, Gender, country, Ethnicity, BMI, CHOL, CR, HDL, HSCRP, LDL, TGL )

  
  browser()
  qinData <- curatedMetagenomicData("QinJ_2012.metaphlan_bugs_list.stool", dryrun=FALSE)
  qinDF <- pData(qinData[[1]])
  browser()
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


