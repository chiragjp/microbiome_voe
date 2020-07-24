library(tidyverse)

# CHANGE FILE NAMES HERE (to be specific for your data)!
metadata_file <- readRDS('2020_T2D_Data/final_metadata.rds')
feature_file <- readRDS('Example_Data/t_metaphlan_abundance_cmd_example_CRC.rds')


# Run the Pipeline
source("prepare_inputs.R")
source("compute_associations.R")
source("compute_metanalysis.R")
source("compute_vibrations.R")