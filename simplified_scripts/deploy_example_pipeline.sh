#!/bin/bash

###### data cleaning

### 1) download and prepare raw data 
#args: outputname of cleaned abundance data, outputname of mapping data between feature names and column number in abundance data 

### 2) prune metadata to remove columns with insufficient metadata or those that will cause the model to fail
#args: metadata downloaded from CMD, phenotype, threshold of metadata presence required in a given column for it to be considered, outputfilename
echo 'Cleaning metadata'
Rscript prepare_inputs.R '2020_T2D_Data/final_metadata.rds' '/n/scratch3/users/e/ea221/karlsson/completed_alignments/find_cags/small_abundance_matrix.tsv' '1_prepared_metadata.rds'

###### meta-analysis
### 3) compute initial associations
#please note that because we are working here with massively subsampled data with a different FDR cutoff, the results will NOT be identical to the manuscript in terms of which features are significantly associated with CRC
#args: pruned metadata, relative abundance data, outputname
 echo 'Computing initial regressions.'
Rscript compute_associations.R '1_prepared_metadata.rds' '2020_T2D_Data/abundance_data.rds' 'Karlsson_associations.rds' CRC
# 
# ### 4) compute meta-analysis across output from prior step
# #args: associations per cohort from prior step, outputname, column to meta-analyze over (study_condition, which refers to the column that holds the case-control encoding)
 echo 'Running meta-analysis.'
Rscript compute_metaanalysis.R CRC_example_associations.rds CRC_example_meta_analysis.rds study_condition 
# 
# ### 5) clean up meta-analytic results by first removing features present in only 1 cohort and then formatting the data in an easy to parse dataframe
# #args: output from meta-analysis, outputname of filtered dataframe
# Rscript filter_single_cohort_sig_feats.R CRC_example_meta_analysis.rds CRC_example_meta_analysis_filtered.rds
# #args: output from meta-analysis filtering step, path to mapping file generated in step 1, cohort type of 'multi' or 'single', outputname of filtered dataframe
Rscript clean_meta_analysis_output.R CRC_example_meta_analysis_filtered.rds metaphlan_mapping_data_example_CRC multi CRC_example_meta_analysis_filtered_cleaned.rds
# 
# ### 6) plot initial meta-analysis results and determine statistically significant features worth vibrating over
# rm -rf example_pipeline_output; mkdir example_pipeline_output
# #this script will generate an example volcano plot in the newly minted example_run_output_plots directory
# #argds: clean meta analysis output, name of output file that will contain features to vibrate over
echo 'Computing vibration of effects.'
Rscript plot_volcano_and_find_vibration_features.R CRC_example_meta_analysis_filtered_cleaned.rds example_CRC_features_to_vibrate.rds
# 
# ###### vibration of effects
# ###7) run vibration of effects analysis over features of interest
# #this step is in essence identical to step 4, except it fits many more models
# #this script has been simplified from that used in the full pipeline
# #the full pipeline had to deal with splitting the millions of features we had to deal with in to small dataframes while tracking the feature number to name mappings
# #seeing as this is a small example meant for low numbers of features (fewer than 10K), this code has been removed from this script 
#args: model dataframe (the merged metadata by cohort file) made in step 2, features to vibrate over, transposed abundance data, outputname, phenotype of interest
Rscript compute_vibrations.R '1_prepared_metadata.rds' 'example_CRC_features_to_vibrate.rds' 't_metaphlan_abundance_cmd_example_CRC.rds example_CRC_vibration_output.rds' 'CRC'  
# 
# ###8) plot vibration volcano plots, regress incidence of different adjusters on effect size
# #args: vibration output, cleaned meta-analysis output, mapping file, phenotype
# echo 'Parsing vibration output.'
Rscript parse_vibration_output.R example_CRC_vibration_output.rds CRC_example_meta_analysis_filtered_cleaned.rds metaphlan_mapping_data_example_CRC CRC
# 
# ###description of output files in the example_pipeline_output folder:
# 
# # 1) A volcano plot for your initial meta-analysis
# # 2) Volcano plots for your vibrations
# # 3) Regressioon output (disease_specific_confounders_estimate_metaphlan_mapping_data_example_CRC.csv) describing the effect of adjusters on model effect sizes
# # 4) Raw and summarized vibration data for all models fit