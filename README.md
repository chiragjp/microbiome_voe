# microbiome_voe

This repository contains the associated code for our manuscript, "XXX." Contents of directories:

###manuscript_scripts:

> Bash and figure plotting scripts used in automating the pipeline for use on a slurm-based compute and/or generating figures.

###pipeline_scripts: 

> Specific scripts modified for our particular use case, parallelized computation of associations across millions of features for a range of microbiome data types (e.g. species, pathways, genes).

###example_scripts:

> A refined example case of our pipeline that will run locally on most machines. 

####to run the example code:

1) Clone the repo:

`git clone https://github.com/chiragjp/microbiome_voe.git`

2) Install miniconda (https://docs.conda.io/en/latest/miniconda.html)

3) Using the YAML file in the example_scripts directory, install the necessary R packages into a conda environment with:

`XXX`

You can alternatively just install the dependencies (listed below) manually. Or using the script example_scripts/install_packages.R, though depending on your system architecture this can occasionally run into build issues that have to be manually debugged.

####R dependencies:

Dependency | Version
-----------|--------
broom|
broom.mixed|
cowplot|
curatedMetagenomicData|
ggplot2|
dplyr|
lme4|
lmerTest|
magrittr|
meta|
metafor|
purrr|
rje|
rlang 
stringr|
tidyverse|
vegan|
