# microbiome_voe

This repository contains the associated code for our manuscript, "XXX." Contents of directories:

### manuscript_scripts:

> Bash and figure plotting scripts used in automating the pipeline for use on a slurm-based compute and/or generating figures.

### pipeline_scripts: 

> Specific scripts modified for our particular use case, parallelized computation of associations across millions of features for a range of microbiome data types (e.g. species, pathways, genes).

### example_scripts:

> A refined example case of our pipeline that will run locally on most machines. If you like, you can adapt this code as a model for your own vibration of effects analyses (this will be most easily done with data from curatedMetagenomicData given the structure of the pipeline.)

#### to run the example code:

1) Clone the repo:

`git clone https://github.com/chiragjp/microbiome_voe.git`

2) Install miniconda (https://docs.conda.io/en/latest/miniconda.html) and cd into the example_scripts directory

3) Using the YAML file in the example_scripts directory, install the necessary R packages into a conda environment with:

```
conda env create --file voe.yml --name voetest 
conda activate voetest
```

4) You can now run the miniature pipeline in the example_scripts directory. 
```
./deploy_example_pipeline.sh
```

On most modern persnal machines it should complete in about 2 minutes tops, generating plots and summary data in a folder called `example_pipeline_output`. If you want more information on each individual step as well as script arguments, you can read the `deploy_example_pipeline.sh` file, which has all of the above as well as a description of the output files you would expect 

If you don't like using Conda, you can alternatively just install the dependencies (listed below) manually. Or using the script example_scripts/install_packages.R, though depending on your system architecture and available dependencies this can occasionally run into build issues that have to be manually debugged.

The R (V3.6.1) dependencies and versions are all listed in example_scripts/voe.yml file.