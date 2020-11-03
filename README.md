# Leveraging vibration of effects to stress-test microbiome associations

This repository contains the associated code for our in-prep manuscript, "Robust, gene-level, metagenomic architectures across 7 diseases yield high-resolution disease indicators
." 

If you want to learn more about vibration of effects (VoE), check out this paper: 

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4555355/

This code is NOT meant to be distributed or deployed as a package, and instead is meant to provide assistance in understanding and reproducing the specific work in our manuscript. If you're interested in running VoE analyses yourself on other observational data, we are currently developing an R package for just this task.

## Contents of directories:

### manuscript_scripts:

> Bash and figure plotting scripts used in automating the pipeline for use on a slurm-based compute and/or generating figures.

### pipeline_scripts: 

> Specific scripts modified for our particular use case, parallelized computation of associations across millions of features for a range of microbiome data types (e.g. species, pathways, genes).

### example_scripts:

> A refined example case of our pipeline that will run locally on linux-based operating systems. It shows the vibration for a few associations between colorectal cancer and gut microbiota, with an emphasis on *F. nucleatum*, a particularly robust (i.e. minimal vibration) correlation. If you like, you can adapt this code as a model for your own vibration of effects analyses (this will be most easily done with data from curatedMetagenomicData given the structure of the pipeline). We are currently working to extend the code here to a more adaptable platform for other datasets.

#### To run the example code:

1) Clone the repo:

`git clone https://github.com/chiragjp/microbiome_voe.git`

2) Install miniconda (https://docs.conda.io/en/latest/miniconda.html) and cd into the example_scripts directory

3) Using the YAML file in the example_scripts directory, install the necessary R packages into and activate a conda environment with:

```
conda env create --file voe.yml --name voetest 
conda activate voetest
```

4) Due to some currently unresolved conda sourcing issues, you currently need to open an interactive session of R (type "R" at the terminal after loading the newly created environment) and install curatedMetagenomicData manually, following any prompts the appear:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager",repos='http://cran.us.r-project.org')

BiocManager::install("curatedMetagenomicData")
```

4) After the install completes, close out of the R environment, returning to your bash command line. You can now run the miniature pipeline in the example_scripts directory. 
```
chmod +x deploy_example_pipeline.sh
./deploy_example_pipeline.sh
```

On most modern personal machines it should complete in about 2 minutes tops (though the install can take longer), generating plots and summary data in a folder called `example_pipeline_output`. If you want more information on each individual step as well as script arguments, you can read the `deploy_example_pipeline.sh` file, which has all of the above as well as a description of the output files you would expect 

If you don't like using Conda, you can alternatively just install R (V4.0.1) the dependencies in the voe.yml manually or using the script example_scripts/install_packages.R, though depending on your system architecture and available dependencies this can occasionally run into build issues that have to be manually debugged.
