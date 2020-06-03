###install packages

install.packages(c("tidyverse","vegan","ggplot2","magrittr","lme4","cowplot","lmerTest","magrittr","meta",",metafor","rje","rlang","broom.mixed"),dependencies=TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("curatedMetagenomicData")

