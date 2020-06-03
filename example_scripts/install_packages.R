###install packages

install.packages("tidyverse",repos='http://cran.us.r-project.org',dependencies=TRUE)
install.packages("vegan",repos='http://cran.us.r-project.org',dependencies=TRUE)
install.packages("ggplot2",repos='http://cran.us.r-project.org',dependencies=TRUE)
install.packages("magrittr",repos='http://cran.us.r-project.org',dependencies=TRUE)
install.packages("lme4",repos='http://cran.us.r-project.org',dependencies=TRUE)
install.packages("cowplot",repos='http://cran.us.r-project.org',dependencies=TRUE)
install.packages("lmerTest",repos='http://cran.us.r-project.org',dependencies=TRUE)
install.packages("meta",repos='http://cran.us.r-project.org',dependencies=TRUE)
install.packages("metafor",repos='http://cran.us.r-project.org',dependencies=TRUE)
install.packages("rje",repos='http://cran.us.r-project.org',dependencies=TRUE)
install.packages("rlang",repos='http://cran.us.r-project.org',dependencies=TRUE)
install.packages("broom.mixed",repos='http://cran.us.r-project.org',dependencies=TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("curatedMetagenomicData")

