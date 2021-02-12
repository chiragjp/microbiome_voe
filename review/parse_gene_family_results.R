#parse gene family validation results for memory efficiency

library(tidyverse)
library(broom)
library(purrr)
library(caret)
library(mixOmics)

select <- dplyr::select
map <- purrr::map
tidy <- broom::tidy

args = commandArgs(trailingOnly=TRUE)

print('Loading raw data...')
inputfile_list = unlist(unname(read.csv(args[[1]],header=FALSE)))
datasetfile = args[[2]]
phenotype = args[[3]]

randomized = strsplit(args[[1]],'_') %>% map_chr(5)
batch = strsplit(args[[1]],'_') %>% map_chr(6)

print('Loading regression results...')
initial_regression_output = list()
for(f in inputfile_list){
	initial_regression_output[[f]] = readRDS(f)
}

print('Merging data...')
initial_regression_output = bind_rows(initial_regression_output)

print('Finding most significant outcomes and running analysis...')
#take 10K  significant results and run through elastic net, random forest, and sparse pls
univariate_output_names = initial_regression_output %>% arrange(p.value,estimate) %>% head(10000) %>% select(feature) %>% unlist %>% unname

saveRDS(univariate_output_names,paste('top_10k_genes','_',args[[1]],sep=''))
