#prep raw abundance data for negative control and method comparisons

library(tidyverse)
library(WGCNA)

chunk <- 10000

data <- as.data.frame(readRDS('metaphlan_abundance_cmd.rds'))
minval = min(unlist(unname(map(data, function(x) min(x[which(x>0)])))))
n <- nrow(data)
r  <- rep(1:ceiling(n/chunk),each=chunk)[1:n]
output <- split(data,r)
for(i in seq(length(output))){
	temp = log(output[[i]] + minval)
	t_output = transposeBigData(as.data.frame(temp),10000) %>% as.data.frame %>% rownames_to_column('sampleID') %>% mutate(dataset_name = unlist(unname(map(sampleID, function(x) strsplit(x,'\\.') %>% map_chr(1)))),sampleID = unlist(unname(map(sampleID, function(x) strsplit(x,':') %>% map_chr(2))))) %>% relocate(dataset_name)
	saveRDS(t_output,paste('t_',as.character(i),'_metaphlan.rds',sep=''))
}

data = as.data.frame(readRDS('pathway_abundance_cmd.rds'))
minval = min(unlist(unname(map(data, function(x) min(x[which(x>0)])))))
n <- nrow(data)
r  <- rep(1:ceiling(n/chunk),each=chunk)[1:n]
output <- split(data,r)
for(i in seq(length(output))){
	temp = log(output[[i]] + minval)
	t_output = transposeBigData(as.data.frame(temp),10000) %>% as.data.frame %>% rownames_to_column('sampleID') %>% mutate(dataset_name = unlist(unname(map(sampleID, function(x) strsplit(x,'\\.') %>% map_chr(1)))),sampleID = unlist(unname(map(sampleID, function(x) strsplit(x,':') %>% map_chr(2))))) %>% relocate(dataset_name)
	saveRDS(t_output,paste('t_',as.character(i),'_pathway.rds',sep=''))
}

data = as.data.frame(readRDS('genefamilies_stool_cmd.rds'))
minval = min(unlist(unname(map(data, function(x) min(x[which(x>0)])))))
n <- nrow(data)
r  <- rep(1:ceiling(n/chunk),each=chunk)[1:n]
output <- split(data,r)
for(i in seq(length(output))){
	temp = log(output[[i]] + minval)
	t_output = transposeBigData(as.data.frame(temp),10000) %>% as.data.frame %>% rownames_to_column('sampleID') %>% mutate(dataset_name = unlist(unname(map(sampleID, function(x) strsplit(x,'\\.') %>% map_chr(1)))),sampleID = unlist(unname(map(sampleID, function(x) strsplit(x,':') %>% map_chr(2))))) %>% relocate(dataset_name)
	saveRDS(t_output,paste('t_',as.character(i),'_genefamilies.rds',sep=''))
}


data = as.data.frame(readRDS('gene_families_overlapping_genes.rds'))
minval = min(unlist(unname(map(data, function(x) min(x[which(x>0)])))))
n <- nrow(data)
r  <- rep(1:ceiling(n/chunk),each=chunk)[1:n]
output <- split(data,r)
for(i in seq(length(output))){
	temp = log(output[[i]] + minval)
	t_output = transposeBigData(as.data.frame(temp),10000) %>% as.data.frame %>% rownames_to_column('sampleID') %>% mutate(dataset_name = unlist(unname(map(sampleID, function(x) strsplit(x,'\\.') %>% map_chr(1)))),sampleID = unlist(unname(map(sampleID, function(x) strsplit(x,':') %>% map_chr(2))))) %>% relocate(dataset_name)
	saveRDS(t_output,paste('t_',as.character(i),'_genefamilies_overlapping.rds',sep=''))
}


