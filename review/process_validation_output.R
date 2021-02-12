library(tidyverse)

gene_val_list = read.csv('genes_for_validation_full.csv')
gene_val_list_ibd = gene_val_list %>% filter(phenotype=='IBD')
gene_val_list_crc = gene_val_list %>% filter(phenotype=='CRC')

ibd = readRDS('ibd_validation.rds') %>% mutate(phenotype='IBD')
crc = readRDS('crc_validation_error_removed.rds') %>% mutate(phenotype='CRC')

ibd = ibd %>% filter(feature!='study_condition') %>%  mutate(id = strsplit(feature,'\\|','') %>% map_chr(2))
crc = crc %>% filter(feature!='study_condition') %>%  mutate(id = strsplit(feature,'\\|','') %>% map_chr(2))

ibd = left_join(gene_val_list_ibd,ibd,by=c('id')) %>% filter(!is.na(p.value)) 
crc = left_join(gene_val_list_crc,crc,by=c('id')) %>% filter(!is.na(p.value)) 

output = bind_rows(ibd,crc)%>% select(phenotype.x,estimate.x,estimate.y,std.error,statistic,p.value,feature) %>% rename(estimate_old = estimate.x,estimate_new = estimate.y, phenotype=phenotype.x)

intersect(ibd %>% filter(p.value<0.05) %>% select(feature_name),crc %>% filter(p.value<0.05) %>% select(feature_name))

temp = ibd %>% filter(p.value<0.05) %>% select(estimate.x,estimate.y) 
temp[temp>0]=1
temp[temp<0]=-1
temp %>% filter(estimate.x==estimate.y) %>% dim

temp =crc%>%filter(p.value<0.05) %>% select(estimate.x,estimate.y) 
temp[temp>0]=1
temp[temp<0]=-1
temp %>% filter(estimate.x==estimate.y) %>% dim