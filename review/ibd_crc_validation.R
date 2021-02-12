# architecture paper validation regressions

library(tidyverse)
library(broom)

options(warn=2)

crc = readRDS('crc_merged.rds')
crc= crc %>% select(-`*`)

genes = colnames(crc %>% select(-sampleID,-study_condition))

crc_complete = list()
for(g in genes){
	print(g)
	tryCatch({
		crc_complete[[g]] = tidy(lm(log(crc[,g]+0.000001) ~ study_condition,data = crc)) %>% filter(term!='(Intercept)') %>% mutate(feature=g)
	},
	error = function(e){})
}

crc_complete = bind_rows(crc_complete)

saveRDS(crc_complete,'crc_validation_error_removed.rds')

ibd = readRDS('ibd_averaged_merged.rds')
ibd = ibd[complete.cases(ibd), ]
 
genes = colnames(ibd %>% select(-subjectID,-study_condition)) 

ibd_complete = list()
for(g in genes){
	print(g)
	tryCatch({
		ibd_complete[[g]] = tidy(lm(log(ibd[,g][[1]]+0.000001) ~ study_condition,data = ibd)) %>% filter(term!='(Intercept)') %>% mutate(feature=g)
	},
	error = function(e){})}

ibd_complete = bind_rows(ibd_complete)

saveRDS(ibd_complete,'ibd_validation_error_removed.rds')



