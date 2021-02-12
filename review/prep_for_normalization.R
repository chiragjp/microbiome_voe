# architecture paper validation regressions

library(tidyverse)
library(broom)

average_feature = function(feature_name,data){
	data = data %>% select(subjectID,study_condition,feature_name) %>% group_by(study_condition,subjectID) %>% summarise(.groups='keep',taxa = mean(get(feature_name))) %>% distinct
	colnames(data)[3]=feature_name
	data = data %>% ungroup %>% arrange(subjectID,study_condition)
	return(data)
}

# load new abundance data

ibd_norm = t(read.csv('ibd_val_normalized',sep='\t') %>% select(-X) %>% column_to_rownames('X0')) %>% as.data.frame %>% rownames_to_column('sampleID')

# in the case of the ibd data, average as necessary to remove repeated measures
ibd_mdat = t(read.csv('ibd_metadata',sep='\t',row.names=1)) %>% as.data.frame %>% rownames_to_column('sampleID')  %>% select(subjectID,sampleID,study_condition) %>% mutate(sampleID = gsub('HMP_2019_ibdmdb_','',sampleID), study_condition = if_else(study_condition=='control',0,1)) 

ibd_merged = left_join(ibd_norm,ibd_mdat) %>% filter(!is.na(subjectID))

output = list()
genes=colnames(ibd_merged %>% select(-study_condition,-subjectID,-sampleID))
count=0
for(taxa in genes[1:5]){
	count=count+1
	print(taxa)
	out = average_feature(taxa,ibd_merged) 
	subs = out %>% select(subjectID,study_condition)
	output[[taxa]] = out %>% select(-subjectID,-study_condition)
}

averaged_data = bind_cols(subs,output)

saveRDS(averaged_data,'ibd_averaged_merged.rds')

crc_norm = t(read.csv('crc_val_normalized',sep='\t',row.names=1)) %>% as.data.frame %>% rownames_to_column('sampleID')


crc_mdat = read.csv('crc_metadata',sep='\t') %>% select(sampleID,study_condition) %>% mutate(sampleID = gsub('-','.',sampleID), study_condition = if_else(study_condition=='control',0,1))


crc_merged = left_join(crc_norm,crc_mdat) 

saveRDS(crc_merged,'crc_merged.rds')

