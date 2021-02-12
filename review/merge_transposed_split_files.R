library(tidyverse)

files=list.files()
files=files[grep('t_',files)]
files=files[grep('rds',files)]
files=files[-grep('modeling',files)]
files=files[grep('genefam',files)]
files=files[-grep('overlap',files)]

data = readRDS(files[[1]])

for(f in files[2:length(files)]){
	print(f)
	data=bind_cols(data,readRDS(f) %>% select(-dataset_name,-sampleID))
}

saveRDS(data,'t_genefamilies_full.rds')

library(tidyverse)

files=list.files()
files=files[grep('t_',files)]
files=files[grep('rds',files)]
files=files[-grep('modeling',files)]
files=files[grep('metaphlan',files)]

data = readRDS(files[[1]])

for(f in files[2:length(files)]){
	print(f)
	data=bind_cols(data,readRDS(f) %>% select(-dataset_name,-sampleID))
}

saveRDS(data,'t_metaphlan_full.rds')

library(tidyverse)

files=list.files()
files=files[grep('t_',files)]
files=files[grep('rds',files)]
files=files[-grep('modeling',files)]
files=files[grep('pathway',files)]

data = readRDS(files[[1]])

for(f in files[2:length(files)]){
	print(f)
	data=bind_cols(data,readRDS(f) %>% select(-dataset_name,-sampleID))
}

saveRDS(data,'t_pathways_full.rds')