#combine association outputs

#script for merging meta analysis outputs

library(tidyverse)

files=list.files()
files=files[grep('meta_analysis.rds',files)]

dataset=c('metaphlan','pathways','genefamilies')

diseases=c('CRC','IBD','ACVD','T1D','T2D','cirrhosis','adenoma','otitis')

for(datatype in dataset){
  print(datatype)
  mapper=read.csv(paste(datatype,'_split_files_temp',sep=''),header=FALSE,stringsAsFactors = FALSE)
  t_ids=mapper$V1
  segment_mapping=data.frame(map_chr(t_ids, ~as.integer(str_split(str_split(., "_")[[1]][[4]],'.rds')[[1]][[1]])),stringsAsFactors = FALSE)
  for(disease in diseases){
    files_sub=files[grep(datatype,files)]
    files_sub=files_sub[grep(disease,files_sub)]
    dataFiles=list()
    for(f in files_sub){
      segment=as.numeric(strsplit(strsplit(f,'_meta_')[[1]][[1]],'_')[[1]][[4]])
      x=readRDS(f)
      names=colnames(x)
      count=0
      for(n in names){
        count=count+1
        base=as.numeric(strsplit(n,'_')[[1]][[2]])
        name=paste('feature',(as.numeric(segment_mapping[segment,])-1)*4000+base,sep='_')
        names[count]=name
      }
      colnames(x)=names
      dataFiles[[f]]=x
    }
    dataFile <- bind_cols(dataFiles)
    print(disease)
    print(dim(dataFile))
    saveRDS(dataFile,paste(datatype,'_',disease,'_meta_assocs_merged.rds',sep=''))
  }
}


