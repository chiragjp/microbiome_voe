library(curatedMetagenomicData)

#use this script to download and merge pathways, gene family, and species abundances
#will require substantial RAM (~150G) and storage space (~100G)

saveRDS(combined_metadata,'cmd_all_metadata.rds')

a=curatedMetagenomicData('*genefamilies_relab*',dryrun=T)
a2=a[-grep('WenC_2017',a)]
a=curatedMetagenomicData(a2,dryrun=F)

count=0
for(i in seq(length(a))){
  print(a[i])
  name=names(a[i])
  i2=mergeData(a[i])
  i2=exprs(i2)
  #i2=i2[rowSums(i2!=0)>(.2*ncol(i2)),]
  saveRDS(i2,paste('gene_families_',name[[1]],'_cmd_full.rds',sep=''))
}


a=curatedMetagenomicData('*pathabundance*',dryrun=T)
a2=a[-grep('WenC_2017',a)]
a=curatedMetagenomicData(a2,dryrun=F)
a=mergeData(a)
a=exprs(a)
saveRDS(a,'pathway_abundance_cmd.rds')

a=curatedMetagenomicData('*metaphlan*',dryrun=T)
a2=a[-grep('WenC_2017',a)]
a=curatedMetagenomicData(a2,dryrun=F)
a=mergeData(a)
a=exprs(a)
saveRDS(a,'metaphlan_abundance_cmd.rds')



