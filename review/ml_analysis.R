#compute genetic architectures with randomized disease vector

###fudge factor


library(tidyverse)
library(broom)
library(purrr)
library(caret)
library(mixOmics)

select <- dplyr::select
map <- purrr::map
tidy <- broom::tidy

deploy_ml_analysis <- function(dataset, univariate_output_names, pheno){
	dataset_sub = dataset %>% select(subjectID,study_condition,dataset_name,study_condition,all_of(univariate_output_names))
	if(length(dataset_sub$subjectID)!=length(unique(dataset_sub$subjectID))){
		dataset_sub = dataset_sub %>% group_by(subjectID,study_condition,dataset_name) %>% summarise(across(everything(), mean)) %>% ungroup()
	}
	dataset_sub = dataset_sub %>% select(-dataset_name,-subjectID)
	#elastic net
	set.seed(42)
	cv_10 = trainControl(method = "cv", number = 10)
	print('fitting elastic net')
	elnet = train(as.factor(study_condition) ~ ., data = dataset_sub, method = "glmnet", trControl = cv_10)
	elnet_coef = data.frame(as.matrix(coef(elnet$finalModel,elnet$finalModel$lambdaOpt))) %>% rename(elastic_net=X1) %>% rownames_to_column('feature') %>% mutate(feature = gsub('`','',feature))
	#sparsePLS
	print('fitting spls')
	spls_out = spls(Y=as.matrix(dataset_sub %>% select(study_condition)),X= as.matrix(dataset_sub %>% select(-study_condition)),ncomp=1)
	spls_coef = as.data.frame(as.matrix(spls_out$loadings$X)) %>% rename(spls=comp1) %>% rownames_to_column('feature')
	#random forest
	print('fitting random forest')
	rf = train(as.factor(study_condition) ~ ., data = dataset_sub, method = "rf", trControl = cv_10)
	rf_imp = data.frame(as.matrix(varImp(rf)$importance)) %>% rename(random_forest=Overall) %>% rownames_to_column('feature')%>% mutate(feature = gsub('`','',feature))
	ml_output=full_join(elnet_coef,rf_imp)
	ml_output=full_join(ml_output,spls_coef) %>% filter(feature!='(Intercept)')
	return(ml_output)
}

args = commandArgs(trailingOnly=TRUE)

print('Loading raw data...')
inputfile_list = unlist(unname(read.csv(args[[1]],header=FALSE)))
datasetfile = args[[2]]
phenotype = args[[3]]

randomized = strsplit(args[[1]],'_') %>% map_chr(6)
batch = strsplit(args[[1]],'_') %>% map_chr(5)

dataset = readRDS(datasetfile) %>% select(dataset_name,sampleID,all_of(readRDS(paste('top_10k_genes_',args[[1]],sep=''))))
dataset$dataset_name =  gsub('gene_families_','',unlist(unname(dataset$dataset_name)))
dataset$sampleID =  unlist(unname(dataset$sampleID))
metadata = readRDS('metaphlan_metadata.rds') %>% select(dataset_name,sampleID,subjectID,study_condition)
dataset = left_join(metadata,dataset) %>% select(-sampleID)
datasets = dataset %>% filter(study_condition ==phenotype) %>% select(dataset_name) %>% distinct %>% unlist %>% unname
dataset = dataset %>% filter(dataset_name %in% datasets, study_condition=='control' | study_condition==phenotype) %>% mutate(study_condition = if_else(study_condition == 'control',0,1))

print('Loading regression results...')
initial_regression_output = list()
for(f in inputfile_list){
	initial_regression_output[[f]] = readRDS(f)
}

print('Merging data...')
initial_regression_output = bind_rows(initial_regression_output)

print('Finding most significant outcomes and running analysis...')
#take 10K  significant results and run through elastic net, random forest, and sparse pls
univariate_output_names = initial_regression_output %>% arrange(p.value,desc(abs(estimate))) %>% head(10000) %>% select(feature) %>% unlist %>% unname
ml_output = deploy_ml_analysis(dataset, univariate_output_names,phenotype)
output = left_join(initial_regression_output,ml_output) %>% mutate(batch_corrected = batch,randomized=randomized)

print('Writing to file...')
saveRDS(output,paste('modeling_comparison_alternative_approaches','_',args[[1]],sep=''))












