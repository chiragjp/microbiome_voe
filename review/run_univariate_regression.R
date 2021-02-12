 ©#compute genetic architectures with randomized disease vector

###fudge factor


library(tidyverse)
library(broom)
library(purrr)
library(caret)

select <- dplyr::select
map <- purrr::map
tidy <- broom::tidy

#univariate associations with optional batch adjustment
univariate_linear_regression <-function(data_table_for_regression,feature_name,pheno,batch){
	#average repeated measures if necessary
	if(length(data_table_for_regression$subjectID)!=length(unique(data_table_for_regression$subjectID))){
		data_table_for_regression = data_table_for_regression %>% group_by(subjectID,study_condition,dataset_name) %>% summarize(feature = get(feature_name)) %>% ungroup()
	}
	else{
		colnames(data_table_for_regression)[3] = 'feature'
	}
	if(batch==TRUE){
		regression_output = broom::tidy(lm(data = data_table_for_regression, I(scale(feature)) ~ study_condition + dataset_name)) %>% filter(term=='study_condition') %>% mutate(feature = feature_name)
		regression_output$term[[1]] = pheno		
	}
	else{
		regression_output = broom::tidy(lm(data = data_table_for_regression, I(scale(feature)) ~ study_condition)) %>% filter(term=='study_condition') %>% mutate(feature = feature_name)
		regression_output$term[[1]] = pheno		
	}
	return(regression_output)
}

args = commandArgs(trailingOnly=TRUE)

#specify input file, disease of interest, whether or not to randomize disease vector, whether nor not to adjust for cohort, and the datatype
inputfile = args[[1]]
phenotype = args[[2]]
batch = as.logical(args[[3]])
randomize = as.logical(args[[4]])

#load data and subset to relevant phenotype
print('Loading and preprocessing data...')
data = readRDS(inputfile)  
data$dataset_name = gsub('gene_families_','',data$dataset_name)
data$dataset_name =  unlist(unname(data$dataset_name))
data$sampleID =  unlist(unname(data$sampleID))
metadata = readRDS('metaphlan_metadata.rds') %>% select(dataset_name,sampleID,subjectID,study_condition)
data = left_join(metadata,data) %>% select(-sampleID)
datasets = data %>% filter(study_condition ==phenotype) %>% select(dataset_name) %>% distinct %>% unlist %>% unname
data_subset = data %>% filter(dataset_name %in% datasets, study_condition=='control' | study_condition==phenotype) %>% mutate(study_condition = if_else(study_condition == 'control',0,1))

#if necessary, randomize the disease column
if(randomize==TRUE){
	print('Randomizing data...')
	sc = data_subset %>% select(study_condition) %>% unlist %>% unname
	sc = sample(sc)
	data_subset$study_condition = sc
}

#compute initial univariate regressions
print('Computing regressions...')
initial_regression_output = map(colnames(data_subset %>% select(-subjectID,-dataset_name,-study_condition)), function(x) univariate_linear_regression(data_subset %>% select(subjectID,dataset_name,x,study_condition),x,phenotype,batch)) %>% bind_rows 

print('Writing to file...')
saveRDS(initial_regression_output,paste('modeling_comparison_univariate','_',phenotype,'_',as.character(batch),'_',as.character(randomize),'_',inputfile,sep=''))
