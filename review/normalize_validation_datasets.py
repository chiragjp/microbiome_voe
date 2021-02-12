import pandas as pd
import os
import sys

locationfile = '../samples_to_keep_ibd'

genedata = pd.read_csv('../gene_information_for_norm',sep='\t',index_col=0,header=None)

with open(locationfile) as f:
	for name in f:
		print(name)
		name=name.rstrip()
		outputloc = '../normalized_crc/' + name + '_normalized.tsv'
		files = os.listdir('.')
		files = [x for x in files if name in x]
		data_files = []
		for fname in files:		
			data_files.append(pd.read_csv(fname,sep='\t',index_col=0,header=None).iloc[:,1])
		data_files = pd.concat(data_files,axis=1).sum(axis=1)
		data = pd.concat([genedata,data_files],axis=1)
		data = data.iloc[:-1,:]
		print('	Normalizing data')
		data[3] = data[0]/data[1]
		output = data[3]/sum(data[3])
		output.name = name
		print('	Writing to file')
		output.to_csv(outputloc,sep='\t',index=False,columns=[name])
		print('	Done')