import pandas as pd
import os
from ete3 import NCBITaxa

ncbi = NCBITaxa()
#ncbi.update_taxonomy_database()

f='for_ete3_overall.csv'
data=pd.read_csv(f)
dicts=pd.DataFrame([[list(ncbi.get_name_translator([x]))[0],ncbi.get_name_translator([x])[x][0]] for x in data.species])
dicts.columns=['species','new_tax_id']
data=pd.merge(data,dicts,right_on='species',left_on='species')
tree = ncbi.get_topology(list(set(data.new_tax_id)))
tree.write(format=1, outfile="ncbi_tree_for_voe_d_%s.nw"%f.split('_')[2].split('.')[0])
data.to_csv('mapping_for_ncbi_tree_overall.csv')
