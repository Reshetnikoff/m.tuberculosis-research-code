# Script to convert output of abess_lr.py to appropriate input make_db_next_abess_iter.py format

# Configutation part #

# Directory with matrix before deduplication
data_dir_pre = ""
# Directory with matrix after deduplication
data_dir = ""
# Output of abess_lr.py
input_dir = "./output_abess"
# Output of this script
output_dir = "./output_abess"

######################
import os 
import sys

import pandas as pd
import numpy as np

drug = sys.argv[1]

thr = "3"

drugs = ['Kanamycin', 'Amikacin', 'Streptomycin', 'Ofloxacin',
         'Moxifloxacin', 'Isoniazid', 'Rifampicin', 'Ethambutol',
         'Pyrazinamide', 'Capreomycin', 'Ethionamide', 'Prothionamide',
         'Ciprofloxacin']

def get_name_by_ind(pos_info, j, drug):
    temp = pos_info.iloc[j]
    gene, pos, ind1, ind2, act = temp['gene'], temp['pos'], temp['ind1'], temp['ind2'], temp['act']
    
    if pos == 'changed':
        return f"{gene}_domain"
    elif pos == 'broken':
        return f"{gene}_broken"
    elif act == 'snp':
        return f"{gene}_{pos}_{ind2}_snp"
    elif ind2 == 'del':
        return f"{gene}_{pos}_{ind1}_{ind2}"
    elif act == 'del':
        return f"{gene}_{pos}_{ind2}_{act}"
    elif act == 'ins':
        return f"{gene}_{pos}_{ind2}_{act}"
    elif gene == 'structure':
        return f"{gene}_{pos}_{structure}"

def find_index(first_pos_info, second_pos_info, ind):
    data = second_pos_info.loc[np.sum(first_pos_info.loc[ind].values == second_pos_info.values, axis=1) == 5]
    if len(data) != 0:
        return int(data.index[0])
    else: return -1

def converter(drug):
    general_output = pd.DataFrame()
    for ind in range(5):
        pos_info = pd.read_csv(f"{data_dir_pre}/{drug}.gt.pos.domains.{thr}.{ind}", sep='\t', index_col=None, header=None, usecols=[0, 1, 2, 3, 4])
        pos_info.columns = ["gene", 'pos', 'ind1', 'ind2', 'act']

        f = open(f"{data_dir}/{drug}.gt.features.domains.{thr}.train.{ind}", 'r')
        features_info = [[int(y) for y in x[:-1].split(', ')] for x in f.readlines()]
        f.close()
       
        data = pd.read_csv(f"{input_dir}/{drug}.domains.3.result.{ind}", sep='\t', index_col=0, header=None, usecols=[0, 1, 2])
        if len(data) == 0: continue
        data = data[data.index == 'gic']
        data = data.astype({1: str, 2: str})        

        print(data)
        indexes =  [int(x[1:])-1 for x in data[1].values[0].split(', ')]
        features = [float(x) for x in data[2].values[0].split(', ')]

        print(indexes)
        print(features)
        features_ind = [features_info[int(i)] for i in indexes]
        
        output = pd.DataFrame()
        output['ind'] = indexes
        output['ext_ind'] = [", ".join([str(i) for i in x]) for x in features_ind]
        output['ext_desc'] = [", ".join([get_name_by_ind(pos_info, int(i), drug) for i in x]) for x in features_ind]
        output['coef'] = features
        general_output = pd.concat([general_output, output], axis=0)
        output.to_csv(f"{output_dir}/{drug}.domains.3.output.{ind}", sep='\t', index=False, header=True)
    general_output['sum'] = 1
    general_output.drop(labels=['ind'], axis=1, inplace=True)
   
    general_output = general_output.groupby(['ext_desc']).agg({'sum': 'sum', 'coef': 'mean'})
    general_output.sort_values(by='sum', ascending=False, inplace=True)
    
    general_output.to_csv(f"{output_dir}/{drug}.domains.3.output", sep='\t', header=True)
converter(drug)
