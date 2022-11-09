# Script for matrix feature deduplication

# Configutation part #

# Input folder with matrices and feature's info
input_dir = ""
# Output folder 
output_dir = ""

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

def view1D(a): # a is array
    a = np.ascontiguousarray(a)
    void_dt = np.dtype((np.void, a.dtype.itemsize * a.shape[1]))
    return a.view(void_dt).ravel()

def group_duplicate_cols(df):
    a = df.values
    sidx = view1D(a.T).argsort()
    b = a[:,sidx]

    m = np.concatenate(([False], (b[:,1:] == b[:,:-1]).all(0), [False] ))
    idx = np.flatnonzero(m[1:] != m[:-1])
    C = df.columns[sidx].tolist()
    return [C[i:j] for i,j in zip(idx[::2],idx[1::2]+1)]

for ind in range(5):
    data_file = f"{input_dir}/{drug}.gt.domains.{thr}.train.{ind}"
    data_pos_file = f"{input_dir}/{drug}.gt.pos.domains.{thr}.{ind}"
    pheno_file = f"{input_dir}/{drug}.phen.domains.{thr}.train.{ind}"

    data = pd.read_csv(data_file, sep=' ', index_col=None, header=None)
    pheno = pd.read_csv(pheno_file, sep=' ', index_col=None, header=None)
    data_pos = pd.read_csv(data_pos_file, sep='\t', index_col=None, header=None,
                      names=['gene', 'pos', 'ind1', 'ind2', 'act', 'add'])

    duplicated_columns = group_duplicate_cols(data)

    new_columns = []
    drop_columns = []
    for column in data.columns:
        flag = 0
        if column in drop_columns:
            continue
        for duplicated in duplicated_columns:
            if column in duplicated:
                new_columns.append(', '.join([str(i) for i in duplicated]))
                drop_columns.extend(duplicated)
                drop_columns.remove(column)
                flag = 1
                break
        if flag == 0:
            new_columns.append(str(column))
    new_data_pos = pd.DataFrame(new_columns)

    data.drop(drop_columns, axis=1, inplace=True)
    data.fillna(0, inplace=True)
    
    data_file_output = f"{output_dir}/{drug}.gt.domains.{thr}.train.{ind}"
    data_pos_file_output = f"{output_dir}/{drug}.gt.features.domains.{thr}.train.{ind}"
    pheno_file_output =  f"{output_dir}/{drug}.phen.domains.{thr}.train.{ind}"
    
    pheno.to_csv(pheno_file_output, sep=' ', index=False, header=False)
    data.to_csv(data_file_output, sep=' ', index=False, header=False)
    new_data_pos.to_csv(data_pos_file_output, sep='\t', index=False, header=False)
        


