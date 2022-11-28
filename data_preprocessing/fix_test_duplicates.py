#  Script for test matrix feature deduplication

# Configutation part #

# Input folder with matrices and feature's info
input_dir = "../db/cv_bess_files"
# Output folder (the same as in fix_duplicates.py)
output_dir = "../db/cv_bess_files.2"

######################

import os
import sys

import pandas as pd
import numpy as np

drug = sys.argv[1]
thr = "3"

for ind in range(5):
    data_file = f"{input_dir}/{drug}.gt.domains.{thr}.test.{ind}"
    data_feature_file = f"{output_dir}/{drug}.gt.features.domains.{thr}.train.{ind}"
    pheno_file = f"{input_dir}/{drug}.phen.domains.{thr}.test.{ind}"
    pheno = pd.read_csv(pheno_file, sep=' ', index_col=None, header=None)
    if os.stat(data_file).st_size == 0: continue
    data = pd.read_csv(data_file, sep=' ', index_col=None, header=None)
    
    f = open(data_feature_file, 'r')
    data_features = f.readlines()
    f.close()
    data_features = [x[:-1].split(', ') for x in data_features]     
    data_features = [[int(x) for x in y] for y in data_features]

    drop_columns = []
    for feature in data_features: 
        if len(feature) > 1:
            temp = feature
            temp.remove(min(temp))
            drop_columns.extend(temp)

    data.drop(drop_columns, axis=1, inplace=True)
    
    pheno_file_output =  f"{output_dir}/{drug}.phen.domains.{thr}.test.{ind}"
    pheno.to_csv(pheno_file_output, sep=' ', index=False, header=False)

    data_file_output = f"{output_dir}/{drug}.gt.domains.{thr}.test.{ind}"
    data.to_csv(data_file_output, sep=' ', index=False, header=False)
        


