# Script for converting file to HHS format 

# Configutation part #

# Directory with matrix after deduplication
input_dir = ""
# Output directory
output_dir = "../db/hhs_files"

######################

import os
import sys

import numpy as np
import pandas as pd

drug = sys.argv[1]
thr = "3"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

for ind in range(5):
    data_test =  pd.read_csv(f"{input_dir}/{drug}.gt.domains.{thr}.test.{ind}", sep=' ', index_col=None, header=None)
    data_train =  pd.read_csv(f"{input_dir}/{drug}.gt.domains.{thr}.train.{ind}", sep=' ', index_col=None, header=None)

    pheno_test = pd.read_csv(f"{input_dir}/{drug}.phen.domains.{thr}.test.{ind}", sep=' ', index_col=None, header=None)
    pheno_train = pd.read_csv(f"{input_dir}/{drug}.phen.domains.{thr}.train.{ind}", sep=' ', index_col=None, header=None)

    data_test = data_test.transpose()
    data_train = data_train.transpose()

    data_train.columns = [f"sam{i}" for i in range(data_train.shape[1])]
    data_test.columns = [f"sam{i}" for i in range(data_train.shape[1], (data_train.shape[1] + data_test.shape[1]))]

    pheno_train.index = [f"sam{int(i)}" for i in range(pheno_train.shape[0])]
    pheno_test.index = [f"sam{int(i)}" for i in range(pheno_train.shape[0], pheno_train.shape[0] + pheno_test.shape[0])]

    pheno_train = pheno_train.astype({0: int})
    pheno_test = pheno_test.astype({0: int})
    
    data_train.index = [f"feature{i}" for i in range(data_train.shape[0])]
    data_test.index = [f"feature{i}" for i in range(data_test.shape[0])]
        
    f = open(f"{output_dir}/{drug}.gt.domains.{thr}.test.{ind}", 'w')
    f.write(",".join(data_test.columns) + "\n")
    for row in list(data_test.index):
        f.write(row + ":" + "".join(data_test.loc[row].to_numpy(dtype=str)) + "\n")
    f.close()

    f = open(f"{output_dir}/{drug}.gt.domains.{thr}.train.{ind}", 'w')
    f.write(",".join(data_train.columns) + "\n")
    for row in list(data_train.index):
        f.write(row + ":" + "".join(data_train.loc[row].to_numpy(dtype=str)) + "\n")
    f.close() 

    pheno_test.to_csv(f"{output_dir}/{drug}.phen.domains.{thr}.test.{ind}", sep=',', index=True, header=False)
    pheno_train.to_csv(f"{output_dir}/{drug}.phen.domains.{thr}.train.{ind}", sep=',', index=True, header=False)
