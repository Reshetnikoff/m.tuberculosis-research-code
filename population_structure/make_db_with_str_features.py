# Script for add TreeBreaker features to the matrix

# Configutation part #

# Directory with matrix files
data_dir = ''
# Directory with split info
split_dir = ''
# Output directory
output_dir = ''
# Output of TreeBreaker features script
feature_dir = ''
# Diretory with pheno
pheno_dir = '../db/pheno'


######################

import os
import sys

import pandas as pd
import numpy as np

drug = sys.argv[1]

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

def func(drug): 
    for ind in range(5):
        f = open(f"{feature_dir}/{drug}.features", 'r')
        feature_data = f.readlines()
        feature_names =  [x.split('\t')[0] for x in feature_data]
        str_feature = [x.split('\t')[1].split(',') for x in feature_data]
        f.close()

        pheno = pd.read_csv(f"{pheno_dir}/{drug}.pheno", sep='\t', index_col=0, header=None)
        split = pd.read_csv(f"{split_dir}/{drug}/fold{ind+1}/train_test_split.tsv", 
                           sep='\t', index_col=0, header=None)

        data_train = pd.read_csv(f"{data_dir}/{drug}.gt.domains.3.train.{ind}", 
                                 sep=' ', index_col=None, header=None)
        data_test = pd.read_csv(f"{data_dir}/{drug}.gt.domains.3.test.{ind}", 
                                 sep=' ', index_col=None, header=None)
        data_pos = pd.read_csv(f"{data_dir}/{drug}.gt.pos.domains.3.{ind}", 
                                 sep='\t', index_col=None, header=None)

        data_train.index=split[split[1] == 'train'].index
        data_test.index=split[split[1] == 'test'].index

        i = 0
        for samples, name in zip(str_feature, feature_names):
            data_train[data_train.columns[-1]+1] = [1 if x in samples else 0 for x in data_train.index]
            data_test[data_test.columns[-1]+1] = [1 if x in samples else 0 for x in data_test.index]
            data_pos.loc[len(data_pos)] = ['structure', i, name, np.NaN, np.NaN]
            i += 1

        data_train.to_csv(f"{output_dir}/{drug}.gt.domains.3.train.{ind}", 
                          sep=' ', index=False, header=False)
        data_test.to_csv(f"{output_dir}/{drug}.gt.domains.3.test.{ind}", 
                          sep=' ', index=False, header=False)
        data_pos.to_csv(f"{output_dir}/{drug}.gt.pos.domains.3.{ind}", 
                                 sep='\t', index=False, header=False)

func(drug)
