# Script for generating hamming distance files

# Configutation part #

# Directory with matrix after deduplication
input_dir = "../db/cv_bess_files.2"
# Output directory with HHS format files
output_dir = "../db/db/hhs_files"

######################

import sys

import gzip

import pandas as pd
import numpy as np

from scipy.spatial.distance import hamming

drug = sys.argv[1]
thr = '3'

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

for ind in range(5):
    data = pd.read_csv(f"{input_dir}/{drug}.gt.domains.{thr}.train.{ind}", sep=' ', index_col=None, header=None)
    n = data.shape[0]    
    
    result = pd.DataFrame(columns=[f"sam{i}" for i in range(n)])
    
    for i, name_row in enumerate([f"sam{i}" for i in range(n)]):
        row = data.loc[i]
        result.loc[name_row] = (row.to_numpy() != data.to_numpy()).sum(axis=1)

    result.to_csv(f"{output_dir}/{drug}.gt.domains.{thr}.train.{ind}.dist", sep=',', index=True, header=True) 
