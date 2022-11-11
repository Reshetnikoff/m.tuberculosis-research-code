import pandas as pd
import numpy as np
import sys
from os.path import exists
import os
from sklearn.model_selection import StratifiedKFold
from joblib import Parallel, delayed

phenofolder='../db/pheno/'

def stratified_split(my_drug, outfolder, num):
    skf = StratifiedKFold(n_splits=5, shuffle=True)
    os.makedirs(outfolder + my_drug + '/' + str(num), exist_ok=True)
    logfile_pattern = outfolder + my_drug + '/splits/' + my_drug + '.' + str(num) + '.'
    print(logfile_pattern[:-1])
    c = 1
    for train_index, test_index in skf.split(X, y):
        X_train, X_test = X[train_index], X[test_index]
        with open(logfile_pattern + str(c) + '_fold_split.txt', 'w') as logfile:
            logfile.write('train\n')
            for sample in X_train:
                logfile.write(sample + '\n')
            logfile.write('test\n')
            for sample in X_test:
                logfile.write(sample + '\n')
        with open(logfile_pattern + 'fold_' + str(c) + '.train_indices', 'w') as logfile_indices:
            for ind in train_index:
                logfile_indices.write(str(ind) + '\n')
        c += 1

if __name__ == "__main__":
    my_drug = sys.argv[1]
    outfolder = sys.argv[2]
    os.makedirs(outfolder + my_drug + '/splits', exist_ok=True)
    X, y = [], []
    with open(phenofolder + my_drug + '.pheno') as f:
        for line in f:
            X.append(line.strip().split('\t')[0])
            y.append(int(line.strip().split('\t')[1]))
    X, y = np.array(X), np.array(y)
    tasks = Parallel(n_jobs=100)(delayed(stratified_split)(my_drug, outfolder, i) for i in range(1, 101))