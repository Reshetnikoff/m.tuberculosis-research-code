# Evualation quality of prediction on logistic regression using Indels & SNPs & broken gene features datasets

# Configuration part #

# Directory with matrix without dataset split and without PFAM domain features
bess_dir = "../db/bess_files"
# Directory with information about data splits
split_dir = f"../db/splits"
# Output directory
output_dir = f"./output/snp_indel"

######################

import os 
import sys

from multiprocessing import Pool

import pandas as pd
import numpy as np

from scipy.sparse import csr_matrix
from sklearn.metrics import roc_auc_score
from sklearn.linear_model import LogisticRegression

drug = sys.argv[1]

output_dir = f"{output_dir}/{drug}"
split_dir = f"{output_dir}/{drug}"

N_ITERS = 100
N_SPLIT = 5
N_THREADS = 12

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

def get_auc_dict(y_test, y_proba):
    auc = 0
    try:
        auc = roc_auc_score(y_test, y_proba)
    except ValueError:
        pass
    return auc

def get_TN(y_test, y_pred):
    TN = 0
    for i in range(len(y_test)):
        if y_test[i] == 0 and y_pred[i] == 0:
            TN += 1
    return TN

def get_TP(y_test, y_pred):
    TP = 0
    for i in range(len(y_test)):
        if y_test[i] == 1 and y_pred[i] == 1:
            TP += 1
    return TP

def get_FP(y_test, y_pred):
    FP = 0
    for i in range(len(y_test)):
        if y_test[i] == 0 and y_pred[i] == 1:
            FP += 1
    return FP

def get_FN(y_test, y_pred):
    FN = 0
    for i in range(len(y_test)):
        if y_test[i] == 1 and y_pred[i] == 0:
            FN += 1
    return FN

if not os.path.exists(f"{output_dir}"):
    os.makedirs(f"{output_dir}")

data = pd.read_csv(f"{bess_dir}/{drug}.gt.3", sep=' ', 
        index_col=None, header=None)
pheno = pd.read_csv(f"{bess_dir}/{drug}.phen", sep=' ', index_col=None, header=None)
pos_info = pd.read_csv(f"{bess_dir}/{drug}.gt.pos.3", sep='\t', 
        index_col=None, header=None, usecols=[0, 1, 2, 3, 4])
pos_info.columns = ['gene', 'pos', 'ind1', 'ind2', 'act']

# Filtration

y = pheno.to_numpy().ravel()


result = pd.DataFrame(columns=['n_iter', 'n_split', 'auc', 'TN', 'TP', 'FP', 'FN'])
k = 0

for i in range(1, N_ITERS + 1):
    for j in range(1, N_SPLIT + 1):
        f = open(f"{split_dir}/{drug}.{i}.fold_{j}.train_indices")
        train_index = [int(x[:-1]) if x[-1] == '\n' else int(x) for x in f.readlines()]
        f.close()

        test_index = [int(x) for x in data.index if int(x) not in train_index]

        X_train = data.iloc[train_index, :]
        X_test = data.iloc[test_index, :]
        y_train = y[train_index]
        y_test = y[test_index]
  
        X_train = X_train.transpose()
        X_train.drop_duplicates(inplace=True)
        X_train = X_train.transpose()

        X_test = X_test[X_train.columns]

        columns = X_train.columns

        X_train = csr_matrix(X_train)
        X_test = csr_matrix(X_test)

        model = LogisticRegression(solver='saga', penalty='l1', max_iter=10000)
        model.fit(X_train, y_train)

        y_proba = model.predict_proba(X_test)[:, 1]
        y_pred = model.predict(X_test) 
        

        coef = pd.DataFrame(data=model.coef_[0],
                                columns=["coef"])

        coef['ext_features'] = [get_name_by_ind(pos_info, int(x), drug) for x in columns]
        coef.sort_values(by=['coef'], ascending=False, inplace=True)
        coef.to_csv(f"{output_dir}/coef_{i}_{j}.csv",
                          sep='\t', header=True, index=True)

        auc = get_auc_dict(y_test, y_proba)
        TP = get_TP(y_test, y_pred)
        FP = get_FP(y_test, y_pred)
        TN = get_TN(y_test, y_pred)
        FN = get_FN(y_test, y_pred)

        result.loc[k] = [i, j, auc, TN, TP, FP, FN]
        k += 1
    if i % 5 == 0:
        result.to_csv(f"{output_dir}/result.csv", sep='\t', 
                index=False, header=True)

result.to_csv(f"{output_dir}/result.csv", sep='\t', 
                index=False, header=True)
