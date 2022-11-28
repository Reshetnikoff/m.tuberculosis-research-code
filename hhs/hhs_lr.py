# Script for applying logistic regression on HHS selected features

# Configutation part #

# Directory with matrix after deduplication
hhs_files_dir = "../db/cv_bess_files.2"
# Directory with matrix before deduplication
hhs_files_dir_bef_dedup = "../db/cv_bess_files"
# Directory with output run_multi_hhs.sh script
input_dir = f"./output_hhs"
# Output directory
output_dir = f"./output_hhs_lr"
# Hyperparameter - minimum frequency of mutations that lead to the resistance phenotype 
f = 1
######################


import os
import sys

import pandas as pd
import numpy as np

from sklearn.linear_model import LogisticRegression

from sklearn.metrics import roc_auc_score

input_dir = f"{input_dir}/{f}"
output_dir = f"{input_dir}/{f}"

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

thr = "3"

drugs = ['Kanamycin', 'Amikacin', 'Streptomycin', 'Ofloxacin',
         'Moxifloxacin', 'Isoniazid', 'Rifampicin', 'Ethambutol',
         'Pyrazinamide', 'Capreomycin', 'Ethionamide', 'Prothionamide',
         'Ciprofloxacin']

def get_auc(y_test, y_proba):
    y_proba = y_proba[:, 1]
    auc = 0
    try:
        auc = roc_auc_score(y_test, y_proba)
    except ValueError:
        pass
    return auc

def ppv_score(y_test, y_pred):
    TP = 0
    FP = 0
    for i in range(len(y_test)):
        if y_test[i] == 1 and y_pred[i] == 1:
            TP += 1
        if y_test[i] == 0 and y_pred[i] == 1:
            FP += 1
    if TP + FP != 0:
        return TP/(TP + FP)
    else:
        return np.NaN

def npv_score(y_test, y_pred):
    TN = 0
    FN = 0
    for i in range(len(y_test)):
        if y_test[i] == 0 and y_pred[i] == 0:
            TN += 1
        if y_test[i] == 1 and y_pred[i] == 0:
            FN += 1
    if TN + FN != 0:
        return TN/(TN + FN)
    else:
        return np.NaN

def sensitivity_score(y_test, y_pred):
    TP = 0
    FN = 0
    for i in range(len(y_test)):
        if y_test[i] == 1 and y_pred[i] == 1:
            TP += 1
        if y_test[i] == 1 and y_pred[i] == 0:
            FN += 1
    if TP + FN != 0:
        return TP/(TP + FN)
    else:
        return np.NaN

def specificity_score(y_test, y_pred):
    TN = 0
    FP = 0
    for i in range(len(y_test)):
        if y_test[i] == 0 and y_pred[i] == 0:
            TN += 1
        if y_test[i] == 0 and y_pred[i] == 1:
            FP += 1
    if TN + FP != 0:
        return TN/(TN + FP)
    else:
        return np.NaN

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

k = 0
result = pd.DataFrame(columns=["drug", "auc", "sensitivity", "specivity", "npv", "ppv"])
for drug in drugs:
    for ind in range(5):
        pos_info = pd.read_csv(f"{hhs_files_dir_bef_dedup}/{drug}.gt.pos.domains.{thr}.{ind}", sep='\t', index_col=None, header=None, usecols=[0, 1, 2, 3, 4])
        pos_info.columns = ["gene", 'pos', 'ind1', 'ind2', 'act']

        f = open(f"{hhs_files_dir}/{drug}.gt.features.domains.{thr}.train.{ind}", 'r')
        features_info = [[int(y) for y in x[:-1].split(', ')] for x in f.readlines()]
        f.close()
        
        sel_output = pd.read_csv(f"{input_dir}/{drug}.gt.domains.{thr}.train.{ind}.result", sep=',',
                                 index_col=0, header=1, skiprows=0)
         
        indeces = [int(x[7:]) for x in sel_output.index]
        features_inds = [features_info[int(x[7:])] for x in sel_output.index]
        features_names = [[get_name_by_ind(pos_info, j, drug) for j in x] for x in features_inds]

        X_train = pd.read_csv(f'{hhs_files_dir}/{drug}.gt.domains.{thr}.train.{ind}', sep=' ',
                      index_col=None, header=None, usecols=indeces)
        y_train = pd.read_csv(f"{hhs_files_dir}/{drug}.phen.domains.{thr}.train.{ind}", sep=' ',
                   index_col=None, header=None).to_numpy().ravel()
        X_test = pd.read_csv(f'{hhs_files_dir}/{drug}.gt.domains.{thr}.test.{ind}', sep=' ',
                      index_col=None, header=None, usecols=indeces)
        y_test = pd.read_csv(f"{hhs_files_dir}/{drug}.phen.domains.{thr}.test.{ind}", sep=' ',
                   index_col=None, header=None).to_numpy().ravel()

        model = LogisticRegression(solver='saga', max_iter=10000)
        model.fit(X_train, y_train)
        y_prob = model.predict_proba(X_test)
        y_pred = model.predict(X_test)

        result.loc[k] = [drug , get_auc(y_test, y_prob), sensitivity_score(y_test, y_pred), specificity_score(y_test, y_pred),
                        npv_score(y_test, y_pred), ppv_score(y_test, y_pred)]
        k += 1

        output = pd.DataFrame()
        output['ind'] = [", ".join([str(y) for y in x]) for x in features_inds]
        output['name'] = [", ".join([str(y) for y in x]) for x in features_names]
        output['coef'] = model.coef_[0]
        output.to_csv(f"{output_dir}/{drug}.gt.domains.{thr}.final.{ind}", sep='\t', index=False, header=True) 
result = result.groupby(by='drug').mean()
result.to_csv(f"{output_dir}/auc.csv", sep='\t', index=True, header=True)
