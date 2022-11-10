# Script for running logistic regression with elastic net regularization 
# Configutation part #

# Directory with matrix after deduplication
data_dir = ""
# Directory with matrix before deduplication
data_dir_pre = ""
# Output directory
output_dir = f'./output_elastic_net'
# Gamma parameter (parameter l1_ration in the sklearn logistic regression)
a = 0.25

######################

import os
import sys

import numpy as np
import pandas as pd

from pycasso import core

from scipy.sparse import csc_matrix

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

drug = sys.argv[1]

thr=3
LAMBDA_START = 100
LAMBDA_STEP = 50

def get_auc(y_test, y_proba):
    y_proba = y_proba[:, 1]
    auc = 0
    try:
        auc = roc_auc_score(y_test, y_proba)
    except ValueError:
        pass
    return auc

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

if not os.path.exists(f"{output_dir}"):
    os.mkdir(f"{output_dir}")

if not os.path.exists(f"{output_dir}/{drug}"):
    os.mkdir(f"{output_dir}/{drug}")

metrics = pd.DataFrame(columns=['l', 'auc', 'sensitivity', 'specificity', 'npv', 'ppv'])
k = 0
lambdas = []

for ind in range(5):
    X_train = pd.read_csv(f'{data_dir}/{drug}.gt.domains.{thr}.train.{ind}', sep=' ',
                      index_col=None, header=None)
    y_train = pd.read_csv(f"{data_dir}/{drug}.phen.domains.{thr}.train.{ind}", sep=' ',
                   index_col=None, header=None).to_numpy().ravel()
    X_test = pd.read_csv(f'{data_dir}/{drug}.gt.domains.{thr}.test.{ind}', sep=' ',
                      index_col=None, header=None)
    y_test = pd.read_csv(f"{data_dir}/{drug}.phen.domains.{thr}.test.{ind}", sep=' ',
                   index_col=None, header=None)
    y_test = y_test.to_numpy()

    pos_info = pd.read_csv(f"{data_dir_pre}/{drug}.gt.pos.domains.{thr}.{ind}", sep='\t', index_col=None, header=None, usecols=[0, 1, 2, 3, 4])
    pos_info.columns = ["gene", 'pos', 'ind1', 'ind2', 'act']

    f = open(f"{data_dir}/{drug}.gt.features.domains.{thr}.train.{ind}", 'r')
    features_info = [[int(y) for y in str(x)[:-1].split(', ')] for x in f.readlines()]
    f.close()
    
    X_train = csc_matrix(X_train.to_numpy())
    X_test = csc_matrix(X_test.to_numpy())
    
    if lambdas == []:
        l = LAMBDA_START
        while True:
            model = LogisticRegression(penalty='elasticnet', C=1/(l*a), 
                                 solver='saga', max_iter=10000, 
                                l1_ratio=a)
            model.fit(X_train, y_train)
            l += LAMBDA_STEP

            if np.all(model.coef_ == 0): break
    
        lambdas = np.linspace(0, l, 22)[1:-1]

    for j, l in enumerate(lambdas):
        model = LogisticRegression(penalty='elasticnet', C=1/(l*a),
                                solver='saga', max_iter=10000,
                                l1_ratio=a)
        model.fit(X_train, y_train)

        beta = model.coef_[0]
        nnz_ind = np.argwhere(beta).transpose()[0]
        if len(nnz_ind) == 0:
            metrics.loc[k] = [l, 0.5, 0, 0, 0, 0]

            k += 1
            continue

        X_train_ind = X_train[:, nnz_ind]
        X_test_ind = X_test[:, nnz_ind]

        model = LogisticRegression(solver='saga', max_iter=10000)
        model.fit(X_train_ind, y_train)
        y_prob = model.predict_proba(X_test_ind)
        y_pred = model.predict(X_test_ind)
        metrics.loc[k] = [l, get_auc(y_test, y_prob), sensitivity_score(y_test, y_pred),
                      specificity_score(y_test, y_pred), npv_score(y_test, y_pred),
                      ppv_score(y_test, y_pred)]
        k += 1

        feature = [features_info[i] for i in nnz_ind]
        feature_name = [" ".join([get_name_by_ind(pos_info, y, drug) for y in x]) for x in feature]

        result = pd.DataFrame()
        result['ind'] = [", ".join([str(i) for i in x]) for x in feature]
        result['name'] = feature_name
        result['coef'] = model.coef_[0]
        result.to_csv(f"{output_dir}/{drug}/{l:.2f}_{a}_{ind}.csv", sep='\t', index=False, header=True)

mean = metrics.groupby(by='l').mean()
mean.columns = ['mean_auc', 'mean_sensitivity', 'mean_specificity', 'mean_npv', 'mean_ppv']
var = metrics.groupby(by='l').var(ddof=1)
var.columns = ['var_auc', 'var_sensitivity', 'var_specificity', 'var_npv', 'var_ppv']

metrics = pd.DataFrame(index=mean.index)
metrics['mean_auc'] = mean['mean_auc']
metrics['var_auc'] = var['var_auc']

metrics['mean_sensitivity'] = mean['mean_sensitivity']
metrics['var_sensitivity'] = var['var_sensitivity']

metrics['mean_specificity'] = mean['mean_specificity']
metrics['var_specificity'] = var['var_specificity']

metrics['mean_npv'] = mean['mean_npv']
metrics['var_npv'] = var['var_npv']

metrics['mean_ppv'] = mean['mean_ppv']
metrics['var_ppv'] = var['var_ppv']

metrics.to_csv(f"{output_dir}/{drug}/auc_{a}.csv", sep='\t', index=True, header=True)
