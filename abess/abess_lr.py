# The script apply logistic regression on non-zero selected features by ABESS

# Configutation part #

# Directory with matrix after deduplication
data_dir = ""
# Directory with output of ABESS algorithm (output of run_abess.R)
result_dir = "./output_abess"
# Output directory: {output_dir}/{exp_name}. Leave exp_name variable is unchanged 
output_dir = f"./output_abess_lr"
exp_name = "lr"

######################

import os

import pandas as pd
import numpy as np

from multiprocessing import Pool

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score

thr = "3"
folds = 5

METHOD = "LR"

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
        return f"{gene}_{pos}_{ind1}"

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

def convert_pd_to_array(x):
    x = str(x).split(", ")
    x = [int(i[1:])-1 for i in x]
    return x


def get_auc_dict(y_test, y_proba):
    auc = 0
    try:
        auc = roc_auc_score(y_test, y_proba)
    except ValueError:
        pass
    return auc

def cv_abbes(drug):
    if not os.path.exists(f"{output_dir}/{exp_name}/{drug}"):
        os.mkdir(f"{output_dir}/{exp_name}/{drug}")

    result = dict()
    result[drug] = dict()
    for method in ['gic']:  
        result[drug][method] = dict()
        result[drug][method]['selection'] = dict()
        
        auc_scores = []
        sensitivity_scores = []
        specivity_scores = []
        npv_scores = []
        ppv_scores = []

        for ind in range(folds):
            if not os.path.exists(f"{result_dir}/{drug}.domains.{thr}.result.{ind}"): continue
            if os.stat(f"{result_dir}/{drug}.domains.3.result.{ind}").st_size == 0: continue
            # if not os.path.exists(f"{data_dir}/{drug}.gt.domains.{thr}.test.{ind}"): continue
            selection = pd.read_csv(f"{result_dir}/{drug}.domains.{thr}.result.{ind}", sep='\t', index_col=0, header=None,
                                    usecols=[0, 1, 2])

            indexes = convert_pd_to_array(selection[selection.index == method][1][0])
            
            X_train = pd.read_csv(f"{data_dir}/{drug}.gt.domains.{thr}.train.{ind}", sep=' ', 
                                   index_col=None, header=None, usecols=indexes)
            X_test = pd.read_csv(f"{data_dir}/{drug}.gt.domains.{thr}.test.{ind}", sep=' ', 
                                  index_col=None, header=None, usecols=indexes)

            y_train = pd.read_csv(f"{data_dir}/{drug}.phen.domains.{thr}.train.{ind}", sep=' ', index_col=None, header=None).to_numpy().ravel()
            y_test = pd.read_csv(f"{data_dir}/{drug}.phen.domains.{thr}.test.{ind}", sep=' ', index_col=None, header=None).to_numpy().ravel()
            
            pos_info = pd.read_csv(f"{data_dir}/{drug}.gt.pos.domains.{thr}.{ind}", sep='\t', index_col=None, header=None, usecols=[0, 1, 2, 3, 4])
            pos_info.columns = ["gene", 'pos', 'ind1', 'ind2', 'act']
            f = open(f"{data_dir}/{drug}.gt.features.domains.{thr}.train.{ind}", 'r')
            features_info = [[int(y) for y in str(x)[:-1].split(', ')] for x in f.readlines()]
            f.close() 
            
            features = [features_info[i] for i in indexes]

            X_train_selection = X_train.loc[:, indexes]
            X_train_selection.columns = [get_name_by_ind(pos_info, i, drug) for i in indexes]

            X_test_selection = X_test.loc[:, indexes]
            X_test_selection.columns = [get_name_by_ind(pos_info, i, drug) for i in indexes]

            X_train_selection = X_train_selection.loc[:, (X_train_selection != 0).any(axis=0)]
            X_test_selection = X_test_selection.loc[:, (X_train_selection != 0).any(axis=0)]

            if METHOD == 'DICT':
                y_pred = [int(x) for x in (X_test != 0).any(axis=1)]
            elif METHOD == 'LR':
                model = LogisticRegression(solver='saga', max_iter=10000)
                model.fit(X_train_selection, y_train)
                y_proba = model.predict_proba(X_test_selection)[:, 1]
                y_pred = model.predict(X_test_selection)
                coef = pd.DataFrame(data=model.coef_[0], 
                                    index=indexes, 
                                    columns=["coef"])
                coef['ext_feature'] = [" ".join([get_name_by_ind(pos_info, i, drug) for i in x]) for x in features]

                coef.to_csv(f"{output_dir}/{exp_name}/{drug}/{method}_{ind}.csv", 
                            sep='\t', header=True, index=True)

            auc_scores.append(get_auc_dict(y_test, y_proba))
            sensitivity_scores.append(sensitivity_score(y_test, y_pred))
            specivity_scores.append(specificity_score(y_test, y_pred))
            npv_scores.append(npv_score(y_test, y_pred))
            ppv_scores.append(ppv_score(y_test, y_pred))

        result[drug][method]['mean_auc'] = np.mean(auc_scores)
        result[drug][method]['var_auc'] = np.var(auc_scores, ddof=1)
        
        result[drug][method]['mean_sensitivity'] = np.mean(sensitivity_scores)
        result[drug][method]['var_sensitivity'] = np.var(sensitivity_scores, ddof=1)

        result[drug][method]['mean_specivity'] = np.mean(specivity_scores)
        result[drug][method]['var_specivity'] = np.var(specivity_scores, ddof=1)

        result[drug][method]['mean_npv'] = np.mean(npv_scores)
        result[drug][method]['var_npv'] = np.var(npv_scores, ddof=1)

        result[drug][method]['mean_ppv'] = np.mean(ppv_scores)
        result[drug][method]['var_ppv'] = np.var(ppv_scores, ddof=1)
    return result

if not os.path.exists(output_dir):
        os.mkdir(output_dir)

if not os.path.exists(f"{output_dir}/{exp_name}"):
    os.mkdir(f"{output_dir}/{exp_name}")

p = Pool(len(drugs))
results = p.map(cv_abbes, drugs)

i = 0
result_pd = pd.DataFrame(columns=['drug', 'method', 'mean_auc', 'var_auc', 'mean_sensitivity', 'var_sensitivity', 
                                  'mean_specivity', 'var_specivity', 'mean_npv', 'var_npv', 'mean_ppv', 'var_ppv'])

for j, result in enumerate(results):
    drug = list(result.keys())[0]
    for method in result[drug].keys():
        result_pd.loc[i] = [drug, method, result[drug][method]['mean_auc'], result[drug][method]['var_auc'],
                            result[drug][method]['mean_sensitivity'], result[drug][method]['var_sensitivity'], 
                            result[drug][method]['mean_specivity'], result[drug][method]['var_specivity'], 
                            result[drug][method]['mean_npv'], result[drug][method]['var_npv'],
                            result[drug][method]['mean_ppv'], result[drug][method]['var_ppv']]
        i += 1

result_pd.to_csv(f"{output_dir}/{exp_name}/auc_result.csv", sep='\t', header=True, index=False)

