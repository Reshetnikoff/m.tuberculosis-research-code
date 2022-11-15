# Script to make matrices for next iteration of ABESS
# The script remove explained resistance samples by logistic regression 
# trained on majorly selected features from previous iteration
# and remove random susceptible samples to preserve the balance of classes

# Configutation part #

# Directory with matrix before deduplication
data_dir_pre = ""
# Directory with matrix after deduplication
data_dir = "2"
# Output of converter.py
selection_dir = "./output_abess"
# Output directory to save next iteration matrices
output_dir = '../db/cv_files_iter_2'
# Directory with stored informations about number of explained samples
result_dir = "./output_iterative_abess/explained_iter_1"

######################

import sys
import os

import pandas as pd
import numpy as np

from sklearn.linear_model import LogisticRegression

drug = sys.argv[1]

if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    
if not os.path.exists(result_dir):
    os.makedirs(result_dir)

drugs = ['Kanamycin', 'Amikacin', 'Streptomycin', 'Ofloxacin',
         'Moxifloxacin', 'Isoniazid', 'Rifampicin', 'Ethambutol',
         'Pyrazinamide', 'Capreomycin', 'Ethionamide', 'Prothionamide',
         'Ciprofloxacin']

def represents_int(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False

def get_ind_by_name(pos_info, name):
    if name.split('_')[-1] == 'domain':
        gene = "_".join(name.split('_')[:-1])
        temp = pos_info[(pos_info['gene'] == gene) & (pos_info['pos'] == 'changed')]
    elif name.split('_')[-1] == 'broken':
        gene = "_".join(name.split('_')[:-1])
        temp = pos_info[(pos_info['gene'] == gene) & (pos_info['pos'] == 'broken')]
    elif name.split('_')[-1] == 'snp':
        gene = "_".join(name.split('_')[:-3])
        pos = name.split('_')[-3]
        ind = name.split('_')[-2]
        temp = pos_info[(pos_info['gene'] == gene) & (pos_info['pos'] == pos) & (pos_info['ind2'] == ind) & (pos_info['act'] == 'snp')]
    elif name.split('_')[-1] == 'del':
        gene = "_".join(name.split('_')[:-3])
        pos = name.split('_')[-3]
        ind = name.split('_')[-2]
        if represents_int(ind):
            temp = pos_info[(pos_info['gene'] == gene) & (pos_info['pos'] == pos) & (pos_info['ind1'] == ind)  & (pos_info['ind2'] == 'del')]
        else:
            temp = pos_info[(pos_info['gene'] == gene) & (pos_info['pos'] == pos) & (pos_info['ind2'] == ind)  & (pos_info['act'] == 'del')]
    elif name.split('_')[-1] == 'ins':
        gene = "_".join(name.split('_')[:-3])
        pos = name.split('_')[-3]
        ind = name.split('_')[-2]
        temp = pos_info[(pos_info['gene'] == gene) & (pos_info['pos'] == pos) & (pos_info['ind2'] == ind)  & (pos_info['act'] == 'del')]
    assert len(temp) <= 1, "In function get_ind_by_name by name given more or less than one ind"
    if len(temp) == 0: return -1
    else: return temp.index[0]


main_selection = pd.read_csv(f"{selection_dir}/{drug}.domains.3.output", sep='\t', index_col=None, header=0)
main_selection = main_selection[main_selection['sum'] >= 3]

result = pd.DataFrame(index=[0, 1, 2, 3, 4], columns=['old', 'new'])

for i in range(5):
    pos_info = pd.read_csv(f"{data_dir_pre}/{drug}.gt.pos.domains.3.{i}", sep='\t', index_col=None, header=None)
    pos_info.columns = ['gene', 'pos', 'ind1', 'ind2', 'act']
    indeces = [get_ind_by_name(pos_info, x.split(', ')[0]) for x in main_selection['ext_desc']]
    while -1 in indeces: indeces.remove(-1)
    X_train = pd.read_csv(f"{data_dir_pre}/{drug}.gt.domains.3.train.{i}", sep=' ', 
                    index_col=None, header=None)
    y_train = pd.read_csv(f"{data_dir_pre}/{drug}.phen.domains.3.train.{i}", sep=' ', 
                    index_col=None, header=None)
    X_train['y'] = y_train

    X_test = pd.read_csv(f"{data_dir_pre}/{drug}.gt.domains.3.test.{i}", sep=' ',
                    index_col=None, header=None)
    y_test = pd.read_csv(f"{data_dir_pre}/{drug}.phen.domains.3.test.{i}", sep=' ',
                    index_col=None, header=None)
    X_test['y'] = y_test

    if len(indeces) != 0:
        model = LogisticRegression(max_iter=10000, solver='saga', penalty='l1')
        model.fit(X_train[indeces], y_train.to_numpy().ravel())
        
        y_pred_train = model.predict(X_train[indeces]).ravel()
        drop_samples_indeces = list(X_train[(X_train['y'] == y_pred_train) & (X_train['y'] == 1)].index)
        new_X_train = X_train[~X_train.index.isin(drop_samples_indeces)]
        new_X_train_1 = new_X_train[new_X_train['y'] == 1]
        new_X_train_0 = new_X_train[new_X_train['y'] == 0]
        new_X_train_0 = new_X_train_0.sample(int(len(new_X_train_0) * len(new_X_train_1) / np.sum(X_train['y'])))
        new_X_train = pd.concat([new_X_train_0, new_X_train_1], axis=0)

        y_pred_test = model.predict(X_test[indeces]).ravel()
        drop_samples_indeces = list(X_test[(X_test['y'] == y_pred_test) & (X_test['y'] == 1)].index)
        new_X_test = X_test[~X_test.index.isin(drop_samples_indeces)]
        new_X_test_1 = new_X_test[new_X_test['y'] == 1]
        new_X_test_0 = new_X_test[new_X_test['y'] == 0]
        new_X_test_0 = new_X_test_0.sample(int(len(new_X_test_0) * len(new_X_test_1) / np.sum(X_test['y'])))
        new_X_test = pd.concat([new_X_test_0, new_X_test_1], axis=0)
        
        result.loc[i] = [np.sum(y_train)[0], np.sum(new_X_train['y'])]
    else:
        result.loc[i] = [np.sum(y_train)[0], np.sum(y_train)[0]]
        new_X_train = X_train
        new_X_test = X_test


    full_indeces = []
    for d in drugs:
        selection = pd.read_csv(f"{selection_dir}/{d}.domains.3.output", sep='\t', index_col=None, header=0)
        selection = selection[selection['sum'] >= 3]
        selection = selection[selection['coef'] >= 0.]
        full_indeces.extend([get_ind_by_name(pos_info, x.split(', ')[0]) for x in selection['ext_desc']])
    full_indeces = list(set(full_indeces))
    while -1 in full_indeces: full_indeces.remove(-1)    

    new_X_train.drop(full_indeces, axis=1, inplace=True)
    new_X_test.drop(full_indeces, axis=1, inplace=True)
    pos_info.drop(full_indeces, axis=0, inplace=True)

    new_X_train = new_X_train.sample(frac=1)
    new_y_train = new_X_train['y']
    new_X_train.drop('y', axis=1, inplace=True)
    new_X_test = new_X_test.sample(frac=1)
    new_y_test = new_X_test['y']
    new_X_test.drop('y', axis=1, inplace=True)
 
    true_indeces = [bool(x) for x in (np.sum(new_X_train, axis=0) >= 3)]
    new_pos_info = pos_info.loc[true_indeces]
    new_X_train = new_X_train.loc[:, true_indeces]
    new_X_test = new_X_test.loc[:, true_indeces]

    new_X_train.to_csv(f"{output_dir}/{drug}.gt.domains.3.train.{i}", sep=' ', index=False, header=False)
    new_y_train.to_csv(f"{output_dir}/{drug}.phen.domains.3.train.{i}", sep=' ', index=False, header=False)
    new_pos_info.to_csv(f"{output_dir}/{drug}.gt.pos.domains.3.{i}", sep='\t', index=False, header=False)
    new_X_test.to_csv(f"{output_dir}/{drug}.gt.domains.3.test.{i}", sep=' ', index=False, header=False)
    new_y_test.to_csv(f"{output_dir}/{drug}.phen.domains.3.test.{i}", sep=' ', index=False, header=False)

result.to_csv(f"{result_dir}/{drug}.csv", sep='\t', index=False, header=True)
