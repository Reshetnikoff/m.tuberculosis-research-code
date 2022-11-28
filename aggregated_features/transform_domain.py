# Script for converting PFAM domain features to different format

# Configutation part #

# Directory with phenotype data
pheno_dir = "../db/pheno"
# Directory with information about data splits
split_dir = f"./db/splits"
# Directory with domain data
root_domain_dir = "./db/domain_files_100"
# Output directory
output_dir = f"./transformed_domain_files"
# List of {gene}_{domain}
pfam_features_file=  "../db/PfamFeatures.csv"

######################

import os
import sys
import time

from multiprocessing import Pool

import pandas as pd
import numpy as np

from scipy.sparse import save_npz, csc_matrix 

start_time = time.time()

drug = sys.argv[1]
split_mode = sys.argv[2]

output_dir = f"{output_dir}/{drug}"
split_dir = f"{output_dir}/{drug}"

if split_mode == "1":
    splits = list(range(1, 17))
elif split_mode == "2":
    splits = list(range(17, 33))
elif split_mode == "3":
    splits = list(range(33, 49))
elif split_mode == "4":
    splits = list(range(49, 65))
elif split_mode == "5":
    splits = list(range(65, 81))
elif split_mode == "6":
    splits = list(range(81, 97))
elif split_mode == "7":
    splits = list(range(97, 101))
elif split_mode == '8':
    splits = [1, 2]


if not os.path.exists(output_dir):
    os.makedirs(output_dir)

def transform_domain(split):
    pheno = pd.read_csv(f"{pheno_dir}/{drug}.pheno", sep='\t', index_col=0, header=None)

    if not os.path.exists(f"{output_dir}/{split}"):
        os.mkdir(f"{output_dir}/{split}")
    
    for fold in range(1, 6):
        start_fold_time = time.time()

        f = open(f"{split_dir}/{drug}.{split}.fold_{fold}.train_indices")
        train_indeces = [int(x[:-1]) if x[-1] == '\n' else int(x) for x in f.readlines()]
        f.close()
        test_indeces = [int(x) for x in range(len(pheno.index)) if x not in train_indeces]

        domain_dir = f"{root_domain_dir}/{drug}/{split}/fold{fold}"

        f = open(pfam_features_file)
        features = [x[:-1] if x[-1] == '\n' else x for x in f.readlines()]
        f.close()

        domain_data = pd.DataFrame(index=features)
        samples_dir = os.listdir(domain_dir)

        for sample_name in pheno.index:
            sample = f"{sample_name}_result.tsv"

            f = open(f"{domain_dir}/{sample}")
            temp = [x.split()[0] for x in f.readlines()  if x.split()[0] != '\n']
            f.close()

            domain_data[sample_name] = 0
            domain_data.loc[temp,sample_name] = 1

        domain_data = domain_data.loc[~(domain_data==0).all(axis=1)]
        features = list(domain_data.index)

        pos_domain_data = pd.DataFrame([[x, 'changed'] for x in features])
        pos_domain_data.to_csv(f"{output_dir}/{split}/{drug}.gt.pos.domains.3.fold.{fold}", sep='\t', index=False, header=False)

        domain_data = domain_data.transpose()        

        output_data_train = csc_matrix(domain_data.iloc[train_indeces])
        output_data_test = csc_matrix(domain_data.iloc[test_indeces])

        save_npz(f"{output_dir}/{split}/{drug}.gt.domains.3.fold.{fold}.train.npz", output_data_train)
        save_npz(f"{output_dir}/{split}/{drug}.gt.domains.3.fold.{fold}.test.npz", output_data_test)

        pheno.iloc[train_indeces].to_csv(f"{output_dir}/{split}/{drug}.pheno.fold.{fold}.train",
                        sep=' ', index=False, header=False)
        pheno.iloc[test_indeces].to_csv(f"{output_dir}/{split}/{drug}.pheno.fold.{fold}.test",
                        sep=' ', index=False, header=False)

        end_fold_time = time.time()
        print(f"Fold {fold} done for {drug}: {end_fold_time - start_fold_time}")


p = Pool(len(splits))
p.map(transform_domain, splits)

end_time = time.time()
print(end_time - start_time)
