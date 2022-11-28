import sys
import pandas as pd
import numpy as np
from os.path import exists
import os
from joblib import Parallel, delayed

#path to drug phenotypes
path_to_pheno = '../db/pheno/'

def make_files_for_drug(drug, outfolder, data_folder, domain_folder, threshold):
    os.makedirs(outfolder + drug, exist_ok=True)
    phenofile = path_to_pheno + drug + '.pheno'
    phenotypes = {line.strip().split('\t')[0]: line.strip().split('\t')[1] 
                  for line in open(phenofile).readlines()}
    if not exists(outfolder + drug + '/' + drug + '.phen'):
        with open(outfolder + drug + '/' + drug + '.phen', 'w') as outf:
            for sample in phenotypes:
                outf.write(phenotypes[sample] + '\n')
    total_count = len(phenotypes)
    org_mut = {}
    mutations_dict = {}
    for organism in phenotypes:
        org_mut[organism] = []
        #annotations
        with open(data_folder + drug + '/' + organism + '_result.tsv') as nonref:
            for line in nonref:
                if ((line[0] != '-') and ((line.rfind('snp') != -1) | (line.rfind('ins') != -1) 
                | (line.rfind('del') != -1) | (line.rfind('broken') != -1) | (line.rfind('changed') != -1))):
                    if (line.rfind('snp') != -1) and (len(line.strip().split('\t')) == 6):
                        mutation = line.strip()[:-4]
                    else:
                        mutation = line.strip()
                    if mutation not in mutations_dict:
                        mutations_dict[mutation] = 1
                    else:
                        mutations_dict[mutation] += 1
                    org_mut[organism].append(mutation)
    for fold in os.listdir(domain_folder + drug):
        os.makedirs(outfolder + drug + '/' + fold, exist_ok=True)
        folder_to_write = outfolder + drug + '/' + fold + '/'
        domain_dict = {}
        domain_counts = {}
        for organism in phenotypes:
            #domains
            domain_dict[organism] = []
            if exists(domain_folder + drug + '/' + fold + '/' + organism + '_result.tsv'):
                with open(domain_folder + drug + '/' + fold + '/' + organism + '_result.tsv') as nonref:
                    for line in nonref:
                        mutation = line.strip()
                        if mutation not in domain_counts:
                            domain_counts[mutation] = 1
                        else:
                            domain_counts[mutation] += 1
                        domain_dict[organism].append(mutation)
        all_positions = []
        idx = {}
        for m in mutations_dict:
            if (mutations_dict[m] >= threshold) and (mutations_dict[m] <= (total_count - threshold)):
                all_positions.append(m)
        for m in domain_counts:
            if (domain_counts[m] >= threshold) and (domain_counts[m] <= (total_count - threshold)):
                all_positions.append(m)
        pos_dict = open(folder_to_write + drug + '.gt.pos.domains.' + str(threshold), 'w')
        for i in range(len(all_positions)):
            idx[all_positions[i]] = i
            pos_dict.write(all_positions[i] + '\n')
        pos_dict.close()
        with open(phenofile, 'r') as pheno:
            outf = open(folder_to_write + drug + '.gt.domains.' + str(threshold), 'w')
            for line in pheno:
                organism, resistance = line.strip().split('\t')
                default_l = ['0']*len(all_positions)
                for p in org_mut[organism]:
                    if p in idx:
                        default_l[idx[p]] = '1'
                for p in domain_dict[organism]:
                    if p in idx:
                        default_l[idx[p]] = '1'
                outf.write(' '.join(default_l) + '\n')
            outf.close()

if __name__ == "__main__":
    drugs = ['Isoniazid', 'Kanamycin', 'Ethambutol', 'Capreomycin', 
             'Ciprofloxacin', 'Moxifloxacin', 'Ofloxacin', 'Pyrazinamide', 
             'Rifampicin', 'Amikacin', 'Streptomycin', 'Prothionamide', 
             'Ethionamide']
    #folder with domain features, should contain separate folders for each fold
    domain_folder = '../db/domain_data'
    #folder with annotations
    data_folder = '../db/annotated_data'
    #mutation frequency threshold
    threshold = int(sys.argv[1])
    #folder for output
    outfolder = f'../db/cv_bess_files/'
    os.makedirs(outfolder, exist_ok=True)
    tasks = Parallel(n_jobs=13)(delayed(make_files_for_drug)(drug, outfolder, data_folder, domain_folder, threshold) for drug in drugs)
