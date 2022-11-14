import sys
import pandas as pd
import numpy as np
from os.path import exists

#path to drug phenotypes
path_to_pheno = '../db/pheno/'

if __name__ == "__main__":
    drug = sys.argv[1]
    #mutation (feature) frequency threshold
    threshold = sys.argv[2]
    #path to annotations
    data_folder = ''
    #if we should include domains
    domains = False
    path_to_domains = ''
    outfolder = f'abess_matrix_general_{drug}_thr_{threshold}_domains_{domains}/'
    mutations_dict = {}
    org_mut = {}
    total_count = 0
    with open(path_to_pheno + drug + '.pheno', 'r') as pheno:
        outf = open(outfolder + drug + '.phen', 'w')
        for line in pheno:
            organism, resistance = line.strip().split('\t')
            total_count += 1
            org_mut[organism] = []
            with open(data_folder + drug + '/' + organism + '_result.tsv') as nonref:
                for line in nonref:
                    if ((line[0] != '-') and ((line.rfind('snp') != -1) | (line.rfind('ins') != -1) 
                    | (line.rfind('del') != -1) | (line.rfind('broken') != -1))):
                        if (line.rfind('snp') != -1) and (len(line.strip().split('\t')) == 6):
                            mutation = line.strip()[:-4]
                        else:
                            mutation = line.strip()
                        if mutation not in mutations_dict:
                            mutations_dict[mutation] = 1
                        else:
                            mutations_dict[mutation] += 1
                        org_mut[organism].append(mutation)
            if domains:
                filename = path_to_domains + drug + '/' + organism + '_result.tsv'
                if exists(filename):
                    with open(filename, 'r') as f:
                        for line in f:
                            temp = line.strip().split('\t')
                            org_mut[organism].append(line.strip())
                            if line.strip() not in mutations_dict:
                                mutations_dict[line.strip()] = 1
                            else:
                                mutations_dict[line.strip()] += 1
            outf.write(resistance + '\n')
        outf.close()
    all_positions = []
    idx = {}
    for m in mutations_dict:
        if (mutations_dict[m] >= threshold) and (mutations_dict[m] <= (total_count - threshold)):
            all_positions.append(m)
    pos_dict = open(outfolder + drug + '.gt.pos.' + str(threshold), 'w')
    for i in range(len(all_positions)):
        idx[all_positions[i]] = i
        pos_dict.write(all_positions[i] + '\n')
    pos_dict.close()
    with open(path_to_pheno + drug + '.pheno', 'r') as pheno:
        outf = open(outfolder + drug + '.gt.' + str(threshold), 'w')
        for line in pheno:
            organism, resistance = line.strip().split('\t')
            default_l = ['0']*len(all_positions)
            for p in org_mut[organism]:
                if p in idx:
                    default_l[idx[p]] = '1'
            outf.write(' '.join(default_l) + '\n')
        outf.close()