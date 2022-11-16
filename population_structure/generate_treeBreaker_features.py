#
# These functions run TreeBreaker for given trees, parse the results and print population_structure features
#
# TreeBreaker should be installed preliinary and added to the PATH

from multiprocessing import Pool
from multiprocessing.pool import ThreadPool
from os import listdir, makedirs
from subprocess import call
from ete3 import Tree

trees_path = '../db/trees/'
phenotypes_path = '../db/pheno/'
tree_breaker_output_path = './tree_breaker_output/'
tree_features_path = './tree_features/'

min_posterior_for_TreeBreaker = 0.5
drugs = ['Rifampicin', 'Isoniazid', 'Ethambutol', 'Pyrazinamide', 'Streptomycin', 'Capreomycin', 'Kanamycin', 'Amikacin', 'Ofloxacin', 'Ciprofloxacin', 'Moxifloxacin', 'Prothionamide', 'Ethionamide']


def run_tree_braker(drug):
    call('treeBreaker ' + trees_path + drug + '.nw ' + phenotypes_path + drug + '.pheno ' + tree_breaker_output_path + drug + '.out', shell=True)
    return 1


def run_tree_breaker_for_all_drugs():
    makedirs(tree_breaker_output_path, exist_ok=True)

    pool = ThreadPool()
    tasks = pool.map(run_tree_braker, (drug for drug in drugs))
    c = 0
    for task in tasks:
        c += task
    print(f'Finished TreeBreaker for {c} drugs')


def gen_features_for_drug_with_names(filename):
    features = []
    drug = filename.split('.')[0]
    sid_to_pheno = {}
    for l in open(phenotypes_path + drug + '.pheno').readlines():
        s = l.strip().split('\t')
        sid_to_pheno[s[0]] = int(s[1])
    with open(tree_breaker_output_path + filename, 'r') as f:
        marked_tree = Tree(f.readlines()[-1], format=1)
        for node in marked_tree.traverse():
            i = node.name.index('{')
            node_id = node.name[:i]
            prob = float(node.name[node.name.rfind('=') + 1: -1])
            if prob > min_posterior_for_TreeBreaker:
                nodes = []
                r = 0
                for n in node.get_leaves():
                    sid = n.name.split('{')[0]
                    if sid_to_pheno[sid] == 1:
                        r += 1
                    nodes.append(sid)
                features.append((node_id, ','.join(nodes), r, len(nodes) - r))
    with open(tree_features_path + drug + '.features', 'w') as f:
        for node_id, node_list, r_num, s_num in features:
            f.write(f'{node_id}\t{node_list}\t{r_num}|{s_num}\n')
    return 1


def gen_features_for_all_trees():
    makedirs(tree_features_path, exist_ok=True)

    pool = Pool()
    tasks = pool.map(gen_features_for_drug_with_names, (fname for fname in listdir(tree_breaker_output_path)))
    sum = 0
    for c in tasks:
        sum += c
    print('Generated features for ' + str(sum) + ' drugs')


if __name__ == '__main__':
    run_tree_breaker_for_all_drugs()
    gen_features_for_all_trees()
