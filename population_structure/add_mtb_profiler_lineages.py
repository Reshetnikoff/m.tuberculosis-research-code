from multiprocessing import Pool
from os import listdir, makedirs
from Bio import SeqIO
from random import shuffle

path_to_ref = '../db/h37rv.fasta'
path_to_dict = '../db/tbdb.barcode.bed'
path_to_dict_converted = './tbdb.barcode.bed.converted'
path_to_snps = '../db/nucl_data/'

path_to_lineages = './lineages.csv'
path_to_tree_breaker_features = './tree_features/'
path_to_tree_breaker_features_with_linages = './tree_features_with_lineages/'
path_to_tree_breaker_features_with_linages_pvalue = './tree_features_with_lineages_with_pvalues/'


def convert_dict():
    ref_seq = str(SeqIO.read(path_to_ref, 'fasta').seq).upper()
    ref_len = len(ref_seq)
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', '-': '-'}
    with open(path_to_dict_converted, 'w') as f:
        for l in open(path_to_dict).readlines()[4:]:
            s = l.strip().split('\t')
            pos = ref_len - int(s[1])
            lineage = s[3]
            if lineage.startswith('lineage'):
                lineage = lineage[7:]
            f.write(f'{lineage} {pos} {ref_seq[pos - 1]}->{complement[s[4]]}\n')
                 

def compute_lineage():
    mut_to_lineage = {}
    ref_to_lineage = []
    for l in open(path_to_dict_converted).readlines():
        s = l.strip().split(' ')
        if s[2][0] == s[2][3]:
            ref_to_lineage.append((s[1], s[0]))
            continue
        mut = s[1] + '\t' + s[2][0] + '\t' + s[2][3]
        mut_to_lineage[mut] = s[0]
        mut = s[1] + '\t' + s[2][3] + '\t' + s[2][0]
        mut_to_lineage[mut] = s[0]        
    with open(path_to_lineages, 'w') as f:
        for fname in listdir(path_to_snps):
            i = fname.find('.variants')
            if i != -1:
                sid = fname[: i]
                f.write(sid + '\t')
                lineages = set()
                mut_poses = set()
                for l in open(path_to_snps + fname).readlines():
                    l = l.strip()
                    s = l.split('\t')
                    lineage = mut_to_lineage.get(l)
                    if lineage is not None:
                        lineages.add(lineage)
                    mut_poses.add(s[0])
                for pos, ln in ref_to_lineage:
                    if pos not in mut_poses:
                        lineages.add(ln)
                max_len = 0              
                longest_linages = set()
                for lineage in lineages:
                    c = 1
                    if lineage.startswith('M.'):
                        if c > max_len:
                            max_len = c
                            longest_linages = set()
                            longest_linages.add(lineage)
                        continue                        
                    for i in range(len(lineage)):
                        if lineage[i] == '.':
                            if lineage[:i] not in lineages:
                                break
                            else:
                                c += 1
                    else:
                        if c > max_len:
                            max_len = c
                            longest_linages = set()
                            longest_linages.add(lineage)
                        elif c == max_len:
                            longest_linages.add(lineage)
                if len(longest_linages) > 0:
                    f.write(';'.join(longest_linages) + '\n')
                else:
                    f.write('-\n')


def add_lineages_to_treebreaker():
    makedirs(path_to_tree_breaker_features_with_linages, exist_ok=True)
    sid_to_linages = {}
    for l in open(path_to_lineages).readlines():
        s = l.strip().split('\t')
        sid_to_linages[s[0]] = s[1].split(';')
    for fname in listdir(path_to_tree_breaker_features):
        with open(path_to_tree_breaker_features_with_linages + fname, 'w') as f:
            for l in open(path_to_tree_breaker_features + fname):
                s = l.strip().split('\t')
                lin_to_c = {}
                for sid in s[1].split(','):
                    lines = sid_to_linages[sid]
                    for line in lines:
                        c = lin_to_c.get(line)
                        if c is None:
                            lin_to_c[line] = 1
                        else:
                            lin_to_c[line] = c + 1
                lineage_str = []
                for lin, c in lin_to_c.items():
                    lineage_str.append(lin + ':' + str(c))
                f.write(s[0] + '\t' + s[2] + '\t' + '|'.join(lineage_str) + '\n')


def match_sid_to_lineages():
    sid_to_lin = {}
    for l in open(path_to_lineages).readlines():
        s = l.strip().split('\t')
        sid_to_lin[s[0]] = s[1]  
    tree_breaker_sids_to_loc = {}
    checked = set()
    for fname in listdir(path_to_tree_breaker_features):
        for l in open(path_to_tree_breaker_features + fname):
            s = l.strip().split('\t')      
            for sid in s[1].split(','):
                if sid in checked:
                    continue
                checked.add(sid)
                loc = sid_to_lin.get(sid)
                if loc is None:
                    continue
                tree_breaker_sids_to_loc[sid] = loc                      
    return tree_breaker_sids_to_loc


def pvalue(fname, sid_to_loc, iter_num):

    def compute_total_freqs(sid_to_loc, node_to_sids):
        sid_to_loc_filtered = {}
        checked = set()
        for _, sids in node_to_sids:
            for sid in sids:
                if sid in checked:
                    continue
                checked.add(sid)
                loc = sid_to_loc.get(sid)
                if loc is None:
                    continue
                sid_to_loc_filtered[sid] = loc
        loc_to_c = {}
        total_sids = 0
        for loc in sid_to_loc_filtered.values():
            c = loc_to_c.get(loc)
            if c is None:                
                loc_to_c[loc] = 1
            else:
                loc_to_c[loc] = c + 1
            total_sids += 1
        loc_freqs = []
        for loc, c in loc_to_c.items():
            loc_freqs.append((loc, c/total_sids))
        return loc_freqs

    def compute_distance(sid_to_loc, node_to_sids, loc_freqs):
        node_id_to_distance = []
        loc_to_c = {}
        for node_id, sids in node_to_sids:
            loc_to_c = {}
            for sid in sids:
                loc = sid_to_loc.get(sid)
                if loc is None:
                    continue
                c = loc_to_c.get(loc)
                if c is None:
                    loc_to_c[loc] = 1
                else:
                    loc_to_c[loc] = c + 1
            total = sum(c for c in loc_to_c.values())
            dist = 0
            for loc, total_freq in loc_freqs:
                c = loc_to_c.get(loc)
                if c is None:
                    dist += total_freq*total_freq
                else:
                    dist += (total_freq - c/total)*(total_freq - c/total)
            node_id_to_distance.append((node_id, dist))
        return node_id_to_distance

    def shuffle_locs(sids, locs):
        shuffle(locs)
        return {sids[i]: locs[i] for i in range(len(sids))}

    node_to_sids = []
    sids = []
    locs = []
    s_to_l = {}
    for l in open(path_to_tree_breaker_features + fname):
        s = l.strip().split('\t')
        sid_list = s[1].split(',')
        node_to_sids.append((s[0], sid_list))
        for sid in sid_list:
            loc = sid_to_loc.get(sid)
            if loc is None:
                continue
            s_to_l[sid] = loc
    for s, l in s_to_l.items():
        sids.append(s)
        locs.append(l)
    loc_freqs = compute_total_freqs(sid_to_loc, node_to_sids)
    node_id_to_stat = compute_distance(sid_to_loc, node_to_sids, loc_freqs)
    pvalues = [0.0 for j in range(len(node_to_sids))]
    for i in range(iter_num):
        fake_node_id_to_stat = compute_distance(shuffle_locs(sids, locs), node_to_sids, loc_freqs)
        for j in range(len(node_to_sids)):
            _, true_stat = node_id_to_stat[j]
            _, fake_stat = fake_node_id_to_stat[j]
            if fake_stat >= true_stat:
                pvalues[j] += 1.0
    with open(path_to_tree_breaker_features_with_linages_pvalue + fname + '.pvalues', 'w') as f:
        for j in range(len(node_to_sids)):
            node_id, true_stat = node_id_to_stat[j]
            f.write(node_id + '\t' + str(true_stat) + '\t' + str(pvalues[j]/iter_num) + '\n')     


def estimate_pvalues(iter_num=100000):

    makedirs(path_to_tree_breaker_features_with_linages_pvalue, exist_ok=True)
    sid_to_loc = match_sid_to_lineages()
    pool = Pool()
    pool.starmap(pvalue, ((fname, sid_to_loc, iter_num) for fname in listdir(path_to_tree_breaker_features)))

        

if __name__ == '__main__':
    convert_dict()
    compute_lineage()
    add_lineages_to_treebreaker()
    estimate_pvalues()
