#
# These functions add location inforation to TreeBreaker features and compute corresponding p-values.
# Finally, location and lineage information is intergated into final tables.
#


from multiprocessing import Pool
from os import listdir, makedirs
from random import shuffle

path_to_info_filtered = '../db/sample.ids'
path_to_ids_with_locations_extended = '../db/locations.tsv'

path_to_lineages = './lineages.csv'
path_to_tree_breaker_features = './tree_features/'
path_to_tree_breaker_features_with_locations = './tree_features_with_countries/'
path_to_tree_breaker_features_with_locations_pvalues = './tree_features_with_pvalues/'
path_to_tree_breaker_features_with_linages = './tree_features_with_lineages/'
path_to_tree_breaker_features_with_linages_pvalue = './tree_features_with_lineages_with_pvalues/'
path_to_tree_breaker_features_tables = './tree_features_tables/'


def match_sid_to_locs(country_only=True):
    err_to_samea = {}
    for l in open(path_to_info_filtered).readlines():
        s = l.strip().split('\t')
        samea = s[0]
        for err in s[1].split(';'):
            err_to_samea[err] = samea
    sid_to_loc = {}
    for l in open(path_to_ids_with_locations_extended).readlines():
        s = l.strip().split('\t')
        if s[1] == '\"missing\"' or s[1] == '-':
            continue
        if country_only:
            loc = s[1].replace('\"', '')#.replace(' ', '_')
            i = loc.find(':')
            if i != -1:
                loc = loc[:i]
        else:
            loc = s[1].replace('\"', '').replace(' ', '_')
        sid_to_loc[s[0]] = loc  
    tree_breaker_sids_to_loc = {}
    checked = set()
    for fname in listdir(path_to_tree_breaker_features):
        for l in open(path_to_tree_breaker_features + fname):
            s = l.strip().split('\t')      
            for sid in s[1].split(','):
                if sid in checked:
                    continue
                checked.add(sid)
                loc = sid_to_loc.get(sid)
                if loc is None:
                    samea = err_to_samea.get(sid)
                    if samea is None:
                        continue
                    else:
                        loc = sid_to_loc.get(samea)
                        if loc is None:
                            continue
                tree_breaker_sids_to_loc[sid] = loc                      
    return tree_breaker_sids_to_loc


def add_locations_to_treebreaker():
    makedirs(path_to_tree_breaker_features_with_locations, exist_ok=True)
    sid_to_loc = match_sid_to_locs()
    for fname in listdir(path_to_tree_breaker_features):
        with open(path_to_tree_breaker_features_with_locations + fname, 'w') as f:
            f.write('Node_ID\tR|S\tLocation:count|...\n')
            for l in open(path_to_tree_breaker_features + fname):
                s = l.strip().split('\t')
                loc_to_c = {}
                for sid in s[1].split(','):
                    loc = sid_to_loc.get(sid)
                    if loc is None:
                        continue
                    c = loc_to_c.get(loc)
                    if c is None:
                        loc_to_c[loc] = 1
                    else:
                        loc_to_c[loc] = c + 1
                loc_to_c_list = [(loc, c) for loc, c in loc_to_c.items()]
                loc_to_c_list.sort(key=lambda x: -x[1])
                loc_str = []
                for loc, c in loc_to_c_list:
                    loc_str.append(loc + ':' + str(c))
                if len(loc_str) > 0:
                    f.write(s[0] + '\t' + s[2] + '\t' + '|'.join(loc_str) + '\n')
                else:
                    f.write(s[0] + '\t' + s[2] + '\t-\n')


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
    with open(path_to_tree_breaker_features_with_locations_pvalues + fname + '.pvalues', 'w') as f:
        for j in range(len(node_to_sids)):
            node_id, true_stat = node_id_to_stat[j]
            f.write(node_id + '\t' + str(true_stat) + '\t' + str(pvalues[j]/iter_num) + '\n') 


def estimate_pvalue(iter_num=100000):

    makedirs(path_to_tree_breaker_features_with_locations_pvalues, exist_ok=True)
    sid_to_loc = match_sid_to_locs()
    pool = Pool()
    pool.starmap(pvalue, ((fname, sid_to_loc, iter_num) for fname in listdir(path_to_tree_breaker_features)))
      

def print_full_table():
    makedirs(path_to_tree_breaker_features_tables, exist_ok=True)
    for fname in listdir(path_to_tree_breaker_features_with_locations):
        if fname.endswith('.features'):
            drug = fname[:fname.index('.')]
            print(drug)
            node_to_locations = []
            node_to_lineages = {}
            node_to_loc_pvalues = {} 
            node_to_lin_pvalues = {}      
            for l in open(path_to_tree_breaker_features_with_locations + fname).readlines()[1:]:
                s = l.split('\t')
                if not s[0].startswith('Node'):
                    continue
                node_to_locations.append((s[0], (s[1] + '\t' + s[2])[:-1]))
            for l in open(path_to_tree_breaker_features_with_locations_pvalues + fname + '.pvalues').readlines():
                s = l.strip().split('\t')
                node_to_loc_pvalues[s[0]] = s[-1]    
            for l in open(path_to_tree_breaker_features_with_linages + fname).readlines():
                s = l.strip().split('\t')
                lineages = s[-1].split('|')
                lineages_filtered = [l for l in lineages if '-' not in l]
                node_to_lineages[s[0]] = '|'.join(lineages_filtered)
            for l in open(path_to_tree_breaker_features_with_linages_pvalue + fname + '.pvalues').readlines():
                s = l.strip().split('\t')
                node_to_lin_pvalues[s[0]] = s[-1]                                
            with open(path_to_tree_breaker_features_tables + fname, 'w') as f:
                f.write('Node_ID\tR|S\tLocation:count|...\tlocation p-value\tlineage:count|...\tlineage p-value\n')
                for node_id, loc in node_to_locations:
                    f.write(node_id + '\t' + loc + '\t' + str(node_to_loc_pvalues[node_id]) + '\t' + 
                    node_to_lineages[node_id] + '\t' + node_to_lin_pvalues[node_id] + '\n')


def merge_tables(big_table_path, small_table_path):
    drug_to_stats = []
    total_features = 0
    for fname in listdir(path_to_tree_breaker_features_tables):
        drug = fname[:fname.index('.')]
        if drug == 'Para-aminosalisylic_acid':
            continue
        lines = open(path_to_tree_breaker_features_tables + fname).readlines()[2:]
        total_features += len(lines)
        drug_to_stats.append((drug, lines))
    with open(big_table_path, 'w') as fbig:
        with open(small_table_path, 'w') as fsmall:
            fbig.write('Drug\tNode_ID\tR\tS\tLocation is known\tlocation p-value\tlineage p-value\n')
            fsmall.write('Drug\tTreeBreaker clades\tAssociated with location\tAssociated with lineage\n')
            for drug, lines in drug_to_stats:
                tree_features = 0
                location_is_significant = 0
                lineage_is_significant = 0
                for l in lines:
                    tree_features += 1
                    s = l.strip().split('\t')
                    rs = s[1]
                    i = rs.index('|')
                    res = [drug]
                    res.append(s[0])
                    res.append(rs[:i])
                    res.append(rs[i + 1:])
                    loc_counts = s[2].split('|')
                    c = 0
                    for lc in loc_counts:
                        i = lc.find(':')
                        if i != -1:
                            c += int(lc[i + 1:])
                    res.append(str(c))
                    loc_pval = float(s[3])*total_features
                    if loc_pval < 1:
                        res.append(f'{loc_pval:.4}')
                        if loc_pval < 0.05:
                            location_is_significant += 1
                    else:
                        res.append('1.0')
                    lin_pval = float(s[5])*total_features
                    if lin_pval < 1:
                        res.append(f'{lin_pval:.4}')
                        if lin_pval < 0.05:
                            lineage_is_significant += 1
                    else:
                        res.append('1.0')
                    fbig.write('\t'.join(res) + '\n')
                fsmall.write(f'{drug}\t{tree_features}\t{location_is_significant}\t{lineage_is_significant}\n')
                print(drug)
                print(f'total features: {tree_features}, loc significant: {location_is_significant}, lineage significant: {lineage_is_significant}')  


def merge_loc_and_lineages(output):
    sid_to_loc = {}
    for l in open(path_to_ids_with_locations_extended).readlines():
        s = l.strip().split('\t')
        if s[1] == '\"missing\"' or s[1] == '-':
            continue
        loc = s[1].replace('\"', '')
        sid_to_loc[s[0]] = loc  
    sid_to_linages = {}
    for l in open(path_to_lineages).readlines():
        s = l.strip().split('\t')
        sid_to_linages[s[0]] = s[1]
    with open(output, 'w') as f:
        f.write('BioSample\tRun ID\tBioProject\tLocation\tLineage\n')
        for l in open(path_to_info_filtered).readlines()[1:]:
            s = l.strip().split('\t')
            samea = s[0]
            loc = sid_to_loc.get(samea)
            if loc is None:
                for err in s[1].split(';'):
                    loc = sid_to_loc.get(err)
                    if loc is not None:
                        break
            lin = sid_to_linages.get(samea)
            if lin is None:
                for err in s[1].split(';'):
                    lin = sid_to_linages.get(err)
                    if lin is not None:
                        break            
            f.write(s[0] + '\t' + s[1] + '\t' + s[2] + '\t')
            if loc is None:
                f.write('-\t')
            else:
                f.write(loc + '\t')
            if lin is None:
                f.write('-\n')
            else:
                f.write(lin + '\n')            


def merge_tree_breaker_features(outpath):
    with open(outpath, 'w') as f:
        f.write('Drug\tNode_ID\tSample_ID\n')
        for fname in listdir(path_to_tree_breaker_features):
            drug = fname[:fname.index('.')]          
            for l in open(path_to_tree_breaker_features + fname).readlines()[1:]:
                if not l.startswith('Node'):
                    continue                
                s = l.strip().split('\t')
                node_id = s[0]
                for sid in s[1].split(','):
                    f.write(drug + '\t' + node_id + '\t' + sid + '\n')


if __name__ == '__main__':
    add_locations_to_treebreaker()
    estimate_pvalue()
    print_full_table()
    merge_tables('./tree_features_supplementary.table(TableS28)', './tree_features_grouped_supplementary.table(TableS26)')
    merge_loc_and_lineages('./id_with_locs_and_lineages.table(TableS27)')
    merge_tree_breaker_features('./tree_breaker.table(TableS25)')