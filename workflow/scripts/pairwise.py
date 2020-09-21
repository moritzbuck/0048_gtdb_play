import sys
import pandas
import json
from tqdm import tqdm
from numpy import abs, log, sqrt
script , bmft_file , dist_json , pairwise_csv, tax_level, traits = sys.argv

bmft = pandas.read_csv(bmft_file, index_col=0)
level = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
radius = 3

radius_rank = [ (i-radius) if (i-radius) > 0 else 0 for i, t in enumerate(level) if t == tax_level][0]
radius_level = level[radius_rank]

with open(dist_json) as handle:
    distances = json.load(handle)
distances = { tuple(l[0:2]): {'phylo_dist' : l[2]} for l in distances if l[0] != l[1] }

for k,v in tqdm(distances.items()):
    for tax_l in level:
        v[tax_l] = bmft.loc[k[0], tax_l] if bmft.loc[k[0], tax_l] == bmft.loc[k[1], tax_l] else 'NA'
    v['genome_ratio'] = bmft.loc[k[0], 'mean_genome_size']/bmft.loc[k[1], 'mean_genome_size']
    v['varfract_ratio'] = bmft.loc[k[0], 'mean_varfract']/bmft.loc[k[1], 'mean_varfract']
    v['cogs_ratio'] = bmft.loc[k[0], 'mean_cogs']/bmft.loc[k[1], 'mean_cogs']
    v['GC_ratio'] = bmft.loc[k[0], 'mean_GC']/bmft.loc[k[1], 'mean_GC']
    v['core_ratio'] = bmft.loc[k[0], 'core_len']/bmft.loc[k[1], 'core_len']
    v['access_ratio'] = bmft.loc[k[0], 'accessory_len']/bmft.loc[k[1], 'accessory_len']
    v['normaccess_ratio'] = (bmft.loc[k[0], 'accessory_len']/log(bmft.loc[k[0], 'nb_estgenomes']))/(bmft.loc[k[1], 'accessory_len']/log(bmft.loc[k[1], 'nb_estgenomes']))
    v['variable_ratio'] = bmft.loc[k[0], 'mean_variable']/bmft.loc[k[1], 'mean_variable']

pandas.DataFrame.from_dict(distances, orient = "index").to_csv(pairwise_csv, index_label=["k","l"])
