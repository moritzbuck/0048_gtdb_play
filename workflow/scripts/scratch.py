from numpy import mean
from tqdm import tqdm
import pandas
import os
from os.path import join as pjoin
import json

motupan_outs = []
for v in tqdm(os.walk("root", followlinks=True)):
    for vv in v[2]:
        if vv.endswith( ".motupan.json"):
            motupan_outs += [pjoin(v[0], vv)]

database = {}
for p in tqdm(motupan_outs):
    with open(p) as handle:
        tt = json.load(handle)
        for k,v in tt.items():
          stats = dict()
          stats['fpr'] = v['fpr']
          stats['recall'] = v['recall']
          stats['lowest_false'] = v['lowest_false']
          stats['core_len'] = 0 if v['core']  is None else len(v['core'])
          stats['nb_genomes'] = v['nb_genomes']
          stats['mean_new_completeness'] = mean([g['new_completness'] for g  in v['genomes']])
          stats['mean_checkm_completeness'] = mean([g['checkm_complet'] for g  in v['genomes']])
          database[k] = stats

pandas.DataFrame.from_dict(database, orient="index").to_csv("/home/moritz/temp/test.csv")
