import sys
import pandas
import json
from tqdm import tqdm
from ete3 import NCBITaxa
from collections import defaultdict
import os
from os.path import join as pjoin
script, ab_initio_db, emapper_db, out_db, tax_level = sys.argv

level = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
prefix = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']

level2pref = {k :v for k,v in zip(level, prefix)}

tax_db = NCBITaxa()

with open(ab_initio_db) as handle:
    stuff = json.load(handle)

genuses = {list(v['cogs'].values())[0][0].split("_silix_")[0] for k, v in stuff.items()}
del stuff

abinit_outs = []
for v in tqdm(os.walk(out_db.split("dbs/")[0] + "root", followlinks = True)):
    for vv in v[2]:
        if vv.endswith(".silix.clusters"):
            abinit_outs += [pjoin(v[0], vv)]

abinit_outs = [ l for l in abinit_outs if os.path.basename(l).split(".")[0] in genuses]

nogs2lowest = {}

nog2ko = defaultdict(lambda : [])
nog2cogcat = defaultdict(lambda : [])
nog2ec = defaultdict(lambda : [])
nog2cazy = defaultdict(lambda : [])

abinit2ko = defaultdict(lambda : [])
abinit2cogcat = defaultdict(lambda : [])
abinit2ec = defaultdict(lambda : [])
abinit2cazy = defaultdict(lambda : [])

for genus in tqdm(abinit_outs):
    with open(genus) as handle:
        gid2abinit = {kk :vv for k,v in json.load(handle).items() for kk,vv in v.items()}


    emapper_outs = []
    for v in os.walk(os.path.dirname(genus), followlinks=True):
        for vv in v[2]:
            if vv.endswith(".emapper") and not "__" in vv:
                emapper_outs += [pjoin(v[0], vv)]

    for f in emapper_outs:
        with open(f) as handle:
            res = json.load(handle)
            gid = os.path.basename(f)[:-8]
            for k,v in res.items():
                if v :
                    nogs = v["eggNOG OGs"]
                    if nogs not in nogs2lowest:
                        taxos = [v.split("@")[1] for v in nogs.split(",")]
                        tax2level = {k : len(v) for k,v in tax_db.get_lineage_translator(list(taxos)).items()}
                        lowest_nog = min([vv.split('@') for vv in nogs.split(",")], key = lambda x : tax2level.get(int(x[1]),1000))[0]
                        nogs2lowest[nogs] = lowest_nog
                    else :
                        lowest_nog = nogs2lowest[nogs]
                    nog2ko[lowest_nog] += [v['KEGG_ko']]
                    nog2cogcat[lowest_nog] += [v['COG Functional cat.']]
                    nog2ec[lowest_nog] += [v['EC']]
                    nog2cazy[lowest_nog] += [v['CAZy']]

                    if k in gid2abinit:
                        abinit2ko[gid2abinit[k]] += [v['KEGG_ko']]
                        abinit2cogcat[gid2abinit[k]] += [v['COG Functional cat.']]
                        abinit2ec[gid2abinit[k]] += [v['EC']]
                        abinit2cazy[gid2abinit[k]] += [v['CAZy']]


out_root = os.path.dirname(out_db)
def write_out_file(dicto, type, trait):
    with open(pjoin(out_root, "{}-{}.json".format(type,trait)), "w") as handle:
        json.dump(dicto, handle, indent=4, sort_keys=True)

annot_dict = defaultdict(lambda : {})
for k, v in tqdm(nog2ko.items()):
    annot_count = len(v)
    cand_counts = {vvv : 0 for vv in v for vvv in vv.split(",")}
    for vv in v:
        for vvv in vv.split(","):
            cand_counts[vvv] += 1
    cand_freq = {f : v/annot_count for f,v in cand_counts.items()}
    if cand_freq.get('', 0) != 1:
        unannt = cand_freq.get('',0)
        annot_dict[k]['KO'] = {}
        annot_dict[k]['KO']['unannotated'] = cand_freq[''] if '' in cand_freq else 0
        if '' in cand_freq:
            del cand_freq['']
        annot_dict[k]['KO']['annotations'] = {f : v/(1-unannt) for f,v in cand_freq.items()}
write_out_file(annot_dict, 'emapper', 'KO')

annot_dict = defaultdict(lambda : {})
for k, v in tqdm(nog2cogcat.items()):
    annot_count = len(v)
    cand_counts = {c : 0 for c in  set("".join(v))}
    cand_counts[''] = sum([vv == '' for vv in v])
    for vv in v:
        for vvv in vv:
            cand_counts[vvv] += 1
    cand_freq = {f : v/annot_count for f,v in cand_counts.items()}
    if cand_freq.get('', 0) != 1:
        unannt = cand_freq.get('',0)
        annot_dict[k]['cogcat'] = {}
        annot_dict[k]['cogcat']['unannotated'] = cand_freq[''] if '' in cand_freq else 0
        if '' in cand_freq:
            del cand_freq['']
        annot_dict[k]['cogcat']['annotations'] = {f : v/(1-unannt) for f,v in cand_freq.items()}
write_out_file(annot_dict, 'emapper', 'cogcat')


annot_dict = defaultdict(lambda : {})
for k, v in tqdm(nog2ec.items()):
    annot_count = len(v)
    cand_counts = {vvv : 0 for vv in v for vvv in vv.split(",")}
    for vv in v:
        for vvv in vv.split(","):
            cand_counts[vvv] += 1
    cand_freq = {f : v/annot_count for f,v in cand_counts.items()}
    if cand_freq.get('', 0) != 1:
        unannt = cand_freq.get('',0)
        annot_dict[k]['EC'] = {}
        annot_dict[k]['EC']['unannotated'] = cand_freq[''] if '' in cand_freq else 0
        if '' in cand_freq:
            del cand_freq['']
        annot_dict[k]['EC']['annotations'] = {f : v/(1-unannt) for f,v in cand_freq.items()}
write_out_file(annot_dict, 'emapper', 'EC')


annot_dict = defaultdict(lambda : {})
for k, v in tqdm(nog2cazy.items()):
    annot_count = len(v)
    cand_counts = {vvv : 0 for vv in v for vvv in vv.split(",")}
    for vv in v:
        for vvv in vv.split(","):
            cand_counts[vvv] += 1
    cand_freq = {f : v/annot_count for f,v in cand_counts.items()}
    if cand_freq.get('', 0) != 1:
        unannt = cand_freq.get('',0)
        annot_dict[k]['CAZy'] = {}
        annot_dict[k]['CAZy']['unannotated'] = cand_freq[''] if '' in cand_freq else 0
        if '' in cand_freq:
            del cand_freq['']
        annot_dict[k]['CAZy']['annotations'] = {f : v/(1-unannt) for f,v in cand_freq.items()}
write_out_file(annot_dict, 'emapper', 'CAZy')


annot_dict = defaultdict(lambda : {})
for k, v in tqdm(abinit2ko.items()):
    annot_count = len(v)
    cand_counts = {vvv : 0 for vv in v for vvv in vv.split(",")}
    for vv in v:
        for vvv in vv.split(","):
            cand_counts[vvv] += 1
    cand_freq = {f : v/annot_count for f,v in cand_counts.items()}
    if cand_freq.get('', 0) != 1:
        unannt = cand_freq.get('',0)
        annot_dict[k]['KO'] = {}
        annot_dict[k]['KO']['unannotated'] = cand_freq[''] if '' in cand_freq else 0
        if '' in cand_freq:
            del cand_freq['']
        annot_dict[k]['KO']['annotations'] = {f : v/(1-unannt) for f,v in cand_freq.items()}
write_out_file(annot_dict, 'ab_initio', 'KO')


annot_dict = defaultdict(lambda : {})
for k, v in tqdm(abinit2cogcat.items()):
    annot_count = len(v)
    cand_counts = {c : 0 for c in  set("".join(v))}
    cand_counts[''] = sum([vv == '' for vv in v])
    for vv in v:
        for vvv in vv:
            cand_counts[vvv] += 1
    cand_freq = {f : v/annot_count for f,v in cand_counts.items()}
    if cand_freq.get('', 0) != 1:
        unannt = cand_freq.get('',0)
        annot_dict[k]['cogcat'] = {}
        annot_dict[k]['cogcat']['unannotated'] = cand_freq[''] if '' in cand_freq else 0
        if '' in cand_freq:
            del cand_freq['']
        annot_dict[k]['cogcat']['annotations'] = {f : v/(1-unannt) for f,v in cand_freq.items()}
write_out_file(annot_dict, 'ab_initio', 'cogcat')

annot_dict = defaultdict(lambda : {})
for k, v in tqdm(abinit2ec.items()):
    annot_count = len(v)
    cand_counts = {vvv : 0 for vv in v for vvv in vv.split(",")}
    for vv in v:
        for vvv in vv.split(","):
            cand_counts[vvv] += 1
    cand_freq = {f : v/annot_count for f,v in cand_counts.items()}
    if cand_freq.get('', 0) != 1:
        unannt = cand_freq.get('',0)
        annot_dict[k]['EC'] = {}
        annot_dict[k]['EC']['unannotated'] = cand_freq[''] if '' in cand_freq else 0
        if '' in cand_freq:
            del cand_freq['']
        annot_dict[k]['EC']['annotations'] = {f : v/(1-unannt) for f,v in cand_freq.items()}
write_out_file(annot_dict, 'ab_initio', 'EC')

annot_dict = defaultdict(lambda : {})
for k, v in tqdm(abinit2cazy.items()):
    annot_count = len(v)
    cand_counts = {vvv : 0 for vv in v for vvv in vv.split(",")}
    for vv in v:
        for vvv in vv.split(","):
            cand_counts[vvv] += 1
    cand_freq = {f : v/annot_count for f,v in cand_counts.items()}
    if cand_freq.get('', 0) != 1:
        unannt = cand_freq.get('',0)
        annot_dict[k]['CAZy'] = {}
        annot_dict[k]['CAZy']['unannotated'] = cand_freq[''] if '' in cand_freq else 0
        if '' in cand_freq:
            del cand_freq['']
        annot_dict[k]['CAZy']['annotations'] = {f : v/(1-unannt) for f,v in cand_freq.items()}
write_out_file(annot_dict, 'ab_initio', 'CAZy')
