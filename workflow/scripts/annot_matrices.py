import sys
import pandas
import json
from tqdm import tqdm
script, db_json, annot_db_json, cogcat_table, tax_level, traits = sys.argv

annot_cutof = 0.9

with open(db_json) as handle:
    db = json.load(handle)

print("Processing all the stuff")

def get_annot(cog, typ):
    if cog in annot_db:
        if annot_db[cog][typ]['unannotated'] < annot_cutof:
            return annot_db[cog][typ]['annotations']

def cog_set_annot(cogs, typ):
    ans = { k : get_annot(k, typ) for k in cogs}
    ans = {k : v for k,v in ans.items() if v}
    kes = {kk : vv for v in ans.values() for kk,vv in v.items()}
    annot2fract = {k : 0 for k in kes} #print(weight")
    for v in ans.values():
        for kk,vv in v.items():
            annot2fract[kk] += vv
    annot2fract = {k : v/len(cogs) for k,v in annot2fract.items()}
    return annot2fract

for typ in ['cogcat', 'KO', 'EC', "CAZy"]:
    print("Processing", typ)
    with open(annot_db_json.replace('cogcat',typ)) as handle:
        annot_db = json.load(handle)
    tables = {}
    for k,v in tqdm(db.items()):
        tables[k + ":accesorry" ] =  cog_set_annot(v['aux_genome'], typ)
        tables[k + ":core" ] =  cog_set_annot(v['core'], typ)
    out_table = pandas.DataFrame.from_dict({k : v for k,v in tables.items() if v}, orient = 'index').fillna(0)
    out_table.to_csv(cogcat_table.replace("cogcat", typ.lower()))
