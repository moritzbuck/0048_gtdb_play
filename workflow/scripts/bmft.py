import sys
import pandas
import json
from tqdm import tqdm
script, db_json, output_csv, tax_level, traits = sys.argv

print("Loading", tax_level, "db for", traits, "traits")

level = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
not_none = lambda l : [ll for ll in l if ll is not None]
mean = lambda l : sum(not_none(l))/len(not_none(l))
def sd(l):
    m = mean(l)
    return mean([(v-m)**2 for v in not_none(l)])**0.5

with open(db_json) as handle:
    db = json.load(handle)

print("Processing all the stuff")

bmft = {}
for k,v in tqdm(db.items()):
    bmft[k] = {}
    bmft[k]['core_len'] = len(v['core'])
    bmft[k]['accessory_len'] = len(v['aux_genome'])
    bmft[k]['nb_singletons'] = len(v['singleton_cogs'])
    bmft[k]['mean_completeness'] = mean([vv['new_completness'] for vv in v['genomes']])
    bmft[k]['mean_variable'] = mean([vv['acc_cogs_count']*100/vv['new_completness'] for vv in v['genomes']])
    bmft[k]['sd_variable'] = sd([vv['acc_cogs_count']*100/vv['new_completness'] for vv in v['genomes']])
    bmft[k]['mean_cogs'] = mean([vv['total_cogs']*100/vv['new_completness'] for vv in v['genomes']])
    bmft[k]['mean_genome_size'] = mean([vv['genome_len']*100/vv['new_completness'] for vv in v['genomes'] if 'genome_len' in vv])
    bmft[k]['mean_GC'] = mean([vv['GC'] for vv in v['genomes'] if 'GC' in vv])
    bmft[k]['mean_varfract'] = mean([vv['varfract'] for vv in v['genomes']])
    bmft[k]['sd_varfract'] = sd([vv['varfract'] for vv in v['genomes']])
    bmft[k]['nb_genomes'] = v['nb_genomes']
    bmft[k]['nb_estgenomes'] = sum([vv['new_completness']/100 for vv in v['genomes']])
    bmft[k]['rep_genome'] = v['rep_genome']
    bmft[k]['tax_level'] =  tax_level
    bmft[k]['traits'] =  traits
    bmft[k]['mOTU_type'] = "gtdbtk"
    if "s__anoxic" in k :
        bmft[k]['mOTU_type'] = "stratfreshdb"
    elif any(["s__anoxic" in vv['name'] for vv in v['genomes']]):
        bmft[k]['mOTU_type'] = "hybrid"
    bmft[k].update({ level[l-1] : ";".join(v['taxonomy'].split(";")[:l])  for l in range(1,8)})

print("Dump all the stuff")

pandas.DataFrame.from_dict(bmft, orient = 'index').to_csv(output_csv, index_label="clade_name")
