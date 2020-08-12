import sys
import json
from tqdm import tqdm
import os
from os.path import join as pjoin
import gzip

script, out_db, tax_level, trait  = sys.argv
level = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
prefix = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']

level2pref = {k :v for k,v in zip(level, prefix)}

motupan_outs = []
for v in tqdm(os.walk(out_db.split("dbs/")[0] + "root"), followlinks=True):
    for vv in v[2]:
        if vv.endswith(trait + ".motupan.json") and vv.startswith(level2pref[tax_level]):
            motupan_outs += [pjoin(v[0], vv)]

def find_rep(mag_dat, method = "checkm", tolerance = 5, max_contamin = 5):
        if method == "checkm":
            magdat = {k : v for k,v in mag_dat.items() if v['contamination'] < max_contamin}
            most_complete = max(magdat.values(), key = lambda a : a['completeness'])['completeness']
            cands = {k : v for k,v in magdat.items() if v['completeness'] > (most_complete - tolerance)}
            winner = min(cands.items(), key = lambda a : a[1]['contamination'])
            return winner[0]

if trait == "ab_initio":
    cogs_files = { pjoin( k.split("/" + level2pref[tax_level])[0], os.path.basename(k.split("/" + level2pref[tax_level])[0]) + ".silix.clusters.gid2cog")  for k in motupan_outs}

    gid2cogs = {}
    for f in tqdm(cogs_files):
        with open(f) as handle:
            gid2cogs.update(json.load(handle))

#cogs_files = { pjoin( k.split("/s__")[0], os.path.basename(k.split("/s__")[0]) + ".emapper")  for k in motupan_outs}


# gid2cogs = {}
# for f in tqdm(motupan_outs):
#     with open(f) as handle:
#         tt = {k : v['cogs'] for k,v in json.load(handle).items()}
#         gid2cogs.update(tt)
#

database = {}
for p in tqdm(motupan_outs):
    with open(p) as handle:
        tt = json.load(handle)
    idd =  list(tt.keys())[0]
    tt[idd]['taxonomy'] = os.path.dirname(p.split("root/")[1]).replace("/",";")
#    tt[idd]['cog_composition'] = {g['name'] : gid2cogs[g['name']] for g in tt[idd]['genomes']}
    core = set(tt[idd]['core'])
    accessory = set(tt[idd]['aux_genome'])
    singles = set(tt[idd]['singleton_cogs'])
    mg_info = {g['name'] : {"completeness" : g['checkm_complet'],"contamination" : g['checkm_contamin']} for g in tt[idd]['genomes']}

    if trait == "ab_initio":
        tt[idd]['cogs'] ={}

    tt[idd]['rep_genome'] = find_rep(mg_info)
    for g in tt[idd]['genomes']:
        gid = g['name']
        if trait == "ab_initio":
            tt[idd]['cogs'][gid] = gid2cogs[gid]
        if os.path.exists(pjoin(os.path.dirname(p), gid,gid + ".fna.gz" )):
            with gzip.open(pjoin(os.path.dirname(p), gid,gid + ".fna.gz" )) as handle:
                ls = [l.decode()[:-1] for l in handle]
            counts = [ (len(l), l.count("G") + l.count("C")) for l in ls if l[0] != ">"]
            g['genome_len'] = sum([l[0] for l in counts])
            g['GC'] = sum([l[1] for l in counts])/g['genome_len']
        g['total_cogs'] = len(tt[idd]['cogs'][gid])
        g['acc_cogs_count'] = len(accessory.intersection(tt[idd]['cogs'][gid]))
        g['single_cogs_count'] = len(singles.intersection(tt[idd]['cogs'][gid]))
        g['varfract'] =  g['acc_cogs_count']/g['total_cogs'] if g['total_cogs'] > 0 else None
        if  g['total_cogs'] == 0 :
            print(gid,"of taxonomy",tt[idd]['taxonomy'], "has no proteins...")
            tt[idd]['contains_broken_proteoms'] = 1
    database.update(tt)

with open(out_db, "w") as handle:
    json.dump(database, handle, indent=4, sort_keys=True)

# broken proteoms :
# d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Yersinia;s__Yersinia pestis
# GCA_002412245.1
# d__Bacteria;p__Actinobacteriota;c__Actinobacteria;o__Mycobacteriales;f__Mycobacteriaceae;g__Mycobacterium;s__Mycobacterium tuberculosis
# GCF_001380255.1
# GCF_001948285.1
# GCF_001372235.1
# GCF_001849425.1
# d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Rhizobiaceae;g__Rhizobium;s__Rhizobium sp900466475
# GCA_900466455.1
# d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Legionellales;f__Legionellaceae;g__Legionella;s__Legionella pneumophila
# GCF_900092465.1
# d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Alcaligenes;s__Alcaligenes faecalis
# GCF_002120045.1
# d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Burkholderia;s__Burkholderia multivorans
# GCF_002981415.1
# d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Xanthomonadales;f__Xanthomonadaceae;g__Xanthomonas;s__Xanthomonas citri
# GCF_003064125.1
# d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Vibrionaceae;g__Vibrio;s__Vibrio parahaemolyticus
# GCA_000877625.1
# d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Streptococcaceae;g__Streptococcus;s__Streptococcus mutans
# GCF_001068415.1
# d__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus_F;s__Lactobacillus_F plantarum
# GCA_002737775.1
# d__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae_G;g__Bacillus_A;s__Bacillus_A thuringiensis
# GCF_002813875.1
# d__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae_G;g__Bacillus_A;s__Bacillus_A thuringiensis_J
# GCF_002148155.1
# d__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;f__Bacillaceae_G;g__Bacillus_A;s__Bacillus_A wiedmannii
# GCF_001044815.1
# d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Pseudomonadales;f__Pseudomonadaceae;g__Pseudomonas;s__Pseudomonas aeruginosa
# GCF_900145395.1
# d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Klebsiella;s__Klebsiella pneumoniae
# GCA_001879345.1
# GCF_001597145.1
# GCF_000822505.1
# d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Salmonella;s__Salmonella enterica
# GCF_001143105.1
# GCF_002106735.1
# d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia flexneri
# GCF_002152115.1
# d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli
# GCF_002468365.1
