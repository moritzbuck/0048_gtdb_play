import sys
from subprocess import call
import os
from os.path import join as pjoin
import json
import shutil
import tempfile
from tqdm import tqdm
from ete3 import NCBITaxa

script, gid_file, emapper_results, motupan_out, checkm_file = sys.argv

clade_folder = os.path.dirname(gid_file)
log_file = gid_file.replace("gids", "log")

tax_db = NCBITaxa()

with open(gid_file) as gids_handle:
    gids = {k.strip() for k in gids_handle}

clade = os.path.basename(gid_file)[:-5]

gid2eggs = {}
for g in tqdm(gids):
    patty = pjoin(clade_folder, g, g + ".emapper")
    with open(patty) as handle:
        gid2eggs[g] = {v["eggNOG OGs"] for k,v in json.load(handle).items() if v}

taxos = {vvv.split("@")[1] for v in gid2eggs.values() for vv in v for vvv in vv.split(",")}
tax2level = {k : len(v) for k,v in tax_db.get_lineage_translator(list(taxos)).items()}
for g in tqdm(gids):
    gid2eggs[g] = list({ min([vv.split('@') for vv in v.split(",")], key = lambda x : tax2level.get(int(x[1]),1000))[0] for v in gid2eggs[g]})


#checking checkm file
completeness_switch = "--checkm " + checkm_file if os.path.exists(checkm_file) else ""

print("executing mOTUpan")
with tempfile.NamedTemporaryFile(mode = "w", suffix = ".gid2cog") as temp:
    json.dump(gid2eggs, temp, indent=4, sort_keys=True)
    temp.flush()
    exec = "mOTUpan.py --output '{output}' {completes} --cog_file '{cogs}' --name '{taxo}' >> '{log}' 2>> '{log}'".format(output = motupan_out, completes = completeness_switch, cogs = temp.name, log = log_file, taxo = "egg_" + clade)
    call(exec, shell = True)
