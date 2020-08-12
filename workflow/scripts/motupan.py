import sys
from subprocess import call
import os
from os.path import join as pjoin
import json
import shutil
import tempfile

script, gid_file, silix_clusts, motupan_out, checkm_file = sys.argv

clade_folder = os.path.dirname(gid_file)
log_file = gid_file.replace("gids", "log")
with open(gid_file) as gids_handle:
    gids = {k.strip() for k in gids_handle}

clade = os.path.basename(gid_file)[:-5]
silix_cogs = silix_clusts + ".gid2cog"
with open(silix_cogs) as handle:
    gid2cog = json.load(handle)

simple_cogs = {k : v for k,v in gid2cog.items() if k in gids}

#checking checkm file
completeness_switch = "--checkm " + checkm_file if os.path.exists(checkm_file) else ""

print("executing mOTUpan")
with tempfile.NamedTemporaryFile(mode = "w", suffix = ".gid2cog") as temp:
    json.dump(simple_cogs, temp, indent=4, sort_keys=True)
    temp.flush()
    exec = "mOTUpan.py --output '{output}' {completes} --cog_file '{cogs}' --name '{taxo}' >> '{log}' 2>> '{log}'".format(output = motupan_out, completes = completeness_switch, cogs = temp.name, log = log_file, taxo = clade)
    call(exec, shell = True)
