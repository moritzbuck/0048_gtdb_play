import sys
from subprocess import call
import os
from os.path import join as pjoin
import gzip
import json
import shutil

script, gid_file, full_proteom, clusters_file, silix_clusts, threads = sys.argv

clade_folder = os.path.dirname(gid_file)
log_file = gid_file.replace("gids", "log")
with open(gid_file) as gids_handle:
    gids = {k.strip() for k in gids_handle}

clade = os.path.basename(gid_file)[:-5]
hits_file = silix_clusts.replace("silix.clusters", "selfhits")



exec = """
echo running selfdiamond
diamond makedb --db {faa} --in {faa} >> {log} 2>&1
diamond blastp --more-sensitive  -e0.001  -p {threads} -f 6 -q {faa} --db {faa} -o {hits} >> {log} 2>&1
echo running silix
silix {faa} {hits} >  {clust_file} 2> {log}
""".format(faa = full_proteom, clust_file = silix_clusts, threads = threads, hits = hits_file, log = log_file)

print("executing silix COG-clustering mini-pipe")

call(exec, shell = True)

print("parsing silix clusters")
with open(silix_clusts) as handle:
    silix_out = { l.strip().split("\t")[1] : l.split("\t")[0] for l in handle.readlines()}

print("parsing preclusters")
with open(clusters_file) as handle:
    derep2clusters = {l.split()[0] : l[:-1].split()[1:] for l in handle}
preclusters = {vv : k for k,v in derep2clusters.items() for vv in v}


proteoms = []
for v in os.walk(clade_folder):
    for vv in v[1]:
        if vv in gids:
            proteoms += [pjoin(v[0], vv, vv + ".faa.gz")]

print("parsing", len(proteoms), "proteims and making COG output")
gid_annot = {}
gid2cogs = {}
for f in proteoms:
    fold = os.path.dirname(f)
    gid = os.path.basename(fold)
    with gzip.open(f) as handle:
        lines = [l.decode() for l in handle.readlines()]
    gid_annot[os.path.basename(f)[:-7]] = {l.split()[0][1:] : clade + "_silix_COG_" + silix_out[preclusters[l.split()[0][1:]]] for l in lines if l.startswith(">")}
    gid2cogs[os.path.basename(f)[:-7]] = list(set(gid_annot[os.path.basename(f)[:-7]].values()))

print("dump stuff")

with open(silix_clusts, "w") as handle:
    json.dump(gid_annot, handle, indent=4, sort_keys=True)

with open(silix_clusts + ".gid2cog", "w") as handle:
    json.dump(gid2cogs, handle, indent=4, sort_keys=True)
