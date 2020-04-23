import sys
from subprocess import call
import os
from os.path import join as pjoin
import gzip
import json
import shutil

script, gid_file, full_proteom, clusters_file, silix_clusts, threads = sys.argv

clade_folder = os.path.dirname(gid_file)
with open(gid_file) as gids_handle:
    gids = {k.strip() for k in gids_handle}

hits_file = silix_clusts.replace("silix.clusters", "selfhits")
    
exec = """
diamond makedb --db {faa} --in {faa}
diamond blastp --more-sensitive  -e0.001  -p {threads} -f 6 -q {faa} --db {faa} -o {hits}
silix <(unpigz -c {faa}) {hits} >  {clust_file}
""".format(faa = full_proteom, clust_file = silix_clusts, threads = threads, hits = hits_file)

call(exec, shell = True)

with open(silix_clusts) as handle:
    silix_clusts = { l.strip().split("\t")[1] : l.split("\t")[0] for l in handle.readlines()}

with open(clusters_file) as handle:
    preclusters = {ll : l.split()[0] for l in handle for ll in l.strip().split()[1:]}

proteoms = []
for v in os.walk(clade_folder):
    for vv in v[1]:
        if vv in gids:
            proteoms += [pjoin(v[0], vv, vv + ".faa.gz")]

for f in proteoms:
    fold = os.path.dirname(f)
    gid = os.path.basename(fold)
    with gzip.open(f) as handle:
        lines = [l.decode() for l in handle.readlines()]
    gid_annot = {l.split()[0][1:] : annotations.get(l.split()[0][1:], None) for l in lines if l.startswith(">")}
    with open(pjoin(fold, gid + ".emapper"), "w") as handle:
        json.dump(gid_annot, handle)