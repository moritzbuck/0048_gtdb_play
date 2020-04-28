import sys
from subprocess import call
import os
from os.path import join as pjoin
import gzip
import json
import shutil

script, gid_file, full_proteom, clusters_file, emap_out, threads = sys.argv

clade_folder = os.path.dirname(gid_file)
with open(gid_file) as gids_handle:
    gids = {k.strip() for k in gids_handle}

exec = "emapper.py  -i {input}  -o {output} --cpu {threads}  -m diamond".format(input = full_proteom, output = emap_out, threads = threads)

call(exec, shell = True)

shutil.move(emap_out + ".emapper.annotations", emap_out)
os.remove(emap_out + ".emapper.seed_orthologs")

with open(emap_out) as handle:
    l = handle.readline()
    l = handle.readline()
    l = handle.readline()
    header = handle.readline()[1:-1].split("\t")
    annotations = {l.split("\t")[0] : {ll[0] : ll[1] for ll in zip(header[1:], l.split("\t")[1:]) } for l in tqdm(handle.readlines())}

with open(clusters_file) as handle:
    preclusters = {ll : l.split()[0] for l in tqdm(handle) for ll in l.strip().split()[1:]}

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
        json.dump(gid_annot, handle, indent=4, sort_keys=True)
