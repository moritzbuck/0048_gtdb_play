import sys
from subprocess import call
import os
from os.path import join as pjoin
import gzip
import json
import shutil

script, gid_file, full_proteom, clusters_file, emap_out, threads = sys.argv

clade_folder = os.path.dirname(gid_file)
log_file = gid_file.replace("gids", "log")

print("annotationg", gid_file)

with open(gid_file) as gids_handle:
    gids = {k.strip() for k in gids_handle}

print("running emapper on preclusters")
    
exec = "emapper.py  -i {input}  -o {output} --cpu {threads}  -m diamond  >> {log} 2>&1".format(input = full_proteom, output = emap_out, threads = threads, log = log_file)

call(exec, shell = True)

shutil.move(emap_out + ".emapper.annotations", emap_out)
os.remove(emap_out + ".emapper.seed_orthologs")

print("parsing emapper results")

with open(emap_out) as handle:
    l = handle.readline()
    l = handle.readline()
    l = handle.readline()
    header = handle.readline()[1:-1].split("\t")
    annotations = {l.split("\t")[0] : {ll[0] : ll[1] for ll in zip(header[1:], l.split("\t")[1:]) } for l in handle.readlines()}

print("parsing preclusters")
    
with open(clusters_file) as handle:
    derep2clusters = {l.split()[0] : l[:-1].split()[1:] for l in handle}
preclusters = {vv : k for k,v in derep2clusters.items() for vv in v} 


proteoms = []
for v in os.walk(clade_folder):
    for vv in v[1]:
        if vv in gids:
            proteoms += [pjoin(v[0], vv, vv + ".faa.gz")]

print("parsing", len(proteoms), "proteoms and attributing annotations")        
            
for f in proteoms:
    fold = os.path.dirname(f)
    gid = os.path.basename(fold)
    with gzip.open(f) as handle:
        lines = [l.decode() for l in handle.readlines()]
    gid_annot = {l.split()[0][1:] : annotations.get(preclusters[l.split()[0][1:]], None) for l in lines if l.startswith(">")}
    with open(pjoin(fold, gid + ".emapper"), "w") as handle:
        json.dump(gid_annot, handle, indent=4, sort_keys=True)
