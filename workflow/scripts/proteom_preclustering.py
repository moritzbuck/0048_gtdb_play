import sys
from subprocess import call
import os
from os.path import join as pjoin
import gzip

script, gid_file, full_proteom, representatives, clusters_file, threads = sys.argv

clade_folder = os.path.dirname(gid_file)
log_file = gid_file.replace("gids", "log")
with open(gid_file) as gids_handle:
    gids = {k.strip() for k in gids_handle}

print("preculstering", gid_file)

proteoms = []
for v in os.walk(clade_folder):
    for vv in v[1]:
        if vv in gids:
            proteoms += [pjoin(v[0], vv, vv + ".faa.gz")]

print("concatenating", len(proteoms), "proteoms")

with gzip.open(full_proteom, "w") as full_p_handle:
    for p in proteoms:
        with gzip.open(p) as prot_handle:
                prots = prot_handle.readlines()
        full_p_handle.writelines(prots)

print("executing cd-hit")

exec = "cd-hit -i {input} -o {output} -c 0.95 -M 0 -T {threads} -d 0 -s 0.95 >> {log} 2>&1".format(input = full_proteom, output = representatives, threads = threads, log = log_file)

call(exec, shell = True)

print("parsing cd-hit")

with open(representatives + ".clstr") as handle:
    clusters = "\n".join(handle.readlines()).split("Cluster ")
    #os.remove(representatives + ".clstr")

clusters = [c.split("\n\n") for c in clusters[1:] if "*" in c]
clusters = [[cc.split(">")[1].split("... ") for cc in c if ">" in cc and cc != ">"] for c in clusters ]
clusters = {[cc[0] for cc in c if cc[1] == "*" or cc[1] == "*\n"][0] : [cc[0] for cc in c] for c in clusters}

print("dumping", len(clusters), "preclusters")

with open(clusters_file, "w") as handle:
    handle.writelines(["\t".join([k] + v) + "\n" for k, v in clusters.items()])
