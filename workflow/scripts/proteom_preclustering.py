import sys
from subprocess import call
from os.path import join as pjoin
import gzip

script, gid_file, full_proteom, representatives, clusters_file, threads = sys.argv

clade_folder = os.path.dirname(gid_file)

with open(gid_file) as gids_handle:
    gids = {k.strip() for k in gids_handle}

proteoms = []
for v in os.walk(clade_folder):
    for vv in v[1]:
        if vv in gids:
            proteoms += [pjoin(v[0], vv, vv + ".faa.gz")]

with gzip.open(full_proteom, "w") as full_p_handle:
    for p in proteoms:
        with gzip.open(p) as prot_handle:
                prots = prot_handle.readlines()
        full_p_handle.writelines(prots)

exec = "cd-hit -i {input} -o {output} -c 0.99 -M 0 -T {threads} -d 0 -s 0.99".format(input = full_proteom, output = representatives, threads = threads)

call(exec, shell = True)

with open(representatives + ".clstr") as handle:
    clusters = "\n".join(handle.readlines()).split("Cluster ")

clusters = [c.split("\n\n")[1:-1] for c in clusters[1:-1]]
clusters = [[cc.split(">")[1].split("... ") for cc in c] for c in clusters]
clusters = {[cc[0] for cc in c if cc[1] == "*"][0] : [cc[0] for cc in c] for c in clusters}

with open(clusters_file, "w") as handle:
    handle.writelines(["\t".join([k] + v) + "\n" for k, v in clusters.items()])