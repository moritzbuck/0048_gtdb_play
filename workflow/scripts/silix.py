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

clade = os.path.basename(gid_file)[:-3]
hits_file = silix_clusts.replace("silix.clusters", "selfhits")

exec = """
unpigz {faa}
diamond makedb --db {faa} --in {faa}
diamond blastp --more-sensitive  -e0.001  -p {threads} -f 6 -q {faa} --db {faa} -o {hits}
silix {faa_minusgz} {hits} >  {clust_file}
pigz {faa}
""".format(faa = full_proteom, faa_minusgz = full_proteom[:-3], clust_file = silix_clusts, threads = threads, hits = hits_file)

call(exec, shell = True)

with open(silix_clusts) as handle:
    silix_out = { l.strip().split("\t")[1] : l.split("\t")[0] for l in handle.readlines()}

with open(clusters_file) as handle:
    preclusters = {ll : l.split()[0] for l in handle for ll in l.strip().split()[1:]}

proteoms = []
for v in os.walk(clade_folder):
    for vv in v[1]:
        if vv in gids:
            proteoms += [pjoin(v[0], vv, vv + ".faa.gz")]

gid_annot = {}
gid2cogs = {}
for f in proteoms:
    fold = os.path.dirname(f)
    gid = os.path.basename(fold)
    with gzip.open(f) as handle:
        lines = [l.decode() for l in handle.readlines()]
    gid_annot[os.path.basename(f)[:-7]] = {l.split()[0][1:] : 'g__TMED189' + "_silix_COG_" + silix_out[l.split()[0][1:]] for l in lines if l.startswith(">")}
    gid2cogs[os.path.basename(f)[:-7]] = list(set(gid_annot[os.path.basename(f)[:-7]].values()))

with open(silix_clusts, "w") as handle:
    json.dump(gid_annot, handle, indent=4, sort_keys=True)

with open(silix_clusts + ".gid2cog", "w") as handle:
    json.dump(gid2cogs, handle, indent=4, sort_keys=True)
