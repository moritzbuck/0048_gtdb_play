import sys
from subprocess import call
import os
from os.path import join as pjoin
import json
import shutil
import tempfile
from ete3 import Tree
from sys import stdout

script, db_json, arc_tree_file, bac_tree_file, distance_json, tax_level, traits = sys.argv
level = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']
prefix = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
radius = 7

radius_level = level[[ (i-radius) if (i-radius) > 0 else 0 for i, t in enumerate(level) if t == tax_level][0]]
radius_prefix = prefix[[ (i-radius) if (i-radius) > 0 else 0 for i, t in enumerate(level) if t == tax_level][0]]

with open(db_json) as handle:
    db = json.load(handle)

bac_dir = pjoin(os.path.dirname(bac_tree_file), "rep_folder/d__Bacteria")
arc_dir = pjoin(os.path.dirname(arc_tree_file), "rep_folder/d__Archaea")

os.makedirs(bac_dir,  exist_ok=True)
os.makedirs(arc_dir,  exist_ok=True)

print("copying representative genomes")
rep2clade = {}
clade2radius = {}
for k, v in db.items():
    rep2clade[v['rep_genome']] = k
    genome_path = pjoin(db_json.split("dbs")[0], "root", v['taxonomy'].replace(";","/"))
    genome_file = pjoin(genome_path,  v['rep_genome'] ,  v['rep_genome'] + ".fna.gz")
    clade2radius[k] = [ ";".join(v['taxonomy'].split(";")[:t]) for t in range(1,7) if v['taxonomy'].split(";")[:t][-1].startswith(radius_prefix)][0]
    if "d__Bacteria" in genome_path:
        shutil.copy(genome_file, bac_dir)
    else :
        shutil.copy(genome_file, arc_dir)

radius2clade = {s : [] for s in set(clade2radius.values())}
for k,v in clade2radius.items():
    radius2clade[v] += [k]



print("Running GTDBtk for making trees")
gtdb_line = """gtdbtk de_novo_wf --cpus 24 --extension gz --genome_dir {genome_dir} {dom_switch} --out_dir ../../data/0039_mOTUlizer_play/analyses/species/{traits}/tree/gtdbtk_{dom} --skip_gtdb_refs --outgroup_taxon None"""

call(gtdb_line.format(dom_switch = "--bacteria", dom="bac", genome_dir=bac_dir, traits = traits), shell = True)
call(gtdb_line.format(dom_switch = "--archaea", dom="arc", genome_dir=arc_dir, traits = traits), shell = True)
#shutil.copy("../../data/0039_mOTUlizer_play/analyses/species/ab_initio/tree/gtdbtk_ab_initio_arc_rep.nwk", "../../data/0039_mOTUlizer_play/analyses/species/ab_initio/tree/gtdbtk_denovo_arc_rep.nwk")
#shutil.copy("../../data/0039_mOTUlizer_play/analyses/species/ab_initio/tree/gtdbtk_ab_initio_bac_rep.nwk", "../../data/0039_mOTUlizer_play/analyses/species/ab_initio/tree/gtdbtk_denovo_bac_rep.nwk")

if os.path.exists(pjoin(os.path.dirname(arc_tree_file), "gtdbtk_arc/infer/intermediate_results/gtdbtk.ar122.unrooted.tree")):
    arc_tree = Tree(pjoin(os.path.dirname(arc_tree_file), "gtdbtk_arc/infer/intermediate_results/gtdbtk.ar122.unrooted.tree"))
    for n in arc_tree.iter_leaves():
        n.name = rep2clade[n.name[:-4]]
else :
    arc_tree = Tree()
if os.path.exists(pjoin(os.path.dirname(arc_tree_file), "gtdbtk_bac/infer/intermediate_results/gtdbtk.bac120.unrooted.tree")):
    bac_tree = Tree(pjoin(os.path.dirname(bac_tree_file), "gtdbtk_bac/infer/intermediate_results/gtdbtk.bac120.unrooted.tree"))
    for n in bac_tree.iter_leaves():
        n.name = rep2clade[n.name[:-4]]
else :
    bac_tree = Tree()

#distances = {(m.name, n.name) : arc_tree.get_distance(m,n) for m in  arc_tree.iter_leaves() for n in  arc_tree.iter_leaves()}
#distances.update({(m.name[:-4], n.name[:-4]) : bac_tree.get_distance(m,n) for m in  tqdm(bac_tree.iter_leaves()) for n in  bac_tree.iter_leaves()})
print("Computing phylo distances!")

distances = {}
for k , v in radius2clade.items():
    print("Distances for ", k )
    tree = bac_tree if k.startswith("d__Bacteria") else arc_tree
    leaves = [l for l in tree.iter_leaves() if l.name in v]
    distances.update({(m.name, n.name) : tree.get_distance(m,n) for m in  leaves for n in  leaves})

print("Dumping stuff!")
arc_tree.write(outfile = arc_tree_file)
bac_tree.write(outfile = bac_tree_file)

with open(distance_json, "w") as handle:
    json.dump([list(k) + [v] for k,v in distances.items()], handle, indent=4, sort_keys=True)
