include : "workflow/rules/clades.smk"
include : "workflow/rules/dbs.smk"
include : "workflow/rules/analyses.smk"

# cd-hit and silix in a bin folder

root = "/home/moritz/data/0039_mOTUlizer_play/root/"
clades = list(zip(*glob_wildcards("/home/moritz/data/0039_mOTUlizer_play/root/d__{dom}/p__{phylum}/c__{clas}/o__{order}/f__{family}/g__{genus}/s__{species}/{genome}/")))

prefixes = ("d__","p__","c__","o__","f__","g__", "s__")

species = {tuple(c[0:7]) for c in clades if "/" not in c[0:7][-1]}

def count_genomes(taxo):
    path = root + "/".join([p + t for p,t in zip(prefixes,taxo)])
    gid_file = path +  "/" + path.split("/")[-1] + ".gids"
    with open(gid_file) as handle:
        return len(handle.readlines())
good_species = [s for s in species if count_genomes(s) > 5 ]


species = list({root + "/".join([p + cc for p,cc in zip(prefixes,c)][0:7])  for c in good_species})
genuses = list({root + "/".join([p + cc for p,cc in zip(prefixes,c)][0:6])  for c in good_species})

rule all :
    input : [g + "/" + os.path.basename(g) +  ".ab_initio.motupan.json" for g in species]
