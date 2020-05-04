include : "workflow/rules/clades.smk"

# cd-hit and silix in a bin folder

root = "/home/moritz/data/gtdb/root/"
clades = list(zip(*glob_wildcards("/home/moritz/data/gtdb/root/d__{dom}/p__{phylum}/c__{clas}/o__{order}/f__{family}/g__{genus}/s__{species}/{genome}/")))

prefixes = ("d__","p__","c__","o__","f__","g__", "s__")

genuses = {root + "/".join([p + cc for p,cc in zip(prefixes,c)][0:6])  for c in clades}

print("#We have ", len(genuses), " for ", len(clades), " genomes")

rule all :
    input : [g + "/" + os.path.basename(g) + ff for g in genuses for ff in [".silix.clusters", ".emapper"]]
