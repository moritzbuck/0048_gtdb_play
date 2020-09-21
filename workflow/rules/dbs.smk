root = "/home/moritz/data/0039_mOTUlizer_play/root/"
db_folder = "/home/moritz/data/0039_mOTUlizer_play/dbs/"

clades = list(zip(*glob_wildcards("/home/moritz/data/0039_mOTUlizer_play/root/d__{dom}/p__{phylum}/c__{clas}/o__{order}/f__{family}/g__{genus}/s__{species}/{genome}/")))

prefixes = ("d__","p__","c__","o__","f__","g__", "s__")

species = {tuple(c[0:7]) for c in clades if "/" not in c[0:7][-1]}

def count_genomes(taxo):
    path = root + "/".join([p + t for p,t in zip(prefixes,taxo)])
    gid_file = path +  "/" + path.split("/")[-1] + ".gids"
    with open(gid_file) as handle:
        return len(handle.readlines())
good_species = [s for s in species if count_genomes(s) > 5 ]


species = list({"{base}/root/" + "/".join([p + cc for p,cc in zip(prefixes,c)][0:7])  for c in good_species})
genuses = list({"{base}/root/" + "/".join([p + cc for p,cc in zip(prefixes,c)][0:6])  for c in good_species})

rule make_species_db:
    input : motupan_outs = [g + "/" + os.path.basename(g) +  ".{traits}.motupan.json" for g in species]
    output : out_db = "{base}/dbs/{tax_level}/{traits}/database.json"
    params : script = "workflow/scripts/species_abinitio_db.py"
    conda : "../envs/python_std.yaml"
    threads : 1
    shell : """
        python3 {params.script} {output.out_db} {wildcards.tax_level} {wildcards.traits}
        """

rule make_annot_db:
    input : ab_initio_db = "{base}/dbs/{tax_level}/ab_initio/database.json",  emapper_db = "{base}/dbs/{tax_level}/emapper/database.json"
    output : out_dbs = ["{base}/dbs/{tax_level}/annotation_databases/" + type + "-" + trait + ".json" for type in ['EC', 'KO', 'cogcat', 'CAZy'] for trait in ['ab_initio', 'emapper']]
    params : script = "workflow/scripts/make_annot_db.py"
    conda : "../envs/python_std.yaml"
    threads : 1
    shell : """
        python3 {params.script} {input.ab_initio_db} {input.emapper_db} {output.out_dbs[0]} {wildcards.tax_level}
        """
