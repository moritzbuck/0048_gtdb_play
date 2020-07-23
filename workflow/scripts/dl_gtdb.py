import pandas
from subprocess import call
from os.path import join as pjoin
import os
import shutil
import pandas
import hashlib
import numpy
from tqdm import tqdm
from pathos.multiprocessing import ProcessingPool as Pool
import json

def make_folder(gid):
    nums = gid.split('_')[-1]
    first = nums[0:3]
    second = nums[3:6]
    third = nums[6:9]
    pat = pjoin(others_folder, "{first}/{second}/{third}/{gid}".format(first=first, second = second, third = third, gid = gid))
    if not os.path.exists(pat):
        os.makedirs(pat, exist_ok=True)
    return pat

def process_gid(gid):
    if gid.startswith("UBA"):
        dl_folder = pjoin(uba_folder, gid)
    else :
        dl_folder = make_folder(gid)

    if metadata.loc[rev_rename.get(gid, gid), "ncbi_translation_table"] != 'none':
        trans_table = "-g " + metadata.loc[rev_rename.get(gid, gid), "ncbi_translation_table"]
    else :
        trans_table = ""

    prodigal_exe = prodigal_string.format(fold = dl_folder, gid = gid, trans_table = trans_table)

    call(prodigal_exe, shell = True)


base = "/home/moritz/data/gtdb"
metadata = pandas.concat([pandas.read_csv(pjoin(base,"bac120_metadata_r95.tsv"), sep="\t", index_col=0, low_memory=False), pandas.read_csv(pjoin(base,"ar122_metadata_r95.tsv"), sep="\t", index_col=0, low_memory=False)])
uba_folder = pjoin(base, "UBAs")
others_folder = pjoin(base, "RS_nd_GBs")

prodigal_string = """unpigz -f -p 1 {fold}/{gid}.fna.gz
prodigal -q -a {fold}/{gid}.faa -d {fold}/{gid}.ffn {trans_table} -f gff -i {fold}/{gid}.fna -o {fold}/{gid}.gff
pigz -f -p 1  {fold}/{gid}.faa {fold}/{gid}.fna {fold}/{gid}.ffn {fold}/{gid}.gff
"""

metadata.index = [f if f.startswith("UBA") else f[3:] for f in metadata.index ]

with open("renamed.json") as handle:
    rename = json.load(handle)
rev_rename = {v :k for k,v in rename.items()}

for t in ["refseq", "genbank"]:
    for k in ['archaea', 'bacteria']:
        for g in tqdm(os.listdir(pjoin(t,k))):
            call("mv {gpat} {new_pat}".format(gpat = pjoin(t,k,g), new_pat = make_folder(g)), shell = True)

there = []
for f in tqdm(metadata.index):
    if not f.startswith("UBA"):
#        f = rename.get(f, f)
        raw_dled = pjoin(make_folder(f), f)
        if os.path.exists(raw_dled):
            ass = [f for f in os.listdir(raw_dled) if f.endswith("fna.gz")]
            if len(ass) > 0:
                shutil.move(pjoin(raw_dled, ass[0]), pjoin(make_folder(f), f + ".fna.gz"))
                shutil.rmtree(raw_dled)
        if os.path.exists(pjoin(make_folder(f), f + ".fna.gz")):
            there += [f]
    else :
        there += [f]

to_genbank = lambda x : x.replace("GCF","GCA")
plus_version = lambda x : x.split(".")[0] + "." + str(int(x.split(".")[1]) +1)
to_refseq = lambda x : x.replace("GCA","GCF")
test_names = lambda x : (x, plus_version(x), to_refseq(x), to_refseq(plus_version(x)), to_genbank(x), to_genbank(plus_version(x)))
test_paths = lambda gid: {x : os.path.exists(make_folder(x) + "/" + x + ".fna.gz") for x in test_names(gid)}

to_dl = [k for k in metadata.index if not k.startswith("UBA") and not os.path.exists(pjoin(make_folder(k), k + ".faa.gz"))]
tt = {k : test_paths(k) for k in to_dl}
rename = {k : [kk for kk,vv in v.items() if vv][0] for k,v in tt.items() if sum(v.values()) == 1}
rename['GCA_900033825.1'] = 'GCA_900033825.3'
rename['GCA_003312795.1'] = 'GCA_003312795.4'
rename['GCA_003312635.3'] = 'GCA_003312635.3'
rename['GCA_003298655.1'] = 'GCA_003298655.3'
rename['GCA_003312775.3'] = 'GCA_003312775.3'
rename['GCA_002850675.2'] = 'GCA_002850675.5'
rename['GCA_003312605.1'] = 'GCA_003312605.3'
rename['GCA_002872495.1'] = 'GCA_002872495.3'

with open(pjoin(base, "genbamk_ids.txt"), "w") as handle :
    handle.writelines([k + "\n" for k in metadata.index if k.startswith("GCA_")] )
with open(pjoin(base, "refsecks_ids.txt"), "w") as handle :
    handle.writelines([k + "\n" for k in metadata.index  if k.startswith("GCF_")] )

call("ncbi-genome-download -p 24  --format fasta -A genbamk_ids.txt -s genbank all")
call("ncbi-genome-download -p 24  --format fasta -A refsecks_ids.txt -s refseq all")

there = [rename.get(k,k) for k in metadata.index if k.startswith("UBA") or os.path.exists(pjoin(make_folder(rename.get(k,k)), rename.get(k,k) + ".fna.gz"))]

there = [k for k in there if not os.path.exists(pjoin(make_folder(k), k + ".faa.gz")) and not os.path.exists(pjoin(uba_folder, k, k + ".faa.gz"))]

pool = Pool(24)
pool.map(process_gid, there)


pats = []
for gid, taxo in tqdm(metadata.gtdb_taxonomy.iteritems()):
    if gid.startswith("UBA"):
        gfold = pjoin(uba_folder, gid)
    else :
        gfold = make_folder(rename.get(gid,gid))
    if os.path.exists(pjoin(gfold, rename.get(gid,gid) + ".faa.gz")):
        pats += [gfold]
        taxo_path = "root/" + taxo.replace(";", "/")
        os.makedirs(taxo_path , exist_ok= True)
        if not os.path.exists(pjoin(taxo_path, rename.get(gid, gid))):
            os.symlink(gfold, pjoin(taxo_path, rename.get(gid, gid)))

def make_gid_files(clade = "root", folder = "root"):
    subclades = [f for f in os.listdir(folder) if os.path.isdir(pjoin(folder, f)) and f[1:3] == "__"]
    if len(subclades) > 0:
        gids = set()
        for f in subclades:
            gids.update(make_gid_files(f, pjoin(folder, f)))
    else :
        gids = {rename.get(f,f) for f in os.listdir(folder) if f in metadata.index or f in rename.values()}

    if os.path.exists(pjoin(folder, clade + ".gids")) :
        with open(pjoin(folder, clade + ".gids")) as handle:
            old_gids = {l.rstrip() for l in handle}
        if old_gids == gids:
            return gids

    with open(pjoin(folder, clade + ".gids"), "w") as handle:
        handle.writelines([gid + "\n" for gid in gids])
    return gids



make_gid_files()

call("find -L root/ -name "*.fna.gz" -exec unpigz -c -p20 {} \; > all_proteoms.faa", shell = True)
