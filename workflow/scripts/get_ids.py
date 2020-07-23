from os.path import join as pjoin
import os
from tqdm import tqdm

with open("/home/moritz/projects/0026_Chlorobi/data/gtdb.txt") as handle:
    ids = [l.strip() for l in handle]

others_folder = "RS_nd_GBs"

def make_folder(gid):
    if gid.startswith("RS") or gid.startswith("GB"):
        gid = gid[3:]
    nums = gid.split('_')[-1]
    first = nums[0:3]
    second = nums[3:6]
    third = nums[6:9]
    pat = pjoin(others_folder, "{first}/{second}/{third}/{gid}".format(first=first, second = second, third = third, gid = gid))
    return pat

for f in tqdm(ids) :
    if os.path.exists(make_folder(f)):
        shutil.copytree(make_folder(f), "chls/" + f)
