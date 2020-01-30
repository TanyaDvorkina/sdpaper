import os
from os import listdir
from os.path import isfile, isdir, join

path = "/Sid/tdvorkina/monomers/new_monomers/cent_reads/"

datasets = ["SRR9087600", "SRR9087599", "SRR9087598", "SRR9087597"]

for d in datasets:
    print(d)
    reads_num = set()
    with open(os.path.join(path, d, "ac", "inferred_monomers.fa"), "r") as fin:
        for ln in fin.readlines():
            if ln.startswith(">"):
                read_name = ln.strip().split("/")[0][1:] 
                reads_num.add(read_name)
    print(len(reads_num))
    with open(os.path.join(path, d, "ac", "true_cen.lst"), "w") as fout:
        fout.write("\n".join(reads_num))