import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio import SearchIO
from Bio.SeqRecord import SeqRecord

from suffix_trees import STree

import numpy as np; np.random.seed(123)
import seaborn as sns; sns.set()
import pandas as pd

import sys
import os
from os import listdir
from os.path import isfile, isdir, join


def load_enum(filename):
    records = list(SeqIO.parse(filename, "fasta"))
    res = {}
    res2 = {}
    for m in records:
        m_name = "_".join(m.name.split("_")[1:])
        m_num = m.name.split("_")[0]
        res[m_name] = m_num
        res2[m_num] = m_name
    return res, res2

def load_fasta(filename, tp = "list"):
    if tp == "map":
        records = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
        for r in records:
            records[r] = records[r].upper() 
    else:
        records = list(SeqIO.parse(filename, "fasta"))
        for i in range(len(records)):
            records[i] = records[i].upper()
    return records

def save_fasta(filename, orfs):
    with open(filename + ".fasta", "w") as output_handle:
        SeqIO.write(orfs, output_handle, "fasta")

def make_record(seq, name, sid, d=""):
    return SeqRecord(seq, id=sid, name=name, description = d)

def simplify(s, monomers):
    return monomers[s]
    #return "/".join([x.split("_")[1] for x in s.split(".")[:2]])

def load_string_decomposition(filename, monomers):
    reads_mapping = {}
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            if len(ln.strip().split("\t")) > 4:
                sseqid, qseqid, sstart, send, idnt  = ln.strip().split("\t")[:5]
                if sseqid not in reads_mapping:
                        reads_mapping[sseqid] = []
                s, e, idnt = int(sstart), int(send), float(idnt)
                rev = False
                if qseqid.endswith("'"):
                    rev = True
                    qseqid = qseqid[:-1]
                qseqid = qseqid.split()[0]
                reads_mapping[sseqid].append({"qid": simplify(qseqid, monomers),"s": s, "e": e, "rev": rev, "idnt": idnt})
            else:
                print(ln)
    new_reads_mapping = {}
    for r in reads_mapping:
        reads_mapping[r] = reads_mapping[r][1:-1]
        all_rev = True
        avg_identity = 0
        for i in range(1, len(reads_mapping[r])):
            if reads_mapping[r][i]["rev"] != reads_mapping[r][i - 1]["rev"]:
                all_rev = False
            avg_identity += reads_mapping[r][i]["idnt"]
        avg_identity /= len(reads_mapping[r])
        if all_rev and avg_identity > 160:
            if reads_mapping[r][0]["rev"]:
                new_reads_mapping[r] = sorted(reads_mapping[r], key=lambda x: (-x["s"], x["e"]))
            else:
                new_reads_mapping[r] = sorted(reads_mapping[r], key=lambda x: (x["s"], -x["e"]))
    return new_reads_mapping

monomers, c_monomers = load_enum("/Sid/tdvorkina/monomers/new_monomers/monomer_clustering/cluter_represent_ed6_enumarated.fasta")
monomers_real = load_fasta("/Sid/tdvorkina/monomers/new_monomers/monomer_clustering/cluter_represent_ed6.fasta", "map")
dec = load_string_decomposition("/Sid/tdvorkina/monomers/new_monomers/cent_reads/SRR9087598/decomposition_cent_reads_raw_ed6.tsv", monomers)
dec2 = load_string_decomposition("/Sid/tdvorkina/monomers/new_monomers/cent_reads/SRR9087597/decomposition_cent_reads_raw_ed6.tsv", monomers)

for d in dec2:
    dec[d] = dec2[d]
print("Reads num", len(dec))

monomer_freq = {}
for m in monomers:
    monomer_freq[monomers[m]] = [0,[], set()] 

for r in dec:
    i = 0
    for h in dec[r]:
        c = h["qid"]
        monomer_freq[c][0] += 1
        monomer_freq[c][1].append(str(min(i, len(dec[r]) - i)))
        monomer_freq[c][2].add(r)
        i += 1

monomer_freq_lst = []
for m in monomer_freq:
    monomer_freq_lst.append([m, monomer_freq[m][0], monomer_freq[m][1], monomer_freq[m][2]])
print(len(monomers), len(monomer_freq_lst))
monomer_freq_lst = sorted(monomer_freq_lst, key = lambda x: -x[1])

for m in monomer_freq_lst[:10]:
    print(m[0], m[1])

new_monomer_enum = []
i = 0
cnt = 0
sm = 0
for m in monomer_freq_lst:
    ini_name = c_monomers[m[0]]
    i += 1
    new_monomer_enum.append(make_record(monomers_real[ini_name].seq, str(i)+"_" + monomers_real[ini_name].name, str(i)+"_" + monomers_real[ini_name].id))
    print(str(i) + "\t" + str(m[1]))
    if m[1] > 0:
        cnt += 1
    sm += m[1]

print("Non zero", cnt)
print("Total", sm)

save_fasta("/Sid/tdvorkina/monomers/new_monomers/monomer_clustering/cluter_represent_ed6_enumarated_order", new_monomer_enum)

fig = plt.figure()
plt.hist([min(x[1], 1000) for x in monomer_freq_lst if x[1] > 0], bins=30)
plt.xticks(np.arange(0, 1001, step=100), ("0", "100", "200", "300", "400", "500", "600", "700", "800", "900", ">1000"))

plt.savefig("./monomers_freq.png")