#!/usr/bin/env python3
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio import SearchIO
from Bio.SeqRecord import SeqRecord

import os
from os import listdir
from os.path import join, isfile
import sys
import argparse

import edlib
import random

random.seed(123)

import re

import numpy as np; np.random.seed(0)
import seaborn as sns; sns.set()
import pandas as pd

ref_hits_file = "/Sid/tdvorkina/monomers/sdpaper/SD/decomposition_cenX_full.tsv"
monomers_file = "../data/DXZ1_inferred_monomers_single.fa"

def load_fasta(filename):
    records = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
    new_records = {}
    for r in records:
        new_records[r] = records[r].upper()
    return new_records

def load_string_decomposition(filename, identity = 0):
    reads_mapping = {}
    filtered = 0
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            sseqid, qseqid, sstart, send, idnt = ln.strip().split("\t")[:5]
            qseqid = qseqid.split()[0]
            if sseqid not in reads_mapping:
                    reads_mapping[sseqid] = []
            s, e, idnt = int(sstart), int(send), float(idnt)
            rev = False
            if qseqid.endswith("'"):
                rev = True
                qseqid = qseqid[:-1]
            if idnt > identity:
                reads_mapping[sseqid].append({"qid": convert_to_abc(qseqid), "s": s, "e": e, "rev": rev, "idnt": idnt})
            else:
                filtered += 1
    for r in reads_mapping:
        reads_mapping[r] = sorted(reads_mapping[r], key=lambda x: (x["s"], -x["e"]))
    print("Unreliable ", filtered)
    return reads_mapping

def convert_to_abc(name):
    res = name.split("_")[0]
    if "fake" in name:
        res = 'M'
    return res

def edist(lst):
    if len(str(lst[0])) == 0:
        return -1, ""
    if len(str(lst[1])) == 0:
        return -1, ""
    result = edlib.align(str(lst[0]), str(lst[1]), mode="NW", task="path")
    return result["editDistance"], result["cigar"]

def aai(ar):
    p1, p2 = str(ar[0]), str(ar[1])
    if p1.endswith("*"):
        p1 = p1[:-1]
    if p2.endswith("*"):
        p2 = p2[:-1]
    ed, cigar = edist([str(p1), str(p2)])
    if ed == -1:
        return 0
    matches = re.findall(r'\d+=', cigar)
    aai = 0.0
    for m in matches:
        aai += int(m[:-1])
    aai /= max(len(p1), len(p2))
    return aai*100

def cnt_identity(hits, filename):
    idnts_x, idnts = [],  []
    start, end = 0, 20000
    s_ind, e_ind = 0, 0
    for r in hits:
        for a in hits[r]:
            idnts.append(a["idnt"])
            idnts_x.append(a["s"])

    fig = plt.figure(figsize =(50,10))
    plt.plot(idnts_x, idnts)
    plt.ylabel("Monomer Identity", fontsize=24)
    plt.yticks(fontsize=24)
    plt.xticks(fontsize=24)
    plt.ylim(20, 100)
    plt.savefig(filename)

def random_stats(hits):
    df = pd.DataFrame(hits)
    fig = plt.figure()
    ax = sns.violinplot(x="qid", y="idnt", data=df, order=["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"])
    plt.savefig("random_mono_d.png")

def draw_pvalue(idnts, filename):
    idnts = sorted(idnts)
    x, y = [], []
    j = 0
    for i in range(60, 71):
        while j < len(idnts) and idnts[j] < i:
            j += 1
        if j < len(idnts):
            print("\t".join([ str(x) for x in [i, j + 1, (len(idnts) - j + 1)/len(idnts), idnts[j]] ]))
        else:
            print("\t".join([ str(x) for x in [i, j, (len(idnts) - j + 1)/len(idnts), idnts[-1]] ]))
        x.append(i)
        y.append( (len(idnts) - j + 1)/len(idnts))
    fig = plt.figure()#figsize =(50,10))
    plt.plot(x, y)
    plt.ylabel("P-value")
    plt.savefig(filename)

ref_hits = load_string_decomposition(ref_hits_file)
cnt_identity(ref_hits, "./test.png")


monomers = load_fasta(monomers_file)

num = 100000
hits = []
idnts = []
for j in range(num):
    seq = ""
    for _ in range(171):
        seq += random.choice(["A","C","G","T"])
    max_indt = 0
    for m in monomers:
        cur_idnt = aai([seq, monomers[m].seq])
        #hits.append({"qid": convert_to_abc(m), "idnt": cur_idnt})
        if cur_idnt > max_indt:
            max_indt = cur_idnt
    idnts.append(max_indt)
print(sum(idnts)/len(idnts))
print(sorted(idnts)[len(idnts)//2])
draw_pvalue(idnts, "./pvalue_monomer_identity.png")



# random_stats(hits)
