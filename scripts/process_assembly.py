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

import numpy as np; np.random.seed(0)
import seaborn as sns; sns.set()
import pandas as pd

ref_file = "/Bmo/miheenko/git/tandemQUAST/t2t7/cenx-t2t7.masked.fa"
ref_hits_file = "/Sid/tdvorkina/monomers/sdpaper/SD/decomposition_cenX_t2t.tsv"

def load_fasta(filename):
    records = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
    new_records = {}
    for r in records:
        new_records[r] = records[r].upper()
    return new_records

def load_string_decomposition(filename, identity = 70):
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

def cnt_identity(hits):
    idnts = []
    for r in hits:
        for a in hits[r]:
            idnts.append(a["idnt"])
    print("Monomers num ", len(idnts))
    print("Median identity ", sorted(idnts)[len(idnts)//2])
    print(sorted(idnts)[:20])

def cnt_gaps(hits):
    gaps = 0
    for r in hits:
        for i in range(len(hits[r])):
            if i > 0:
                if hits[r][i]["s"] - hits[r][i-1]["e"] > 1:
                    print(hits[r][i-1]["qid"], hits[r][i]["qid"], hits[r][i]["s"] - hits[r][i-1]["e"])
                    gaps += 1
    print(gaps)

def assembly_stats(hits):
    for r in hits:
        df = pd.DataFrame(hits[r])
        fig = plt.figure()
        ax = sns.violinplot(x="qid", y="idnt", data=df, order=["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"])
        plt.savefig("ref_mono_d.png")

ref = load_fasta(ref_file)
ref_hits = load_string_decomposition(ref_hits_file)
cnt_identity(ref_hits)
cnt_gaps(ref_hits)
assembly_stats(ref_hits)
