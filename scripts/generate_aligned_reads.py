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

path = "/Sid/tdvorkina/monomers/sdpaper"

def load_fasta(filename):
    records = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
    return records

def load_good_reads():
    trash = ["86722ef4"]
    filename = os.path.join(path, "data_extended", "gooood_reads_t2t7.txt")
    res = []
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            if len(ln) > 1:
                good_read = True
                for r in trash:
                    if ln.startswith(r):
                        good_read = False
                        break
                if good_read:
                    res.append(ln.strip())
    print("Reliable reads " + str(len(res)))
    return res

def replacer(r):
    chars = ["_", "=", "|", ":"]
    r = r.lower()
    for c in chars:
        r = r.replace(c, "-")
    r = r.replace("--", "-")
    return r

def load_reads_positions():
    good_reads = load_good_reads()
    tm_res = os.path.join(path, "data_extended/cenx-t2t7_alignment.bed")
    res = {}
    with open(tm_res, "r") as fin:
        for ln in fin.readlines():
            lst = ln.strip().split()
            name = lst[3]
            is_good = False
            for r in good_reads:
                if name.startswith(r):
                    is_good = True
                    break
            if is_good:
                ref_start, ref_end, read_start, read_end = [int(lst[1]), int(lst[2]), int(lst[4]), int(lst[5])]
                if read_start > read_end:
                    res[name] = [read_end, read_start]
                else:
                    res[name] = [read_start, read_end]
    return res

def save_fasta(filename, orfs):
    with open(filename + ".fasta", "w") as output_handle:
        SeqIO.write(orfs, output_handle, "fasta")

def make_record(seq, name, sid, d=""):
    return SeqRecord(seq, id=sid, name=name, description = d)

reads = load_fasta("/Sid/tdvorkina/monomers/sdpaper/data_extended/centromeric_reads.fasta")
reads_pos = load_reads_positions()

aligned_reads = []
for r in reads:
    r_r = replacer(r)
    if r_r in reads_pos:
        start, end = reads_pos[r_r]
        aligned_reads.append(make_record(reads[r].seq[start: end + 1], r, reads[r].id))

save_fasta("./aligned_reads", aligned_reads) 


