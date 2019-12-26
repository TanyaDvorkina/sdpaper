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
import re

import numpy as np; np.random.seed(0)
import seaborn as sns; sns.set()
import pandas as pd


path = "/Sid/tdvorkina/monomers/sdpaper"
data_map = {"SD":{"read_hits": os.path.join(path, "SD", "decomposition_reads.tsv"), "ref_hits": os.path.join(path, "SD", "decomposition_cenX_t2t.tsv")},
            "AC_hmmer": {"read_hits": os.path.join(path, "AC_hmmer", "decomposition_reads.tsv"), "ref_hits": os.path.join(path, "SD", "decomposition_cenX_t2t.tsv")}}
reads_file = os.path.join(path, "data_extended/centromeric_reads.fasta")

monomers_file = "../data/DXZ1_inferred_monomers_single.fa" 

tm_res =  os.path.join(path, "data_extended/cenx-t2t7_alignment.bed")
ref_file =  os.path.join(path, "data_extended/cenx-t2t7.fa")

def replacer(r):
    chars = ["_", "=", "|", ":"]
    r = r.lower()
    for c in chars:
        r = r.replace(c, "-")
    r = r.replace("--", "-")
    return r

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

def load_fasta(filename):
    records = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
    new_records = {}
    for r in records:
        new_records[replacer(r)] = records[r].upper()
    return new_records


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

def analyze_ncrf(filename, reads_aligned):
    total_idnt = 0
    total_gaps = 0
    total_monomers = 0
    total_length = 0
    covered_length = 0
    with open(filename, "r") as fin:
        read_name = None
        for ln in fin.readlines():
            lst = ln.strip().split()
            if len(lst) > 4:
                read_name, rlen, alen, c1, c2, aln_seq1 = replacer(lst[0]), int(lst[1]), int(lst[2][:-2]), int(lst[3].split("-")[0]), int(lst[3].split("-")[1]), lst[4]
                aln_seq1 = aln_seq1.upper()
            elif len(lst) > 3:
                name, rplen, aln_seq2 = lst[0], int(lst[1][:-2]), lst[3]
                aln_seq2 = aln_seq2.upper()
                if read_name in reads_aligned:
                    start, end = reads_aligned[read_name]["s"], reads_aligned[read_name]["e"]
                    print(read_name, c1, c2, start, end)
                    if start > c2 or end < c1:
                        total_length += end - start
                        continue
                    if c1 < start:
                        pos_aln, pos_target = 0, 0
                        while c1 + pos_target < start:
                            if aln_seq1[pos_aln] != "-":
                                pos_target += 1
                            pos_aln += 1
                        aln_seq1 = aln_seq1[pos_aln:]
                        aln_seq2 = aln_seq2[pos_aln:]

                    if c2 > end:
                        pos_aln, pos_target = 0, 0
                        while c2 - pos_target > end + 1:
                            if aln_seq1[len(aln_seq1) - pos_aln - 1] != "-":
                                pos_target += 1
                            pos_aln += 1
                        aln_seq1 = aln_seq1[:-pos_aln]
                        aln_seq2 = aln_seq2[:-pos_aln]
                    aln_seq1 = aln_seq1.replace("-","")
                    aln_seq2 = aln_seq2.replace("-","")
                # if name.endswith("-"):
                #     aln_seq2 = str(Seq(aln_seq2).reverse_complement())
                    idd = aai([aln_seq1, aln_seq2])
                    print(len(aln_seq1), len(aln_seq2), idd)
                    covered_length += len(aln_seq1)
                    total_length += end - start
                    total_idnt += idd
                    total_monomers += 1

    print("Covered ", 100*covered_length/total_length, covered_length, total_length)
    print("Avg identity", total_idnt/total_monomers)
    print("Gaps %", 100*total_gaps/total_monomers)

def load_string_decomposition(filename, identity = 69):
    reads_mapping = {}
    filtered = 0
    good_reads = load_good_reads()
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            sseqid, qseqid, sstart, send, idnt = ln.strip().split("\t")[:5]
            qseqid = qseqid.split()[0]
            sseqid = replacer(sseqid)
            is_good = False
            for c in good_reads:
                if sseqid.startswith(c):
                    is_good = True
                    break
            if is_good:
                if sseqid not in reads_mapping:
                        reads_mapping[sseqid] = []
                s, e, idnt = int(sstart), int(send), float(idnt)
                rev = False
                if qseqid.endswith("'"):
                    rev = True
                    qseqid = qseqid[:-1]
                if idnt > identity:
                    reads_mapping[sseqid].append({"qid": qseqid, "s": s, "e": e, "rev": rev, "idnt": idnt})
                else:
                    filtered += 1
    print("Unreliable ", filtered)
    for r in reads_mapping:
        reads_mapping[r] = sorted(reads_mapping[r], key=lambda x: (x["s"], -x["e"]))
    return reads_mapping

def decomposition_filtering(reads_mapping):
    new_mapping = {}
    bad_reads_num = 0
    for r in reads_mapping:
        cnt_rev = 0
        for i in range(len(reads_mapping[r])):
            h1 = reads_mapping[r][i]
            if h1["rev"]:
                cnt_rev += 1
        if cnt_rev >= 36 or len(reads_mapping[r]) - cnt_rev >= 36:
            if cnt_rev == 0 or cnt_rev == len(reads_mapping[r]):
                new_mapping[r] = reads_mapping[r]
            elif cnt_rev < 0.1*len(reads_mapping[r]):
                new_mapping[r] = []
                for i in range(len(reads_mapping[r])):
                    h1 = reads_mapping[r][i]
                    if not h1["rev"]:
                        new_mapping[r].append(h1)
            elif cnt_rev >= 0.9*len(reads_mapping[r]):
                new_mapping[r] = []
                for i in range(len(reads_mapping[r])):
                    h1 = reads_mapping[r][i]
                    if h1["rev"]:
                        new_mapping[r].append(h1)
            else:
                bad_reads_num += 1
        else:
            bad_reads_num += 1

    print("Bad reads " + str(bad_reads_num))
    for r in new_mapping:
        if len(new_mapping[r]) > 0 and new_mapping[r][0]["rev"]:
            new_mapping[r] = sorted(new_mapping[r], key=lambda x: (-x["e"], x["s"]))
        else:
            new_mapping[r] = sorted(new_mapping[r], key=lambda x: (x["s"], -x["e"]))
    return new_mapping

def find_ends_indexes(hits, start, end, name):
    i = 0
    while i < len(hits) and hits[i]["s"] < start - 50:
        i += 1
    if i >= len(hits):
        print("Left border violation " + str(name) + " " + str(start) + "," + str(end))
        exit(-1)

    l, r = i, i
    while r < len(hits) and hits[r]["e"] < end + 50:
        r += 1

    if r > len(hits):
        print("Right border violation " + str(name) + " " + str(start) + "," + str(end) )
        exit(-1)
    return l,r - 1

def covered_aligned(hits, reads, pos):
    total_idnt = 0
    total_gaps = 0
    total_monomers = 0
    total_length = 0
    covered_length = 0
    for read in hits:
        if read in pos:
            l, r = find_ends_indexes(hits[read], pos[read]["s"], pos[read]["e"], read)
            total_length += pos[read]["e"] - pos[read]["s"]
            sum_idnt = 0
            sum_gap = 0
            for i in range(l, r + 1):
                sum_idnt += hits[read][i]["idnt"]
                covered_length += hits[read][i]["e"] - hits[read][i]["s"] 
                if i != r:
                    sum_gap += int( max(0, hits[read][i + 1]["s"] - hits[read][i]["e"]) > 10)
            total_idnt += sum_idnt
            total_gaps += sum_gap
            total_monomers += r - l + 1
    print("Covered ", 100*covered_length/total_length, covered_length, total_length)
    print("Avg identity", total_idnt/total_monomers)
    print("Gaps %", 100*total_gaps/total_monomers)

def convert_to_abc(name):
    res = name.split("_")[0]
    if "fake" in name:
        res = 'M'
    return res


reads = load_fasta(reads_file)
monomers = load_fasta(monomers_file)
dxz1_len = 0
for m in monomers:
    dxz1_len += len(monomers[m].seq)

good_reads = load_good_reads()
reads_aligned = {}
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
            read_start, read_end = [int(lst[4]), int(lst[5])]
            orient = "+"
            if read_start > read_end:
                read_start, read_end = read_end, read_start
                orient = "-"
            reads_aligned[name] = {"s": read_start, "e": read_end, "o": orient}


print("NCRF")
analyze_ncrf("/Sid/tdvorkina/monomers/NCRF/decomposition_reads_ncrf.ncrf", reads_aligned)

print("SD")
print(data_map["SD"]["read_hits"])
mappings = load_string_decomposition(data_map["SD"]["read_hits"])
covered_aligned(mappings, reads, reads_aligned)


print("AC")
mappings = load_string_decomposition(data_map["AC_hmmer"]["read_hits"])
covered_aligned(mappings, reads, reads_aligned)