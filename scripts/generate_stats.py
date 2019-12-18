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


prefix = "AC_hmmer"
path = "/Sid/tdvorkina/monomers/sdpaper"
data_map = {"SD":{"read_hits": os.path.join(path, "SD", "decomposition_reads.tsv"), "ref_hits": os.path.join(path, "SD", "decomposition_cenx_flye_polished.tsv")},
            "AC_hmmer": {"read_hits": os.path.join(path, prefix, "decomposition_reads.tsv"), "ref_hits": os.path.join(path, "SD", "decomposition_cenx_flye_polished.tsv")}}

tm_res = "/Bmo/miheenko/git/tandemQUAST/cenx_flye_polished/cenx-flye-polished_alignment.bed"
reads_file = os.path.join(path, "data/centromeric_reads.fasta")
ref_file = "/Bmo/miheenko/git/tandemQUAST/assemblies/cenx_flye_polished.fasta"
monomers_file = "/home/tdvorkina/projects/centroFlye/longreads_decomposer/git/sdpaper/data/DXZ1_inferred_monomers_single.fa" 

def load_good_reads():
    trash = []
    filename = os.path.join(path, "data", "gooood_reads.txt")
    res = []
    with open(filename, "r") as fin:
        for ln in fin.readlines():
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

def convert_to_abc(name):
    res = name.split("_")[0]
    if "fake" in name:
        res = 'M'
    return res

def load_fasta(filename):
    records = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
    new_records = {}
    for r in records:
        new_records[replacer(r)] = records[r].upper()
    return new_records

def load_string_decomposition(filename, identity = 55):
    reads_mapping = {}
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            sseqid, qseqid, sstart, send, idnt = ln.strip().split("\t")[:5]
            qseqid = qseqid.split()[0]
            if replacer(sseqid) not in reads_mapping:
                    reads_mapping[replacer(sseqid)] = []
            s, e, idnt = int(sstart), int(send), float(idnt)
            rev = False
            if qseqid.endswith("'"):
                rev = True
                qseqid = qseqid[:-1]
            if idnt > identity:
                reads_mapping[replacer(sseqid)].append({"qid": qseqid, "s": s, "e": e, "rev": rev, "idnt": idnt})
    for r in reads_mapping:
        reads_mapping[r] = sorted(reads_mapping[r], key=lambda x: (x["s"], -x["e"]))
    return reads_mapping

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

def extract_monomer_string(hits, start, end, orientation, name, avg_monomer_len):
    l, r = find_ends_indexes(hits, start, end, name)
    if l > r:
        print("No alignments!")
        return ""
    print(start, end, hits[l]["s"], hits[r]["e"], l, r)
    res = ""
    num = 0
    for i in range(l, r + 1):
        c = convert_to_abc(hits[i]["qid"])
        if orientation == "-":
            if hits[i]["rev"]:
                c = c
            else:
                c = c + "'"
        if c.endswith("'"):
            c = c[:-1].lower()
            print("wrong orientation")
        res += c
        num += 1
        if i != r:
            dist = max(0, hits[i + 1]["s"] - hits[i]["e"])
            q, rr = divmod(dist, avg_monomer_len)
            m_num = int(q + int(2 * rr >= avg_monomer_len))
            if rr >= avg_monomer_len:
                print(rr)
            for _ in range(m_num):
                res += "?"
                num += 1
    ind = 0
    if orientation == "-":
        res = res[::-1]
    print(num, num/12)
    return res

def remove_corner_errors(niceAlign, stats_sq):
    if niceAlign["query_aligned"][0] == "-":
        i = 0
        while niceAlign["query_aligned"][i] == "-":
            if niceAlign["query_aligned"][i] == "?":
                stats_sq["-"]["?"] -= 1
            else:
                stats_sq["-"]["m"] -= 1
            i += 1

    if niceAlign["query_aligned"][-1] == "-":
        i = -1
        while niceAlign["query_aligned"][i] == "-":
            if niceAlign["query_aligned"][i] == "?":
                stats_sq["-"]["?"] -= 1
            else:
                stats_sq["-"]["m"] -= 1
            i -= 1

    if niceAlign["target_aligned"][0] == "-":
        i = 0
        while niceAlign["target_aligned"][i] == "-":
            if niceAlign["target_aligned"][i] == "?":
                stats_sq["?"]["-"] -= 1
            else:
                stats_sq["m"]["-"] -= 1
            i += 1

    if niceAlign["target_aligned"][-1] == "-":
        i = -1
        while niceAlign["target_aligned"][i] == "-":
            if niceAlign["target_aligned"][i] == "?":
                stats_sq["?"]["-"] -= 1
            else:
                stats_sq["m"]["-"] -= 1
            i -= 1
    return stats_sq

def cnt_edist(lst):
    if len(str(lst[0])) == 0:
        return -1
    if len(str(lst[1])) == 0:
        return -1
    result = edlib.align(str(lst[0]), str(lst[1]), mode="NW", task="path")
    niceAlign = edlib.getNiceAlignment(result, str(lst[0]), str(lst[1]))
    print(niceAlign["query_aligned"])
    print(niceAlign["matched_aligned"])
    print(niceAlign["target_aligned"])
    print("")
    return niceAlign

def get_monomer_alignment(ref_hits, read_hits, ref_start, ref_end, read_start, read_end, name, orientation, avg_monomer_len):
    ref_monomer_str = extract_monomer_string(ref_hits, ref_start, ref_end, "+", name, avg_monomer_len)
    read_monomer_str = extract_monomer_string(read_hits, read_start, read_end, orientation, name, avg_monomer_len)
    print("Ref sz = " + str(len(ref_monomer_str)) + " Approx " + str((ref_end - ref_start)//avg_monomer_len))
    print("Read sz = " + str(len(read_monomer_str)) + " Approx " + str((read_end - read_start)//avg_monomer_len))
    if "'" in ref_monomer_str or "'" in read_monomer_str:
        print(ref_monomer_str)
        print(read_monomer_str)
        exit(-1)
    l, r = find_ends_indexes(read_hits, read_start, read_end, name)
    if (r - l + 1) > 0:
        avg_monomer_len_c = 0
        avg_monomer_idnt = 0
        for i in range(l, r + 1):
            avg_monomer_len_c += read_hits[i]["e"] - read_hits[i]["s"]
            avg_monomer_idnt += read_hits[i]["idnt"]
        print(avg_monomer_len_c//(r- l +1), avg_monomer_idnt//(r - l +1))
    else:
        print("WAT")
        print(l,r)
        exit(-1)

    niceAlign = cnt_edist([read_monomer_str, ref_monomer_str])
    stats_sq = {"?": {"?": 0, "m": 0, "mm": 0, "-": 0}, "m": {"?": 0, "m": 0, "mm": 0, "-": 0}, "-": {"?": 0, "m": 0, "mm": 0, "-": 0}}
    if niceAlign == -1:
        for c in ref_monomer_str:
            stats["mm_pairs"].append([c, "?"])
        exit(-1)
        return -1, -1

    for i in range(len(niceAlign["query_aligned"])):
        if niceAlign["query_aligned"][i] == niceAlign["target_aligned"][i]:
            if niceAlign["query_aligned"][i] == "?":
                stats_sq["?"]["?"] += 1
            else:
                stats_sq["m"]["m"] += 1
        elif niceAlign["query_aligned"][i] == "-":
            if niceAlign["target_aligned"][i] == "?":
                stats_sq["-"]["?"] += 1
            else:
                stats_sq["-"]["m"] += 1
        elif niceAlign["target_aligned"][i] == "-":
            if niceAlign["query_aligned"][i] == "?":
                stats_sq["?"]["-"] += 1
            else:
                stats_sq["m"]["-"] += 1
        else:
            if niceAlign["target_aligned"][i] == "?":
                stats_sq["m"]["?"] += 1
            elif niceAlign["query_aligned"][i] == "?":
                stats_sq["?"]["m"] += 1
            else:
                stats_sq["m"]["mm"] += 1

    stats_sq = remove_corner_errors(niceAlign, stats_sq)
    return stats_sq

def get_avg_len(seqs):
    res = 0
    for m in seqs:
        res += len(seqs[m].seq)
    res //= len(seqs)
    return res

ref = load_fasta(ref_file)

reads = load_fasta(reads_file)
total_len = 0
for r in reads:
    total_len += len(reads[r].seq)

monomers = load_fasta(monomers_file)
avg_monomer_len = get_avg_len(monomers)

reads_hits = load_string_decomposition(data_map[prefix]["read_hits"], 65)
ref_hits = load_string_decomposition(data_map[prefix]["ref_hits"], 70)
ref_name = list(ref_hits.keys())[0]

overall_stats = {"?": {"?": 0, "m": 0, "mm": 0, "-": 0}, "m": {"?": 0, "m": 0, "mm": 0, "-": 0}, "-": {"?": 0, "m": 0, "mm": 0, "-": 0}}
good_reads = load_good_reads()
aligned_len = 0
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
            orient = "+"
            if read_start > read_end:
                read_start, read_end = read_end, read_start
                orient = "-"
            aligned_len += read_end - read_start
            print(name, orient)
            stats = get_monomer_alignment(ref_hits[ref_name], reads_hits[name],\
                                    ref_start, ref_end, read_start, read_end, name, orient, avg_monomer_len)
            print(stats)
            for p1 in stats:
                for p2 in stats[p1]:
                    overall_stats[p1][p2] += stats[p1][p2]
            print("")

print("Aligned len ", aligned_len, " Total len ", total_len, " %", aligned_len*100/total_len)
total_syms = 0
for p1 in overall_stats:
    for p2 in overall_stats[p1]:
        total_syms += overall_stats[p1][p2]
        print(p1, p2, overall_stats[p1][p2])
errors = total_syms - overall_stats["?"]["?"] - overall_stats["m"]["m"]
print("Percent error ", errors*100/total_syms)