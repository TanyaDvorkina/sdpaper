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

import edlib
import re

import numpy as np; np.random.seed(0)
import seaborn as sns; sns.set()
import pandas as pd

prefix = "AC_hmmer"
path = "/Sid/tdvorkina/monomers/sdpaper"
reads_file = os.path.join(path, "data_extended/centromeric_reads.fasta")
monomers_file = "../data/DXZ1_inferred_monomers_single.fa"
full = True

def load_alns(filename):
    res = []
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            read, m, m_s, m_e, m_i, tm, tm_s, tm_e, tm_i = ln.strip().split("\t")
            res.append({"r": read, "m": m, "m_s": int(m_s), "m_e": int(m_e), "m_i": float(m_i),
                                   "tm": tm, "tm_s": int(tm_s), "tm_e": int(tm_e), "tm_i": float(tm_i)})
    return res

def load_alns_map(filename):
    res = {}
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            read, m, m_s, m_e, m_i, tm, tm_s, tm_e, tm_i = ln.strip().split("\t")
            if read not in res:
                res[read] = []
            strand = "+"
            if int(m_s) > int(m_e):
                m_s, m_e = m_e, m_s
                strand = "-"
            res[read].append({"r": read,  "m": m, "m_s": int(m_s), "m_e": int(m_e), "m_i": float(m_i), "strand": strand,
                         "tm": tm, "tm_s": int(tm_s), "tm_e": int(tm_e), "tm_i": float(tm_i)})
    for r in res:
        res[r] = sorted(res[r], key=lambda x: (x["m_s"], -x["m_e"]))
    return res

def replacer(r):
    chars = ["_", "=", "|", ":"]
    r = r.lower()
    for c in chars:
        r = r.replace(c, "-")
    r = r.replace("--", "-")
    return r

def load_fasta(filename, rep = False):
    records = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
    new_records = {}
    if rep:
        for r in records:
            new_records[replacer(r)] = records[r].upper()
    else:
        for r in records:
            new_records[r] = records[r].upper()
    return new_records

def draw_heatmap(alns, filename):
    mp = {}
    alf = [" ", "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "?"]
    for a in alf[1:]:
        mp[a] = {}
        for aa in alf[1:]:
            mp[a][aa] = 0
    for a in alns:
        if a["m"] != "-" and a["tm"] != "-" and not a["m"].islower():
            if a["tm"] != a["m"]:
                mp[a["tm"]][a["m"]] += 1
    dtf = pd.DataFrame(mp)
    print(dtf)
    ax = sns.heatmap(dtf, annot=True, fmt="d", annot_kws={'size':10}, vmin=0, vmax=200)
    plt.savefig(filename)


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

def convert_to_abc(name):
    res = name.split("_")[0]
    if "fake" in name:
        res = 'M'
    return res

def get_color(aln, reads, monomers):
    true_m = aln["tm"]
    s, e = aln["m_s"], aln["m_e"]
    seq = reads[aln["r"]].seq[s: e+1]
    if aln["strand"] == "-":
        seq = seq.reverse_complement()
    idnts = []
    for m in monomers:
        idnt = aai([seq, monomers[m].seq])
        idnts.append([idnt, convert_to_abc(m)])
    idnts = sorted(idnts, key=lambda x: -x[0])
    print(true_m, idnts[0], idnts[1])
    if idnts[0][1] == true_m:
        return "b", idnts[0][0], idnts[1][0]
    elif idnts[1][1] == true_m:
        return "g", idnts[0][0], idnts[1][0]
    else:
        return "r", idnts[0][0], idnts[1][0]

def process_unpredicted_ac(alns_ac, alns_sd, reads, monomers):
    res = {"r": [[], []], "b": [[], []], "g": [[], []]}
    covered, total = 0, 0
    alns = []
    for r in alns_ac:
        j = 0
        for i in range(len(alns_ac[r])):
            if alns_ac[r][i]["m"] == "?":
                ac_start, ac_end = alns_ac[r][i]["m_s"], alns_ac[r][i]["m_e"]
                while j < len(alns_sd[r]) and alns_sd[r][j]["m_e"] < ac_start + 50:
                    j += 1
                if j < len(alns_sd[r]):
                    if min(alns_sd[r][j]["m_e"], ac_end) - max(alns_sd[r][j]["m_s"], ac_start) > 100 \
                        and alns_sd[r][j]["m"] != "?" and alns_sd[r][j]["m"] != "-":
                        covered += 1
                        alns.append(alns_sd[r][j])
                if j < len(alns_sd[r]):
                    if alns_sd[r][j]["tm"] != "?":
                        total += 1
                else:
                    total += 1

    print("Covered ", covered, "Total", total, "%", covered*100/total)
    for a in alns:
        color, x, y = get_color(a, reads, monomers)
        print(color, x, y)
        res[color][0].append(x)
        res[color][1].append(y)
    print("r", len(res["r"][0]), len(res["r"][0])/covered)
    print("g", len(res["g"][0]), len(res["g"][0])/covered)
    print("b", len(res["b"][0]), len(res["b"][0])/covered)
    return res

def draw_2d(unpredicted_monomers, filename):
    fig = plt.figure()

    x, y = 60, 100
    sns.kdeplot(unpredicted_monomers["b"][0], unpredicted_monomers["b"][1], cmap="Blues", shade=False, shade_lowest=False)
    sns.kdeplot(unpredicted_monomers["g"][0], unpredicted_monomers["g"][1], cmap="Greens", shade=False, shade_lowest=False)
    #sns.kdeplot(unpredicted_monomers["r"][0], unpredicted_monomers["r"][1], cmap="Reds", shade=False, shade_lowest=False)
    plt.plot(unpredicted_monomers["r"][0], unpredicted_monomers["r"][1], 'ro')
    plt.plot([x, y], [x, y], '--')
    plt.ylim(x, y)
    plt.xlim(x, y)
    plt.xlabel('Best monomer identity')
    plt.ylabel('Second best monomer identity')
    plt.savefig(filename)

def draw_dif_hist(unmonomers, filename):
    fig = plt.figure()
    sns.distplot([unmonomers["b"][0][i] - unmonomers["b"][1][i] for i in range(len(unmonomers["b"][0]))], norm_hist = False, kde=False, hist=True, color="b")
    sns.distplot([unmonomers["g"][0][i] - unmonomers["g"][1][i] for i in range(len(unmonomers["g"][0]))], norm_hist = False, kde=False, hist=True, color="g")
    sns.distplot([unmonomers["r"][0][i] - unmonomers["r"][1][i] for i in range(len(unmonomers["r"][0]))], norm_hist = False, kde=False, hist=True, color="r")
    plt.xlabel('Best monomer - second best')
    plt.savefig(filename)

if full:
    alns_ac = load_alns_map(os.path.join(path, "AC_hmmer", "monomer_alns.tsv")) 
    alns_sd = load_alns_map(os.path.join(path, "SD", "monomer_alns.tsv")) 
    reads = load_fasta(reads_file, True)
    monomers = load_fasta(monomers_file)
    unpredicted_monomers = process_unpredicted_ac(alns_ac, alns_sd, reads, monomers)
    draw_2d(unpredicted_monomers, "./unpredicted_2D.png")
    draw_dif_hist(unpredicted_monomers, "./unpredicted_hist.png")
else:
    alns = load_alns(os.path.join(path, prefix, "monomer_alns.tsv"))
    draw_heatmap(alns, os.path.join(path, prefix, "heatmap_mm.png"))

