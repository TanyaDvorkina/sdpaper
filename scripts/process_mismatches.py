#!/usr/bin/env python3
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os
from os import listdir
from os.path import join, isfile
import sys

import numpy as np; np.random.seed(0)
import seaborn as sns; sns.set()
import pandas as pd

prefix = "AC_hmmer"
path = "/Sid/tdvorkina/monomers/sdpaper"

full = False

def load_alns(filename):
    res = []
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            read, m, m_s, m_e, m_i, tm, tm_s, tm_e, tm_i = ln.strip().split("\t")
            res.append({"r": read, "m": m, "m_s": int(m_s), "m_e": int(m_e), "m_i": float(m_i),
                                   "tm": tm, "tm_s": int(tm_s), "tm_e": int(tm_e), "tm_i": float(tm_i)})
    return res

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

def identify_monomers(seq, monomers):
    pass

def process_unpredicted(alns, reads, monomers):
    res = []
    for a in alns:
        if a["m"] == "?" and a["tm"] != "?":
            continue
    return res

def draw_unpredicted(unpredicted_monomers):
    pass

alns = load_alns(os.path.join(path, prefix, "monomer_alns.tsv"))
draw_heatmap(alns, os.path.join(path, prefix, "heatmap_mm.png"))

if full:
    reads = load_fasta()
    monomers = load_fasta()
    unpredicted_monomers = process_unpredicted(alns, reads, monomers)