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

import numpy as np; np.random.seed(0)
import seaborn as sns; sns.set()
import pandas as pd

import edlib
import re

prefix = "SD"
path = "/Sid/tdvorkina/monomers/sdpaper"
data_map = {"SD":{"read_hits": os.path.join(path, "SD", "decomposition_reads.tsv"), "ref_hits": os.path.join(path, "SD", "decomposition_cenX_t2t.tsv")},
            "AC_hmmer": {"read_hits": os.path.join(path, prefix, "decomposition_reads.tsv"), "ref_hits": os.path.join(path, "SD", "decomposition_cenX_t2t.tsv")}}
reads_file = os.path.join(path, "data/centromeric_reads.fasta")

ref_file = "/Bmo/miheenko/git/tandemQUAST/t2t7/cenx-t2t7.masked.fa"
monomers_file = "/home/tdvorkina/projects/centroFlye/longreads_decomposer/git/sdpaper/data/DXZ1_inferred_monomers_single.fa"

abnormal_hors = {"FK": ["ABCDEFGHIJFGHIJKL", "ABCDEFGHIJKGHIJKL"]} #, "KL": ["ABCDEFGHIJKA", "ABCDEFGHIJLA"]}
abnormal_pos = {"FK": -6, "KL": -2}

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

def load_fasta(filename, replace = False):
    records = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
    new_records = {}
    for r in records:
        if replace:
            new_records[replacer(r)] = records[r].upper()
        else:
            new_records[r] = records[r].upper()
    return new_records

def load_alns(filename):
    res_read = []
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            read, m, m_s, m_e, m_i, tm, tm_s, tm_e, tm_i = ln.strip().split("\t")
            res_read.append({"r": read, "m": m, "m_s": int(m_s), "m_e": int(m_e), "m_i": float(m_i)})
    return res_read

def load_string_decomposition(filename, identity = 70):
    reads_mapping = []
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            sseqid, qseqid, sstart, send, idnt = ln.strip().split("\t")[:5]
            qseqid = qseqid.split()[0]
            s, e, idnt = int(sstart), int(send), float(idnt)
            if idnt > identity:
                reads_mapping.append({"r": replacer(sseqid), "m": convert_to_abc(qseqid), "m_s": s, "m_e": e, "m_i": idnt})
    reads_mapping = sorted(reads_mapping, key=lambda x: (x["m_s"], -x["m_e"]))
    return reads_mapping

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

def get_abnormal_hor_monomer(alns, seqs, hor_tp):
    prev_read = ""
    cur_seq = ""
    res = []
    for i in range(len(alns)):
        a = alns[i]
        if a["r"] != prev_read:
            cur_seq = ""
            prev_read = a["r"]
        cur_seq += a["m"]
        if cur_seq[-len(abnormal_hors[hor_tp][0]):] in abnormal_hors[hor_tp]:
            print(cur_seq[-len(abnormal_hors[hor_tp][0]):])
            s, e = alns[i + abnormal_pos[hor_tp]]["m_s"], alns[i + abnormal_pos[hor_tp]]["m_e"]
            if s < e:
                print(s, e, str(seqs[a["r"]].seq[s:e]))
                res.append({"hor": cur_seq[abnormal_pos[hor_tp]], "seq": seqs[a["r"]].seq[s:e]})
            else:
                print(s, e, str(seqs[a["r"]].seq[e:s].reverse_complement()))
                res.append({"hor": cur_seq[abnormal_pos[hor_tp]], "seq": seqs[a["r"]].seq[e:s].reverse_complement()})
    print("Abnormal num ", len(res))
    return res

def calculate_identities_mon(alns, monomers, hor_tp):
    ms = list(hor_tp)
    res = {}
    for c in ms:
        res[c] = []
    for a in alns:
        for m in monomers:
            for c in list(hor_tp):
                if m.startswith(c):
                    s1 = monomers[m].seq
                    s2 = a["seq"]
                    res[c].append(aai([s1, s2]))
            s1 = monomers[m].seq
            s2 = a["seq"]
        #     print(m, aai([s1, s2]))
        # print("")
    s1, s2 = None, None
    mons = list(hor_tp)
    for m in monomers:
        if m.startswith(mons[0]):
            s1 = monomers[m].seq
        if m.startswith(mons[1]):
            s2 = monomers[m].seq
    res[mons[0]].append(aai([s1, s2]))
    res[mons[1]].append(aai([s1, s2]))
    return res

def calculate_identities(alns):
    res = {}
    for i in range(len(alns)):
        res[str(i + 1)] = {}
        for j in range(len(alns)):
            res[str(i + 1)][str(j + 1)] = aai([alns[i]["seq"], alns[j]["seq"]])
    return res


def draw_plot(idnts, mons, filename):
    fig = plt.figure()
    m_lst = list(mons)
    plt.plot(idnts[m_lst[0]], idnts[m_lst[1]],'bo', alpha=0.5)
    plt.xlabel("Identity " + m_lst[0])
    plt.ylabel("Identity " + m_lst[1])
    plt.xlim(65, 100)
    plt.ylim(65, 100)
    plt.savefig(filename)


def draw_heatmap(idnts, idnts_m, filename):
    mons = list(idnts_m.keys())
    for i in range(len(idnts_m[mons[0]]) - 1):
        c = mons[0]
        idnts[str(i + 1)][c] = idnts_m[c][i] 
    idnts[mons[1]] = {}
    for i in range(len(idnts_m[mons[1]]) - 1):
        c = mons[1]
        idnts[c][str(i + 1)] = idnts_m[c][i]
    idnts[mons[1]][mons[0]] = idnts_m[mons[1]][-1]
    for m in idnts:
        for mm in idnts[m]:
            idnts[m][mm] = int(idnts[m][mm])
    fig = plt.figure()
    dtf = pd.DataFrame(idnts)
    cols = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"]
    cols.append(mons[1])
    dtf = dtf.ix[:, cols]
    rows = cols[:-1]
    rows.append(mons[0])
    dtf = dtf.reindex(rows)
    print(dtf)
    ax = sns.heatmap(dtf, annot=True, fmt="d")
    plt.savefig(filename)


ref = load_fasta(ref_file, True)
ref_name = list(ref.keys())[0]
alns_read = load_alns(os.path.join(path, prefix, "monomer_alns.tsv"))
alns_ref = load_string_decomposition(data_map[prefix]["ref_hits"])
reads = load_fasta(reads_file, True)
monomers = load_fasta(monomers_file)

for h in abnormal_hors:
    print(h)
    print("Calculating reads stats")
    nuc_seqs = get_abnormal_hor_monomer(alns_read, reads, h)
    idnts = calculate_identities_mon(nuc_seqs, monomers, h)
    draw_plot(idnts, h, os.path.join(path, prefix, "reads_" + h + "_plot.png"))

    print("Calculating assembly stats")
    nuc_seqs = get_abnormal_hor_monomer(alns_ref, ref, h)
    idnts_m = calculate_identities_mon(nuc_seqs, monomers, h)
    idnts = calculate_identities(nuc_seqs)
    print(idnts_m)
    draw_heatmap(idnts, idnts_m, os.path.join(path, prefix, "heatmap_" + h + ".png")) 