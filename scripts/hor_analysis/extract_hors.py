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

import re

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

def simplify(s, monomers):
    return monomers[s]
    #return "/".join([x.split("_")[1] for x in s.split(".")[:2]])

def load_string_decomposition(filename, monomers):
    reads_mapping = {}
    reads_length = {}
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            if ln.startswith("SRR9087598.153892 ") or ln.startswith("SRR9087597.98571 "):
                break
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
    print("Initial Read num", len(reads_mapping))
    for r in reads_mapping:
        ln = int(r.split()[2][len("length="):])
        reads_length[r] = ln
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

    total_len = 0
    for r in new_reads_mapping:
        total_len += reads_length[r]
    print("Num", len(new_reads_mapping), "Length", total_len)
    return new_reads_mapping

def extract_hors_old(dec):
    hors_set = {}
    for r in dec:
        cur_hors = []
        cur_lx_min = None
        is_first = True
        for h in dec[r]:
            if h["qid"] not in cur_hors:
                cur_hors.append(h["qid"])
                if cur_lx_min == None or cur_hors[cur_lx_min] > h["qid"]:
                    cur_lx_min = len(cur_hors) - 1
            else:
                cur_hors = cur_hors[cur_lx_min:] + cur_hors[:cur_lx_min]
                s = "_".join(cur_hors)
                if not is_first:
                    if s not in hors_set:
                        hors_set[s] = 0
                    hors_set[s] += 1
                is_first = False
                cur_hors = [h["qid"]]
                cur_lx_min = 0
    hors_lst = []
    for h in hors_set:
        hors_lst.append([h, hors_set[h]])

    hors_lst = sorted(hors_lst, key=lambda x: -x[1])
    return hors_lst

def cyclic_shift(s):
    min_s = s
    lst = s[1:-1].split("_")
    for c in range(1, len(lst)):
        cur_s = "_" + "_".join(lst[c:] + lst[:c]) + "_"
        if cur_s < min_s:
            min_s = cur_s
    return min_s

def extract_hors_new(dec):
    hors_set = {}
    hors_set_plus = {}
    for r in dec:
        c_hors_set = {}
        c_hors_set_plus = {}
        r_s = "_".join([x["qid"] for x in dec[r]])
        r_s = "_" + r_s + "_"
        for k in range(2, 20):
            for i in range(len(dec[r]) - k):
                cur_s = [dec[r][j]["qid"] for j in range(i, i + k)]
                cur_ss = "_".join(cur_s)
                dcur_s = "_" + cur_ss + "_" + cur_ss + "_"
                if "_" + cur_ss + "_" not in c_hors_set:
                    cnt = len(re.findall('(?=' + dcur_s + ')', r_s))
                    if cnt > 0:
                        c_hors_set["_" + cur_ss + "_"] = cnt
                        c_hors_set_plus["_" + cur_ss + "_"] = cnt + 1
        if "_37_33_32_34_38_40_25_29_35_36_30_31_27_28_53_" in c_hors_set or "_47_37_33_32_34_38_40_25_29_35_36_30_31_27_28_" in c_hors_set:
            print(r)
            print(c_hors_set.keys())
        for it in c_hors_set:
            if it not in hors_set:
                hors_set[it] = 0
                hors_set_plus[it] = 0
            hors_set[it] += c_hors_set[it]
            hors_set_plus[it] += c_hors_set_plus[it]
    exit(-1)
    hors_lst = []
    for h in hors_set:
        hors_lst.append([h, hors_set[h], hors_set_plus[h]])

    hors_lst = sorted(hors_lst, key=lambda x: -x[1])

    new_hors_lst = []
    cntt_new = 0
    for i in range(len(hors_lst)):
        h = hors_lst[i][0]
        toadd = True
        min_h = cyclic_shift(h)
        for j in range(i):
            cur_h = cyclic_shift(hors_lst[j][0])
            s = cur_h
            for k in range(5):
                if min_h == s:
                    toadd = False
                    break
                s += cur_h[1:]
            if not toadd:
                break
        if toadd:
            new_hors_lst.append(hors_lst[i])
            # if len(new_hors_lst) > 3:
            #     exit(-1)
    i = 1
    for h in new_hors_lst:
        nm = len(h[0][1:-1].split("_"))
        print("\t".join([str(i),str(nm),h[0], str(h[1]), str(h[2])]))
        i += 1
    return new_hors_lst

def find_set(ind, parent):
    if ind == parent[ind]:
        return ind
    parent[ind] = find_set(parent[ind], parent)
    return parent[ind]

def union_sets(a, b, parent, rank):
    a = find_set(a, parent)
    b = find_set(b, parent)
    if a != b:
        if rank[a] < rank[b]:
            a, b = b, a
        parent[b] = a
        if rank[a] == rank[b]:
            rank[a] += 1

def divide_into_clusters(hors, th):
    clusters = []
    parent = [i for i in range(len(hors))]
    rank = [0 for _ in range(len(hors))]
    for i in range(len(hors)):
        for j in range(i + 1, len(hors)):
            s1 = set(hors[i][0][1:-1].split("_"))
            s2 = set(hors[j][0][1:-1].split("_"))
            inter = len(s1 & s2)
            if inter >= th:
                union_sets(i, j, parent, rank)

    clusters_id = {}
    for i in range(len(hors)):
        ind = find_set(i, parent)
        if ind not in clusters_id:
            clusters_id[ind] = []
        clusters_id[ind].append(hors[i])

    for cl in clusters_id:
        isgood = False
        for h in clusters_id[cl]:
            if h[1] > 20:
                isgood = True
                break
        if isgood:
            clusters.append(clusters_id[cl])

    print(len(clusters))
    for cl in clusters:
        for h in sorted(cl, key=lambda x: -x[1]):
            print(h[0], h[1])
        print("")

    return clusters

def load_lst():
    filename = "known_hors.tsv"
    res = []
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            #print(ln)
            name, seq = ln.strip().split()
            res.append([seq, name])
    return res

known_hors_seq = load_lst()

monomers, c_monomers = load_enum("/Sid/tdvorkina/monomers/new_monomers/monomer_clustering/cluter_represent_ed6_enumarated_order.fasta")
dec = load_string_decomposition("/Sid/tdvorkina/monomers/new_monomers/cent_reads/SRR9087598/decomposition_cent_reads_raw_ed6.tsv", monomers)
dec2 = load_string_decomposition("/Sid/tdvorkina/monomers/new_monomers/cent_reads/SRR9087597/decomposition_cent_reads_raw_ed6.tsv", monomers)

for d in dec2:
    dec[d] = dec2[d]
print("Reads num", len(dec))

hors_lst = extract_hors_new(dec)

res = {}
for kh in known_hors_seq:
    res[kh[1]] = [kh[0]]

unused_hors = []
for h in hors_lst:
    used = False
    for kh in known_hors_seq:
        h_s = set(h[0][1:-1].split("_"))
        kh_s = set(kh[0][1:-1].split("_"))
        if len(h_s & kh_s) > max(2, 0.5*min(len(h_s), len(kh_s))):
            res[kh[1]].append(h)
            used = True
    if not used:
        unused_hors.append(h)

for k in res:
    print(k, res[k][0])
    for it in res[k][1:10]:
        ln = "\t".join(it[0][1:-1].split("_"))
        nm = len(it[0][1:-1].split("_"))
        print(" ", ln, nm, it[1])

print("")
for h in unused_hors:
    if h[1] > 10:
        print(h[0] +"\t" + str(h[1]))




