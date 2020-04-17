from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio import SearchIO
from Bio.SeqRecord import SeqRecord

import numpy as np; np.random.seed(123)
import seaborn as sns; sns.set()
import pandas as pd

import sys
import os
from os import listdir
from os.path import isfile, isdir, join

import re
import edlib


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



def load_dec(filename, monomers_seq):
    reads_mapping = {}
    ss = ""
    monomers = set()
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            sseqid, qseqid, sstart, send, idnt  = ln.strip().split("\t")[:5]
            sseqid = sseqid.split()[0]
            if sseqid not in reads_mapping:
                    reads_mapping[sseqid] = []
            s, e, idnt = int(sstart), int(send), float(idnt)
            rev = False
            if qseqid.endswith("'"):
                rev = True
                qseqid = qseqid[:-1]
            qseqid = qseqid.split()[0]
            if idnt >= 75:
                monomers.add(qseqid)
                reads_mapping[sseqid].append({"qid": qseqid, "s": s, "e": e, "rev": rev, "idnt": idnt})
            else:
                reads_mapping[sseqid].append({"qid": "?","s": s, "e": e, "rev": rev, "idnt": idnt})
    print(len(monomers))

    monomers_lst = [x.id for x in monomers_seq]
    monomers_mp = {}
    monomers_mp_r = {}
    for i in range(len(monomers_lst)):
        qid = monomers_lst[i]
        new_qid = "m" + str(i + 1)
        print("\t".join([qid, new_qid]))
        monomers_mp[qid] = new_qid
        monomers_mp_r[new_qid] = qid
        monomers_mp_r[new_qid + "'"] = qid + "'"
    monomers_mp["?"] = "?"
    monomers_mp_r["?"] = "?"
    for c in reads_mapping:
        for j in range(len(reads_mapping[c])):
            if reads_mapping[c][j]["rev"]:
                reads_mapping[c][j]["qid"] = monomers_mp[reads_mapping[c][j]["qid"]] + "'"
            else:
                reads_mapping[c][j]["qid"] = monomers_mp[reads_mapping[c][j]["qid"]]
    new_reads_mapping = {}
    for r in reads_mapping:
        cur_mapping = []
        inside_good = False
        left, right = -1, -1
        for i in range(len(reads_mapping[r])):
            if left == -1 and reads_mapping[r][i]["idnt"] > 95:
                left = i
            if reads_mapping[r][i]["idnt"] > 95:
                right = i
        print(reads_mapping[r][left]["s"], reads_mapping[r][right]["e"])
        #print(reads_mapping[r][left])
        if len(reads_mapping[r][left: right + 1]) > 36:
            new_reads_mapping[r] = sorted(reads_mapping[r][left: right + 1], key = lambda x: x["s"])
    return new_reads_mapping, monomers_mp_r

def cyclic_shift(s):
    min_s = s
    lst = s[1:-1].split("_")
    for c in range(1, len(lst)):
        cur_s = "_" + "_".join(lst[c:] + lst[:c]) + "_"
        if cur_s < min_s:
            min_s = cur_s
    return min_s

def is_equal(s1, s2):
    #return s1 == s2
    for i in range(len(s1)):
        if s1[i:] + s1[:i] == s2:
            return True
    return False

def cyclic_shift2(s1):
    min_s = s1
    for i in range(len(s1)):
        if s1[i:] + s1[:i] < min_s:
            min_s = s1[i:] + s1[:i]
    return min_s

def annotate(seq):
    annotation = []
    prev, cnt, start = '', 0, 0
    for i in range(len(seq)):
        c = seq[i][0]
        if c == prev:
            cnt += 1
        else:
            if prev != '':
                annotation.append([prev, cnt, {"start":start, "end": seq[i-1][2], "hor_full": prev}])
            prev = c
            cnt = 1
            start = seq[i][1]
            idnt = seq[i][4]
    annotation.append([prev, cnt])
    return annotation

def collapse_annotation(annotation):
    new_annotation = []
    prev, cnt, start = "", 0, 0
    for i in range(len(annotation)):
        c, c_cnt = annotation[i][:2]
        if c == prev:
            cnt += c_cnt
        else:
            if prev != "":
                new_annotation.append([prev, cnt, {"s": start, "e": annotation[i-1][2]["e"]}])
            prev = c 
            cnt = c_cnt
            start = annotation[i][2]["s"]
    i = len(annotation) - 1
    new_annotation.append([prev, cnt, {"s": start, "e": annotation[i-1][2]["e"]}])
    return new_annotation

def find_potential_hors(annotation, annotation_seq, min_hor_len, max_hor_len, hors):
    potential_hors = {}
    annotation_str = "_".join(annotation_seq)
    for i in range(len(annotation)):
        end_ind = i
        len_subseq, subseq = 0, []
        len_monomer_subseq = 0
        while end_ind < len(annotation) and len_subseq < max_hor_len and annotation[end_ind][0] != "?":
            len_subseq += annotation[end_ind][1]
            if annotation[end_ind][0] in hors:
                len_monomer_subseq += len(hors[annotation[end_ind][0]].split("_"))*annotation[end_ind][1]
            else:
                len_monomer_subseq += annotation[end_ind][1]
            subseq.append(annotation_seq[end_ind])
            end_ind += 1
            if len_monomer_subseq < max_hor_len:
                subseq_str = "_".join(subseq)
                if "?" in subseq_str:
                    print(subseq_str)
                    exit(-1)
                if subseq_str not in potential_hors:
                    annotation_new_lst = annotation_str.split(subseq_str)
                    annotation_new_str = annotation_str.replace(subseq_str, "X")
                    annotation_new_str = annotation_new_str.replace("X_X", "X")

                    new_set_size = len(annotation_new_str.split("_"))
                    potential_hors[subseq_str] = {"set_size": new_set_size, "cnt": len(annotation_new_lst) - 1}
    return potential_hors

def build_full_hor(new_hor, hors, name):
    new_hor_represent = []
    for c in new_hor.split("_"):
        cc, cnt = c.split("[")[0], int(c.split("[")[1][:-1])
        for j in range(cnt):
            if cc in hors:
                new_hor_represent.append(hors[cc])
            else:
                new_hor_represent.append(cc + "[1]")
    hors[name] = "_".join(new_hor_represent)
    return hors

def update_annotation(annotation, annotation_seq, new_hor, name):
    new_annotation = []
    k = len(new_hor.split("_"))
    i = 0
    while i < len(annotation):
        if i + k < len(annotation):
            cur_seq = "_".join(annotation_seq[i:i+k])
        else:
            cur_seq = ""
        if cur_seq == new_hor:
            new_annotation.append([name, 1, {"s": annotation[i][2]["s"], "e": annotation[i + k - 1][2]["e"]} ])
            i += k
        else:
            new_annotation.append(annotation[i])
            i += 1
    return new_annotation
    
def update_hor_alignment(seq, name, hor, m_hor):
    linear_hor = []
    len_hor = len(m_hor.split("_"))
    for m in hor.split("_"):
        c, cnt = m.split("[")[0], int(m.split("[")[1][:-1])
        for _ in range(cnt):
            linear_hor.append(c)
    linear_hor_str = "_".join(linear_hor)
    new_seq = []
    i = 0
    while i < len(seq):
        if i + len(linear_hor) < len(seq):
            cur_seq = "_".join([x["qid"] for x in seq[i: i + len(linear_hor)]])
        else:
            cur_seq = ""
        if cur_seq == linear_hor_str:
            idnt = 0
            for p in seq[i: i + len(linear_hor)]:
                idnt += p["idnt"]
            new_seq.append({"qid": name, "len": len_hor, "s": seq[i]["s"], "e": seq[i + len(linear_hor)]["e"], "idnt": idnt/len(linear_hor)})
            i += len(linear_hor)
        else:
            new_seq.append(seq[i])
            i += 1
    return new_seq

def build_hor_annotation(dec, max_hor_len_, monomers_mp, filename):
    min_idnt, min_cnt, min_weight, min_hor_len, max_hor_len, mult = 75, 5, 5, 2, max_hor_len_, 2
    annotation = []
    seq = []
    for i in range(len(dec)):
        if dec[i]["idnt"] < min_idnt:
            annotation.append(["?", 1, {"s": dec[i]["s"], "e": dec[i]["e"]}])
            seq.append({"qid": "?", "len": 1, "s": dec[i]["s"], "e": dec[i]["e"], "idnt": dec[i]["idnt"]})
        else:
            annotation.append([dec[i]["qid"], 1, {"s": dec[i]["s"], "e": dec[i]["e"]}])
            seq.append({"qid": dec[i]["qid"], "len": 1, "s": dec[i]["s"], "e": dec[i]["e"], "idnt": dec[i]["idnt"]})

    annotation = collapse_annotation(annotation)
    hors = {}
    hors_lst = []
    h_cnt = 0
    while True:
        annotation_seq = []
        for a in annotation:
            annotation_seq.append(a[0] + "[" + str(a[1]) + "]")
        if h_cnt == 0:
            print("_".join(annotation_seq).replace("[1]", ""))
        set_size = len(annotation)
        #print(set_size)
        potential_hors = find_potential_hors(annotation, annotation_seq, min_hor_len, max_hor_len, hors)
        potential_hors_lst = []
        for h in potential_hors:
            if potential_hors[h]["cnt"] > min_cnt:
                potential_hors_lst.append([h, potential_hors[h]])
        potential_hors_lst = sorted(potential_hors_lst, key=lambda x: x[1]["set_size"])
        if len(potential_hors_lst) == 0 or set_size - potential_hors_lst[0][1]["set_size"] < min_weight:
            break
        #print(potential_hors_lst[:20])
        h_cnt += 1
        name = "h" + str(h_cnt)
        hors= build_full_hor(potential_hors_lst[0][0], hors, name)
        hors_lst.append([name, potential_hors_lst[0][0]])
        print("\t".join([name, potential_hors_lst[0][0].replace("[1]",""), hors[name].replace("[1]",""), str(len(hors[name].split("_"))), str(potential_hors_lst[0][1]["cnt"]), \
                                                     str(set_size - potential_hors_lst[0][1]["set_size"]), str(potential_hors_lst[0][1]["set_size"]) ]), flush=True)
        annotation = update_annotation(annotation, annotation_seq, potential_hors_lst[0][0], name)
        annotation = collapse_annotation(annotation)

    annotation_seq = []
    for a in annotation:
        if a[1] > 1:
            annotation_seq.append(a[0] + "[" + str(a[1]) + "]")
        else:
            annotation_seq.append(a[0])
    print("_".join(annotation_seq))

    for n in hors_lst:
        seq = update_hor_alignment(seq, n[0], n[1], hors[n[0]])

    with open(filename, "w") as fout:
        for c in seq:
            if c["qid"] in monomers_mp:
                fout.write("\t".join([monomers_mp[c["qid"]], str(c["len"]), "{0:.2f}".format(c["idnt"]), str(c["s"]), str(c["e"])]) + "\n") 
            else:
                fout.write("\t".join([c["qid"], str(c["len"]), "{0:.2f}".format(c["idnt"]), str(c["s"]), str(c["e"])]) + "\n") 


def create_default_valid():
    valid_mp = {}
    hor = "ABCDEFGHIJLKMNOPQR"
    cnt = 0
    for i in range(len(hor) - 1):
        valid_mp[hor[i]] = {}
        for j in range(i + 1, len(hor)):
            valid_mp[hor[i]][hor[j]] = j - i
    return valid_mp

def build_name(j, i):
    if j < 26:
        mono = chr(ord("A") + j) + str(i)
    else:
        mono = "Z" + str(i) + str(j-26)
    return mono

def create_numbered_valid():
    valid_mp = {}
    for i in range(5):
        for j in range(30):
            mono = build_name(j, i)
            valid_mp[mono] = {}
            for k in range(j + 1, 30):
                mono2 = build_name(k, i)
                valid_mp[mono][mono2] = k - j
    return valid_mp

def create_cen7_valid():
    valid_mp = {}
    hors = ["ABCDEFG", "HIJLKMNOPQRSTUV"]
    for hor in hors:
        for i in range(len(hor) - 1):
            valid_mp[hor[i]] = {}
            for j in range(i + 1, len(hor)):
                valid_mp[hor[i]][hor[j]] = j - i
    return valid_mp

ref = load_fasta(sys.argv[1], "map")
monomers = load_fasta(sys.argv[2])
dec, monomers_mp =load_dec(sys.argv[3], monomers)
filename = sys.argv[4]
max_hor_len = int(sys.argv[5])
name = list(ref.keys())[0]

build_hor_annotation(dec[name], max_hor_len, monomers_mp, filename)
