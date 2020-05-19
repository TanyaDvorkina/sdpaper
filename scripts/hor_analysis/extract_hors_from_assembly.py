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

def make_record(seq, name, sid, d=""):
    return SeqRecord(seq, id=sid, name=name, description = d)

def load_dec(filename, monomers_seq, ref):
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

    new_ref = {}
    new_reads_mapping = {}
    for r in reads_mapping:
        num_rev = 0
        for mp in reads_mapping[r]:
            num_rev += int(mp["rev"])
        if num_rev*2 > len(reads_mapping[r]):
            print(r)
            name, start, end = ref[r].id.split("_")
            start, end = int(start), int(end)
            new_reads_mapping[r + "'"] = []
            for mp in reads_mapping[r]:
                new_reads_mapping[r + "'"].append({"qid": mp["qid"], "s":end - mp["e"], "e": end - mp["s"], "rev": not mp["rev"], "idnt": mp["idnt"]})
            new_reads_mapping[r + "'"] = new_reads_mapping[r + "'"][::-1]
            new_ref[r + "'"] = make_record(ref[r].seq.reverse_complement(), ref[r].id + "'", ref[r].name + "'")
        else:
            new_reads_mapping[r] = reads_mapping[r][:]
            new_ref[r] = ref[r]

    reads_mapping = {}
    for r in new_reads_mapping:
        reads_mapping[r] = new_reads_mapping[r][:]
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
    cnt = 0
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
            cnt += right -left+ 1
    print(cnt)
    return new_ref, new_reads_mapping, monomers_mp_r

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

def find_potential_hors(annotation, annotation_seq, min_hor_len, max_hor_len, hors, potential_hors_all, set_size):
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
    for h in potential_hors_all:
        if h not in potential_hors:
            potential_hors_all[h]["set_size"] += len(annotation)
        else:
            potential_hors_all[h]["set_size"] += potential_hors[h]["set_size"]
            potential_hors_all[h]["cnt"] += potential_hors[h]["cnt"]

    for h in potential_hors:
        if h not in potential_hors_all:
            potential_hors_all[h] = {}
            potential_hors_all[h]["set_size"] = potential_hors[h]["set_size"] + set_size
            potential_hors_all[h]["cnt"] = potential_hors[h]["cnt"]
    return potential_hors_all

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
            new_seq.append({"qid": name, "len": len_hor, "s": seq[i]["s"], "e": seq[i + len(linear_hor) - 1]["e"], "idnt": idnt/len(linear_hor)})
            i += len(linear_hor)
        else:
            new_seq.append(seq[i])
            i += 1
    return new_seq

def known_hors_annotation(reads_annotation, known_hors, hors, hors_lst, h_cnt, min_cnt):
    hors_log = []
    for kh in known_hors:
        cnt = 0
        set_size, new_set_size = 0, 0
        for r in reads_annotation:
            annotation = reads_annotation[r]
            annotation_seq = []
            for a in annotation:
                annotation_seq.append(a[0] + "[" + str(a[1]) + "]")
            annotation_str = "_".join(annotation_seq)
            set_size += len(annotation_seq)
            annotation_new_lst = annotation_str.split(kh)
            annotation_new_str = annotation_str.replace(kh, "X")
            annotation_new_str = annotation_new_str.replace("X_X", "X")
            new_set_size += len(annotation_new_str.split("_"))
            cnt += len(annotation_new_lst) - 1
            #print(kh, new_set_size, cnt)
        if cnt == 0:
            continue
        h_cnt += 1
        name = "h" + str(h_cnt)
        hors= build_full_hor(kh, hors, name)
        hors_lst.append([name, kh])
        hors_log.append([name, kh.replace("[1]",""), hors[name].replace("[1]",""), str(len(hors[name].split("_"))), str(cnt), \
                                                     str(set_size - new_set_size), str(new_set_size) ])
        print("\t".join([name, kh.replace("[1]",""), hors[name].replace("[1]",""), str(len(hors[name].split("_"))), str(cnt), \
                                                     str(set_size - new_set_size), str(new_set_size) ]), flush=True)
        for r in reads_annotation:
            annotation = reads_annotation[r]
            annotation_seq = []
            for a in annotation:
                annotation_seq.append(a[0] + "[" + str(a[1]) + "]")
            annotation = update_annotation(annotation, annotation_seq, kh, name)
            reads_annotation[r] = collapse_annotation(annotation)

    return reads_annotation, hors, hors_lst, h_cnt, hors_log

def cluster_hors(hors, hors_lst, monomers_mp, known_hors):
    clustered_hors = {}
    processed = []
    for h in hors_lst:
        m_h = set(hors[h[0]].replace("[1]", "").split("_"))
        related_hors = []
        related_monomers = []
        for ph in processed:
            m_ph = set(hors[ph].replace("[1]", "").split("_") )
            if len(m_ph & m_h) > 0:
                related_hors.append(ph)
                related_monomers.append(m_ph)
                m_h = m_h - m_ph
        print(h[0], hors[h[0]], related_hors, related_monomers)
        if len(related_hors) == 0:
            clustered_hors[h[0]] = monomers_mp[hors[h[0]].replace("[1]", "").split("_")[0]].split(".")[0]
            print(clustered_hors[h[0]])
            processed.append(h[0])
            continue
        for j in range(len(known_hors)):
            m_ph = set(known_hors[j].replace("[1]", "").split("_") )
            if len(m_ph & m_h) > 0:
                related_hors.append("kh" + str(j))
                related_monomers.append(m_ph)
                m_h = m_h - m_ph
        substrs = []
        for m in hors[h[0]].replace("[1]", "").split("_"):
            rc = m[-1] == "'"
            if len(substrs) != 0 and (not rc and int(substrs[-1][-1][1:]) + 1 == int(m[1:]) or rc and  int(substrs[-1][-1][1:-1]) + 1 == int(m[1:-1])):
                add = True
                for ms in related_monomers:
                    if m in ms and substrs[-1][-1] not in ms:
                        add = False
                        break
                if add:
                    substrs[-1].append(m)
                else:
                    substrs.append([m])
            else:
                substrs.append([m])
        # for s in substrs:
        #     print("_".join(s))
        name = []
        cur_hor = -1
        for s in substrs:
            prev_hor = cur_hor
            for i in range(len(related_hors)):
                if s[0] in related_monomers[i]:
                    cur_hor = i
                    break
            rc = s[0][-1] == "'"
            monomer_name = monomers_mp[s[0][:-1]] if rc else monomers_mp[s[0]]
            #print(s[0], cur_hor, related_monomers)
            cur_hor_name = monomer_name.split(".")[0]
            if prev_hor != cur_hor:
                name.append(cur_hor_name)
            is_full = False
            if (not related_hors[i].startswith("k")) and hors[related_hors[i]].replace("[1]", "") == "_".join(s):
                is_full = True
            if is_full:
                name.append("F")
            else:
                if rc:
                    first, last = monomer_name.split(".")[1] + "'", monomers_mp[s[-1][:-1]].split(".")[1] + "'"
                else:
                    first, last = monomer_name.split(".")[1], monomers_mp[s[-1]].split(".")[1]
                if first == last:
                    name.append(first)
                else:
                    name.append(first + "-" + last)
        print("_".join(name))
        clustered_hors[h[0]] = "_".join(name)
    return clustered_hors



def build_hor_annotation(reads_dec, max_hor_len_, monomers_mp, filename, known_hors = []):
    min_idnt, min_cnt, min_weight, min_hor_len, max_hor_len, mult = 75, 5, 5, 2, max_hor_len_, 2
    annotation = {}
    seq = {}
    for read in reads_dec:
        dec = reads_dec[read]
        seq[read] = []
        annotation[read] = []
        for i in range(len(dec)):
            if dec[i]["idnt"] < min_idnt:
                annotation[read].append(["?", 1, {"s": dec[i]["s"], "e": dec[i]["e"]}])
                seq[read].append({"qid": "?", "len": 1, "s": dec[i]["s"], "e": dec[i]["e"], "idnt": dec[i]["idnt"]})
            else:
                annotation[read].append([dec[i]["qid"], 1, {"s": dec[i]["s"], "e": dec[i]["e"]}])
                seq[read].append({"qid": dec[i]["qid"], "len": 1, "s": dec[i]["s"], "e": dec[i]["e"], "idnt": dec[i]["idnt"]})

    hors = {}
    hors_lst = []
    h_cnt = 0
    for r in annotation:
        annotation[r] = collapse_annotation(annotation[r])

    hors_log = []
    if len(known_hors) > 0:
        annotation, hors, hors_lst, h_cnt, hors_log = known_hors_annotation(annotation, known_hors, hors, hors_lst, h_cnt, min_cnt)

    while True:
        potential_hors = {}
        set_size = 0
        for r in annotation:
            annotation_seq = []
            for a in annotation[r]:
                annotation_seq.append(a[0] + "[" + str(a[1]) + "]")
            if h_cnt == 0:
                print("_".join(annotation_seq).replace("[1]", ""))
            #print(set_size)
            potential_hors = find_potential_hors(annotation[r], annotation_seq, min_hor_len, max_hor_len, hors, potential_hors, set_size)
            set_size += len(annotation[r])

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
        hors_log.append([name, potential_hors_lst[0][0].replace("[1]",""), hors[name].replace("[1]",""), str(len(hors[name].split("_"))), str(potential_hors_lst[0][1]["cnt"]), \
                                                     str(set_size - potential_hors_lst[0][1]["set_size"]), str(potential_hors_lst[0][1]["set_size"]) ])
        print("\t".join(hors_log[-1]), flush=True)
        for r in annotation:
            annotation_seq = []
            for a in annotation[r]:
                annotation_seq.append(a[0] + "[" + str(a[1]) + "]")
            annotation[r] = update_annotation(annotation[r], annotation_seq, potential_hors_lst[0][0], name)
            annotation[r] = collapse_annotation(annotation[r])

    for r in annotation:
        annotation_seq = []
        for a in annotation[r]:
            if a[1] > 1:
                annotation_seq.append(a[0] + "[" + str(a[1]) + "]")
            else:
                annotation_seq.append(a[0])
        print(r)
        print("_".join(annotation_seq))

    for n in hors_lst:
        for r in seq:
            seq[r] = update_hor_alignment(seq[r], n[0], n[1], hors[n[0]])

    clustered_hors = cluster_hors(hors, hors_lst, monomers_mp, known_hors)
    for h in hors_log:
        h = [clustered_hors[h[0]]] + h
        print("\t".join(h), flush=True)
    prev = 0
    with open(filename, "w") as fout:
        for r in seq:
            for c in seq[r]:
                if c["qid"] in monomers_mp:
                    fout.write("\t".join([r, monomers_mp[c["qid"]], str(c["len"]), "{0:.2f}".format(c["idnt"]), str(c["s"]), str(c["e"]), str(c["e"] - c["s"] + 1), str(c["s"] - prev)]) + "\n")
                else:
                    fout.write("\t".join([r, clustered_hors[c["qid"]], str(c["len"]), "{0:.2f}".format(c["idnt"]), str(c["s"]), str(c["e"]), str(c["e"] - c["s"] + 1),  str(c["s"] - prev)]) + "\n")
                prev = c["e"]

def build_known_hors(lst):
    known_hors = []
    for c in lst:
        s, e = c[0], c[1]
        kh = []
        kh_rc = []
        for i in range(s, e + 1):
            kh.append("m" + str(i) + "[1]")
            kh_rc.append("m" + str(i) + "'[1]")
        known_hors.append("_".join(kh[::-1]))
        known_hors.append("_".join(kh))
        known_hors.append("_".join(kh_rc))
        known_hors.append("_".join(kh_rc[::-1]))
    return known_hors

ref = load_fasta(sys.argv[1], "map")
monomers = load_fasta(sys.argv[2])
ref, dec, monomers_mp =load_dec(sys.argv[3], monomers, ref)
filename = sys.argv[4]
max_hor_len = int(sys.argv[5])
name = list(ref.keys())[0]

lst2 = [[1,4], [5,14], [15, 24]]
lst3 = [[1, 17], [55, 64]]
lst7 = [[1, 6], [7, 22]]
lst8 = [[1, 11]]
lst10 = [[35, 42], [53, 60], [184, 201]]
lst12 = [[1, 8], [17, 34]]
lst16 = [[1,10], [11, 22], [23, 33]]
lst19 = [[1, 6], [7, 23], [24, 36], [37, 51], [52, 83], [84, 99], [5, 6]]
lst20 = [[1, 16], [17, 24], [25, 35], [36, 46], [47, 54]]

lst1 = [[1,6], [7,17], [18,28],[29, 37],[38, 46], [47, 50], [51, 53], [54,55]]
lst4 = [[1,19], [20, 32]]
lst5 = [[1,6], [7, 22], [23, 25], [26, 27], [28, 39],[40, 49], [50, 62], [63, 78], [79, 93], [94, 111], [112, 132]]
lst6 = [[1,18]]
lst9 = [[1,7]]
lst11 = [[1,5], [6, 17], [18, 24], [25, 40]]
lst17 = [[1,14], [15, 28], [29, 44], [45, 54]]
lst18 = [[1,12], [13, 22], [23, 32], [33, 42], [43, 47]]
lstX = [[1,12]]
known_hors = build_known_hors(lst18)
print(known_hors)
build_hor_annotation(dec, max_hor_len, monomers_mp, filename, known_hors)

