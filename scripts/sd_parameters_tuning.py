#!/usr/bin/env python3
import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

import os

import numpy as np; np.random.seed(0)
import seaborn as sns; sns.set()
import pandas as pd

from sklearn import preprocessing
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report
from sklearn.decomposition import PCA

import joblib


def convert_to_abc(name):
    res = name.split("_")[0]
    if "fake" in name:
        res = 'M'
    return res

def load_reads_positions():
    path = "/Sid/tdvorkina/monomers/sdpaper"
    tm_res = os.path.join(path, "data_extended/cenx-t2t7_alignment.bed")
    res = {}
    with open(tm_res, "r") as fin:
        for ln in fin.readlines():
            lst = ln.strip().split()
            name = lst[3]
            is_good = False
            ref_start, ref_end, read_start, read_end = [int(lst[1]), int(lst[2]), int(lst[4]), int(lst[5])]
            if read_start > read_end:
                res[name] = [read_end, read_start]
            else:
                res[name] = [read_start, read_end]
    return res

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

def replacer(r):
    chars = ["_", "=", "|", ":"]
    r = r.lower()
    for c in chars:
        r = r.replace(c, "-")
    r = r.replace("--", "-")
    return r

def filter_decomposition(reads_mapping, positions):
    res = {}
    for read in reads_mapping:
        r_read = replacer(read)
        if r_read in positions:
            l, r = find_ends_indexes(reads_mapping[read], positions[r_read][0], positions[r_read][1], read)
            res[read] = reads_mapping[read][l: r + 1]
        else:
            print("Hmm", read)
            exit()
    return res

def load_string_decomposition(filename, identity = 0):
    reads_mapping = {}
    reads_mapping_homo = {}
    filtered = 0
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            sseqid, qseqid, sstart, send, idnt, q2seqid, idnt2, qseqid_homo, idnt_homo, q2seqid_homo, idnt2_homo  = ln.strip().split("\t")[:11]
            qseqid = qseqid.split()[0]
            q2seqid = q2seqid.split()[0]
            qseqid_homo = qseqid_homo.split()[0]
            q2seqid_homo = q2seqid_homo.split()[0]
            if sseqid not in reads_mapping:
                    reads_mapping[sseqid] = []
                    reads_mapping_homo[sseqid] = []
            s, e, idnt, idnt2, idnt_homo, idnt2_homo = int(sstart), int(send), float(idnt), float(idnt2), float(idnt_homo), float(idnt2_homo)
            rev = False
            if qseqid.endswith("'"):
                rev = True
                qseqid = qseqid[:-1]
            if idnt > identity:
                reads_mapping[sseqid].append({"qid": convert_to_abc(qseqid), "qid2": convert_to_abc(q2seqid), \
                                              "s": s, "e": e, "rev": rev, "idnt": idnt, "idnt_diff": idnt - idnt2, \
                                              "qid_homo": convert_to_abc(qseqid_homo), "qid2_homo": convert_to_abc(q2seqid_homo), \
                                              "idnt_homo": idnt_homo, "idnt_diff_homo": idnt_homo - idnt2_homo})
            else:
                filtered += 1

    window = 2
    for r in reads_mapping:
        reads_mapping[r] = reads_mapping[r][1:-1]
        for i in range(len(reads_mapping[r])):
            avg_idnt = 0
            for j in range(max(0, i - window), min(i + window, len(reads_mapping[r]))):
                if i != j:
                    avg_idnt += reads_mapping[r][j]["idnt"]
            reads_mapping[r][i]["avg_idnt"] = avg_idnt/2*window
    return reads_mapping


def train_logreg(decomp, prefix):
    df_cenX = pd.DataFrame(decomp["cenX"])
    df_nocenX = pd.DataFrame(decomp["nocenX"])
    df_cenX["type"] = [1 for _ in range(len(df_cenX))]
    df_nocenX["type"] = [0 for _ in range(len(df_nocenX))]
    df = df_cenX.append(df_nocenX, ignore_index = True)
    X = pd.concat([df["idnt"], df["idnt_diff"]], axis=1, keys = ["idnt", "idnt_diff"])
    X_scaled = preprocessing.scale(X)
    y = df["type"]
    print("Fitting..")
    clf = LogisticRegression(random_state=0).fit(X_scaled, y)

    print("Saving..")
    joblib.dump(clf, prefix + "logreg_model.sav")

    print("Predicting..")
    y_pred = clf.predict(X_scaled)
    print(classification_report(y,y_pred))

    y_lst, y_pred_lst = list(y), list(y_pred)
    tp = len(y[(y == y_pred)])
    fp = len(y[(y == 0) & (y_pred == 1)])
    fn = len(y[(y == 1) & (y_pred == 0)])
    precision = tp/(tp + fp)
    recall = tp/(tp + fn)
    print(tp, fp, fn, precision, recall, 2*precision*recall/(precision + recall))

    fig = plt.figure()
    plt.plot(df[y == 1].idnt, df[y == 1].idnt_diff, 'b^', df[y == 0].idnt, df[y == 0].idnt_diff, 'rs', alpha=0.2)
    plt.savefig(prefix + "original_classes.png")

    fig = plt.figure()
    plt.plot(df[y_pred == 1].idnt, df[y_pred == 1].idnt_diff, 'b^', df[y_pred == 0].idnt, df[y_pred == 0].idnt_diff, 'rs', alpha=0.2)
    plt.savefig(prefix + "logreg_classes.png")

print("ONT")
cenX_decomp = load_string_decomposition("/Sid/tdvorkina/monomers/new_monomers/params/ont_true_reads_decomposition.tsv")
positions = load_reads_positions()
cenX_decomp = filter_decomposition(cenX_decomp, positions)
nocenX_decomp = load_string_decomposition("/Sid/tdvorkina/monomers/new_monomers/params/ont_false_reads_decomposition.tsv")
decomp = {"cenX":[], "nocenX": []}
for r in cenX_decomp:
    for c in cenX_decomp[r]:
        decomp["cenX"].append(c)
for r in nocenX_decomp:
    for c in nocenX_decomp[r]:
        decomp["nocenX"].append(c)

train_logreg(decomp, "ont_")

print("Hifi")
decomp = load_string_decomposition("/Sid/tdvorkina/monomers/new_monomers/params/hifi_seq_decomposition.tsv")
train_logreg(decomp, "hifi_")


