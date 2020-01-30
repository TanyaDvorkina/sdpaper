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

from random import sample
import random
random.seed(123)

import joblib


TRAIN_SZ = 0.8

def convert_to_abc(name):
    res = name.split("_")[0]
    if "fake" in name:
        res = 'M'
    return res

def load_reads_positions():
    path = "/Sid/tdvorkina/monomers/sd_ismb_submission"
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

def load_string_decomposition(filename, train, identity = 0):
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
    train_set = set(sample(sorted(list(reads_mapping.keys())), int(train*len(reads_mapping))))
    window = 2
    for r in reads_mapping:
        reads_mapping[r] = reads_mapping[r][1:-1]
        for i in range(len(reads_mapping[r])):
            avg_idnt = 0
            for j in range(max(0, i - window), min(i + window, len(reads_mapping[r]))):
                if i != j:
                    avg_idnt += reads_mapping[r][j]["idnt"]
            reads_mapping[r][i]["avg_idnt"] = avg_idnt/2*window
            if r in train_set:
                reads_mapping[r][i]["train"] = 1
            else:
                reads_mapping[r][i]["train"] = 0
    return reads_mapping


def train_logreg(decomp, prefix):
    df_cenX = pd.DataFrame(decomp["cenX"])
    df_nocenX = pd.DataFrame(decomp["nocenX"])
    df_cenX["type"] = [1 for _ in range(len(df_cenX))]
    df_nocenX["type"] = [0 for _ in range(len(df_nocenX))]
    df = df_cenX.append(df_nocenX, ignore_index = True)
    df_train = df[df.train == 1]
    df_test = df[df.train == 0]
    X = pd.concat([df_train["idnt"], df_train["idnt_diff"]], axis=1, keys = ["idnt", "idnt_diff"]) 
    y = df_train["type"]
    print(len(y))
    print("Fitting..")
    clf = LogisticRegression(random_state=0).fit(X, y)

    print("Saving..")
    joblib.dump(clf, prefix + "logreg_model.sav")

    print("Predicting..")
    X_test = pd.concat([df_test["idnt"], df_test["idnt_diff"]], axis=1, keys = ["idnt", "idnt_diff"]) 
    y_pred = clf.predict(X_test)
    y_test = df_test["type"]
    print(len(y_pred), len(y_test))
    print(classification_report(y_test,y_pred))

    tp = len(y_test[(y_test == y_pred) & (y_test == 1)])
    tn = len(y_test[(y_test == y_pred) & (y_test == 0)])
    fp = len(y_test[(y_test == 0) & (y_pred == 1)])
    fn = len(y_test[(y_test == 1) & (y_pred == 0)])
    precision = tp/(tp + fp)
    recall = tp/(tp + fn)
    print(tp, tn, fp, fn, precision, recall, 2*precision*recall/(precision + recall))
    print("FPR", fp/(fp + tn), "FNR", fn/(fn +tp))

    fig = plt.figure()
    plt.plot(df_test[y_test == 1].idnt, df_test[y_test == 1].idnt_diff, 'b^', df_test[y_test == 0].idnt, df_test[y_test == 0].idnt_diff, 'rs', alpha=0.2)
    plt.savefig(prefix + "original_classes.png")

    fig = plt.figure()
    plt.plot(df_test[y_pred == 1].idnt, df_test[y_pred == 1].idnt_diff, 'b^', df_test[y_pred == 0].idnt, df_test[y_pred == 0].idnt_diff, 'rs', alpha=0.2)
    plt.savefig(prefix + "logreg_classes.png")

def train_logreg_1d(decomp, prefix):
    df_cenX = pd.DataFrame(decomp["cenX"])
    df_nocenX = pd.DataFrame(decomp["nocenX"])
    df_cenX["type"] = [1 for _ in range(len(df_cenX))]
    df_nocenX["type"] = [0 for _ in range(len(df_nocenX))]
    df = df_cenX.append(df_nocenX, ignore_index = True)
    df_train = df[df.train == 1]
    df_test = df[df.train == 0]
    X = pd.concat([df_train["idnt"]], axis=1, keys = ["idnt"]) 
    y = df_train["type"]
    print("Fitting..")
    clf = LogisticRegression(random_state=0).fit(X, y)

    print("Saving..")
    joblib.dump(clf, prefix + "logreg_model1d.sav")

    print("Predicting..")
    X_test = pd.concat([df_test["idnt"]], axis=1, keys = ["idnt"]) 
    y_pred = clf.predict(X_test)
    y_test = df_test["type"]
    print(len(y_pred), len(y_test))
    print(classification_report(y_test,y_pred))

    tp = len(y_test[(y_test == y_pred) & (y_test == 1)])
    tn = len(y_test[(y_test == y_pred) & (y_test == 0)])
    fp = len(y_test[(y_test == 0) & (y_pred == 1)])
    fn = len(y_test[(y_test == 1) & (y_pred == 0)])
    precision = tp/(tp + fp)
    recall = tp/(tp + fn)
    print(tp, tn, fp, fn, precision, recall, 2*precision*recall/(precision + recall))
    print("FPR", fp/(fp + tn), "FNR", fn/(fn +tp))

    fig = plt.figure()
    plt.plot(df_test[y_test == 1].idnt, df_test[y_test == 1].idnt_diff, 'b^', df_test[y_test == 0].idnt, df_test[y_test == 0].idnt_diff, 'rs', alpha=0.2)
    plt.savefig(prefix + "original_classes.png")

    fig = plt.figure()
    plt.plot(df_test[y_pred == 1].idnt, df_test[y_pred == 1].idnt_diff, 'b^', df_test[y_pred == 0].idnt, df_test[y_pred == 0].idnt_diff, 'rs', alpha=0.2)
    plt.savefig(prefix + "logreg_classes.png")

print("ONT")
cenX_decomp = load_string_decomposition("/Sid/tdvorkina/monomers/new_monomers/params/ont_true_reads_decomposition.tsv", train = TRAIN_SZ)
positions = load_reads_positions()
cenX_decomp = filter_decomposition(cenX_decomp, positions)
nocenX_decomp = load_string_decomposition("/Sid/tdvorkina/monomers/new_monomers/params/ont_false_reads_decomposition.tsv", train = TRAIN_SZ)
decomp = {"cenX":[], "nocenX": []}
for r in cenX_decomp:
    for c in cenX_decomp[r]:
        decomp["cenX"].append(c)
for r in nocenX_decomp:
    for c in nocenX_decomp[r]:
        decomp["nocenX"].append(c)

train_logreg(decomp, "ont_")
train_logreg_1d(decomp, "ont_")

print("SimONT")
cenX_decomp = load_string_decomposition("/Sid/tdvorkina/monomers/sd_ismb_submission/data_extended/chrX_simulated/cenX/decomposition_cenX_aligned_reads.tsv", train = TRAIN_SZ)
nocenX_decomp = load_string_decomposition("/Sid/tdvorkina/monomers/sd_ismb_submission/data_extended/chrX_simulated/noCenX/decomposition_noCenX_aligned_reads.tsv", train = TRAIN_SZ)
decomp = {"cenX":[], "nocenX": []}
for r in cenX_decomp:
    for c in cenX_decomp[r]:
        decomp["cenX"].append(c)
for r in nocenX_decomp:
    for c in nocenX_decomp[r]:
        decomp["nocenX"].append(c)

train_logreg(decomp, "ont_")
train_logreg_1d(decomp, "ont_")

# print("Hifi")
# decomp = load_string_decomposition("/Sid/tdvorkina/monomers/new_monomers/params/hifi_seq_decomposition.tsv")
# train_logreg(decomp, "hifi_")


