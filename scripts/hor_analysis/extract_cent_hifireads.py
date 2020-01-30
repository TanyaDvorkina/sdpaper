import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from gzip import open as gzopen 

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio import SearchIO
from Bio.SeqRecord import SeqRecord


import numpy as np; np.random.seed(0)
import pandas as pd

def load_aligned_reads(filename):
    res = {}
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            q_name, q_len, q_start, q_end, q_or, t_name, t_len, t_start, t_end, match, aln_len = ln.split("\t")[:11]
            if int(q_end) - int(q_start) > 0.9*int(q_len): # and int(match)/int(aln_len) > 0.8:
                if q_name not in res:
                    res[q_name] = []
                res[q_name].append(t_name)
    new_res = set()
    for q in res:
        if len(res[q]) == 1:
            new_res.add(q)
    return new_res

def load_fasta(filename):
    records = list(SeqIO.parse(filename, "fasta"))
    for i in range(len(records)):
        records[i] = records[i].upper()
    return records

def save_fasta(filename, orfs):
    with open(filename + ".fasta", "w") as output_handle:
        SeqIO.write(orfs, output_handle, "fasta")

datasets = ["SRR9087600", "SRR9087599", "SRR9087598", "SRR9087597"]

for d in datasets:
    print(d)
    aligned_reads = load_aligned_reads("/Sid/tdvorkina/monomers/new_monomers/full_hg/assembly/flye_" + d + ".paf")
    print("Loaded")
    reads = load_fasta("/Sid/tdvorkina/monomers/new_monomers/" + d + "_1.fasta")
    print("Read")

    cent_reads = [r for r in reads if r.name not in aligned_reads]
    print(len(reads), len(cent_reads))
    print("Total reads number", len(reads))
    print("Total reads length", sum([len(r.seq) for r in reads]))
    print("Total reads number 10-12kbp", len([len(r.seq) for r in reads if len(r.seq) > 10000 and len(r.seq) < 12001]))
    print("Total reads length 10-12kbp", sum([len(r.seq) for r in reads if len(r.seq) > 10000 and len(r.seq) < 12001]))
    save_fasta("/Sid/tdvorkina/monomers/new_monomers/cent_reads/" + d + "/cent_reads", cent_reads)

    lens = [len(r.seq) for r in cent_reads]
    total_len = sum(lens)
    print("Total cent reads number", len(lens))
    print("Total cent reads length", total_len)

    mp_th = {10: 0, 11:0, 12: 0, 13: 0, 14: 0, 15:0}
    mp_th_len = {10: 0, 11:0, 12: 0, 13: 0, 14: 0, 15:0}
    for r in cent_reads:
        for k in mp_th:
            if len(r.seq) > k*1000:
                mp_th[k] += 1
                mp_th_len[k] += len(r.seq)

    for k in sorted(list(mp_th.keys())):
        print(k, mp_th[k], mp_th_len[k], mp_th_len[k]/total_len)

    plt.hist(lens, bins=50)
    plt.savefig("/Sid/tdvorkina/monomers/new_monomers/cent_reads/" + d + "/readlen_hist.png")