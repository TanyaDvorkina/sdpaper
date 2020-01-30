import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

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


def load_fasta(filename):
    records = list(SeqIO.parse(filename, "fasta"))
    for i in range(len(records)):
        records[i] = records[i].upper()
    return records

flye_assembly = load_fasta("/Sid/tdvorkina/monomers/new_monomers/full_hg/assembly/flye.contigs.fasta")
lens = [len(r.seq) for r in flye_assembly]
print(len(lens), sum(lens))

exit(-1)
path = "/Sid/tdvorkina/monomers/new_monomers/cent_reads/"

datasets = ["SRR9087600", "SRR9087599", "SRR9087598", "SRR9087597"]

for d in datasets:
    print(d)
    reads = load_fasta(os.path.join(path, d, "ac_cent_reads.fasta"))
    lens = [len(r.seq) for r in reads]
    print(len(lens), sum(lens))
