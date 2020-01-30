from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio import SearchIO
from Bio.SeqRecord import SeqRecord

import random
random.seed(123)


def load_fasta(filename):
    records = SeqIO.to_dict(SeqIO.parse(filename, "fasta"))
    return records

def save_fasta(filename, orfs):
    with open(filename + ".fasta", "w") as output_handle:
        SeqIO.write(orfs, output_handle, "fasta")

def load_lst(filename):
    reads = []
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            reads.append(ln.strip().split()[-1])
    return reads

def load_chrX_alns(filename, n):
    reads = []
    total_len = 0
    reads_set = {}
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            q_name, q_len, q_start, q_end = ln.strip().split("\t")[:4]
            matches, total_len2 = map(int, ln.strip().split("\t")[9: 11])
            if int(q_end) - int(q_start) > 0.9*int(q_len) and int(q_len) > 20000 and matches/total_len2 > 0.8:
                reads_set[q_name] = int(q_len)
                total_len += int(q_len)
                #print(q_name)
    print(total_len)
    aligned_len = 0
    false_read_names = random.sample(list(reads_set.keys()), n)
    for name in false_read_names:
        aligned_len += reads_set[name]
    print(aligned_len)
    for name in false_read_names:
        print(name)
    return reads

cenX_names = load_lst("/Sid/tdvorkina/monomers/new_monomers/params/cenX_reads.lst")
cenX_reads = load_fasta("/Sid/tdvorkina/monomers/chrX/centromeric_reads.fasta")

num = 1000
true_read_names = random.sample(cenX_names, num)
true_reads = []
true_reads_len = 0
for name in true_read_names:
    name_prefix = name[:10]
    for r in cenX_reads:
        if r.startswith(name_prefix):
            new_name = r
            break
    if len(cenX_reads[new_name]) > 10000:
        true_reads.append(cenX_reads[new_name])
        true_reads_len += len(cenX_reads[new_name])

save_fasta("/Sid/tdvorkina/monomers/new_monomers/params/ont_true_reads", true_reads)
print("CenX reads length ", true_reads_len)

chrX_reads = load_fasta("/Sid/tdvorkina/monomers/new_monomers/rel3.fasta")
chrX_names = load_chrX_alns("/Sid/tdvorkina/monomers/new_monomers/output_rel3.paf", num)
print(len(chrX_names))


