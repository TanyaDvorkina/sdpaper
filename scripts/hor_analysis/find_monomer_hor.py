import os
from os import listdir
from os.path import isfile, isdir, join

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio import SearchIO
from Bio.SeqRecord import SeqRecord

from joblib import Parallel, delayed


import edlib
import re

def load_enum(filename1, filename2):
    records = list(SeqIO.parse(filename1, "fasta"))
    records_new = list(SeqIO.parse(filename2, "fasta"))
    res_ini = {}
    for m in records:
        m_name = "_".join(m.name.split("_")[1:])
        m_num = m.name.split("_")[0]
        res_ini[m_name] = m.name

    res = {}
    for m in records_new:
        m_name = "_".join(m.name.split("_")[1:])
        m_num = m.name.split("_")[0]
        res[res_ini[m_name]] = m_num    
    return res

def simplify(s, monomers):
    return monomers[s]

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

def edist_hw(lst):
    if len(str(lst[0])) == 0:
        return -1, ""
    if len(str(lst[1])) == 0:
        return -1, ""
    result = edlib.align(str(lst[0]), str(lst[1]), mode="HW", k=500)
    return result["editDistance"]

def edist_nw(lst):
    if len(str(lst[0])) == 0:
        return -1, ""
    if len(str(lst[1])) == 0:
        return -1, ""
    result = edlib.align(str(lst[0]), str(lst[1]), mode="NW", k=500)
    return result["editDistance"]

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
    total_length = 0 #max(len(p1), len(p2))
    n = 0
    for c in cigar:
        if c.isdigit():
            n = n*10 + int(c)
        else:
            total_length += n
            n = 0
    matches = re.findall(r'\d+=', cigar)
    aai = 0.0
    for m in matches:
        aai += int(m[:-1])
    aai /= total_length
    return aai*100

def save_fasta(filename, orfs):
    with open(filename + ".fasta", "w") as output_handle:
        SeqIO.write(orfs, output_handle, "fasta")

def make_record(seq, name, sid, d=""):
    return SeqRecord(seq, id=sid, name=name, description = d)


def load_string_decomposition(filename, monomers, monomers_real, reads):
    reads_mapping = {}
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            if len(ln.strip().split("\t")) > 4:
                sseqid, qseqid, sstart, send, idnt  = ln.strip().split("\t")[:5]
                s, e, idnt = int(sstart), int(send), float(idnt)
                rev = False
                if qseqid.endswith("'"):
                    rev = True
                    qseqid = qseqid[:-1]
                qseqid = qseqid.split()[0]
                sseqid = sseqid.split()[0]
                if sseqid not in reads_mapping:
                    reads_mapping[sseqid] = []
                if rev:
                    ed = edist_nw([monomers_real[qseqid].seq, reads[sseqid][s:e].reverse_complement().seq])
                    idnt = aai([monomers_real[qseqid].seq, reads[sseqid][s:e].reverse_complement().seq])
                    reads_mapping[sseqid].append({"qid": simplify(qseqid, monomers), \
                                                 "s": s, "e": e, "rev": rev, "ed": ed, "idnt": idnt, "seq": reads[sseqid][s:e].reverse_complement().seq})
                else:    
                    ed = edist_nw([monomers_real[qseqid].seq, reads[sseqid][s:e].seq])
                    idnt = aai([monomers_real[qseqid].seq, reads[sseqid][s:e].seq])
                    reads_mapping[sseqid].append({"qid": simplify(qseqid, monomers), \
                                                 "s": s, "e": e, "rev": rev, "ed": ed, "idnt": idnt, "seq": reads[sseqid][s:e].seq })
    for r in reads_mapping:
        reads_mapping[r] = sorted(reads_mapping[r], key=lambda x: (x["s"], -x["e"]))
    return reads_mapping


monomers_enum = load_enum("/Sid/tdvorkina/monomers/new_monomers/monomer_clustering/cluter_represent_ed6_enumarated.fasta", \
                                            "/Sid/tdvorkina/monomers/new_monomers/monomer_clustering/cluter_represent_ed6_enumarated_order.fasta")
hors = load_fasta("/Sid/tdvorkina/monomers/new_monomers/known_hors/known_hors_doubled.fasta", "map")
monomers = load_fasta("/Sid/tdvorkina/monomers/new_monomers/monomer_clustering/cluter_represent_ed6_enumarated.fasta", "map")

dec = load_string_decomposition("/Sid/tdvorkina/monomers/new_monomers/known_hors/decomposition_known_hors_doubled.tsv", monomers_enum, monomers, hors)
print(len(dec))
valid = []
for r in sorted(list(dec.keys()), key=lambda x: x[1:].split("Z")[0]):
    cur_monomer_lst, cur_monomer_ed, cur_monomer_idnt = [], [], []
    cur_valid_lst = []
    for it in dec[r]:
        cur_monomer_lst.append("{:^3}".format(it["qid"]))
        cur_monomer_ed.append("{:^3}".format(it["ed"]))
        cur_monomer_idnt.append("{:^3.0f}".format(it["idnt"]))
        if it["ed"] < 26:
            cur_valid_lst.append(it["qid"])
        else:
            cur_valid_lst.append("?")
    cur_valid_lst = cur_valid_lst[1:len(cur_valid_lst)//2 + 1]
    cur_monomer_lst = cur_monomer_lst[1:len(cur_monomer_lst)//2 + 1]
    cur_monomer_idnt = cur_monomer_idnt[1:len(cur_monomer_idnt)//2 + 1]
    cur_monomer_ed = cur_monomer_ed[1:len(cur_monomer_ed)//2 + 1]
    print(r)
    valid.append(" ".join([r, "_" + "_".join(cur_valid_lst) + "_"]))
    print("\t".join(cur_monomer_lst))
    print("\t".join(cur_monomer_idnt))
    #print("\t".join(cur_monomer_ed))
    print(cur_valid_lst)  
    print("")  

print("\n".join(valid))