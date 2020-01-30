from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio import SearchIO
from Bio.SeqRecord import SeqRecord


import re
import edlib

IDENTITY = 85

EDTHRESHOLD = 3

def edist(lst):
    if len(str(lst[0])) == 0:
        return -1, ""
    if len(str(lst[1])) == 0:
        return -1, ""
    result = edlib.align(str(lst[0]), str(lst[1]), mode="NW", task="path")
    return result["editDistance"], result["cigar"]

def edist_hw(lst):
    if len(str(lst[0])) == 0:
        return -1, ""
    if len(str(lst[1])) == 0:
        return -1, ""
    result = edlib.align(str(lst[0]), str(lst[1]), mode="HW")
    return result["editDistance"]

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

def load_fasta(filename):
    records = list(SeqIO.parse(filename, "fasta"))
    for i in range(len(records)):
        records[i] = records[i].upper()
    return records

def save_unique(monomer):
    unique_m = []
    for i in range(len(monomers)):
        add = True
        for j in range(i):
            if monomers[i].seq == monomers[j].seq:
                add = False
                break
        if add:
            unique_m.append(monomers[i])

    save_fasta("/Sid/tdvorkina/monomers/MigaKH.HigherOrderRptMon_unique", unique_m)

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

def divide_into_clusters(monomers, th):
    clusters = []
    parent = [i for i in range(len(monomers))]
    rank = [0 for _ in range(len(monomers))]
    for i in range(len(monomers)):
        for j in range(i + 1, len(monomers)):
            #idnt = aai([monomers[i].seq, monomers[j].seq])
            ed, _ = edist([monomers[i].seq, monomers[j].seq])
            if ed <= th:
                union_sets(i, j, parent, rank)

    clusters_id = {}
    for i in range(len(monomers)):
        ind = find_set(i, parent)
        if ind not in clusters_id:
            clusters_id[ind] = []
        clusters_id[ind].append(monomers[i])

    for cl in clusters_id:
        clusters.append(clusters_id[cl])


    mn_all = -1
    mx_sz = 0
    cl_i = 0
    for it in clusters:
        mn = -1
        cl_i += 1
        mn_mx, mn_mx_it = 1000, None
        if len(it) > mx_sz:
            mx_sz = len(it)
        if len(it) > 1:
            print(cl_i)
            for i in range(len(it)):
                print(" ", it[i].name)
                mx = -1
                for j in range(i + 1, len(it)):
                    ed, _ = edist([it[i].seq, it[j].seq])
                    if ed > mn:
                        mn = ed
                    if ed > mx:
                        mx = ed
                if mx < mn_mx:
                    mn_mx, mn_mx_it = mx, it[i].name
            #print(mn_mx_it)
            if mn > mn_all:
                mn_all = mn
        # else:
        #     print(it[0].name)

    # for i in range(len(monomers)):
    #     s = monomers[i].name
    #     cl_id = parent[i]
    #     for m in clusters_id[cl_id]:
    #         if m.name !=  monomers[i].name:
    #             s += " " + m.name
    #     print(s) 
    print("Threshold=", th, len(clusters), mn_all, mx_sz)

    return clusters


monomers = load_fasta("/Sid/tdvorkina/monomers/MigaKH.HigherOrderRptMon.fa")
print(len(monomers))
for th in range(6, 7):
    divide_into_clusters(monomers, th)
