import sys
import edlib

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio import SearchIO
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalwCommandline

import argparse
import os, shutil

def edist(lst):
    if len(str(lst[0])) == 0:
        return 100500
    if len(str(lst[1])) == 0:
        return 100500
    result = edlib.align(str(lst[0]), str(lst[1]), mode="NW")
    return result["editDistance"]

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

def save_fasta(filename, orfs):
    with open(filename, "w") as output_handle:
        SeqIO.write(orfs, output_handle, "fasta")


def load_nomono_regions(filename, reads, th_border = 97, th_mono = 70):
    reads_regions = {}
    with open(filename, "r") as fin:
        for ln in fin.readlines():
            sseqid, qseqid, sstart, send, idnt  = ln.strip().split("\t")[:5]
            sseqid = sseqid.split()[0]
            if sseqid not in reads_regions:
                reads_regions[sseqid] = []
            s, e, idnt = int(sstart), int(send), float(idnt)
            reads_regions[sseqid].append({"seq": reads[sseqid].seq[s: e + 1], "idnt": idnt, "s": s, "e": e})

    window = 5
    for r in reads_regions:
        l, e = -1, -1
        idnt_sum = 0
        for j in range(len(reads_regions[r])):
            if j - window >= 0:
                idnt_sum -= reads_regions[r][j - window]["idnt"]
            idnt_sum += reads_regions[r][j]["idnt"]
            if l == -1 and idnt_sum/window >= th_border:
                l = j
            if idnt_sum/window >= th_border:
                e = j
        if e -l > 20:
            #print(reads_regions[r][l]["s"], reads_regions[r][e]["e"])
            reads_regions[r] = reads_regions[r][l:e + 1]
        else:
            reads_regions[r] = [] 

    collapsed_nonmono_regions = {}
    for r in reads_regions:
        if len(reads_regions[r]) > 0:
            in_region = False
            collapsed_nonmono_regions[r] = []
            for i in range(len(reads_regions[r])):
                if reads_regions[r][i]["idnt"] < th_mono:
                    if in_region:
                        collapsed_nonmono_regions[r][-1].append(reads_regions[r][i])
                    else:
                        in_region = True
                        #print("new region")
                        collapsed_nonmono_regions[r].append([reads_regions[r][i]])
                    #print(reads_regions[r][i]["s"], reads_regions[r][i]["e"])
                else:
                    in_region = False
    new_collapsed_nonmono_regions = {}
    for r in collapsed_nonmono_regions:
        if len(collapsed_nonmono_regions[r]) > 0:
            new_collapsed_nonmono_regions[r] = []
            for region in collapsed_nonmono_regions[r]:
                new_collapsed_nonmono_regions[r].append({ "seq": "".join([str(x["seq"]) for x in region]), "s": region[0]["s"], "e": region[-1]["e"] })
    return new_collapsed_nonmono_regions

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

def collapse_homo(sequence):
    new_sequence = sequence[0]
    for i in range(1, len(sequence)):
        if sequence[i] != sequence[i-1]:
            new_sequence += sequence[i]
    return Seq(new_sequence)

def cluster_by_outer_ed(sequences, th=5):
    clusters = []
    parent = [i for i in range(len(sequences))]
    rank = [0 for _ in range(len(sequences))]
    for i in range(len(sequences)):
        for j in range(i + 1, len(sequences)):
            h_sequence_i, h_sequence_j = collapse_homo(sequences[i]), collapse_homo(sequences[j])
            ed = edist([h_sequence_i, h_sequence_j])
            min_ed = min(ed, edist([h_sequence_i.reverse_complement(), h_sequence_j]))
            if min_ed*100/min(len(h_sequence_j), len(h_sequence_i)) <= th:
                union_sets(i, j, parent, rank)

    clusters_id = {}
    for i in range(len(sequences)):
        ind = find_set(i, parent)
        ed = edist([sequences[i], sequences[ind]])
        min_ed = min(ed, edist([sequences[i].reverse_complement(), sequences[ind]]))
        if ind not in clusters_id:
            clusters_id[ind] = []
        if min_ed != ed:
            clusters_id[ind].append(sequences[i].reverse_complement())
        else:
            clusters_id[ind].append(sequences[i])

    for cl in clusters_id:
        clusters.append(clusters_id[cl])
        print(len(clusters), len(clusters[-1]), len(clusters[-1][0]))
    return clusters

def save_clusters(clusters, prefix):
    names = []
    for i in range(len(clusters)):
        sequences = []
        for j in range(len(clusters[i])):
            sequences.append(make_record(clusters[i][j], "cl" + str(i) + "it" + str(j), "cl" + str(i) + "it" + str(j)))
        num = len(sequences)
        if len(sequences) == 1:
            sequences.append(sequences[0])
        save_fasta(prefix + "/cl" +  str(i) + "_" + str(num) + ".fasta", sequences)
        names.append(prefix + "/cl" +  str(i) + "_" + str(num) + ".fasta")
    return names

def construct_consensus(clusters):
    clustalw_exe = r"/home/tdvorkina/soft/clustalo-1.2.4-Ubuntu-x86_64"
    for cl in clusters:
        cline = ClustalwCommandline(clustalw_exe, infile=cl, outfile=cl[:-len("fasta")] + "clu")
        stdout, stderr = cline()

    seq_consensus = []
    for cl in clusters:
        filename = cl[:-len("fasta")] + "clu"
        total_alns = []
        with open(filename, "r") as fin:
            alns = []
            for ln in fin.readlines():
                ln = ln.replace("  ", " ")
                seq = "*"
                if len(ln.split()) >= 2:
                    seq = ln.strip().split()[1]
                if seq[0] in {"A","C","G", "T", "-"}:
                    alns.append(seq)
                elif len(alns) > 0:
                    if len(total_alns) == 0 or len(total_alns) == len(alns):
                        if len(total_alns) == 0:
                            total_alns = ["" for _ in range(len(alns))]
                        for i in range(len(alns)):
                            total_alns[i] += alns[i]
                        alns = []
                    else:
                        print("Something went wrong")
                        exit(-1)
        s_consensus = ""
        for i in range(len(total_alns[0])):
            score = {"A": 0, "C": 0, "G": 0, "T": 0, "-": 0}
            for j in range(len(total_alns)):
                score[total_alns[j][i]] += 1
            max_score = "A"
            for s in score:
                if score[s] > score[max_score]:
                    max_score = s
            if max_score != "-":
                s_consensus += max_score
        #print(cl, s_consensus)
        name = cl.split("/")[-1].split(".")[0]
        seq_consensus.append(make_record(Seq(s_consensus), name + "_" + str(len(s_consensus)), name + "_" + str(len(s_consensus))))
    return seq_consensus

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Searches for non-monomeric regions')
    parser.add_argument('sequences', help='fasta-file with long reads or genomic sequences')
    parser.add_argument('decomposition', help='tsv-file with sequences decomposition')
    parser.add_argument('output', help='output directory')

    args = parser.parse_args()

    prefix = args.output
    if os.path.isdir(prefix):
        shutil.rmtree(prefix)
    os.mkdir(prefix)
    reads = load_fasta(args.sequences, "map")
    regions = load_nomono_regions(args.decomposition, reads)

    cnt = 0
    seqs = []
    for r in regions:
        for region in regions[r]:
            if region["e"] -  region["s"] > 200:
                print(r, len(regions[r]), region["s"],  region["e"], region["e"] -  region["s"])
                #print(region["seq"])
                seqs.append(Seq(region["seq"]))
                cnt += 1

    print("Found", cnt, "potential non-monomeric regions in", len(regions), "reads")
    print("Clustering regions..")
    clusters = cluster_by_outer_ed(seqs)
    print("Saving clusters to", prefix)
    filenames = save_clusters(clusters, prefix)
    print("Constructing consensus..")
    consensus = construct_consensus(filenames)
    print("Saving potential non-monomeric regions to ", prefix + "/consensus_regions.fasta")
    consensus = sorted(consensus, key = lambda x: -int(x.id.split("_")[1]))
    save_fasta(prefix + "/consensus_regions.fasta", consensus)






