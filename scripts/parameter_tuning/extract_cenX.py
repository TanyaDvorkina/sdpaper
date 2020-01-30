from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from Bio import SearchIO
from Bio.SeqRecord import SeqRecord

import random
random.seed(123)

start, end = 57861200, 60968784 
            #57826352, 60935032
start_before = 57809997
end_after = 60945933
#
start_extended, end_extended = 57809997, 60945933
#2773652-2779726 - LINE

def load_fasta(filename):
    records = list(SeqIO.parse(filename, "fasta"))
    return records

def save_fasta(filename, orfs):
    with open(filename + ".fasta", "w") as output_handle:
        SeqIO.write(orfs, output_handle, "fasta")

def make_record(seq, name, sid, d=""):
    return SeqRecord(seq, id=sid, name=name, description = d)


chrX = load_fasta("/Sid/tdvorkina/monomers/new_monomers/chm13.chrX_v0.7.fasta")
noCenX_1 = chrX[0].seq[:start_extended]
noCenX_2 = chrX[0].seq[end_extended:]
lst = [make_record(noCenX_1, "chrX1", "chrX1"), make_record(noCenX_2, "chrX2", "chrX2")]
save_fasta("/Sid/tdvorkina/monomers/new_monomers/chm13.chrX_v0.7_noCenX_new", lst)

# sz = 1000000
# pos = random.randrange(0, start - sz - 10000)
# false_subseq = chrX[0].seq[pos: pos + sz]
# pos = random.randrange(start + 10000, start + 2500000 - sz)
# true_subseq = chrX[0].seq[pos: pos + sz]

# save_fasta("/Sid/tdvorkina/monomers/new_monomers/params/hifi_true_seq", make_record(true_subseq, "cenX", "cenX"))
# save_fasta("/Sid/tdvorkina/monomers/new_monomers/params/hifi_false_seq", make_record(false_subseq, "nocenX", "nocenX"))