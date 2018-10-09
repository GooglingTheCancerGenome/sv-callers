# This script removes a set of sequences from a FASTA file

# Firse argument: FASTA file including sequences to remove
# Second argument: list of FASTA IDs to remove

from Bio import SeqIO
import sys
import gzip

header_set = set(line.strip() for line in open(sys.argv[2]))
print('Sequences to remove: %s' % header_set)

with gzip.open(sys.argv[1], "rt") as gzffile:
    ffile = SeqIO.parse(gzffile, "fasta")
    for seq_record in ffile:
        if seq_record.name not in header_set:
            print(seq_record.format("fasta"))
