import sys
import getopt
import math
import numpy as np
import re

from Bio import SeqIO

argv = sys.argv[1:]
try:
    opts, args = getopt.getopt(argv, "hs:f:")
except getopt.GetoptError:
    print 'promoter_id.py -s <upstream_seq_len>[200] -f <input_gff_file>'
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print "#--- Upstream region extract ---#\n"
        print "Usage:"
        print 'promoter_id.py -s <upstream_seq_len>[200] -f <input_gff_file> \n'
        print "Goals:"
        print "Extract the coordinates upstream of genes distance -s"
        print "\n"

        sys.exit()
    elif opt in ("-s"):
        upstream_seq_len = arg
    elif opt in ("-f"):
        gff_file = arg
print "Input file:", gff_file
print "Upstream sequence length:", upstream_seq_len

fasta_sequences = SeqIO.parse(open(input_file),'fasta')


outfile =
with open(output_file) as out_file:
    for fasta in fasta_sequences:
        name, sequence = fasta.id, fasta.seq.tostring()
        new_sequence = some_function(sequence)
        write_fasta(out_file)
