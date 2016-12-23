import sys
import getopt
import math
import numpy as np
import re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

coordinates_file =""
seq_fasta=""
fasta_name=""
output_file =""
argv = sys.argv[1:]
try:
    opts, args = getopt.getopt(argv, "hn:c:f:")
except getopt.GetoptError:
    print 'promoter_id.py -n <fasta_header_name>[input_fasta_file_name] -c <coordinates_file> -o <output_file> -f <input_fasta_file>'
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print "#--- Upstream Sequence Extract ---#\n"
        print "Usage:"
        print 'promoter_id.py -n <fasta_header_name>[input_fasta_file_name] -c <coordinates_file> -o <output> -f <input_fasta_file>\n'
        print "Goals:"
        print "Extract the the sequence upstream of genes according to cordinates file"
        print "Output fasta file of upstream seqs with correct frame 5' to 3'"
        print "\n"

        sys.exit()
    elif opt in ("-c"):
        coordinates_file = arg
    elif opt in ("-f"):
        seq_fasta = arg
    elif opt in ("-n"):
        fasta_name = arg
#print "Input Sequence file:", seq_fasta
#print "Upstream coordinates_file:", coordinates_file
fasta_name = seq_fasta.split("/")[-1][0:-6]
#print "Fasta_output_header:", fasta_name
seq=""
# Import seq from fasta file
for seq_record in SeqIO.parse(seq_fasta,'fasta'):
    seq=seq_record.seq

coordinates_file_fh = open(coordinates_file, 'r')

upstream_seq_len= re.findall('_([0-9]+).tsv',coordinates_file.split("/")[-1])

# skip headers
next(coordinates_file_fh)
next(coordinates_file_fh)
next(coordinates_file_fh)

gene_list_fasta = []
for line in coordinates_file_fh:
    #strip new line char
    line = line.strip('\n')
    line = line.replace(" ", "")
    # split tabs
    line_split = line.split('\t')

    if line_split[1] == "-":
        fasta_out = '>'+ fasta_name + "_" + line_split[0] + "_" + "Range:" + str(line_split[3])+"-"+str(line_split[2])+ "_Upstream_len:" + str(upstream_seq_len[0])
        print fasta_out
        print seq[int(line_split[3]):int(line_split[2])].reverse_complement()
    else:
        fasta_out = '>'+ fasta_name + "_" + line_split[0] + "_" + "Range:" + str(line_split[2])+"-"+str(line_split[3])+ "_Upstream_len:" + str(upstream_seq_len[0])
        print fasta_out
        print seq[int(line_split[2]):int(line_split[3])]
