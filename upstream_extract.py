import sys
import getopt
import math
import numpy as np
import re

# First get all upstream seqeuneces from gff file
upstream_seq_len = 200 #bp
gff_file = ""
#
# Promoter Identification
#

## Three types of promoter sequences have been identified in eukaryotic DNA:
# - TATA box, the most common, is prevalent in rapidly transcribed genes.
# - Initiator promoters infrequently are found in some genes.
# - CpG islands are characteristic of transcribed genes.

## Promoter-proximal elements occur within 200 base pairs of the start site.
# Several such elements, containing up to 20 base pairs, may help regulate a particular gene.
# Enhancers, which are usually 100-200 base pairs in length, contain multiple 8 to 20bp control elements.
# They may be located from 200 base pairs to tens of kilobases upstream or downstream from a promoter, within an intron,
# or downstream from the final exon of a gene.
# Promoter-proximal elements and enhancers often are cell-type specific, functioning only in specific differentiated cell types.
# Source: https://www.ncbi.nlm.nih.gov/books/NBK21745/

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

gff_f_h = open(gff_file, 'r')

# skip headers
next(gff_f_h)
next(gff_f_h)

gene_list = []
gene_list_w_cord =[]
gene_list_w_cord.append(["Gene", "Sense","Start_cord","Upstream (" + str(upstream_seq_len) + ")"])
for line in gff_f_h:
    #strip new line char
    line = line.strip('\n')
    # remove spaces
    line= line.replace(" ", "_")
    # split tabs
    line_split = line.split('\t')
    if len(line_split) > 1:

        if line_split[2] == 'CDS':
            product = re.findall(';product=(.+);protein_id', line)
            if product not in gene_list:
                upstream_val =0
                run_me='upstream_val= int(line_split[3]) %s int(upstream_seq_len)' % line_split[6]
                exec(run_me)
                gene_list.append(product)
                gene_list_w_cord.append([product[0],line_split[6],line_split[3],upstream_val])
            else:
                continue
    else:
        continue
gff_f_h.close

# output tsv of gene list w cords to terminal 
for i in gene_list_w_cord:
    print i[0],'\t', i[1],'\t', i[2],'\t',i[3]
