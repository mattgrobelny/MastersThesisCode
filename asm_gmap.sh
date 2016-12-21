#
# Gmap alignment of assembled contigs to raw reads
#
# Requires:
# Gmap

# List of 454 newbler file locations
newbler_454_asm_contigs_file_list="/home/mgrobelny/Scripts/github/Thesis_code/454_opt_newbler/Mid1/Mid1_FINAL_3/,
/home/mgrobelny/Scripts/github/Thesis_code/454_opt_newbler/Mid2/Mid2_FINAL_11/,
/home/mgrobelny/Scripts/github/Thesis_code/454_opt_newbler/Mid3/Mid_3_134/,
/home/mgrobelny/Scripts/github/Thesis_code/454_opt_newbler/Mid4/Mid_4_46/:
/home/mgrobelny/Scripts/github/Thesis_code/454_opt_newbler/Mid5/Mid_5_45/
/home/mgrobelny/Scripts/github/Thesis_code/454_opt_newbler/Mid6/Mid2_FINAL_21/"

# newbler contig.fasta name
newbler_454_contig_file_name="454AllContigs.fna"

for contig_fasta in $newbler_454_asm_contigs_file_list
do

  directory="/home/mgrobelny/Scripts/github/Thesis_code/454_gmap"
  genome="/home/a-m/ib501_stud12/shell/Hw7/mmu/Mus_musculus.GRCm38.dna.toplevel.fa"
  database_dir="/home/a-m/ib501_stud12/shell/Hw7/mmu/Mus_musculus.GRCmgsnap38"
  gmap_build -D $directory -d Mus_musculus.GRCm38 $genome

  # alaign reads to reference
  reads="/home/classroom/ib501/alignment/mmu_700_tumor_cells.fq.gz"
  splicesites="Mus_musculus.GRCm38.splicesites.iit"
  output="/home/a-m/ib501_stud12/shell/Hw7/mmu/Mus_musculus.GRCm38.out"
  gsnapl -D $directory -d Mus_musculus.GRCm38 -t 8 -o $output -A sam --gunzip -s $splicesites $reads
