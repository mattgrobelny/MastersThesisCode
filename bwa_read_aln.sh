#!/usr/bin/bash -x

mid=1
alignment_directory="/home/mgrobelny/Scripts/github/Thesis_code/454_bwa_aln/read_aln/"

##########################################################################################

contigs_454=('/home/mgrobelny/Scripts/github/Thesis_code/454_opt_newbler/Mid1/Mid1_FINAL_3/454AllContigs.fna'
'/home/mgrobelny/Scripts/github/Thesis_code/454_opt_newbler/Mid2/Mid2_FINAL_11/454AllContigs.fna'
'/home/mgrobelny/Scripts/github/Thesis_code/454_opt_newbler/Mid3/Mid_3_134/mid3_454AllContigs.fna'
'/home/mgrobelny/Scripts/github/Thesis_code/454_opt_newbler/Mid4/Mid_4_46/454AllContigs.fna'
'/home/mgrobelny/Scripts/github/Thesis_code/454_opt_newbler/Mid5/Mid_5_45/454AllContigs.fna'
'/home/mgrobelny/Scripts/github/Thesis_code/454_opt_newbler/Mid6/Mid2_FINAL_21/454AllContigs.fna')

read_list_454=('/home/mgrobelny/Scripts/github/Thesis_code/454_opt_spades/mid_MID1/mid_MID1_AT_TE_VF_QC.fastq'
'/home/mgrobelny/Scripts/github/Thesis_code/454_opt_spades/mid_MID2/mid_MID2_AT_TE_VF_QC.fastq'
'/home/mgrobelny/Scripts/github/Thesis_code/454_opt_spades/mid_MID3/mid_MID3_AT_TE_VF_QC.fastq'
'/home/mgrobelny/Scripts/github/Thesis_code/454_opt_spades/mid_MID4/mid_MID4_AT_TE_VF_QC.fastq'
'/home/mgrobelny/Scripts/github/Thesis_code/454_opt_spades/mid_MID5/mid_MID5_AT_TE_VF_QC.fastq'
'/home/mgrobelny/Scripts/github/Thesis_code/454_opt_spades/mid_MID6/mid_MID6_AT_TE_VF_QC.fastq')
########################################
# Check newbler assembly
########################################

# # Map reads to contigs
for i in {0..5}
do
  file="mid_bwa_aln_NEWBLER_"
  aln='_alignment.sai'
  sam='_alignment.sam'
  file_name_aln=$alignment_directory$file$mid$aln
  bwa index ${contigs_454[$i]}
  bwa aln -t 8 -B 15 ${contigs_454[$i]} ${read_list_454[$i]} > $file_name_aln

  file_name_sam=$alignment_directory$file$mid$sam
  bwa samse ${contigs_454[$i]} $file_name_aln ${read_list_454[$i]} > $file_name_sam
  let "mid++"
done

sam_list="mid_bwa_aln_NEWBLER_1_alignment.sam
mid_bwa_aln_NEWBLER_2_alignment.sam
mid_bwa_aln_NEWBLER_3_alignment.sam
mid_bwa_aln_NEWBLER_4_alignment.sam
mid_bwa_aln_NEWBLER_5_alignment.sam
mid_bwa_aln_NEWBLER_6_alignment.sam"

# pull out stats from aligments
stats_file='bwa_aln_stats_NEWBLER.tsv'
echo "" > $alignment_directory$stats_file
for sam_output in $sam_list
do
  echo "File: $sam_output" >> $alignment_directory$stats_file
  samtools stats $alignment_directory$sam_output| grep ^SN | cut -f 2- >> $alignment_directory$stats_file
  echo "#######################################################################" >> $alignment_directory$stats_file

done
##########################################################################################

########################################
# Check Spades assembly
########################################

contigs_Spades=('/home/mgrobelny/Scripts/github/Thesis_code/454_opt_spades/mid_MID1/Spades_output_trimmed_careful_mid_MID1/contigs.fasta'
'/home/mgrobelny/Scripts/github/Thesis_code/454_opt_spades/mid_MID2/Spades_output_trimmed_careful_mid_MID2/contigs.fasta'
'/home/mgrobelny/Scripts/github/Thesis_code/454_opt_spades/mid_MID3/Spades_output_trimmed_careful_mid_MID3/contigs.fasta'
'/home/mgrobelny/Scripts/github/Thesis_code/454_opt_spades/mid_MID4/Spades_output_trimmed_careful_mid_MID4/contigs.fasta'
'/home/mgrobelny/Scripts/github/Thesis_code/454_opt_spades/mid_MID5/Spades_output_trimmed_careful_mid_MID5/contigs.fasta'
'/home/mgrobelny/Scripts/github/Thesis_code/454_opt_spades/mid_MID6/Spades_output_trimmed_careful_mid_MID6/contigs.fasta')

mid=1

# # Map reads to contigs
for i in {0..5}
do
  file="mid_bwa_aln_Spades_"
  aln='_alignment.sai'
  sam='_alignment.sam'
  file_name_aln=$alignment_directory$file$mid$aln
  bwa index ${contigs_Spades[$i]}
  bwa aln -t 8 -B 15 ${contigs_Spades[$i]} ${read_list_454[$i]} > $file_name_aln

  file_name_sam=$alignment_directory$file$mid$sam
  bwa samse ${contigs_Spades[$i]} $file_name_aln ${read_list_454[$i]} > $file_name_sam
  let "mid++"
done

sam_list="mid_bwa_aln_Spades_1_alignment.sam
mid_bwa_aln_Spades_2_alignment.sam
mid_bwa_aln_Spades_3_alignment.sam
mid_bwa_aln_Spades_4_alignment.sam
mid_bwa_aln_Spades_5_alignment.sam
mid_bwa_aln_Spades_6_alignment.sam"

# pull out stats from aligments
stats_file='bwa_aln_stats_SPADES.tsv'
echo "" > $alignment_directory$stats_file
for sam_output in $sam_list
do
  echo "File: $sam_output" >> $alignment_directory$stats_file
  samtools stats $alignment_directory$sam_output| grep ^SN | cut -f 2- >> $alignment_directory$stats_file
  echo "#######################################################################" >> $alignment_directory$stats_file

done
