#!/usr/bin/bash -x

# best asssembled cotigs
# contigs_454="/home/mgrobelny/Scripts/github/Thesis_code/454_opt_newbler/Mid1/Mid1_FINAL_3/454AllContigs.fna
# /home/mgrobelny/Scripts/github/Thesis_code/454_opt_newbler/Mid2/Mid2_FINAL_11/454AllContigs.fna
# /home/mgrobelny/Scripts/github/Thesis_code/454_opt_newbler/Mid3/Mid_3_134/mid3_454AllContigs.fna
# /home/mgrobelny/Scripts/github/Thesis_code/454_opt_newbler/Mid4/Mid_4_46/454AllContigs.fna
# /home/mgrobelny/Scripts/github/Thesis_code/454_opt_newbler/Mid5/Mid_5_45/454AllContigs.fna
# /home/mgrobelny/Scripts/github/Thesis_code/454_opt_newbler/Mid6/Mid2_FINAL_21/454AllContigs.fna"

db_location="/home/mgrobelny/Scripts/github/Thesis_code/454_gmap/"
mid=1

# # Make dbs for each assembly
# for contig_file in $contigs_454
# do
#   db="mid"
#   db_name=$db$mid
#   db_out=$db_location
#   mkdir $db_out
#   gmap_build -D $db_out -d $db_name  $contig_file
#   let "mid++"
# done

mid_db_list=('mid1'
'mid2'
'mid3'
'mid4'
'mid5'
'mid6')

alignment_directory="/home/mgrobelny/Scripts/github/Thesis_code/454_gmap/read_aln/"


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

# # Map reads to contigs
# for i in {0..5}
# do
#   gmap -d $db_location${mid_db_list[i]} -t 8 -o $alignment_directory${mid_db_list[i]} -A sam ${read_list_454[i]}
# done


# # Make dbs for each assembly
# for contig_file in $contigs_454
# do
#   db="mid"
#   db_name=$db$mid
#   db_out=$db_location
#   #mkdir $db_out
#   bwa index $contig_file
# done
# #

# # Map reads to contigs
for i in {0..5}
do
  file="mid_bwa_aln_"
  aln='_alignment.sai'
  sam='_alignment.sam'
  file_name_aln=$alignment_directory$file$mid$aln
  bwa index ${contigs_454[$i]}
  bwa aln -t 8 -B 15 ${contigs_454[$i]} ${read_list_454[$i]} > $file_name_aln

  file_name_sam=$alignment_directory$file$mid$sam
  bwa samse ${contigs_454[$i]} $file_name_aln ${read_list_454[$i]} > $file_name_sam
  let "mid++"
done

sam_list="mid_bwa_aln_1_alignment.sam
mid_bwa_aln_2_alignment.sam
mid_bwa_aln_3_alignment.sam
mid_bwa_aln_4_alignment.sam
mid_bwa_aln_5_alignment.sam
mid_bwa_aln_6_alignment.sam"

stats_file='gsnap_aln_stats.tsv'
echo "" > $alignment_directory$stats_file
for sam_output in $sam_list
do
  echo "File: $sam_output" >> $alignment_directory$stats_file
  samtools stats $alignment_directory$sam_output| grep ^SN | cut -f 2- >> $alignment_directory$stats_file
  echo "#######################################################################" >> $alignment_directory$stats_file

done
