#!/usr/bin/bash -x
# requires seq_crumbs
# requires SPades

files="mid_MID1.sff"
# mid_MID2.sff
# mid_MID3.sff
# mid_MID4.sff
# mid_MID5.sff
# mid_MID6.sff"


# location of inpurt mid files
mid_file_dir='/home/mgrobelny/Scripts/github/Thesis_code/454_raw_reads/'

# location of outputs
output_dir='/home/mgrobelny/Scripts/github/Thesis_code/454_opt_spades/'

# name of stats collection file
file_out_text='mid1_6_cleaning_data.txt'



time=$(date)
echo $time > $output_dir$file_out_text
echo "Cleaning of:" >> $output_dir$file_out_text
echo $files >> $output_dir$file_out_text

for file in $files
do
	echo "	"
	working_file=$mid_file_dir$file
	echo "working on $file"

	# extract file name string
	fastq=".fastq"

	filename=${file::-4}


	# extract sff to fastq
	echo "extracting sff"
	echo "	" >> $output_dir$file_out_text
	sff_extract $working_file > $output_dir$filename$fastq


	# trim adapters
	echo "Trimming adapters..."
	adapter_trim="_AT"
	filename_out1=$output_dir$filename$adapter_trim
	adapter_array=($(tagcleaner -predict -fastq $output_dir$filename$fastq | sed -E "1d; s/tag.\t([ATGCN]+)\t.+\t.+/\1\n/"))
	tagcleaner -tag3 ${adapter_array[0]} -tag5 ${adapter_array[1]} -out_format 3 -out $filename_out1 -fastq $output_dir$filename$fastq

	num_nuc_end2=$(sed -n 2~4p $filename$adapter_trim$fastq | tr -d '\n' | wc -m)
	num_nuc_start2=$(sed -n 2~4p $filename$fastq | tr -d '\n' | wc -m)
	delta2=$(bc <<< "scale = 2; (1-($num_nuc_end2/$num_nuc_start2))*100")
	echo "Tags used for $file:" >> $output_dir$file_out_text
	echo "-tag3 ${adapter_array[0]}" >> $output_dir$file_out_text
	echo "-tag5 ${adapter_array[1]}" >> $output_dir$file_out_text

	echo "Percent of nucleotides trimmed as tags/adapters $delta2 %" >> $output_dir$file_out_text
	echo "	" >> $output_dir$file_out_text

	# trim edges
	left_clip=30
	right_clip=20
	trim_edge="_TE"
	echo "Clipping edges..."
	echo "Left clip: $left_clip"
	echo "Right clip: $right_clip"
	filename_out2=$output_dir$filename$adapter_trim$trim_edge$fastq
	trim_edges -l $left_clip -r $right_clip -o $filename_out2 $filename_out1$fastq

	num_nuc_start3=$(sed -n 2~4p $filename_out1$fastq | tr -d '\n' | wc -m)
	num_nuc_end3=$(sed -n 2~4p $filename_out2 | tr -d '\n' | wc -m)
	delta3=$(bc <<< "scale = 2; (1-($num_nuc_end3/$num_nuc_start3))*100")
	echo "Percent of nucleotides trimmed using edge clipping: ~ $delta3 %" >> $output_dir$file_out_text
	echo "	" >> $output_dir$file_out_text

	# filter by blast
	echo "filtering by blast "
	vector_trim="_VF"
	filtered_out_file="_filtered_out_seqs.fastq"
	filter_id=98
	echo "Filtering with $filter_id % id" >> $output_dir$file_out_text
	filter_by_blast -b vectors_454.fasta -e $output_dir$filename$filtered_out_file -s $filter_id -o $output_dir$filename$adapter_trim$trim_edge$vector_trim$fastq $filename_out2

	# Write how much was filtered
	echo "Filtering by blast stats for: $filename.sff" >> $output_dir$file_out_text
	start_num=$(grep -c "^@I" $filename_out2)
	echo "Starting num of reads: $start_num" >> $output_dir$file_out_text
	after_filter=$(grep -c "^@I" $filename$adapter_trim$trim_edge$vector_trim$fastq)
	echo "After filtering: $after_filter" >> $output_dir$file_out_text
	delta=$(bc <<< "scale = 2; $start_num-$after_filter")
	echo "Num of reads filtered out: $delta" >> $output_dir$file_out_text
	echo "	" >> $output_dir$file_out_text

	# trim by quality
	echo "Trimming by quality..."
	quality_thresh=15
	quality_clip="_QC"
	echo "Quality threshold: $quality_thresh" >> $output_dir$file_out_text
	filename_out3=$output_dir$filename$adapter_trim$trim_edge$vector_trim$quality_clip$fastq
	trim_quality -q $quality_thresh -o $filename_out3 $output_dir$filename$adapter_trim$trim_edge$vector_trim$fastq

	num_nuc_start4=$(sed -n 2~4p $output_dir$filename$adapter_trim$trim_edge$vector_trim$fastq | tr -d '\n' | wc -m)
	num_nuc_end4=$(sed -n 2~4p $filename_out3 | tr -d '\n' | wc -m)
	delta4=$(bc <<< "scale = 2; (1-($num_nuc_end4/$num_nuc_start4))*100")
	echo "Percent of nucleotides trimmed using quality clipping: $delta4 %" >> $output_dir$file_out_text

	# assembly with Spades
	echo "assembling with Spades"
	out_file="/home/mgrobelny/Scripts/github/Thesis_code/454_opt_spades/Spades_output_trimmed_careful_"
	if (("$filename" == "mid_MID6"));
	then
	spades.py --trusted-contigs ~/Desktop/Spades_run/trust.fasta -m 10 -t 16 --careful --s1 $filename_out3 -o $out_file$filename

	else
	spades.py -m 10 -t 16 --careful --s1 $filename_out3 -o $out_file$filename
	fi

	# Pull out assemby info
	#cd $out_file$filename
	contigs_file="/home/mgrobelny/Data/mids_all/cleaned_mids/Spades_output_trimmed_careful_$filename/contigs.fasta"
	contig_stats=$(cat $contigs_file |grep ">" | sed -E 's/>NODE_([0-9]+)_length_([0-9]+)_cov_([0-9]+)/\1\t\2\t\3/')
	num_contigs=$(cat $contigs_file |grep ">" | sed -E 's/>NODE_([0-9]+)_length_([0-9]+)_cov_([0-9]+)/\1\t\2\t\3/'| cut -f 1 | tail -n 1)
	largest_contig=$( cat $contigs_file |grep ">" | sed -E 's/>NODE_([0-9]+)_length_([0-9]+)_cov_([0-9]+)/\1\t\2\t\3/'| cut -f 2 | head -n 1)

	#write contig stats to file
	echo "	" >> $output_dir$file_out_text
	echo "Contig stats for $filename_out3 assembly" >> $output_dir$file_out_text
	echo "Total contigs: $num_contigs" >> $output_dir$file_out_text
	echo "Largest contig: $largest_contig" >> $output_dir$file_out_text
	echo "	" >> $output_dir$file_out_text
	echo "#-------------------------------------------------------------------------------#" >> $output_dir$file_out_text

	assembly_data="_Asm_Data.tsv"
	echo "Contig_num	Contig_len	Contig_cov" > $output_dir$filename$adapter_trim$trim_edge$vector_trim$quality_clip$assembly_data
	echo "$contig_stats" >> $output_dir$filename$adapter_trim$trim_edge$vector_trim$quality_clip$assembly_data
	echo "done with $file"
done


# # Kmer genome estimation using jellyfish and ge
#
# #TYPICAL SESSION
#       # count k-mers (see jellyfish documentation for options) gzip -dc
# reads1.fastq.gz reads2.fastq.gz | jellyfish count -m 31 -o fastq.counts -C -s 10000000000 -U 500 -t 30 /dev/fd/0
#
#       # generate a histogram
# jellyfish histo fastq.counts_0 > fastq.counts_0.histo
#
#       # generate a pdf graph of the histogram
# jellyplot.pl fastq.counts_0.histo
#
#       # look at fastq.counts_0.histo.pdf and identify the approximate peak
#
#       # use find_valleys.pl to help pinpoint the actual peak
# find_valleys.pl fastq.counts_0.histo
#
#       # estimate the size and coverage
# estimate_genome_size.pl --kmer=31 --peak=42 --fastq=reads1.fastq.gz reads2.fastq.gz
