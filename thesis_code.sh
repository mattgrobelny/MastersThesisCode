#
# Thesis Code Master Run Script
#
# Matt Grobelny

################################################################################
# Assembly of 454 reads with Newbler and SPades

# If true run 454 newbler assembly optimization (Newbler)
calculate_opi_454_asm=0

# Run spades assember with 454 data (SPades)
calculate_opi_spades_asm=1

################################################################################
# Assembly of Sanger reads with PhredPhrap

# Run Sanger data assembly (PhredPhrap)
run_phredPhrap=0

################################################################################
# Check how reads map back to assembly with BWA for the best Newbler assembly
# and SPADES assemblies for each mid

# Run ALL reads alignment
run_BWA_aln=0

################################################################################
# Run all upstream seq extract

#

################################################################################
if (($calculate_opi_454_asm == "1"))
	then
	# 454 Data - Optimize 454 assembly
	mid_file_dir="/home/mgrobelny/Scripts/github/Thesis_code/454_raw_reads/"

	reads_454_file_list="mid_MID1.sff
	mid_MID2.sff
	mid_MID3.sff
	mid_MID4.sff
	mid_MID5.sff
	mid_MID6.sff"

	output="_output"
	for file in $reads_454_file_list
	do
	wk_output=${file::-4}$output
	wk_file=$mid_file_dir$file
	echo "perl runAssembly_PS.pl 15 45 5 15 50 5 95 99 1 $wk_output ../bothtrimfiles.fasta ../$wk_file TRUE TRUE"

	done

################################################################################
elif (($calculate_opi_spades_asm == "1"))
	then
		/home/mgrobelny/Scripts/github/Thesis_code/Spades_asm.sh

# elif (($run_phredPhrap == "1"))
# 	then
# 		# some dir for 7k and 10k sanger data so phred and phrap can run
#

################################################################################
elif (($run_BWA_aln == "1"))
	then
		/home/mgrobelny/Scripts/github/Thesis_code/bwa_read_aln.sh

################################################################################

else
		echo "done"
fi
