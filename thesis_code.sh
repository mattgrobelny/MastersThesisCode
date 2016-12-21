#
# Thesis Code
#
# Matt Grobelny

# If true run 454 newbler assembly optimization
calculate_opi_454_asm=0
ls

# Run spades assember with 454 data
calculate_opi_spades_asm=0

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
go
	done
elif (($calculate_opi_spades_asm == "1"))
	then
		/home/mgrobelny/Scripts/github/Thesis_code/Spades_asm.sh
else
		echo "done"
fi
