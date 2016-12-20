#!/usr/bin/perl
#use strict;
use warnings;
use Getopt::Long;

#Trigger time count
BEGIN { our $start_run = time(); }

#required parameters:
my $help=$ARGV[0];

print "Hi, here is what you gave me to use:\n @ARGV\n\n";

if (( $help eq "--help")||($help eq "-h")){
  die print "NEWBLER PARAMETERS SWEEP v1.0 \n
	Looks like you need some help... \n
	\n
Requires:\n
	GS denovo Assembler\n
	.sff file\n
	vectortrimfiles.fasta\n
	\n
  Internal Defaults for runAssembly:\n
  -ss 1  --> seed step parameter \n\n

	-sl 10 --> seed length parameter 6 to 15\n\n

	-sc 10 --> seed count parameter\n\n

  -force --> will overwrite files of the same name\n\n

  -cpu 0 --> use all cpus\n\n

  -nobig --> do not output big files (ace, bam ...)\n\n

To run this program you MUST write parameter in this order: \n
	\n
	-readlength_minimum -readlength_maximum -readlength_step -min_minoverlap -max_minoverlap -step_minoverlap \n
	same line continued... \n
	-min_min_id -max_min_id -step_min_id -Projectname(will also be foldername) -vectortrimfiles -sff_file) \n
  same line continued... \n
  -PrintGraphs (TRUE or FALSE) -PrintRecommendedParameters (TRUE or FALSE)\n

##### Example:\n
	perl runAssembly_PS.pl 15 45 5 15 50 5 95 99 1 Project1 ../bothtrimfiles.fasta ../mid_MID1.sff TRUE TRUE\n

Some parameters where hard coded in the version... sorry if you dont see the parameter you want to iterate throught...\n

See README file for more info\n
	"};

# Done with Intro
##################################################################################################################################################

# Set Default Parameters

my $min_min_readlength= '15'; #minimum min read length, minimum min is 15
my $max_min_readlength= '45'; #maximum min read length
my $step_min_readlength='5'; #step size for min read length

##########################SEED PARAMETERS####################################
#my $min_seed= @ARGV[4]; #minimum seed length , minimum min is 15
#my $max_seed= @ARGV[5]; #maximum min read length
#my $step_seed=@ARGV[6]; #step size for min read length

####################MIN OVERLAP LENGHT PARAMETERS#############################
my $min_minoverlap= '15'; #minimum min over lap length, minimum min is 15
my $max_minoverlap= '50'; #maximum min over lap length
my $step_minoverlap='5'; #step size for interation through minlength range

####################MIN OVERLAP LENGHT PARAMETERS#############################
my $min_min_id= '95'; #minimum min id
my $max_min_id= '99'; #maximum min id
my $step_min_id='1'; #step size for interation through minimum identity


####################GOBAL ITERATION PARAMETERS#############################
#setting up global variables to pass to system
my $qsub = ('qsub -I -q classroom');
my $module_load= ('module load 454');
my $runAssembly= ('runAssembly');


#setting global variables for loops
my $readlength   =$min_min_readlength;
my $overlap      =$min_minoverlap;
my $id           =$min_min_id;
my $folder_name  =''; #add empty string to make string as type
my $folder_it =0;
my $folder_step=1;
my $runAssembly_paras=0;
my $vectorTrimfiles='';
my $sff_file='';
my $printGraphs= 'TRUE' ;
my $printrecommendedParameters= 'TRUE';

#Get options
GetOptions (

  #Optional options
  "readlength_Min:i" => \$min_min_readlength,
  "readlength_Max:i" => \$max_min_readlength,
  "readlength_Step:i" => \$step_min_readlength,
  "minoverlap_Min:i" => \$min_minoverlap,
  "minoverlap_Max:i" => \$max_minoverlap,
  "minoverlap_Step:i" => \$step_minoverlap,
  "min_id_Min:i" => \$min_min_id,
  "min_id_Max:i" => \$max_min_id,
  "min_id_Step:i" => \$step_min_id,
  "printGraphs:s" => \$printGraphs,
  "printrecommendedParameters:s" => \$printrecommendedParameters,

  #Required options
  "projectName=s" => \$folder_name,
  "vectorTrimfiles=s" => \$vectorTrimfiles,
  "sff_file=s" => \$sff_file);

$folder_name="" . $folder_name;
# to implemet make sure to change all ARGV


#########################
#prep working directory
system("mkdir $folder_name");
chdir "./$folder_name/";

my $current_dir= `pwd`;
chomp $current_dir;

#prep final output file seprate each data point by comma
my $score_file_name= $current_dir. "/" . $folder_name . "_newbler_scores.txt";
open(my $fh, '>', "$score_file_name");
print {$fh} "Folder_name,Readlength,Overlap,Id,numberOfContigs,numberOfBases,avgContigSize,N50ContigSize,largestContigSize\n";
close ($fh);

# Accessory code
#########################################################################
#print("running $qsub \n");
#system($qsub);

#print("running $module_load \n");
#system($module_load);

#Subroutines
##################################################################################################################################################
sub run_Assembly {
#put together intial argumets to pass to runAssembly
my @paravals_internal = @_;

#default arguments to run
my @args_def= ('-o',"$paravals_internal[0]",
  '-force',
  #'-m',
  '-cpu','0',
  '-nobig',
  '-vs',$vectorTrimfiles,
  '-vt',$vectorTrimfiles,  # hard code.....
  $sff_file
  ); # hard code.....

# relabel out from para_combo_gen to paravals


#set up parameters to pass into system
my @para_args=(
	'-minlen',"$paravals_internal[1]",
	#Minimum length of reads to use in assembly Default: 50 Min: 15

	'-ss','1',
	#Set seed step parameter must be >=1

	'-sl','10',
	#Set seed length parameter 6 to 15

	'-sc','10 ',
	#Set seed count parameter must be >=1

	'-ml',"$paravals_internal[2]",
	#Set minimum overlap must be >=1

	'-mi',"$paravals_internal[3]",
	#Set minimum overlap identify 0 to 99
);
print"running $runAssembly with:\n
@para_args and @args_def \n\n";
system($runAssembly,@para_args,@args_def);

#convert array of parameters to string
#@string_para_args= joing(",",@para_args);

#output string paras and folder name
@runAssembly_paras=($paravals_internal[0],@para_args);

return @runAssembly_paras;
}
##########################################################
# Collect quality scores,  takes in a newblermetric files and outputs a
# 5 metric array
### numberOfContigs   = ######;
### numberOfBases     = #########;

### avgContigSize     = ####;
### N50ContigSize     = ####;
### largestContigSize = #####;

sub col_scores{

#initilize collection vals
my $numberOfContigs =0;
my $numberOfBases   =0;
my $avgContigSize   =0;
my $N50ContigSize   =0;
my $largestContigSize =0;
my $file_name="454NewblerMetrics.txt";

#orient directory
chdir "$current_dir/$folder_name_for_it/";
my @scores=("","","","","");

#whole metrics file variable
my $whole_file=0;
#open 454NewblerMetrics.txt
open (my $FILE,'<', $file_name ) or return @scores;

#import complete metrics files for processing
while (my $line = <$FILE>){
$whole_file=$whole_file . $line;
}
close ($FILE);

$whole_file=~ /numberOfContigs   = (.+)\;/ ;
$numberOfContigs  =$1;
$whole_file=~ /numberOfBases     = (.+)\;/ ;
$numberOfBases    =$1;
$whole_file=~ /avgContigSize     = (.+)\;/ ;
$avgContigSize    =$1;
$whole_file=~ /N50ContigSize     = (.+)\;/ ;
$N50ContigSize    =$1;
$whole_file=~ /largestContigSize = (.+)\;/ ;
$largestContigSize=$1;

chdir "$current_dir";

@scores = ($numberOfContigs,$numberOfBases,$avgContigSize,$N50ContigSize,$largestContigSize);
return @scores;
}

#Done with subs
##################################################################################################################################################

#Forloop
##################################################################################################################################################
#For loop for iterating through each parameter
for ($readlength=$max_min_readlength; $readlength>=$min_min_readlength;$readlength=$readlength-$step_min_readlength) {

	for ($overlap=$max_minoverlap;$overlap>=$min_minoverlap;$overlap=$overlap-$step_minoverlap) {

		for ($id=$max_min_id; $id>=$min_min_id ; $id=$id-$step_min_id) {

  		#set up folder names
  			$folder_name_for_it =$folder_name ."_" . $folder_it;
  			print"Folder name: $folder_name_for_it\n";

  		#organize iterated parametes
    		my @paravals=($folder_name_for_it,$readlength,$overlap,$id);
    		print"Current parameters outside of runAssembly: @paravals\n\n";

  	 	#run Newbler with paramters --> ouputs @runAssembly_paras
  		  run_Assembly(@paravals);

  		    my @scores =col_scores($folder_name_for_it);
          push @paravals, @scores;

  		#split into comma seperated text
  		  my $Final_output=0;
  		    $Final_output=join(",", @paravals);

  		      print"Here are the results for this run: @paravals\n";

  		#save masterlist data to file
  		  open(my $fh, '>>', "$score_file_name");

  		    print {$fh} "$Final_output\n";

  		      close ($fh);

  		        print"DONE with $folder_name_for_it\n \n starting new run...\n\n";
  		#  add 1 to $folder_it interation \n";
      $folder_it= $folder_it+1
    }
  }
}

#Done with loop
print "Done with loop\n";
##################################################################################################################################################

#Data Analysis
##################################################################################################################################################
chdir "$current_dir";


if (( $printGraphs eq "TRUE")||($printrecommendedParameters  eq "TRUE")){
  system("Rscript ../score_analysis.R $printGraphs $printrecommendedParameters ./ $folder_name")
}

##################################################################################################################################################
## print out time at the end of script
my $end_run = time();
my $run_time = $end_run - our $start_run;
die print "Job took $run_time seconds\n";
