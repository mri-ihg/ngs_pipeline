#!/usr/bin/perl

use strict;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);
umask(0002);

#include Utilities.pm
my $prog_path = dirname(abs_path($0));
require $prog_path."/Utilities.pm";

my $command		  = "";
my $outprefix     = "";
my $help          = 0;
my $settings      = "hg19";
my $infile1       = "";
my $infile2       = "";
my $annotationfile = "";
my $threads       = -1;
my $logfile  	  = "";
my $loglevel 	  = "INFO";
my $doDelete	  = 0;
my @files2Delete;
my $starProgramPath = "";
my $starFusionProgramPath = "";
my $starIndex = "";
my $starlogfile = "";
my $doSingleEnd = 0;
my $readGroupSample = "";
my $readGroupLib    = "Lib1";
my $platform        = "Illumina";
my $bam = 0;
my $readGroup = "";
my $tmpArgument ="";

my $params = Utilities::getParams();

my $helptext      = 
"
-f	<infile reads 1>	
-F	<infile reads 2>	
-o	<output prefix>
-b 					if infile is in bam-format
-l	<library name>		the name of the library (default: $readGroupLib)
-s	<read group sample>	read group sample name
--lf	<log file>		the log file for the pipeline logging (default: pipeline.log)
--ll	<log level>		the log level for the pipeline logging; available options: ERROR, INFO, DEBUG; (default: $loglevel)
--se	<settings>		(default: $settings)
-t	<threads>		(default: take it from SGE environment variable NSLOTS)
-nd				do not delete the temporary files
-h				show the help text\n
";

GetOptions(
"b" => \$bam,
"o=s" => \$outprefix,
"f=s" => \$infile1, 
"F=s" => \$infile2,
"t=s" => \$threads,
"h" => \$help,
"rg=s" => \$readGroup,
"lf=s" => \$logfile,
"ll=s" => \$loglevel,
"se=s" => \$settings,
"s=s" => \$readGroupSample,
"l=s" => \$readGroupLib,
"nd" => \$doDelete,
"ta=s" => \$tmpArgument
);


if($threads == -1){
	$threads = 1;
	if($ENV{NSLOTS}){		#get number of slots given by SGE
		$threads = $ENV{NSLOTS};
		
	}
}

$starProgramPath = $params->{programs}->{star}->{path};
$starFusionProgramPath = $params->{programs}->{starfusion}->{path};
$starIndex = $params->{settings}->{$settings}->{starindex};
$starlogfile = $params->{programs}->{star}->{command}->{logfile};
my $bedtools      = $params->{programs}->{bedtools}->{path};
my $sammaxmem     = $params->{programs}->{samtools}->{maxmem};
my $samtools      = $params->{programs}->{samtools}->{path};		

if ($tmpArgument) {
	$starIndex = $tmpArgument;
}


#init the logger
Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();

if ($help == 1) {
	print $helptext;
	exit(1);
}
if ($starIndex eq "") {
	$logger->error("Path to STAR index files missing");
	exit(1);
}
if ($starProgramPath eq "") {
	$logger->error("STAR program path missing");
	exit(1);
}
if ($infile1 eq "") {
	$logger->error("Infile 1 missing");
	exit(1);
}
if ($infile2 eq "") {
	$logger->info("File with read pairs (\"-F <infile read 2>\") not specified! Mapping will be performed in single end mode");
	if (!$doSingleEnd) {
		$doSingleEnd = 1;
	}
}
if ($outprefix eq "") {
	$logger->error("Output prefix missing");
	exit(1);
}

my $outDir = dirname($outprefix);

my @infiles1 = split(",", $infile1);
my @infiles2 = split(",", $infile2);
if ($infile2 ne "" && scalar @infiles1 != scalar @infiles2) {
	$logger->error("$infile1 and $infile2 have different number of input files");
	exit(1);
}

if($bam){
	my $bedCommand = "";
	my $if1 = "";
	my $if2 = "";
	my $fqFile = "";
	for (my $i=0; $i<@infiles1; $i++) {
		$bedCommand = "$bedtools/bedtools bamtofastq -i $infiles1[$i]"; # -fq r1.fastq -fq2 r2.fastq "
		$fqFile = $outprefix . "_" . $i . "_R1.fastq";
		$if1 .= $fqFile . ",";
		push(@files2Delete, $fqFile);
		$bedCommand .= " -fq $fqFile";
		if($infile2 ne ""){
			$fqFile = $outprefix . "_" . $i . "_R2.fastq";
			$bedCommand .= " -fq2 $fqFile";
			$if2 .= $fqFile . ",";
			push(@files2Delete, $fqFile);
		}
		&Utilities::executeCommand($bedCommand, "Converting BAM to FASTQ", $logger);
	}
	$if1 =~ s/,$//;
	$if2 =~ s/,$//;
	$infile1 = $if1;
	$infile2 = $if2;
}

$starlogfile = $outprefix . "." . $starlogfile;
#build readgroup
my @rg = split(",", $readGroup);
$readGroup = "";
for (my $i=0; $i<@rg; $i++) {
	$readGroup .= "ID:" . $rg[$i] . " SM:" . $readGroupSample . " LB:Lib1 PL:Illumina , ";
}
$readGroup =~ s/ , $//;

#prepare the call of the STAR alignment
$command = &prepareSTARCommand();

#start STAR alignment
&Utilities::executeCommand($command,"Start STAR alignment with command:\n$command", $logger);

#delete STAR tmp folders
$command="";
$command="rm -r $outprefix/$readGroupSample"."_STARgenome $outprefix/$readGroupSample"."_STARpass1 $outprefix/$readGroupSample"."_STARtmp"; 
&Utilities::executeCommand($command,"Removing STAR tmp folders:\n$command", $logger);

#sort bam file
$command = "$samtools sort -\@ $threads -m $sammaxmem  $outprefix"."Aligned.out.bam $outprefix"."Aligned.out.sort";
&Utilities::executeCommand($command,"Sort BAM file:\n$command", $logger);


###################################
##delete the temporary files if set
###################################
sub deleteTempFiles() {
	foreach my $delFile (@files2Delete) {
		&Utilities::executeCommand("rm $delFile", "Deleting $delFile", $logger);
	}
}

###########################################
##function to prepare the commands to call in this script
###########################################
sub prepareSTARCommand() {
	my $c = "";
	$c = "$starProgramPath";
	$c .= " --runThreadN $threads";
	$c .= " --genomeDir $starIndex";
	$c .= " --readFilesIn $infile1";
	if ($infile2 ne "") {
		$c .= " $infile2";
	}
	#set output folder
	$c .= " --outFileNamePrefix $outprefix";
	$c .= " --readFilesCommand ";
	if ($infile1 =~ m/\.gz$/) {
		$c .= "zcat";
	} else {
		$c .= "cat";
	}
	$c .= " --outSAMtype BAM Unsorted"; # SortedByCoordinate";
	$c .= " --chimSegmentMin 20";
	$c .= " --quantMode GeneCounts";
	$c .= " --twopassMode Basic";
	$c .= " --outSAMunmapped Within";
	$c .= " --limitOutSJcollapsed 10000000";
	$c .= " --limitIObufferSize 500000000";
	$c .= " --outSAMattrRGline $readGroup";
	$c .= " --outSAMattributes All";
	$c .= " > $starlogfile 2>&1";
	
	return $c;
}
