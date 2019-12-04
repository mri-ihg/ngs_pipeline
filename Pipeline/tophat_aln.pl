#!/usr/bin/perl

############################################################
## 15.07.2015 (Thomas Schwarzmayr): created
############################################################

use strict;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);
umask(002);

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
my $tophat2ProgramPath = "";
my $bowtieIndex = "";
my $bowtie2Index = "";
my $tophat2logfile = "";
my $doSingleEnd = 0;
my $readGroupSample = "";
my $readGroupLib    = "Lib1";
my $platform        = "Illumina";
my $readGroupId = "";
my $bam = 0;
my $readGroup = "";
my $isStrandedRNA = 0;

my $params = Utilities::getParams();

my $helptext      = 
"
-f	<infile reads 1>	
-F	<infile reads 2>	
-o	<output prefix>
-b 					if infile is in bam-format
-s	<read group sample>	read group sample name
--lf	<log file>		the log file for the pipeline logging (default: pipeline.log)
--ll	<log level>		the log level for the pipeline logging; available options: ERROR, INFO, DEBUG; (default: $loglevel)
--se	<settings>		(default: $settings)
--rg	<readgroup>
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
"rg=s" => \$readGroupId,
"lf=s" => \$logfile,
"ll=s" => \$loglevel,
"se=s" => \$settings,
"s=s" => \$readGroupSample,
"sr" => \$isStrandedRNA,
"nd" => \$doDelete
);


if($threads == -1){
	$threads = 1;
	if($ENV{NSLOTS}){		#get number of slots given by SGE
		$threads = $ENV{NSLOTS};
		
	}
}

$tophat2ProgramPath = $params->{programs}->{tophat2}->{path};
$bowtieIndex = $params->{settings}->{$settings}->{bowtieindex};
$bowtie2Index = $params->{settings}->{$settings}->{bowtie2index};
$tophat2logfile = $params->{programs}->{tophat2}->{command}->{logfile};
$annotationfile = $params->{settings}->{$settings}->{gemannotation}; #file of the gtf annotation file
my $bedtools      = $params->{programs}->{bedtools}->{path};

my $sammaxmem     = $params->{programs}->{samtools}->{maxmem};
my $samtools      = $params->{programs}->{samtools}->{path};	

#init the logger
Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();

if ($help == 1) {
	print $helptext;
	exit(1);
}
if ($bowtieIndex eq "") {
	$logger->error("Path to bowtie index files missing");
	exit(1);
}
if ($tophat2ProgramPath eq "") {
	$logger->error("tophat2 program path missing");
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
if (-d $outDir) {
	&Utilities::executeCommand("rm -rf $outDir", "Deleting existing output directory $outDir", $logger);
}
&Utilities::executeCommand("mkdir $outDir", "Creating output folder $outDir", $logger);

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

$tophat2logfile = $outprefix . "." . $tophat2logfile;

#prepare the call of the tophat2 alignment
$command = &prepareTophatCommand();

#start tophat2 alignment
&Utilities::executeCommand($command,"Start tophat2 alignment with command:\n$command", $logger);






#sort bam file
$command ="cat /PATHTO/hg19.dict > $outDir/"."accepted_hits.header.sam";
&Utilities::executeCommand($command,"Sort BAM file:\n$command", $logger);
$command ="samtools view -H $outDir/"."accepted_hits.bam | grep \"\@RG\" >> $outDir/"."accepted_hits.header.sam";
&Utilities::executeCommand($command,"Sort BAM file:\n$command", $logger);
$command="samtools view -H $outDir/"."accepted_hits.bam | grep \"\@PG\" >> $outDir/"."accepted_hits.header.sam";
&Utilities::executeCommand($command,"Sort BAM file:\n$command", $logger);
$command = "(cat $outDir/"."accepted_hits.header.sam; samtools view $outDir/"."accepted_hits.bam) | $samtools view -Sb - | $samtools sort -\@ $threads -m $sammaxmem - $outDir/"."accepted_hits.sort";
&Utilities::executeCommand($command,"Sort BAM file:\n$command", $logger);


#the last thing to do!
if (!$doDelete) {
	&deleteTempFiles();
}


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
sub prepareTophatCommand() {
	my $c = "";
	$c = "$tophat2ProgramPath";
	$c .= " -p $threads";
	if ($isStrandedRNA) {
		$c .= " --library-type fr-firststrand";	#for Illumina stranded RNA samples"
	} else {
		$c .= " --library-type fr-unstranded";
	}
	$c .= " --rg-id $readGroupId";
	$c .= " --rg-sample $readGroupSample";
	$c .= " --rg-library $readGroupSample"."_$readGroupLib";
	$c .= " --rg-description NULL";
	$c .= " --rg-platform-unit 1";
	$c .= " --rg-center IHG";
	my @t = localtime;
	$t[5] += 1900;
	$t[4]++;
	$c .= sprintf(" --rg-date %04d-%02d-%02d", @t[5,4,3]);
	$c .= " --rg-platform $platform";
	$c .= " -G $annotationfile";
	#set output folder
	$c .= " -o $outDir";
	
	$c .= " $bowtie2Index";
	if ($isStrandedRNA) {
		if ($infile2 ne "") {
			$c .= " $infile2";
		}
		$c .= " $infile1";
	} else {
		$c .= " $infile1";
		if ($infile2 ne "") {
			$c .= " $infile2";
		}
	}

	$c .= " > $tophat2logfile 2>&1";
	
	return $c;
}
