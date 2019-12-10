#!/usr/bin/perl

## calculate RNA-SeQC metrics for RNA sample

use strict;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);
umask(002);

#include Utilities.pm
my $prog_path = dirname(abs_path($0));
require $prog_path."/Utilities.pm";

my $settings      = "hg19_test";
my $logfile  	  = "./pipeline.log";
my $loglevel 	  = "INFO";
my $input		  = "";
my $help		  = "";
my $outprefix	  = "";
my $featureFile   = "";
my $sampleId	  = "";
my $command       = "";
my $noDelFolder   = 0;
my $tmpArgument ="";

my $helptext      = 
"Perform RNA-SeQC analysis for rna seq bam files

-i	<infile>	bam file; output of the mapping program
-o	<output prefix>
-s	<sampleId>
-nodel				do not delete the RNA-SeQC folder in the samples project directory if it already exists (default: FALSE)
-lf	<log file>		the log file for the pipeline logging (default: pipeline.log)
-ll	<log level>		the log level for the pipeline logging; available options: ERROR, INFO, DEBUG; (default: $loglevel)
-se	<settings>		(default: $settings)
-h				show help text\n
";

GetOptions(
"h" => \$help,
"i=s" => \$input,
"lf=s" => \$logfile,
"nodel" => \$noDelFolder,
"ll=s" => \$loglevel,
"s=s" => \$sampleId,
"se=s" => \$settings,
"o=s" => \$outprefix,
"ta=s" => \$tmpArgument
);

if ($help == 1) {
	print $helptext;
	exit(1);
}

if ($input eq "" || $outprefix eq "" || $sampleId eq "") {
	print $helptext;
	exit(1);
}

my $params = Utilities::getParams();
my $java           = $params->{programs}->{java7}->{path};
my $rnaseqcProg	   = $params->{programs}->{rnaseqc}->{path};
my $rnaseqcOut     = $params->{programs}->{rnaseqc}->{outputfolder};
my $ref            = $params->{settings}->{$settings}->{reference};
my $annotationfile = $params->{settings}->{$settings}->{gemannotation};
my $rrnalist	   = $params->{settings}->{$settings}->{rrnalist};

my $dir = dirname($outprefix);
my $sampleFile = "$dir/sample.file";
my $outputDir = "$dir/$rnaseqcOut/";

#init the logger
Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();

#check if outputdir exists
if (-d $outputDir || -e $outputDir) {
	$logger->info("Output directory $outputDir already exists");
	if (!$noDelFolder) {
		$logger->info("Argument \"-nodel\" not specified! Deleting folder $outputDir");
		$command = "rm -r $outputDir";
		&Utilities::executeCommand($command, "Perform deletion process for $outputDir", $logger);
		$command = "mkdir $outputDir";
		&Utilities::executeCommand($command, "Perform mkdir for $outputDir", $logger)
	}
} else {
	$command = "mkdir $outputDir";
	&Utilities::executeCommand($command, "Perform mkdir for $outputDir", $logger)
}
$command = "";

$logger->info("creating sample file for RNA-SeQC");
open OUT, ">$sampleFile";
##RNA-SeQC defined formatting
print OUT "Sample ID\tBam File\tNotes\n";
print OUT "$sampleId\t$input\tRNASeQCfor$sampleId";
close OUT;

$command = "$java -XX:ParallelGCThreads=1 -jar $rnaseqcProg -s $sampleFile -o $outputDir -r $ref -t $annotationfile -rRNA $rrnalist > $outputDir/RNASeQC.log 2>&1";
if (&Utilities::executeCommand($command, "Running: RNA-SeQC for sample $sampleId (file: $input)", $logger)) {
	exit(11);
}