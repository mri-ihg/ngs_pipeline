#!/usr/bin/perl

############################################################
## 10.09.2013 (Thomas Schwarzmayr): created
## 21.07.2014 (Thomas Schwarzmayr): bugfix of the getTranscriptLengthes() algorithm
############################################################
## calculates the FPKM value for a given sample

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

my $transcriptId  = ""; 
my $transcriptCount = 0;
my $fpkmValue = 0;
my $featureFile = "";

my $helptext      = 
"Calculate the FPKM count for a sample based on the numbers priour calculated by htseq-count

-i	<infile>	count file; output of the htseqCount.pl script
-o	<output prefix>
-r  <feature file>	optional; the gtf file with the features in 
-lf	<log file>		the log file for the pipeline logging (default: pipeline.log)
-ll	<log level>		the log level for the pipeline logging; available options: ERROR, INFO, DEBUG; (default: $loglevel)
-se	<settings>		(default: $settings)
-h				show help text\n
";

GetOptions(
"h" => \$help,
"i=s" => \$input,
"lf=s" => \$logfile,
"ll=s" => \$loglevel,
"se=s" => \$settings,
"o=s" => \$outprefix,
"r=s" => \$featureFile
);

if ($help == 1) {
	print $helptext;
	exit(1);
}

my $params = Utilities::getParams();
my $samtools      = $params->{programs}->{samtools}->{path};
my $sammaxmem     = $params->{programs}->{samtools}->{maxmem};

my $htseqCount	  = $params->{programs}->{htseq}->{count};
if (!$featureFile) {
	$featureFile	  = $params->{settings}->{$settings}->{gemannotation};	#use same gtf file as for the alignment!
}

#init the logger
Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();

$logger->info("Start calculating transcript lengthes from gtf file $featureFile");
my $start_time = time();
my %transcriptLength = &getTranscriptLengthes($featureFile);
my $end_time = time();
$logger->info("Finished calculation in " . ($end_time - $start_time) . " seconds");

$logger->info("Start calculating FPKM values");
$start_time = time();
open TOTAL, "awk 'BEGIN{SUM=0}{SUM+=\$2}END{print SUM}' $input | " or exit $logger->error("Couldn't open awk 'BEGIN{SUM=0}{SUM+=$2}END{print SUM}' $input | ");
my $totalMapped = <TOTAL>;
chomp $totalMapped;
close TOTAL;
open IN, "$input" or exit $logger->error("Couldn't open file $input");
open OUT, ">$outprefix.fpkm" or exit $logger->error("Couldn't open file for writing: $outprefix.fpkm");

my $i=0;
while (<IN>) {
	$i++;
	($transcriptId, $transcriptCount) = split("\t");
	chomp $transcriptCount;
	next if ($transcriptId =~ m/^$/);
	if ($transcriptLength{$transcriptId} > 0) {
		$fpkmValue = (10**9 * $transcriptCount) / ($totalMapped * $transcriptLength{$transcriptId});
	} else {
		$logger->info("Feature $transcriptId has length 0 or IDs do not match");
	}
	print OUT "$transcriptId\t$fpkmValue\n";
	$fpkmValue = 0;
}
close IN;
close OUT;
$end_time = time();
$logger->info("Successfully finished calculating FPKM values in " . ($end_time - $start_time) . " seconds");


#calc transcript lengthes from gtf file
sub getTranscriptLengthes {
	my $ff = shift; #featureFile
	my %gl = {};
	my $currentTranscriptId = "";
	my $currentGeneId = "";
	my $currentTranscriptLength = 0;
	
	my $inProcessTranscriptId = "NOTDEFINEDRIGHTNOW";
	my $inProcessGeneId = "";
	
	my $i = 0;
	open REF, qq{cat $featureFile | grep exon | sort -t\$'\\t' -k9 | } or exit $logger->error("Couldn't open file cat $featureFile | sort -t\$'\\t' -k9 | ");
	while (<REF>) {
		$i++;
		my ($chr, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split("\t");
		next unless ($type eq "exon");
		$attributes =~ /.*gene_id "([^"]+)"/;
		$currentGeneId = $1;
		$attributes =~ /.*transcript_id "([^"]+)"/;
		$currentTranscriptId = $1;
		if ($inProcessTranscriptId eq "NOTDEFINEDRIGHTNOW") {
			$inProcessTranscriptId = $currentTranscriptId;
		}
		if ($currentTranscriptId eq $inProcessTranscriptId) {
			$currentTranscriptLength += $end - $start + 1;
		} else {
			$inProcessTranscriptId = $currentTranscriptId;
			if (!exists $gl{$inProcessGeneId} || $gl{$inProcessGeneId} < $currentTranscriptLength) {
				$gl{$inProcessGeneId} = $currentTranscriptLength;
			}
			
			$currentTranscriptLength = $end - $start + 1;
			if ($inProcessGeneId ne $currentGeneId) {
				$inProcessGeneId = $currentGeneId;
			}
		}
	}
	if (!exists $gl{$inProcessGeneId} || $gl{$inProcessGeneId} < $currentTranscriptLength) {
		$gl{$inProcessGeneId} = $currentTranscriptLength;
	}
	close REF;
	return %gl;
}