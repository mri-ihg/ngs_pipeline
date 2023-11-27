#!/usr/bin/perl

use strict;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);
use DateTime;
umask(002);

my $prog_path = dirname(abs_path($0));
require $prog_path."/Utilities.pm";

my $help     = 0;
my $dbh = "";
my $settings = "hg19";
my $query = "";
my $params = Utilities::getParams();
my $logfile  = "pipeline.log";
my $loglevel = "INFO";
my @caseSampleIds;
my @controlSampleIds;
my $bamFileLocation = "";
my $caseCountFile = "";
my $controlCountFile = "";
my $outputFolder = "";
my $experiment = "";
my $countFile = "";
my $conditionFile = "";
my $Rexecutable = "R";
my $fileFolder = "";
my $noCreateHeatMap = 1;
my $deSeq = "";
my $fileExtension = "htseqcounts";
my $bamFileName = "merged.bam";
my $noFCshrinkage = 0;
my $helptext      = 
"
This scripts starts the RNA-seq pipeline

-ca	<case sample id's>
-co	<control sample id's>
-o	<output directory>
-e	<experiment>
-lf	<log file>	default: pipeline.log
-ll	<log level>	ERROR,INFO,DEBUG; default: INFO
-h	this help text\n";

GetOptions(
"o=s" => \$outputFolder,
"e=s" => \$experiment,
"ca=s" => \@caseSampleIds,
"co=s" => \@controlSampleIds,
"lf=s" => \$logfile,
"nofcs" => \$noFCshrinkage,
"se=s" => \$settings,
"ll=s" => \$loglevel,
"h" => \$help);

#create outputfolder if it doesn't exist
unless (-d $outputFolder) {
	my $ret = system("mkdir $outputFolder");
	if ($ret) {
		print STDERR "CMD died with return value: $ret\n";
		exit(11);
	}
}

Utilities::initLogger($outputFolder."/".$logfile,$loglevel);
my $logger = Utilities::getLogger();

#connect to database
my $db     = $params->{coredb}->{database};
my $host   = $params->{coredb}->{host};
my $port   = $params->{coredb}->{port};
my $user   = $params->{coredb}->{user};
my $passwd = $params->{coredb}->{password};

unless ( $dbh = DBI->connect( "DBI:mysql:database=$db;host=$host;port=$port", $user, $passwd ) ) {
	DBI->connect( "DBI:mysql:database=$db;host=$host;port=$port", $user, $passwd ) || die print "$DBI::errstr\n";
}

my $fileMerger = $params->{programs}->{rnaseq}->{htseqfilemerger};

if ($experiment eq "") {
	$experiment = "NONAME";
}

$countFile = $experiment . ".all.counts";
$conditionFile = $countFile . ".conditions";

if (!($caseCountFile && $controlCountFile)) {
	$caseCountFile = &getFileLocations(\@caseSampleIds, "case", $fileFolder);
	$controlCountFile = &getFileLocations(\@controlSampleIds, "control", $fileFolder);
}

#merge the count files with htseqOutputFileMerger.pl; kind of a safety-feature and easier for the later R script
my $command = "perl $prog_path/$fileMerger -as $caseCountFile -us $controlCountFile -cfn $countFile -condfn $conditionFile -o $outputFolder -fe $fileExtension";
if (&Utilities::executeCommand($command, "Running: $fileMerger to merge count files in $caseCountFile and $controlCountFile", $logger)) {
	exit(11);
}

##DESeq2 analysis
my $deseq2OutputFolder = "$outputFolder/DESeq2";
$command = "mkdir -p $deseq2OutputFolder";
if (&Utilities::executeCommand($command, "making output directory for DESeq2", $logger)) {
	exit(11);
}
$deSeq = $params->{programs}->{rnaseq}->{RScripts}->{DESeq2};
$command = "$Rexecutable --no-restore --no-save --args $outputFolder/$countFile $outputFolder/$conditionFile $deseq2OutputFolder $experiment fpkmFile $noCreateHeatMap $noFCshrinkage < $prog_path/$deSeq > $deseq2OutputFolder/DESeq2.Rlog 2>&1";
if (&Utilities::executeCommand($command, "Running: DESeq2 R script for differential expression analysis", $logger)) {
	exit(11);
}

##Gene Ontology Analysis - goseq
#if (-f "$deseq2OutputFolder/session.RData" ) {
#	my $goseqOutputFolder = "$deseq2OutputFolder/goseq";
#	$command = "mkdir -p $goseqOutputFolder";
#	if (&Utilities::executeCommand($command, "making output directory for goseq", $logger)) {
#		exit(11);
#	}
#	my $goseq = $params->{programs}->{rnaseq}->{RScripts}->{goseq};
#	my $command = "$Rexecutable --no-restore --no-save --args $goseqOutputFolder/ $deseq2OutputFolder/session.RData $settings < $prog_path/$goseq > $goseqOutputFolder/goseq.Rlog 2>&1";
#	
#	if (&Utilities::executeCommand($command, "Running: GOseq script for Gene Ontology analysis", $logger)) {
#		exit(11);
#	}
#} else {
#	print "Neccessary file 'session.RData' missing in DESeq2 output folder\n";
#	exit(11);
#}
##PATHWAY analysis
#my $goseqOutputFolder = "$deseq2OutputFolder/pathway";
#$command = "mkdir -p $goseqOutputFolder";
#if (&Utilities::executeCommand($command, "making output directory for pathway enrichment analysis", $logger)) {
#	exit(11);
#}
#my $gage = $params->{programs}->{rnaseq}->{RScripts}->{gage};
#my $command = "$Rexecutable --no-restore --no-save --args $goseqOutputFolder/ $settings $deseq2OutputFolder/$experiment.csv < $prog_path/$gage > $goseqOutputFolder/gage.Rlog 2>&1";
#
#if (&Utilities::executeCommand($command, "Running: GAGE script for Pathway enrichment analysis", $logger)) {
#	exit(11);
#}



#######
## subs:
########
sub getFileLocations {
	my $sampleIdsRef = shift;
	my $name = shift;
	my $ff = shift;
	my $bamFileLocation = "";
	my $resultFileLocation = "$outputFolder/$name.tmp";
	
	open OUT, ">$resultFileLocation";
	if ($ff eq "") {
		foreach my $sId (@$sampleIdsRef) {
			my $c = "find -L /data/isilon/seq/analysis/exomehg19/S*/ -maxdepth 1 -type d -name $sId 2> /dev/null";			#maxdepth 2, otherwise it'd be possible that backupfolders like /old/sId/ could be found instead of the actual one
			$logger->info("CMD: $c");
			open FIND, "$c | " or exit $logger->error("Error opening $c |");
			my $loc = <FIND>;
			chomp($loc);
			close FIND;
			if ($loc eq "") {
				$logger->error("Couldn't find bam file on system: $sId");
				exit(11);
			} else {
				$loc .= "/RNAout/paired-endout/$bamFileName";
				print OUT "$sId\t$loc\n";
			}
		}
	}
	close OUT;
	return $resultFileLocation;
}
