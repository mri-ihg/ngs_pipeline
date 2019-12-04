#!/usr/bin/perl

use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use DBI;
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);
use Data::Dumper;
use List::MoreUtils qw/ uniq /;

my $prog_path = dirname( abs_path($0) );
require $prog_path . "/Utilities.pm";

# database
my $dbh           = "";
my $sql           = "";
my $sth           = "";
my $logfile       = "pipeline.log";
my $loglevel      = "INFO";
my $sampleName	  = "";
my $settings      = "hg19_test";
my $help          = 0;
my $folder		  = "";
my $mapper        = "gem"; #mapper that mapped the reads
my $fileName      = ""; #specify the file the metrics are based on
my $noDelete      = 0;
my $error         = 0;

#hash for mapping the header to the db attribute names since columns in RNA-SeQC result file can be unsorted
my %mappingHash = (	"End 2 Mapping Rate",						"end2MappingRate",
					"Chimeric Pairs",							"chimericPairs",
					"Intragenic Rate",							"intragenicRate",
					"Num. Gaps",								"numGaps",
					"Exonic Rate",								"exonicRate",
					"Mapping Rate",								"mappingRate",
					"5' Norm",									"norm5p",
					"Genes Detected",							"genesDetected",
					"Unique Rate of Mapped",					"uniqueRateOfMapped",
					"3' Norm",									"norm3p",
					"Read Length",								"readLength",
					"Mean Per Base Cov.",						"meanPerBaseCov",
					"End 1 Mismatch Rate",						"end1MismatchRate",
					"Fragment Length StdDev",					"fragmentLengthStdDev",
					"Estimated Library Size",					"estimatedLibrarySize",
					"Mapped",									"mapped",
					"Intergenic Rate", 							"intergenicRate",
					"Total Purity Filtered Reads Sequenced",	"totalPurityFilteredReadsSequenced",
					"rRNA",										"rRNA",
					"Failed Vendor QC Check",					"failedVendorQCCheck",
					"Mean CV",									"meanCV",
					"Transcripts Detected",						"transcriptsDetected",
					"Mapped Pairs",								"mappedPairs",
					"Cumul. Gap Length",						"cumulGapLength",
					"Gap %",									"gapPercent",
					"Unpaired Reads",							"unpairedReads",
					"Intronic Rate",							"intronicRate",
					"Mapped Unique Rate of Total",				"mappedUniqueRateofTotal",
					"Expression Profiling Efficiency",			"expressionProfilingEfficiency",
					"Mapped Unique",							"mappedUnique",
					"End 2 Mismatch Rate",						"end2MismatchRate",
					"End 2 Antisense",							"end2Antisense",
					"Alternative Aligments",					"alternativeAligments",
					"End 2 Sense",								"end2Sense",
					"Fragment Length Mean",						"fragmentLengthMean",
					"End 1 Antisense",							"end1Antisense",
					"Split Reads",								"splitReads",
					"Base Mismatch Rate",						"baseMismatchRate",
					"End 1 Sense",								"end1Sense",
					"End 1 % Sense",							"end1PercentSense",
					"rRNA rate",								"rRNARate",
					"End 1 Mapping Rate",						"end1MappingRate",
					"No. Covered 5'",							"noCovered5p",
					"Duplication Rate of Mapped",				"duplicationRateOfMapped",
					"End 2 % Sense",							"end2PercentSense"
				);
					
GetOptions(
	"m=s"  => \$mapper,
	"n=s"  => \$fileName,
	"s=s"  => \$sampleName,
	"f=s"  => \$folder,
	"nodel"=> \$noDelete,
	"se=s" => \$settings,
	"lf=s" => \$logfile,
	"ll=s" => \$loglevel,
	"h"    => \$help
);

Utilities::initLogger( $logfile, $loglevel );
my $logger = Utilities::getLogger();

if ( $help == 1 ) {
	print "\nInsert the results of RNA-SeQC into the DB
	
-m	<mapper>	name of the mapper that aligned the reads (e.g. gem, star, tophat) (default: $mapper)
-n	<file name>	name of the file the RNA-SeQC results are based on (i.e. merged.rmdup.bam or merged.bam) (default: $fileName)
-f	<folder>	the folder where to find the RNA-SeQC output folder is located
-s	<sample name>
-nodel			do not delete all entries for this sample before inserting (default: FALSE) 
-se	<settings>
-lf	log file; default: pipeline.log
-ll	log level: ERROR,INFO,DEBUG; default: INFO
-h	this help\n";
	exit(1);
}


if ( $sampleName eq "" ) {
	$logger->error("no sampleName specified");
	exit(1);
}
if ( $folder eq "" ) {
	$logger->error("no folder specified");
	exit(1);
}
if ( $settings eq "" ) {
	$logger->error("no settings given");
	exit(1);
}
if (! -d $folder ) {
	$logger->error("$folder cannot be found!");
	exit(1);
}

my $params         = Utilities::getParams();
my $coreDB		   = $params->{coredb}->{database};
my $sampleTable    = $params->{coredb}->{sampletable};
my $insertTable    = $params->{coredb}->{rnaseqctable};
my $rnaseqcOut     = $params->{programs}->{rnaseqc}->{outputfolder};
my $rnaseqcOutFile = $params->{programs}->{rnaseqc}->{outputfilename};

my $resultDir  = "$folder/$rnaseqcOut/";
my $resultFile = "$resultDir/$rnaseqcOutFile";

$dbh = Utilities::connectCoreDB($settings);

#read the RNA-SeQC result file
open(IN, "$resultFile") || exit $logger->error("Error opening $resultFile");
my $headerLine = <IN>;
if ($headerLine !~ m/^Sample/) {
	$logger->error("RNA-SeQC result file ($resultFile) seems to be corrupt! Please check...");
	exit(2);
}
my $resultLine = <IN>;
close IN;

chomp $headerLine;
chomp $resultLine;
my @header = split('\t', $headerLine);
my @result = split('\t', $resultLine);

if ((scalar @header) != (scalar @result)) {
	$logger->error("Lines of the RNA-SeQC result file do not have the same number of columns");
	exit(2);
}

#check if sample id is valid
$sql = qq{ select s.idsample 
			from $coreDB.$sampleTable s 
			inner join solexa.sample2library sl on sl.idsample=s.idsample 
			inner join solexa.library l on l.lid=sl.lid 
			inner join solexa.libtype lt on lt.ltid=l.libtype 
			where lt.ltlibtype='RNA'
			and name = '$sampleName' };
			
$sth = $dbh->prepare($sql) || $logger->error("Can't prepare statement: $DBI::errstr");
$sth->execute() || $logger->error($DBI::errstr);
my $idSample = $sth->fetchrow_array();
if ( $sth->rows < 1 ) {
	$logger->error("sample name $sampleName not found in $coreDB.$sampleTable - exiting");
	exit(1);
}

#build insert statement
$fileName = basename($fileName);
my $attributes = "(idsample, mapper, file, ";
my $values = "($idSample, '$mapper', '$fileName', ";
for (my $i=0; $i<@header; $i++) {
	next if ($header[$i] eq "Sample" || $header[$i] eq "Note");
	$attributes .= "$mappingHash{$header[$i]}, ";
	$values .= "$result[$i], ";
}
#delete last comma and substitute with )
$attributes =~ s/,\s$/)/;
$values =~ s/,\s$/)/;
my $insertStatement = qq{ INSERT INTO $insertTable $attributes VALUES $values };

#delete from rnastats table
if (!$noDelete) {
	&deleteSampleEntries($idSample, $coreDB . "." . $insertTable, $mapper,$fileName);
}

$logger->debug("Insert Statement: $insertStatement");
$sth = $dbh->prepare($insertStatement) || $logger->error("Can't prepare statement: $DBI::errstr");
$sth->execute() or $error=1;

if ($error) {
	$logger->error("One or more errors happend during inserts! Maybe you should delete all inserted lines for sample $sampleName with db-sampleid $idSample in table $coreDB.$insertTable and rerun insertion");
}

###########################
# deleteSampleEntries: delete all lines in the table $insertDB for a specific sample
###########################
sub deleteSampleEntries {
	my $idsample = shift;
	my $insertDB = shift;
	my $aligner  = shift; 
	my $fn		 = shift;
	
	$sql = qq{delete from $insertDB where idsample=$idsample and mapper='$aligner' and (file='$fn' or file='')};
	$sth = $dbh->prepare($sql) || $logger->error("Can't prepare statement: $DBI::errstr");
	$sth->execute() || $logger->error($DBI::errstr);
}
