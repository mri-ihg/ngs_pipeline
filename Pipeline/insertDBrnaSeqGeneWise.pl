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

my $delete        = 0;
my $countFile    = "";
my $fpkmFile    = "";
my $sampleName	= "";
my $error = 0;

my $line      = "";
my $settings  = "";
my $help      = 0;

GetOptions(
	"c=s"  => \$countFile,
	"f=s"  => \$fpkmFile,
	"s=s"  => \$sampleName,
	"d"    => \$delete,
	"se=s" => \$settings,
	"lf=s" => \$logfile,
	"ll=s" => \$loglevel,
	"h"    => \$help,
);

my $params         = Utilities::getParams();

my $rnaDB          = $params->{settings}->{$settings}->{rnadb}->{database};
my $geneBasedTable = $params->{settings}->{$settings}->{rnadb}->{genebasedtable};
my $geneTable	   = $params->{settings}->{$settings}->{exomedb}->{genetable};
my $geneDB        = $params->{settings}->{$settings}->{exomedb}->{database};
my $coreDB		   = $params->{coredb}->{database};
my $sampleTable    = $params->{coredb}->{sampletable};

Utilities::initLogger( $logfile, $loglevel );
my $logger = Utilities::getLogger();

if ( $help == 1 ) {
	print "-c	<count file>
-f	<fpkm file> 
-d	delete all entries for this sample before inserting
-s	<sample name>
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
if ( $countFile eq "" ) {
	$logger->error("no count file given");
	exit(1);
}
if ( $fpkmFile eq "" ) {
	$logger->error("no fpkm file given");
	exit(1);
}
if ( $settings eq "" ) {
	$logger->error("no settings given");
	exit(1);
}

$dbh = Utilities::connectRnaDB($settings);

#check if sample id is valid
$sql = qq{ select idsample from $coreDB.$sampleTable where name = '$sampleName' };
$sth = $dbh->prepare($sql) || $logger->error("Can't prepare statement: $DBI::errstr");
$sth->execute() || $logger->error($DBI::errstr);
my $idSample = $sth->fetchrow_array();
if ( $sth->rows < 1 ) {
	$logger->error("sample name $sampleName not found in $coreDB.$sampleTable - exiting");
	exit(1);
}

#delete from snvsample table
if ($delete) {
	&deleteSampleEntries($idSample, $rnaDB . "." . $geneBasedTable);
}

#insert counts
open(IN, $countFile) || exit $logger->error("Error opening $countFile");
my @countLine = "";
while ( <IN>) {
	chomp $_;
	@countLine = split('\t', $_);
	&dbInsert($idSample, \@countLine, $rnaDB . "." . $geneBasedTable, "readcount");
}
close IN;

#insert fpkm values
open(IN, $fpkmFile) || exit $logger->error("Error opening $fpkmFile");
my @fpkmLine = "";
while( <IN> ) {
	chomp $_;
	@fpkmLine = split('\t', $_);
	&dbInsert($idSample, \@fpkmLine, $rnaDB . "." . $geneBasedTable, "fpkm");
}
close IN;

if ($error) {
	$logger->error("One or more errors happend during inserts! Maybe you should delete all inserted lines for sample $sampleName with db-sampleid $idSample in table $rnaDB.$geneBasedTable and rerun insertion");
	#&deleteSampleEntries($idSample, $rnaDB . "." . $geneBasedTable);
}


###########################
# dbInsert: fill databases#
###########################
sub dbInsert {
	my $idsample = shift;
	my $il = shift;
	my @insertLine = @$il;
	my $insertDB = shift;
	my $attribute = shift;
		
	#get the correct geneid from db
	$sql = qq{select idgene from $geneDB.$geneTable where genesymbol = '$insertLine[0]'};
	$sth = $dbh->prepare($sql) || $logger->error("Can't prepare statement: $DBI::errstr");
	$sth->execute() || $logger->error($DBI::errstr);
	
	if ($sth->rows != 1) {
		$logger->info("Error retrieving gene id for gene symbol $insertLine[0]! " . $sth->rows . " rows found in DB");
	} else {
		my ($idgene) = $sth->fetchrow_array();
		
		#check if a db entry for geneid<->sampleid already exists
		$sql = qq{select idsample, idgene from $insertDB where idsample=$idsample and idgene=$idgene};
		$sth = $dbh->prepare($sql) || $logger->error("Can't prepare statement: $DBI::errstr");
		$sth->execute() || $logger->error($DBI::errstr);
		if($sth->rows > 0) {
			$sql = qq{update $insertDB set $attribute=$insertLine[1] where idsample=$idsample and idgene=$idgene};
		} else {
			$sql = qq{insert into $insertDB (idsample, idgene, $attribute) values($idsample, $idgene, $insertLine[1])};
		}
		$logger->debug($sql);
		$sth = $dbh->prepare($sql) || $logger->error("Can't prepare statement: $DBI::errstr");
		$sth->execute() or $error=1;
	}
}

###########################
# deleteSampleEntries: delete all lines in the table $rnaCountDB for a specific sample
###########################
sub deleteSampleEntries {
	my $idsample = shift;
	my $insertDB = shift;
	
	$sql = qq{delete from $insertDB where idsample=$idsample};
	$sth = $dbh->prepare($sql) || $logger->error("Can't prepare statement: $DBI::errstr");
	$sth->execute() || $logger->error($DBI::errstr);
}
