#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use DBI;
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);
use Pod::Usage;

my $prog_path = dirname(abs_path($0));
require $prog_path."/Utilities.pm";

my $sql	       = "";
my $sth        = "";
my $infile   = "";
my $broad    = 0;
my $help     = 0;
my $man		 = 0;
my $logfile  = "SCREEN";
my $loglevel = "INFO";

my ($sample,$dup,$reads,$mapped,$mappedPer,$seq,$uncov,$x1,$x4,$x8,$x20,$avgCov,$std,$medianCov,$mstd,$autosomalCov,$onBait);
my $idsample = 0;

GetOptions(
"i=s"  => \$infile,
"b"    => \$broad,
"lf=s" => \$logfile,
"ll=s" => \$loglevel, 
"h"    => \$help,
"man"  => \$man);

pod2usage( {-exitval => 0  ,-verbose => 1} ) if $help;
pod2usage( {-exitval => 0  ,-verbose => 2} ) if $man;
pod2usage( {-exitval => 1  ,-verbose => 1} ) if $infile eq "";

Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();
my $params = Utilities::getParams();

my $database       = $params->{coredb}->{database};
my $statTable      = $params->{coredb}->{chipseqstatstable};
my $sampleTable    = $params->{coredb}->{sampletable};
my $dbh = Utilities::connectCoreDB();

open(IN, "$infile") || exit $logger->error("Cannot open $infile!");

#print "#ChIP sample\tInput sample\tUniquely mapped reads\tPeaks\tNRF(>0.8)\tFRiP(>0.01)\tPeaks overlapping coding regions\tReplicate sample\tPeaks overlapping replicate sample\n";

while (<IN>)
{
	chomp;
	$_ =~ s/\%//g;
	my ($sample,$inputsample,$reads,$peaks,$nrf,$frip,$coding,$repsample,$repratio,$filename,$numReads,$estFragLen,$corr_estFragLen,$phantomPeak,$corr_phantomPeak,$argmin_corr,$rmin_corr,$normSCC,$relSCC, $qualityTag) = split("\t",$_);
	
	next if $sample =~ /^#/;	#skip header

	$logger->debug("$sample,$inputsample,$reads,$peaks,$nrf,$frip,$coding,$repsample,$repratio,$normSCC,$relSCC, $qualityTag");
	&insert($sample,$inputsample,$reads,$peaks,$nrf,$frip,$coding,$repsample,$repratio,$normSCC,$relSCC, $qualityTag);
}

my $elapsed = (time - $^T)/60;
$logger->info("finished in $elapsed minutes");

#############
#subroutines#
#############
sub getIdsample {
	my $patId = shift;
	my $sql = "";
	my $sth = "";
	my $id = "";
	
	#check if patient id is valid
	$sql = qq{select idsample from $database.$sampleTable where name = '$patId' };
	$sth = $dbh->prepare( $sql ) || $logger->error("Can't prepare statement: $DBI::errstr");
	$sth->execute();
	($id)=$sth->fetchrow_array();
	if(!defined($id))
	{
		print "patient name $patId not found in $database.$sampleTable - exiting\n"; 
		exit(1);
	}
	return $id;
}

##############################################################################################################
sub insert {
	my ($sample,$inputsample,$reads,$peaks,$nrf,$frip,$coding,$repsample,$repratio,$nsc,$rsc, $qualitytag) = @_;
	my $sql = "";
	my $sth = "";

	my $idsample = &getIdsample($sample);
	my $idinput  = "NULL";
	$idinput     = &getIdsample($inputsample) if defined $inputsample && $inputsample ne ""; 
	my $idrep    = "NULL";
	$idrep       = &getIdsample($repsample) if defined $repsample && $repsample ne ""; 
	$repratio    = "NULL" unless defined $repsample && $repsample ne "";
	
	$sql = qq{insert into $database.$statTable ( idsample, idsamplecontrol,uniquereads,peaks,nrf,frip,codingpeaks,idsampleoverlap,overlap,broad, nsc, rsc, qualitytag) 
		      values($idsample,$idinput,$reads,$peaks,$nrf,$frip,$coding,$idrep,$repratio,$broad, $nsc, $rsc, $qualitytag) 
	          ON DUPLICATE KEY UPDATE uniquereads=$reads, peaks=$peaks, nrf=$nrf, frip=$frip, codingpeaks=$coding, idsampleoverlap=$idrep, overlap=$repratio, broad=$broad, nsc=$nsc, rsc=$rsc, qualitytag=$qualitytag};
	$logger->debug($sql);
	$sth = $dbh->prepare( $sql ) || $logger->error("Can't prepare statement: $DBI::errstr");
	$sth->execute() || $logger->error("Can't execute statement: $DBI::errstr");
}


=head1 NAME

insertChIPstats.pl

=head1 SYNOPSIS

insertChIPstats.pl -i ChIP.stats.tsv

=head1 DESCRIPTION

This script inserts the tab separated ChIP-seq statistics file produced by calcChIPstats.pl into the database. 


=head1 OPTIONS

 -i	<ChIP.stats.tsv> output file of the calcChIPstats.pl script; REQUIRED
 -b	"--broad" option was specified for macs2
 -lf	log file; default: print to screen
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page

=head1 AUTHOR

Thomas Wieland

=cut