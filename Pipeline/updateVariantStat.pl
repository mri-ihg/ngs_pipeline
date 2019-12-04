#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use DBI;
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);
use List::MoreUtils qw/ uniq /;
use Pod::Usage;
use Data::Dumper;

my $prog_path = dirname( abs_path($0) );
require $prog_path . "/Utilities.pm";

# database
my $dbh           = "";
my $sql           = "";
my $sth           = "";
my $logfile       = "SCREEN";
my $loglevel      = "INFO";

my $settings  = "hg19_wholegenome";
my $help      = 0;
my $man       = 0;

my $sampleName = "";
my $idsample   = -1;
my $idpipeline = -1;

my $numSNVs       = 0;
my $numIndels     = 0;
my $numPindel     = 0;
my $numExomeDepth = 0;
my $snvGTqual     = 0;
my $snvDepth      = 0;

GetOptions(
	"s=s"  => \$sampleName,
	"se=s" => \$settings,
	"pipe=s"=>\$idpipeline,
	"lf=s" => \$logfile,
	"ll=s" => \$loglevel,
	"h"    => \$help,
	"man"  => \$man
);

pod2usage( { -exitval => 0, -verbose => 1 } ) if $help;
pod2usage( { -exitval => 0, -verbose => 2 } ) if $man;
pod2usage( { -exitval => 1, -verbose => 1 } ) if $sampleName eq "";

my $params = Utilities::getParams();

my $snvTable                  = $params->{settings}->{$settings}->{exomedb}->{snvtable};
my $snvsampleTable            = $params->{settings}->{$settings}->{exomedb}->{snvsampletable};
my $variantStatTable          = $params->{settings}->{$settings}->{exomedb}->{variantstattable};

my $coredb                    = $params->{coredb}->{database};
my $sampleTable               = $params->{coredb}->{sampletable};

Utilities::initLogger( $logfile, $loglevel );
my $logger = Utilities::getLogger();
$dbh = Utilities::connectExomeDB($settings);


## PREPARE DB-STATEMENTS ##
#checkSampleIds 
my $checkSampleId_sql = qq{select idsample from $coredb.$sampleTable where name = ? };
my $checkSampleId_sth = $dbh->prepare($checkSampleId_sql) || $logger->error("Can't prepare statement: $DBI::errstr");
##############################

$checkSampleId_sth->execute($sampleName) || $logger->error($DBI::errstr);
$idsample = $checkSampleId_sth->fetchrow_array();

if ($idsample <= 0) {
	$logger->error("Could not find sample with name \"$sampleName\" in $coredb.$sampleTable! Exiting ...");
	exit 9;
} else {
	$logger->debug("Variantstat:\tidsample: ".$idsample."\tname: ".$sampleName."\n");
	## PREPARE DB-STATEMENTS ##
	#numSNVs
	my $numSNVs_sql = qq{select count(idsnvsample) from ( select * from $snvsampleTable where idsample=? ) ss inner join $snvTable v on v.idsnv=ss.idsnv where (caller='gatk' or caller='samtools') and v.class='snp'};
	my $numSNVs_sth = $dbh->prepare($numSNVs_sql) || $logger->error("Can't prepare statement: $DBI::errstr");
	#numIndels
	my $numIndels_sql = qq{select count(idsnvsample) from ( select * from $snvsampleTable where idsample=? ) ss inner join $snvTable v on v.idsnv=ss.idsnv where (caller='gatk' or caller='samtools') and v.class='indel'};
	my $numIndels_sth = $dbh->prepare($numIndels_sql) || $logger->error("Can't prepare statement: $DBI::errstr");
	#numPindel
	my $numPindel_sql = qq{select count(idsnvsample) from ( select * from $snvsampleTable where idsample=? ) ss where caller='pindel'};
	my $numPindel_sth = $dbh->prepare($numPindel_sql) || $logger->error("Can't prepare statement: $DBI::errstr");
	#numExomeDepth
	my $numExomedepth_sql = qq{select count(idsnvsample) from ( select * from $snvsampleTable where idsample=? ) ss where caller='exomedepth'};
	my $numExomedepth_sth = $dbh->prepare($numExomedepth_sql) || $logger->error("Can't prepare statement: $DBI::errstr");
	#snvGTqual
	my $snvGTqual_sql = qq{select avg(gtqual) from ( select * from $snvsampleTable where idsample=? ) ss inner join $snvTable v on v.idsnv=ss.idsnv where (caller='gatk' or caller='samtools') and v.class='snp'};
	my $snvGTqual_sth = $dbh->prepare($snvGTqual_sql) || $logger->error("Can't prepare statement: $DBI::errstr");
	#snvDepth
	my $snvDepth_sql = qq{select avg(coverage) from ( select * from $snvsampleTable where idsample=? ) ss inner join $snvTable v on v.idsnv=ss.idsnv where (caller='gatk' or caller='samtools') and v.class='snp'};
	my $snvDepth_sth = $dbh->prepare($snvDepth_sql) || $logger->error("Can't prepare statement: $DBI::errstr");
	#####
	#Insert
	my $insert_sql = qq{insert into $variantStatTable (idsample,snv,indel,pindel,exomedepth,snvgtqual, snvdepth) values (?,?,?,?,?,?,?) ON DUPLICATE KEY UPDATE snv=?, indel=?, pindel=?, exomedepth=?, snvgtqual=?, snvdepth=?};
	my $insert_sth = $dbh->prepare($insert_sql) || $logger->error("Can't prepare statement: $DBI::errstr");	
	#####
	#Update Pipeline record
	my $update_sql = qq{update $coredb.pipeline set  numofvars=?, numofsvs=?, numofcnvs=? where idpipeline=?};
	my $update_sth = $dbh->prepare($update_sql) || $logger->error("Can't prepare statement: $DBI::errstr");
	##############################
	
	##
	## CALCULATE NUMBERS ##
	$logger->info("Stats for sample: $sampleName");
	
	#numSNVs
	$numSNVs_sth->execute($idsample) || $logger->error($DBI::errstr);
	$numSNVs = $numSNVs_sth->fetchrow_array();
		$logger->info("numsnv: $numSNVs\n");
	#numIndels
	$numIndels_sth->execute($idsample) || $logger->error($DBI::errstr);
	$numIndels = $numIndels_sth->fetchrow_array();
		$logger->info("numindels: $numIndels\n");
	#numPindel
	$numPindel_sth->execute($idsample) || $logger->error($DBI::errstr);
	$numPindel = $numPindel_sth->fetchrow_array();
		$logger->info("numpindel: $numPindel\n");
	#numExomeDepth
	$numExomedepth_sth->execute($idsample) || $logger->error($DBI::errstr);
	$numExomeDepth = $numExomedepth_sth->fetchrow_array();
		$logger->info("numed: $numExomeDepth\n");
	#snvGTqual
	$snvGTqual_sth->execute($idsample) || $logger->error($DBI::errstr);
	$snvGTqual = $snvGTqual_sth->fetchrow_array();
	$snvGTqual = "NULL" if (!defined($snvGTqual));
		$logger->info("snvgtqual: $snvGTqual\n");
	#snvDepth
	$snvDepth_sth->execute($idsample) || $logger->error($DBI::errstr);
	$snvDepth = $snvDepth_sth->fetchrow_array();
	$snvDepth = "NULL" if (!defined($snvDepth));
		$logger->info("snvdepth: $snvDepth\n");
	##############################
	
	## INSERT ##
	$insert_sth->execute($idsample,$numSNVs,$numIndels,$numPindel,$numExomeDepth,$snvGTqual,$snvDepth,$numSNVs,$numIndels,$numPindel,$numExomeDepth,$snvGTqual,$snvDepth) || $logger->error($DBI::errstr);
	############
	
	## UPDATE PIPELINE RECORD ##
	$update_sth->execute($numSNVs+$numIndels, $numPindel, $numExomeDepth, $idpipeline) || $logger->error($DBI::errstr);
	############
}	
	
=head1 NAME

updateVariantStat.pl

=head1 SYNOPSIS

updateVariantStat.pl -id <sample> 

=head1 DESCRIPTION

This script updates the variantstat table

=head1 OPTIONS

 -s		<sample>; name of the sample; REQUIRED
 -se	name of the settings in the current.config.xml file that holds path to reference genome, 
 	to the annotation file and to possible additional annotation files; use default settings if nothing is given
 -pipe	id pipeline to update values
 -lf	log file; default: print to screen
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page

=head1 AUTHOR

Riccardo Berutti

=cut

