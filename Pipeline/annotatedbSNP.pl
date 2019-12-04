#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use DBI;
use POSIX;
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);
use Pod::Usage;


my $prog_path = dirname(abs_path($0));
require $prog_path."/Utilities.pm";


my $infile    = "";
my $outfile   = "";
my $settings  = "default";
my $logfile  = "SCREEN";
my $loglevel = "INFO";

#program flags
my $help       = 0;
my $man		   = 0;


GetOptions(
"i=s"  => \$infile, 
"o=s"  => \$outfile, 
"se=s" => \$settings,
"lf=s" => \$logfile,
"ll=s" => \$loglevel,
"man"  => \$man,
"h"    => \$help,);


pod2usage( { -exitval => 0, -verbose => 1 } ) if $help;
pod2usage( { -exitval => 0, -verbose => 2 } ) if $man;
pod2usage( { -exitval => 1, -verbose => 1 } ) if $infile eq "";




Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();



my $params	  = Utilities::getParams();
my $tabix     = $params->{programs}->{tabix}->{path};
my $snpsift   = $params->{programs}->{SnpSift}->{path};
my $java      = $params->{programs}->{java}->{path};
my $dbsnpfile = $params->{settings}->{$settings}->{dbsnpfile};


unless(-e "$snpsift"){
	$logger->error("Can't find $snpsift!");
	exit(1);
}

if($outfile eq ""){
	$outfile = $infile;
	$outfile =~ s/vcf$/dbSNP.vcf/;
}

my $outdir = dirname($outfile);

my $filtereddbSNPs = $infile;
$filtereddbSNPs =~ s/vcf$/filtered_dbSNPs.vcf/;

my $bedFile        = $infile;
$bedFile		   =~ s/vcf$/bed/;

#since the dbSNP file is usually quite big (>50M lines), filter it for snps that are intersecting
#with the snps from $infile
$logger->info("Filtering dbSNP file for variants in $infile...");
my $command = "awk '{if(\$1 !~ /^#/){print \$1\"\\t\"(\$2-1)\"\\t\"\$2 } }' $infile > $bedFile";
$logger->debug("Command: ".$command);
system($command);

$command = "$tabix -h -B $dbsnpfile $bedFile > $filtereddbSNPs";
$logger->debug("Command: ".$command);
system($command);

#the filtered set of dbSNP variants can then be used to do the actual annotation
$command = "$java -jar $snpsift annotate  $filtereddbSNPs $infile > $outfile 2> $outdir/snpsift.log";
$logger->info("Annotating variants with dbSNP variants...");
$logger->debug("Command: ".$command);
system($command);


my $elapsed = ( time - $^T ) / 60;
$logger->info("finished in $elapsed minutes");



=head1 NAME

annotatedbSNP.pl

=head1 SYNOPSIS

annotatedbSNP.pl -i variants.vcf

=head1 DESCRIPTION

This script annotates a VCF file with known SNPs from dbSNP. It needs a VCF file containing all dbSNP variants. Also requires SnpSift
and intersectBed.

=head1 OPTIONS

 -i	<infile.vcf> file containing the variants to annotate; REQUIRED
 -o	<outfile.vcf> output file; default: infile.dbSNP.vcf
 -se	name of the settings in the current.config.xml file that holds path to reference genome, 
 	to the annotation file and to possible additional annotation files; use default settings if nothing is given
 -lf	log file; default: print to screen
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page
 

=head1 AUTHOR

Thomas Wieland

=cut

