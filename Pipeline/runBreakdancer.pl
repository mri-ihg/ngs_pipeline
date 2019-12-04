#!/usr/bin/perl

use strict;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;
use Cwd qw(abs_path);
#use Log::Log4perl qw(get_logger :levels);
use IO::File;
use Pod::Usage;
umask(002);

my $prog_path = dirname( abs_path($0) );
require $prog_path."/Utilities.pm";


my $bamfile = "";
my $outdir  = "";
my $help	= 0;
my $man     = 0;
my $logfile  	  = "SCREEN";
my $loglevel 	  = "INFO";
my $chrom   = "";
my $isArrayJob = 0;
my $prefix     = "";
my $settings   = "default";



GetOptions(
"b=s"  => \$bamfile,
"o=s"  => \$outdir,
"c=s"  => \$chrom,
"aj"   => \$isArrayJob,
"se=s"   => \$settings,
"lf=s" => \$logfile,
"ll=s" => \$loglevel,
"man"  => \$man,
"h"	   => \$help
);

pod2usage( {-exitval => 0  ,-verbose => 1} ) if $help;
pod2usage( {-exitval => 0  ,-verbose => 2} ) if $man;
pod2usage( {-exitval => 1  ,-verbose => 1} ) if $bamfile eq "";

Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();



my $params  = Utilities::getParams();
my $breakdancer = $params->{programs}->{breakdancer}->{path};
my $bam2cfg     = $params->{programs}->{breakdancer}->{config};


if($isArrayJob &&  $ENV{SGE_TASK_ID}){
	
	
	open BED, $params->{settings}->{$settings}->{normalchromosomes} || exit $logger->error("Can't open $params->{settings}->{$settings}->{normalchromosomes}!");
	$. = 0;
	my $line;
	do { $line = <BED> } until $. == $ENV{SGE_TASK_ID} || eof;				#read the line corresponding to this job of the array from the sequence dictionary to extract the chromosome to call
	close BED;
	
	chomp $line;
	my ($chr,$start,$end) = split("\t",$line);
	$chrom  = $chr;
	$prefix = $chr.".";
	
	
}

$outdir = dirname($bamfile) if $outdir eq "";

chdir $outdir;

#generate config file
my $cfgFile = $prefix."breakdancer.cfg";
my $command = "perl $bam2cfg $bamfile > $cfgFile";
$logger->info("Running breakdancer bam2cfg...");
$logger->debug("breakdancer command: $command");
system($command);

#run breakdancer
my $outfile = $cfgFile;
$outfile    =~ s/cfg$/out/;
$chrom = " -o $chrom " if $chrom ne "";
$command = "$breakdancer $chrom $cfgFile > $outfile";
$logger->info("Running breakdancer-max...");
$logger->debug("breakdancer command: $command");
system($command);

#convert breakdancer output to VCF
my $vcffile = $outfile;
$vcffile    =~ s/out$/vcf/;
$command = "perl $prog_path/breakdancer2VCF.pl -i $outfile -o $vcffile -lf $logfile -ll $loglevel";
$logger->info("Converting breakdancer2VCF ...");
$logger->debug("breakdancer command: $command");
system($command);

=head1 NAME

runBreakdancer.pl

=head1 SYNOPSIS

 runBreakdancer.pl -b merged.rmdup.bam

=head1 DESCRIPTION

This script is a wrapper script to run breakdancer on WG datasets. 

=head1 OPTIONS

 -b	<bamfile> to run breakdancer on; required
 -o	outdir; default: directory where bamfile lies in
 -c	chromosome to run breakdancer on; default: all chromosomes
 -aj	is SGE Array job; if this flag is chosen, the script runs as a part of a SGE array job
	this means, that the environmental variable SGE_TASK_ID holds the number of this job
	The script will process only the chromosome in line SGE_TASK_ID of the "normal chromosomes"
	file specified in the config file.
 -se	name of the settings in the current.config.xml file that holds path to reference genome, 
 	to the annotation file and to possible additional annotation files; use default settings if nothing is given
 -lf	<log file>; default: print to screen
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page

=head1 AUTHOR

Thomas Wieland

=cut