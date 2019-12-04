#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;
use Getopt::Long qw(:config no_ignore_case);
use Cwd qw(abs_path);
use File::Basename;
use DateTime;
use Pod::Usage;
umask(002);

my $prog_path = dirname( abs_path($0) );
require $prog_path."/Utilities.pm";

my $infile     = "";
my $outfile    = "";
my $targetfile = "";
my $help	   = 0;
my $man		   = 0;
my $logfile    = "pipeline.log";
my $loglevel   = "INFO";
my $fraction   = "";
my $minV	   = 0;
my $minU	   = 0;
my $getMinimal = 0;
my $settings   = "";




	GetOptions(
		"i=s"  => \$infile,
		"o=s"  => \$outfile,
		"t=s"  => \$targetfile,
		"f=s"  => \$fraction,
		"v"	   => \$minV,
		"u"    => \$minU,
		"m"    => \$getMinimal,
		"h"    => \$help,
		"se=s" => \$settings,
		"lf=s" => \$logfile,
		"ll=s" => \$loglevel,
		"man"  => \$man
	);

pod2usage( {-exitval => 0  ,-verbose => 1} ) if $help;
pod2usage( {-exitval => 0  ,-verbose => 2} ) if $man;
pod2usage( {-exitval => 1  ,-verbose => 1} ) if $infile eq "" || $outfile eq "";

Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();

my $params  = Utilities::getParams();
my $bedtools= $params->{programs}->{bedtools}->{path};

my $command; 

unlink($outfile);
if($targetfile ne ""){
	my $options = "";
	$options   .= " -f $fraction -r" if $fraction ne "";
	$options   .= " -v" if $minV;
	$options   .= " -u" if $minU;
	# $targetfile .= " " . $params->{settings}->{$settings}->{miRNAregion} . " -u" if ($settings eq "hg19_plus"); #for now, keep miRNA variants only if hg19_plus is used
	my $catInfile = "cat $infile";
	$catInfile = "perl $prog_path/getMinimalIndels.pl -i $infile" if $getMinimal;
	
	$command = "(grep \"#\" $infile; $catInfile | $bedtools/intersectBed -a stdin -b $targetfile -wa $options) > $outfile";
	$logger->info("Filtering VCF file for target region...");
	$logger->debug("VCF filter command: $command");
}else{
	$command = "ln -s ".basename($infile)." $outfile";
	$logger->info("Creating VCF file link...");
	$logger->debug("VCF filter command: $command");
}


system($command);



=head1 NAME

filterVCFforRegion.pl

=head1 SYNOPSIS

 filterVCFforRegion.pl -i all.vcf -o filtered.vcf -t target.bed

=head1 DESCRIPTION

This script filters a VCF file for a given region using BEDtools.

=head1 OPTIONS

 -i	<infile.vcf>; REQUIRED
 -o	<outfile.vcf>; REQUIRED
 -t <target.bed>; if no region file is given, a link will be created
 -f	<fraction>; fraction of variants that have to overlap in order to be included
 -v include variants that are NOT in the target file
 -u	include variants only once
 -m	minimize indels in -i before checking the overlap (usefull to compare e.g. Pindel to GATK VCFs)
 -se <setting> settings 
 -lf	<log file>; default: pipeline.log
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page

=head1 AUTHOR

Thomas Wieland

=cut
