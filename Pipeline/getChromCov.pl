#!/usr/bin/perl

use strict;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);
use Pod::Usage;

#include Utilities.pm
my $prog_path = dirname(abs_path($0));
require $prog_path."/Utilities.pm";







my $help          = 0;
my $settings      = "";
my $infile        = "";
my $outfile       = "";
my $logfile  	  = "pipeline.log";
my $loglevel 	  = "INFO";
my $man				= 0;




GetOptions(
	"i=s"  => \$infile, 
	"o=s"  => \$outfile,
	"man"  => \$man, 
	"h"    => \$help, 
	"lf=s" => \$logfile,
	"ll=s" => \$loglevel, 
	"se=s" => \$settings);
	
pod2usage( {-exitval => 0  ,-verbose => 1} ) if $help;
pod2usage( {-exitval => 0  ,-verbose => 2} ) if $man;
pod2usage( {-exitval => 1  ,-verbose => 1} ) if $settings eq "" || $infile eq "";


Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();
my $params = Utilities::getParams();

	
my $samtools = $params->{programs}->{samtools}->{path};		
my $ref      = $params->{settings}->{$settings}->{reference};


if($outfile ne ""){
	close STDOUT;
	if($outfile =~ /\.gz$/){
		open(STDOUT,"| bgzip -c > $outfile") or die "Can't open $outfile for writing!\n";
	}else{
		open(STDOUT,">",$outfile) or die "Can't open $outfile for writing!\n";
	}
	
}

open FAI, $ref.".fai" or exit $logger->error("Can't open $ref.fai!");

while(<FAI>){
	my ($chrom,$length,@rest) = split;
	my $depth = `$samtools depth -r $chrom:0-$length $infile | awk '{sum+=\$3}END{print sum/NR}'`;
	print $chrom."\t".$depth;
}
close FAI;
close STDOUT;




=head1 NAME

 getChromCov.pl

=head1 SYNOPSIS

 getChromCov.pl -i merged.rmdup.bam -o cov_per_chrom.txt -se hg19

=head1 DESCRIPTION

This script calculates the average coverage of each chromosome from a BAM file.
It takes the list and length of chromosomes from the fasta index file of the given reference genome

=head1 OPTIONS

 -i	<infile.bam>; REQUIRED
 -o	<outfile.txt>; default: STDOUT
 -se	settings; REQUIRED
 -lf	log file; default: pipeline.log
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -man	show man page
 -h	this help

=head1 AUTHOR

Thomas Wieland

=cut

