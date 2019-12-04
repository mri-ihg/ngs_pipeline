#!/usr/bin/perl

use strict;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);
use Pod::Usage;
umask(002);


my $prog_path = dirname(abs_path($0));
require $prog_path."/Utilities.pm";
require $prog_path."/CollectMetrics.pm";

my $outfile    = "";
my $help       = 0;
my $man		   = 0;
my $params     = Utilities::getParams();
my $logfile    = "pipeline.log";
my $loglevel   = "INFO";
my $sample     = "";
my $infile     = "";


GetOptions(
"i=s"  => \$infile,
"o=s"  => \$outfile, 
"s=s"  => \$sample,
"lf=s" => \$logfile,
"ll=s" => \$loglevel,
"man"  => \$man,
"h"    => \$help);


pod2usage( {-exitval => 0  ,-verbose => 1} ) if $help;
pod2usage( {-exitval => 0  ,-verbose => 2} ) if $man;
pod2usage( {-exitval => 1  ,-verbose => 1} ) if $outfile eq "" || $sample eq "" || $infile eq "";

Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();

my @infiles;
if(-d $infile){
	@infiles = glob("$infile/*.sort.bam");
	if(@infiles == 0){
		$logger->error("Directory without *.sort.bam files given!");
		exit(1);
	}
}else{
	push(@infiles,$infile);
}

my ($insertSize,$dummy,%dummy2) = CollectMetrics::getInsertSize($infiles[0],dirname($outfile));
open(CFG,">$outfile") or exit $logger->error("Can't open $outfile for writing!");
foreach(@infiles){
	print CFG "$_\t$insertSize\t$sample\n";
}
close CFG;


=head1 NAME

createPindelCfg.pl

=head1 SYNOPSIS

 createPindelCfg.pl -i merged.rmdup.bam -o pindel.cfg -s SAMPLENAME

=head1 DESCRIPTION

This is script creates a config file for pindel. It writes the BAM file name and sample name into the file and
gets the insert size using CollectMetrics::getInsertSize

=head1 OPTIONS

 -i	merged.rmdup.bam; if a directory is given: use all *.sort.bam files. Only one file is then used to get the insert size; required
 -o	<outfile>; required
 -s	samplename; required
 -lf	<log file>; default: pipeline.log
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page

=head1 AUTHOR

Thomas Wieland

=cut