#!/usr/bin/perl

use strict;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;
use File::Copy;
use Cwd qw(abs_path);
use Data::Dumper;
use IO::File;
use Pod::Usage;
umask(002);

my $prog_path = dirname( abs_path($0) );
require $prog_path."/Utilities.pm";


my $bamfile = "";
my $outdir  = "";
my $help	= 0;
my $man     = 0;
my $logfile  	  = "pipeline.log";
my $loglevel 	  = "INFO";
my $settings  = "";
my $threads       = -1;



GetOptions(
"b=s"  => \$bamfile,
"o=s"  => \$outdir,
"t=s"  => \$threads,
"se=s" => \$settings,
"lf=s" => \$logfile,
"ll=s" => \$loglevel,
"man"  => \$man,
"h"	   => \$help
);

pod2usage( {-exitval => 0  ,-verbose => 1} ) if $help;
pod2usage( {-exitval => 0  ,-verbose => 2} ) if $man;
pod2usage( {-exitval => 1  ,-verbose => 1} ) if $bamfile eq "" || $outdir eq "" || $settings eq "" ;

Utilities::initLogger($logfile,$loglevel);

my $logger = Utilities::getLogger();

$outdir = dirname($bamfile) if $outdir eq "";

my $params 		  = Utilities::getParams();
my $manta  		  = $params->{programs}->{manta}->{path};
my $normalChrs  	  = $params->{settings}->{$settings}->{normalchromosomes};
my $ref		 	  = $params->{settings}->{$settings}->{reference};

my $discordantBAM = $bamfile;
$discordantBAM    =~ s/bam$/discordant/;

my $splitBAM      = $bamfile;
$splitBAM         =~ s/bam$/splitters/;

if($threads == -1){
	$threads = 1;
	if($ENV{NSLOTS}){		#get number of slots given by SGE
		$threads = $ENV{NSLOTS};
		
	}
}


#creating config file
my $command = "$manta --bam=$bamfile --referenceFasta=$ref --runDir=$outdir";

#adding normal chromosomes as region
open CHROMS, $normalChrs or exit $logger->error("Can't open normal chromosome BED file: $normalChrs!");
while(<CHROMS>){
	chomp;
	my @columns = split;
	$command .= " --region=$columns[0]";
}

$logger->info("Creating manta config...");
$logger->debug("Command: ".$command);
system($command);


#running manta
$command = "$outdir/runWorkflow.py -j $threads -m local";
$logger->info("Running manta...");
$logger->debug("Command: ".$command);
system($command);





=head1 NAME

runWhamg.pl

=head1 SYNOPSIS

 runWhamg.pl -b merged.rmdup.bam -o mantadir/ -se hg19_plus

=head1 DESCRIPTION

This script is a wrapper script to run Whamg
Only tested for whole genome data. It calls variants only on chromosomes 1-22,X,Y and M.

=head1 OPTIONS

 -b	<bamfile> to run cnvnator on; required
 -o	outdir; required
 -se settings; required
 -t	<threads>, default: take number of slots from SGE environment variable NSLOTS or 1 if no SGE job
 -lf	<log file>; default: pipeline.log
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page

=head1 AUTHOR

Thomas Wieland

=cut
