#!/usr/bin/perl

use strict;
use Getopt::Long;
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);

my $prog_path = dirname( abs_path($0) );
require $prog_path . "/Utilities.pm";


my $bamfile    = "";
my $outfile    = "";
my $help       = 0;
my $logfile    = "pipeline.log";
my $loglevel   = "INFO";
my $settings   = "";
my $genotype   = "";
my $site = 0;


GetOptions(
	"b=s"  => \$bamfile,
	"o=s"  => \$outfile,
	"g=s"  => \$genotype,
	"s"	   => \$site,
	"se=s" => \$settings,
	"lf=s" => \$logfile,
	"ll=s" => \$loglevel,
	"h"    => \$help
);

if($help || $outfile eq ""){
	print "
This script checks for possible contaminations in a bam file using
the script \"verifyBamID\" (http://genome.sph.umich.edu/wiki/VerifyBamID).
It takes a HapMap file with population specific allele frequencies as
an input from the settings file.

-o	<outfile>; required
-b	<bamfile> to check; default: <outdir>/merged.rmdup.bam
-g	<genotypefile> genotype file in VCF format
-se	<settings> to take a frequency file from if -g is not specified
-lf	log file; default: pipeline.log
-ll	log level: ERROR,INFO,DEBUG; default: INFO
-h	print this help
";
	exit(0);
}

Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();

my $outdir = dirname($outfile);

if($bamfile eq ""){
	$bamfile = $outdir."/merged.rmdup.bam";
}

if(!(-e $bamfile)){
	$logger->error("$bamfile not found!");
	exit(1);
}



my $params	    = Utilities::getParams();
die $logger->error("ExAC file is not specified in settings $settings!") unless $params->{settings}->{$settings}->{exac};
my $freqFile    = $params->{settings}->{$settings}->{exac};
my $ref         = $params->{settings}->{$settings}->{reference};
my $verifyBamID = $params->{programs}->{verifybamid}->{path};

my $command  = "$verifyBamID --precise --maxDepth 999 --out $outfile --bam $bamfile ";
if($genotype ne ""){
	$command .= "--vcf $genotype";
	if ($site) { $command .= " --site" };
}else{
	$command .= "--vcf $freqFile --site";
}
$logger->info ("Running verifyBamID...");
$logger->debug("Command: ".$command);
system($command);
