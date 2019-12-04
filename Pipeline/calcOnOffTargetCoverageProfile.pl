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
use Bio::DB::Sam;
use List::Util qw(sum);
use Pod::Usage;

my $prog_path = dirname( abs_path($0) );
require $prog_path . "/Utilities.pm";

	# Tool paths
	my $params  		= Utilities::getParams();
	my $rscript 		= $params->{programs}->{Rscript}->{path};
	my $samtools 		= $params->{programs}->{samtools}->{path};
	my $bedtools		= $params->{programs}->{bedtools}->{path}."/bedtools";
	
	# Binner
	my $binner = "perl " . $prog_path . "/calcOnOffTargetCoverageProfileBinner.pl";
	# PlotR
	my $plotR  = $prog_path . "/calcOnOffTargetCoverageProfilePlot.R";


my $samplename = "Sample";
my $bam      = "";
my $logfile  = "pipeline.log";
my $loglevel = "INFO";
my $outstats = "";
	my $outstats_ontarget 	= "";
	my $outstats_offtarget  = "";
my $region   = "chr21";
my $plotonly = 0;
my $outpic	 = "";
my $factor	 = -1;

#program flags
my $help = 0;
my $settings = 0;


GetOptions(
	"b=s"  => \$bam,
	"se=s" => \$settings,
	"o=s"  => \$outstats,
	"op=s" => \$outpic,
	"n=s"  => \$samplename,
	"r=s"  => \$region,
	"P"	   => \$plotonly,
	"lf=s" => \$logfile,
	"ll=s" => \$loglevel,
	"h"    => \$help,
);

Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();

if ( $help == 1 || $bam eq "" || $settings eq "") {
	print
	"\n
	Calculates the depth profile for the samples

	Required:
		-b	<bam file> to check
		-se	settings

	Optional:
		-n  <samplename> sample name
		-o	<out> stats with the depth plot
		    will be set out.ontarget out.offtarget 
		-op	<outpic> name of the .png file for the combo plot
		-P  regenerate only the plot
		-lf	log file; default: pipeline.log
		-ll	log level: ERROR,INFO,DEBUG; default: INFO
		-h	this help
	\n\n";
	exit(1);
}

# Outdir, ref and coding region
my $outdir 			= dirname($bam);
my $ref             = $params->{settings}->{$settings}->{reference};
my $codingregion	= $params->{settings}->{$settings}->{codingregions};

if($outstats eq ""){
	$outstats = $outdir."/"."coverageprofile";
}
	$outstats_ontarget 	= $outstats."_ontarget";
	$outstats_offtarget = $outstats."_offtarget";

if($outpic eq ""){
	$outpic = $outdir."/"."coverageprofile_onofftarget.png";
}

if( !-e $bam){
	$logger->error("BAM file not found: $bam");
	exit(1);
}

# Filterflags good reads
my $regionOption = ( Utilities::validateChrPosFull($region) ? "-r $region" : "" ); #
my $filterflags = "--rf 3 --ff 3852"; 
my $BLOCKAGEOPTION = "";	# IE " | head -n 1000000"  # Debug

# Return status
my $ret = 0;

# ON TARGET (INTERSECT -u option bedtools)
my $cmd  = "$samtools mpileup $filterflags -f $ref $regionOption $bam $BLOCKAGEOPTION | awk '{print \$1\"\t\"\$2\"\t\"\$2+1\"\t\"\$4}' | $bedtools intersect -u -a - -b $codingregion | awk '{print \$4}' | $binner > $outstats_ontarget";
	if ( $plotonly eq "0" )
	{
		$ret = system($cmd);
		if ( $ret ne "0" )
		{
			$logger->error("Pileup for ontarget depth calculation failed");
			exit(1);
		}
	}

# OFF TARGET (REVERSE INTERSECT -v option bedtools)
my $cmd2 = "$samtools mpileup $filterflags -f $ref $regionOption $bam $BLOCKAGEOPTION | awk '{print \$1\"\t\"\$2\"\t\"\$2+1\"\t\"\$4}' | $bedtools intersect -v -a - -b $codingregion | awk '{print \$4}' | $binner > $outstats_offtarget";
	if ( $plotonly eq "0" )
	{
		$ret = system($cmd2);
		if ( $ret ne "0" )
		{
			$logger->error("Pileup for offtarget depth calculation failed");
			exit(1);
		}
	}

# Generate Plot
# 	Rscript doublePlotXY.R "PLOTTITLE" XMIN XMAX "XAXISLABEL" YMIN YMAX "YAXISLABEL" DATASERIESXY1.dat DATASERIESXY2.dat OUTFILE.png
my $cmd3 = "$rscript $plotR \"$samplename - On/Off Target Coverage Profile\" 0 100 \"Coverage X\" 0 1 \"Fraction covered\" $outstats_ontarget $outstats_offtarget $outpic";
	$ret = system($cmd3);
	if ( $ret ne "0" )
	{
		$logger->error("Plot generation failed");
		exit(1);
	} 

exit(0);



=head1 NAME

calcOnOffTargetCoverageProfile.pl

=head1 SYNOPSIS

 calcOnOffTargetCoverageProfile.pl -se <settings> -b <bamfile> [-n <samplename> ] [-o <outstats_basename> ] [ -op <plot.png> ]

=head1 DESCRIPTION

This script generates coverage profile curves differentiated for on and off target regions. On and off
target regions will be selected from settings.

=head1 OPTIONS

-b	<bam file> to check
-se	settings

Optional:
-n  <samplename> sample name
-o	<outstats> stats basename
-op	<outpic> name of the .png file for the plot
-P  regenerate just the plot
-lf	log file; default: pipeline.log
-ll	log level: ERROR,INFO,DEBUG; default: INFO
-h	this help\n\n";
	exit(1);
=head1 AUTHOR

Riccardo Berutti

=cut
