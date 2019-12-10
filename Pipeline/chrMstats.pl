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

	#Tool paths
	my $params  = Utilities::getParams();
	my $rscript = $params->{programs}->{Rscript}->{path};
	my $samtools = $params->{programs}->{samtools}->{path};
	
	my $mtcoverageR =$prog_path."/"."chrMstats_depthplot.R";

my $samplename = "mtDNA_Sample";
my $bam      = "";
my $logfile  = "pipeline.log";
my $loglevel = "INFO";
my $outstats = "";
my $outpic	 = "";
my $factor	 = -1;

#program flags
my $help = 0;
my $settings = 0;


GetOptions(
	"b=s"  => \$bam,
	"se=s" => \$settings,
	"o=s"  => \$outstats,
	"op=s"  => \$outpic,
	"n=s"  => \$samplename,
	"lf=s" => \$logfile,
	"ll=s" => \$loglevel,
	"h"    => \$help,
);

Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();

if ( $help == 1 || $bam eq "" || $settings eq "") {
	print
	"\n
	Calculates the coverage of the mtDNA and produces a depth plot with the genes

	Required:
		-b	<bam file> to check
		-se	settings

	Optional:
		-n  <samplename> sample name
		-o	<outstats> stats with the depth per base
		-op	<outpic> name of the .png file for the plot
		-lf	log file; default: pipeline.log
		-ll	log level: ERROR,INFO,DEBUG; default: INFO
		-h	this help
	\n\n";
	exit(1);
}

my $outdir = dirname($bam);
my $ref      = $params->{settings}->{$settings}->{reference};


if($outstats eq ""){
	$outstats = $outdir."/"."chrM.depth.csv";
}

if($outpic eq ""){
	$outpic = $outdir."/"."chrM.depth.png";
}

if( !-e $bam){
	$logger->error("BAM file not found: $bam");
	exit(1);
}

my $filterflags = "--rf 3 --ff 3852 -A"; 
my $mtregion = "chrM:1-16569";

my $fai=$ref.".fai";
if ( ! -e $fai )
{
	$logger->error("Fasta \"$ref\" not indexed");
	exit(1);
}

# TODO Insert here some check with the reference rRCRS

# Reference command chain to use:
#       $samtools mpileup --rf 3 --ff 3852 /PATHTO/S0100/$i/mtDNAout/paired-endout/merged.rmdup.bam  -r $mtRegion | 
#			awk '{print $2" "$4}' | 
#			$mtCoveragePlot "$i - mtDNA Coverage" "Position" "Read Depth" $i.cov.jpg ; 


my $cmd = "$samtools mpileup $filterflags $bam -r $mtregion | awk '{print \$2\" \"\$4}' > $outstats;";
my $ret = system($cmd);

if ( $ret ne "0" )
{
	$logger->error("Pileup for depth calculation failed");
	exit(1);
}


my $cmd2 = "cat $outstats | $rscript $mtcoverageR \"$samplename\" \"Position\" \"Read Depth\" $outpic";
$ret = system($cmd2);

if ( $ret ne "0" )
{
	$logger->error("chrM depth plot generation failed");
	exit(1);
}


# TODO Insert here Haplogroup calculation with haplogrep.

exit(0);



=head1 NAME

chrMstats.pl

=head1 SYNOPSIS

 chrMstats.pl -se <settings> -b <bamfile> [-n <samplename>] [-o <depthfile> ] [ -op <depthplot.png> ]

=head1 DESCRIPTION

This script generates depth profile for mtDNA and its plot with annotated mt-genes. This script
is suited only for the rCRS revised Cambridge Reference Sequence for mtDNA, embedded into hg38

=head1 OPTIONS

-b	<bam file> to check
-se	settings

Optional:
-n  <samplename> sample name
-o	<outstats> stats with the depth per base
-op	<outpic> name of the .png file for the plot
-lf	log file; default: pipeline.log
-ll	log level: ERROR,INFO,DEBUG; default: INFO
-h	this help\n\n";
	exit(1);
=head1 AUTHOR

Riccardo Berutti

=cut
