#!/usr/bin/perl 

use strict;
use Getopt::Long;
use DBI;
use warnings;
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);

my $prog_path = dirname(abs_path($0));
	require $prog_path."/Utilities.pm";
	require $prog_path."/CollectMetrics.pm";

my $sample   = "";
my $logfile  = "pipeline.log";
my $loglevel = "INFO";
my $help	 = 0;
my $bamfile  = "";
my $libtype  = ""; 
my $settings = "";

my $stats_run = 0;
my $ret = 0;

# Options:
# 
my $align    	= 0;
my $basedistrib = 0;
my $gc		 	= 0;
my $isz			= 0;
my $qualbycycle = 0;
my $mumetrics	= 0;
my $qualdistr   = 0;
my $wgs			= 0;
my $libcomp		= 0;
my $test		= 0;
my $getISZ		= 0;
my $getLibcomp  = 0;
my $getLibcompN = 0;
my $getQFract   = 0;
my $getQ30Fract = 0;

my $getAlignmentMetric="";
my $getBaseDistribution=0;
my $getGCMetric = "";
my $getMeanBaseQual = 0;
my $getMeanBaseQualLastN = 0;
my $getMeanBaseQualFull = 0;
my $getRmDupMetric = 0;

my $helptext = "
	collectMetrics.pl
	
	Runs Picard tools to collect metrics for the given bam file
	
	  Parameters:
		-b	 </path/to/file.bam> - bam file for which the metrics are calculated; the directory of this file is used as output directory; required
		-se  <settings>
	  Metrics to collect: to run it in a parallel fashion, run multiple times each with a different statistics tool
		-a	 calculate alignment metrics
		-bd  calculate base distribution by cycle
		-g	 calculate gc metrics
		-i   calculate insert size distribution
		-qcy calculate mean quality by cycle 
		-mm  calculate multiple metrics (ie pct Q20/Q30)
		-qd  calculate quality score distribution
		-w   calculate wgs metrics
		-e	 estimate library complexity
	  Logging:
		-lf	 log file; default: pipeline.log
		-ll	 log level: ERROR,INFO,DEBUG; default: INFO
	  Help:	
		-h	 this help\n";

GetOptions(
	"b=s"	=> \$bamfile,
	"se=s"  => \$settings,
	"a"		=> \$align,
	"bd"	=> \$basedistrib,
	"g"		=> \$gc,
	"i"		=> \$isz,
	"qcy"   => \$qualbycycle,
	"mm"	=> \$mumetrics,
	"qd"    => \$qualdistr,
	"w"		=> \$wgs,
	"e"		=> \$libcomp,
	"t"		=> \$test,
	"getISZ"=> \$getISZ,
	"getLC" => \$getLibcomp,
	"getLCN" => \$getLibcompN,
	"getQ30"=> \$getQ30Fract,
	"getQ"  => \$getQFract,
	"get_a=s" => \$getAlignmentMetric,
	"get_bd" => \$getBaseDistribution,
	"get_g=s" => \$getGCMetric,
	"get_qcy" => \$getMeanBaseQual,
	"get_qcyn=i" => \$getMeanBaseQualLastN,
	"get_qcyf"=> \$getMeanBaseQualFull,
	"get_dup"=> \$getRmDupMetric,
	"lf=s" 	=> \$logfile,
	"ll=s" 	=> \$loglevel,
	"h"	   	=> \$help
);

if ( $help == 1 ||  $bamfile eq "") {
	print $helptext;
	exit(1);
}

my $params = Utilities::getParams();

Utilities::initLogger( $logfile, $loglevel );
my $logger = Utilities::getLogger();

my $outdir = dirname($bamfile);

# Select which statistical tool to run. 
# for uniformity settings is passed to all the functions - even when it's not strictly needed

#CollectAlignmentSummaryMetrics
if($align){
	$logger->info("Running Picard CollectAlignmentSummaryMetrics...");
	$ret = CollectMetrics::calcAlignmentMetrics($bamfile,$outdir,$settings);
	$logger->info("Finished CollectAlignmentSummaryMetrics!");
	$stats_run++
}

#CollectBaseDistributionByCycle
if($basedistrib){
	$logger->info("Running Picard CollectBaseDistributionByCycle...");
	$ret = CollectMetrics::calcBaseDistributionByCycle($bamfile,$outdir,$settings);
	$logger->info("Finished CollectBaseDistributionByCycle!");
	$stats_run++
}

#CollectGcBiasMetrics
if($gc){
	$logger->info("Running Picard CollectGcBiasMetrics...");
	$ret = CollectMetrics::calcGCMetrics($bamfile,$outdir,$settings);
	$logger->info("Finished CollectGcBiasMetrics!");
	$stats_run++
}

#CollectInsertSizeMetrics
if($isz){
	$logger->info("Running Picard CollectInsertSizeMetrics...");
	$ret = CollectMetrics::calcInsertSizeMetrics($bamfile,$outdir,$settings);
	$logger->info("Finished CollectInsertSizeMetrics!");
	$stats_run++
}

#MeanQualityByCycle
if($qualbycycle){
	$logger->info("Running Picard MeanQualityByCycle...");
	$ret = CollectMetrics::calcMeanQualityByCycle($bamfile,$outdir,$settings);
	$logger->info("Finished MeanQualityByCycle!".withErr($ret));
	$stats_run++
}

#CollectMultipleMetrics
if($mumetrics){
	$logger->info("Running Picard CollectMultipleMetrics...");
	$ret = CollectMetrics::calcMultipleMetrics($bamfile,$outdir,$settings);
	$logger->info("Finished CollectMultipleMetrics!");
	$stats_run++
}

#QualityScoreDistribution
if($qualdistr){
	$logger->info("Running Picard QualityScoreDistribution...");
	$ret = CollectMetrics::calcQualityScoreDistribution($bamfile,$outdir,$settings);
	$logger->info("Finished QualityScoreDistribution!");
	$stats_run++
}

#CollectWgsMetrics
if ($wgs)
{
	$logger->info("Running Picard CollectWgsMetrics...");
	$ret = CollectMetrics::calcWgsMetrics($bamfile,$outdir,$settings);
	$logger->info("Finished CollectWgsMetrics!");
	$stats_run++
}

#EstimateLibraryComplexity
if($libcomp){
	$logger->info("Running Picard EstimateLibraryComplexity...");
	$ret = CollectMetrics::calcEstimateLibraryComplexity($bamfile,$outdir);
	$logger->info("Finished EstimateLibraryComplexity!");
	$stats_run++
}

# Check that at least one tool was selected or return a warning
if ( $stats_run == 0 && $test != 1 ){
	$logger->info("collectMetrics.pl run with 0 statistics selected!");
}

#Test routine
if ($test){
	my @isz=CollectMetrics::getInsertSize($bamfile, $outdir);
	my $libs=join($isz[2],",");
	print( "InsertSize:\t\t".$isz[0]."+/-".$isz[1]." - Lib:\t".$libs."\n");
	my $libCplx=CollectMetrics::estimateLibraryComplexity($bamfile, $outdir);
	print( "LibComplex:\t\t".$libCplx."\n" );
}

# These following function most for test purposes or for third party apps requiring a sample statistic

# Get Parameters [ISZ]
if ($getISZ){
	my @iszItems=CollectMetrics::getInsertSize($bamfile, $outdir);
	print ($iszItems[0]."\t".$iszItems[1]."\n");
}

# Get Parameters [LibComplexity]
if ($getLibcomp){
	my $libcomplexity=CollectMetrics::estimateLibraryComplexity($bamfile, $outdir);
	print ($libcomplexity."\n");
}

# Get Parameters [LibComplexityN] Normalised to genome size
if ($getLibcompN){
	my $libcomplexityN=CollectMetrics::estimateLibraryComplexityNormalised($bamfile, $outdir, $settings);
	print ($libcomplexityN."\n");
}

# Get Parameters [Q102030]
if ($getQFract){
	my @qfract=CollectMetrics::getQFraction($bamfile, $outdir);
	print ($qfract[0]."\t".$qfract[1]."\t".$qfract[2]."\n");
}

# Get Parameters [Q30Fraction]
if ($getQ30Fract){
	my $q30fract=CollectMetrics::getQ30Fraction($bamfile, $outdir);
	print ($q30fract."\n");
}

# Get AlignmentMetric [PAR]
if ($getAlignmentMetric ne "" && $getAlignmentMetric ne "ALL"){
	my $am=CollectMetrics::getAlignmentMetric($bamfile, $outdir, $settings, $getAlignmentMetric);
	print ($getAlignmentMetric."\t".$am."\n");
}

# Get AlignmentMetric [PAR]
if ($getAlignmentMetric eq "ALL"){
	my %data=CollectMetrics::getAlignmentMetric($bamfile, $outdir, $settings);
	foreach my $key (keys %data)
	{
		if ( defined $data{$key}){
			print $key."\t".$data{$key}."\n";
		}
	}
}

# Get Base Distribution
if ( $getBaseDistribution ){
	my %data=CollectMetrics::getBaseDistribution($bamfile, $outdir, $settings);
	foreach my $key (sort keys %data)
	{
		if ( defined $data{$key}){
			print $key."\t".$data{$key}."\n";
		}
	}
	
	my $bc1=CollectMetrics::getBaseDistributionBadCycles($bamfile, $outdir, $settings, 1);
	my $bc2=CollectMetrics::getBaseDistributionBadCycles($bamfile, $outdir, $settings, 2);
	my $bc3=CollectMetrics::getBaseDistributionBadCycles($bamfile, $outdir, $settings, 3);
	my $bc4=CollectMetrics::getBaseDistributionBadCycles($bamfile, $outdir, $settings, 4);
	
	print "Cycles with imbalance >n sigmas\n";
	print "1\t2\t3\t4\n";
	print $bc1."\t".$bc2."\t".$bc3."\t".$bc4."\n";
}

# Get GCMetric [PAR]
if ($getGCMetric ne "" && $getGCMetric ne "ALL"){
	my $am=CollectMetrics::getGCMetric($bamfile, $outdir, $settings, $getGCMetric);
	print ($getGCMetric."\t".$am."\n");
}

# Get GCMetric [PAR]
if ($getGCMetric eq "ALL"){
	my %data=CollectMetrics::getGCMetric($bamfile, $outdir, $settings);
	foreach my $key (keys %data)
	{
		if ( defined $data{$key}){
			print $key."\t".$data{$key}."\n";
		}
	}
}

# Get BQ 
if ($getMeanBaseQual){
	my %data=CollectMetrics::getMeanQuality($bamfile, $outdir, $settings);
	print ($data{"Q"}."\t".$data{"QSD"}."\n");
}

# Get BQ LastN 
if ($getMeanBaseQualLastN){
	my %data=CollectMetrics::getMeanQualityLastN($bamfile, $outdir, $settings, $getMeanBaseQualLastN);
	print ($data{"Q"}."\t".$data{"QSD"}."\n");
}

# Get BQFull
if ($getMeanBaseQualFull){
	my %data=CollectMetrics::getMeanQuality($bamfile, $outdir, $settings);
	print ($data{"Q"}."\t".$data{"QSD"}."\n");
	my $bc1=CollectMetrics::getMeanQualityBadCycles($bamfile, $outdir, $settings, 1);
	my $bc2=CollectMetrics::getMeanQualityBadCycles($bamfile, $outdir, $settings, 2);
	my $bc3=CollectMetrics::getMeanQualityBadCycles($bamfile, $outdir, $settings, 3);
	my $bc4=CollectMetrics::getMeanQualityBadCycles($bamfile, $outdir, $settings, 4);
	
	print "Cycles with Qdiff >n sigmas from average\n";
	print "1\t2\t3\t4\n";
	print $bc1."\t".$bc2."\t".$bc3."\t".$bc4."\n";	
}

# Get Rmdup
if ($getRmDupMetric){
	if ( CollectMetrics::getDuplicates($bamfile, $outdir, $settings) ne "NULL")
	{
		my %data=CollectMetrics::getDuplicates($bamfile, $outdir, $settings);
		print ("DUP: ".$data{"DUPLICATION_RATE"}."  OPTDUP:".$data{"OPTICAL_DUPLICATION_RATE"}."\n");
	}
	else
	{
		print ("DUPS NOT REMOVED\n");
	}
}

# Make errstring
sub withErr
{
   my $ret=shift;
 return ( ($ret != "0") ? " WITH ERRORS" : "") 
}
