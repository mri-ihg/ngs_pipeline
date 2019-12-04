#!/usr/bin/perl

use strict;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);
use IO::File;
umask(002);


my $prog_path = dirname( abs_path($0) );
require $prog_path."/Utilities.pm";


my $outdir  = ".";
my $help	= 0;
my $settings= "";
my $logfile  	  = "pipeline.log";
my $loglevel 	  = "INFO";
my $threads = -1;
my $isArrayJob	= 0;
my $cfgFile  = "";
my $region = "ALL";
my $he = -1;
my $ho = -1;

GetOptions(
"o=s" => \$outdir,
"r=s" => \$region,
"se=s"=> \$settings,
"c=s" => \$cfgFile,
"t=s" => \$threads,
"aj"   => \$isArrayJob,
"lf=s" => \$logfile,
"ll=s" => \$loglevel,
"he=f" => \$he,
"ho=f" => \$ho,

"h"	  => \$help
);

if($settings eq "" || $cfgFile eq "" || $help){
	print "
Usage:	runPindel.pl -c pindel.cfg -se <settings>

-c	<pindel.cfg> pindel config file; required
-se	<settings>; required
-o	<outdir> output directory; default: .
-t	<number_of_threads>; overwrites the number of threads given by SGE in environmental variable NSLOTS
-aj	is SGE Array job; if this flag is chosen, the script runs as a part of a SGE array job
	this means, that the environmental variable SGE_TASK_ID holds the number of this job
	The script will process only the chromosome in line SGE_TASK_ID of the dict file of the reference genome.
	It will use the chromsome name as prefix for all generated output files. -t will be set to 1

Parameters forwarded to Pindel:
-r  <region> chromosme to run Pindel on; default: ALL
-he The proportion of reads to call het (Pindel default 0.2)
-ho The proportion of reads to call hom (Pindel default 0.8)

-lf	<log file>; default: pipeline.log
-ll	log level: ERROR,INFO,DEBUG; default: INFO
-h	show this help
";
exit(0);
}



Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();





my $params  = Utilities::getParams();
my $pindel  = $params->{programs}->{pindel}->{path};
my $pindel2vcf = $params->{programs}->{pindel}->{pindel2vcf};
my $ref     = $params->{settings}->{$settings}->{reference};
my $refname = $params->{settings}->{$settings}->{refname};
my $refdate = $params->{settings}->{$settings}->{refdate};
my $vcftools= $params->{programs}->{vcftools}->{path};
my $vcfconcat= $params->{programs}->{vcftools}->{concat};
my $vcflib  = $params->{programs}->{vcftools}->{lib};
my $bedtools= $params->{programs}->{bedtools}->{path};


#get parameters
my $a       = $params->{programs}->{pindel}->{pindelparams}->{A};
my $pr		= $params->{programs}->{pindel}->{pindel2vcfparams}->{pr};
my $is		= $params->{programs}->{pindel}->{pindel2vcfparams}->{is};


my $prefix = "";

if($isArrayJob &&  $ENV{SGE_TASK_ID}){
	
	$threads = 1;
	
	open BED, $params->{settings}->{$settings}->{normalchromosomes} || exit $logger->error("Can't open $params->{settings}->{$settings}->{normalchromosomes}!");
	$. = 0;
	my $line;
	do { $line = <BED> } until $. == $ENV{SGE_TASK_ID} || eof;				#read the line corresponding to this job of the array from the sequence dictionary to extract the chromosome to call
	close BED;
	
	chomp $line;
	my ($chr,$start,$end) = split("\t",$line);
	$region = $chr;
	$prefix = $chr.".";
}
elsif ( $region ne "ALL")
{
	$prefix = $region.".";	
}

if($threads == -1){
	$threads = 1;
	if($ENV{NSLOTS}){		#get number of slots given by SGE
		$threads = $ENV{NSLOTS};
		
	}
}

my $hehooptions="";
if ( $he != -1 )
{
		$hehooptions=" -he $he ";
}

if ( $ho != -1 )
{
		$hehooptions=$hehooptions." -ho $ho ";
}

$logger->debug("region: $region");
$logger->debug("he: $he");
$logger->debug("ho: $ho");


chdir $outdir;

#run pindel
my $command = "$pindel -f $ref -i $cfgFile -c $region -o $outdir/".$prefix."pindel -T $threads -A $a -L $outdir/".$prefix."pindel.log";
$logger->info("Running pindel...");
$logger->debug("pindel command: $command");
system($command);


#convert pindel to vcf																				      
$command = "$pindel2vcf -G -is $is -p $outdir/".$prefix."pindel_D -R $refname -d $refdate -r $ref -co 255 $hehooptions";
$logger->info("Converting pindel file to vcf...");
$logger->debug("pindel2vcf command: $command");
system($command);

$command = "$pindel2vcf -G -is $is -p $outdir/".$prefix."pindel_SI -R $refname -d $refdate -r $ref -co 255 $hehooptions";
$logger->info("Converting pindel file to vcf...");
$logger->debug("pindel2vcf command: $command");
system($command);

$command = "$pindel2vcf -G -is $is -p $outdir/".$prefix."pindel_INV -R $refname -d $refdate -r $ref -co 2 $hehooptions";
$logger->info("Converting pindel file to vcf...");
$logger->debug("pindel2vcf command: $command");
system($command);


$command ="(awk '{if(\$1 ~ /^#/){print}}' $outdir/".$prefix."pindel_D.vcf;(awk '{if(\$1 !~ /^#/ && \$8 !~ /SVTYPE=RPL/){print}}' $outdir/".$prefix."pindel_D.vcf;awk '{if(\$1 !~ /^#/ && \$8 !~ /SVTYPE=RPL/){print}}' $outdir/".$prefix."pindel_SI.vcf;awk '{if(\$1 !~ /^#/ && \$8 !~ /SVTYPE=RPL/){print}}' $outdir/".$prefix."pindel_INV.vcf) ) | perl $prog_path/vcfsorter.pl -se $settings -ll $loglevel -lf $logfile > $outdir/".$prefix."all.pindel.vcf";
$logger->debug("pindel concat & filter command: $command");
system($command);


#system("awk '{if(\$1 ~ /^#/ || \$8 !~ /SVTYPE=RPL/){print} }' $outdir/".$prefix."all.pindel.vcf > $outdir/".$prefix."tmp.vcf");
#system("mv $outdir/".$prefix."tmp.vcf $outdir/".$prefix."all.pindel.vcf");
