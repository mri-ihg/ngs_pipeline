#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;
use Getopt::Long qw(:config no_ignore_case);
use Cwd qw(abs_path);
use File::Basename;
use Pod::Usage;
use Scalar::Util;
umask(002);

my $start_time = time();

my $prog_path = dirname( abs_path($0) );
require $prog_path."/Utilities.pm";

my $chipfile        = "";
my $peakfile        = "";
my $inputsample     = "";
my $replicatesample = "";
my $replicatepeaks  = "";
my $window		    = 0;
my $outfile         = "";
my $settings        = "";
my $sppFile         = "";
my $help	        = 0;
my $man		        = 0;
my $logfile         = "SCREEN";
my $loglevel        = "INFO";


GetOptions(
	"c=s"  => \$chipfile,
	"p=s"  => \$peakfile,
	"is=s" => \$inputsample,
	"rs=s" => \$replicatesample,
	"rp=s" => \$replicatepeaks,
	"s=s"  => \$sppFile,
	"w=s"  => \$window,
	"o=s"  => \$outfile,
	"se=s" => \$settings,
	"h"    => \$help,
	"lf=s" => \$logfile,
	"ll=s" => \$loglevel,
	"man"  => \$man
);

pod2usage( {-exitval => 0  ,-verbose => 1} ) if $help;
pod2usage( {-exitval => 0  ,-verbose => 2} ) if $man;
pod2usage( {-exitval => 1  ,-verbose => 1} ) if $chipfile eq "" || $settings eq "" || $peakfile eq "" || (($replicatesample ne "" && $replicatepeaks eq "") || ($replicatesample eq "" && $replicatepeaks ne ""));

Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();

my $params   = Utilities::getParams();
my $bedtools = $params->{programs}->{bedtools}->{path};
my $samtools = $params->{programs}->{samtools}->{path};
my $coding   = $params->{settings}->{$settings}->{codingregions};


#open output file, if required
open(OUT,">",$outfile) ||  exit $logger->info("Can't open $outfile for writing!");

#get sample name from BAM file
open BAM, "$samtools view -H $chipfile |" || exit $logger->error("Could not open $samtools view -H $chipfile");
my $sample = "";
while(<BAM>){
	if($_ =~ /^\@/){
		if($_ =~ /^\@RG/){
			my @columns = split();
			foreach my $field(@columns){
				if($field =~ /SM:(.+)/){
					$sample = $1;
					last;
				}
			}
		}
	}
}
close BAM;

#calculate statistics
print OUT "#ChIP sample\tInput sample\tUniquely mapped reads\tPeaks\tNRF(>0.8)\tFRiP(>0.01)\tPeaks overlapping coding regions\tReplicate sample\tPeaks overlapping replicate sample\n";
	
$logger->info("Calculating NRF...");
open NRF,"$samtools view $chipfile | cut -f 3,4  | uniq -c | awk '{sum+=\$1}END{print sum\"\t\"(NR/sum)}' |" or exit $logger->error("Can't open $chipfile!");
my $tmp = <NRF>;
chomp $tmp;
close NRF;
my ($reads,$nrf) = split("\t",$tmp);

$logger->info("Calculating Number of Peaks...");
open PEAKS, "wc -l $peakfile | cut -d \" \" -f 1 |" or exit $logger->error("Can't open $peakfile!");
my $peaks = <PEAKS>;
chomp $peaks;
close PEAKS;
	
$logger->info("Calculating FRiP...");
open FRIP, "$bedtools/intersectBed -abam $chipfile -b $peakfile -u -wa | wc -l |" or exit $logger->error("Can't open $peakfile or $chipfile!");
my $frip = <FRIP>;
chomp $frip;
$frip = $frip/$reads;
close FRIP;
	
$logger->info("Calculating proportion of peeks in coding regions...");
open CREG, "$bedtools/intersectBed -a $peakfile -b $coding -u  | wc -l |" or exit $logger->error("Can't open $peakfile or $coding!");
my $creg = <CREG>;
chomp $creg;
$creg = $creg/$peaks;
close CREG;

my $rep = "";
if($replicatepeaks ne ""){
	$logger->info("Calculating proportion of peeks overlaping replicate...");
	open REP, "$bedtools/windowBed -a $peakfile -b $replicatepeaks -u -w $window | wc -l |" or exit $logger->error("Can't open $peakfile or $replicatepeaks!");
	$rep = <REP>;
	chomp $rep;
	$rep = $rep/$peaks;
	close REP;
}

print OUT "$sample\t$inputsample\t$reads\t$peaks\t$nrf\t$frip\t$creg\t$replicatesample\t$rep\t";

## read spp outputfile
open(SPP, "$sppFile") or exit $logger->error("Can't open $sppFile");
while (<SPP>) {
	chomp;
	$_ =~ s/\%//g;
	
	#format:Filename<tab>numReads<tab>estFragLen<tab>corr_estFragLen<tab>PhantomPeak<tab>corr_phantomPeak<tab>argmin_corr<tab>min_corr<tab>Normalized SCC (NSC)<tab>Relative SCC (RSC)<tab>QualityTag
	my ($filename,$numReads,$estFragLen,$corr_estFragLen,$phantomPeak,$corr_phantomPeak,$argmin_corr,$rmin_corr,$normSCC,$relSCC, $qualityTag) = split("\t",$_);
	print OUT "$filename\t$numReads\t$estFragLen\t$corr_estFragLen\t$phantomPeak\t$corr_phantomPeak\t$argmin_corr\t$rmin_corr\t$normSCC\t$relSCC\t$qualityTag\n";
}
close SPP;
close OUT;

my $end_time = time();
$logger->info("calcChIPstats.pl finished in " . &Utilities::seconds_to_ddhhmmss($end_time - $start_time) . " (ddd:hh:mm:ss)");



=head1 NAME

calcChIPstats.pl

=head1 SYNOPSIS

calcChIPstats.pl -c chip.bam -p peaks.bed -is INPUTSAMPLE -rs REPLICATESAMPLE -rp replicate_peaks.bed
 
=head1 DESCRIPTION

This script calculates quality metrics for a ChIP-seq experiment:
 -) The number of uniquely mapped reads
 -) The number of called peaks
 -) The nonredundant fraction (NRF) of reads (1); is roughly 1-duplicate rate
 -) The fraction of reads in called peaks (FRiP) (1)
 -) The fraction of peaks overlapping coding regions
 -) OPTIONAL: The fraction of peaks overlapping with peaks of a technical/biological replicate experiment
 
(1) Landt et al. 2012, ChIP-seq guidelines and practices of the ENCODE and modENCODE consortia. 

=head1 OPTIONS

 -c	<chip.bam>; BAM file of the ChIP experiment. Best choice is a file with duplicate reads (i.e. merged.bam); REQUIRED
 -p	<peaks.bed> BED file containing the called peaks; REQUIRED
 -o	<outfile.tsv> output file name; default: STDOUT
 -is	name of the input sample, if applicable 
 -rs	name of the technical/biological replicate sample; REQUIRED if -rp is defined
 -rp	<replicate_peaks.bed> BED file containing the called peaks of the technical/biological replicate sample; REQUIRED if -rs is defined
 -w	window specifying the maximum distance of a peak to a replicate-peak to be counted as "overlapping"; default: 0 --> only really overlapping peaks are counted
 -se	settings; REQUIRED
 -lf	<log file>; default: SCREEN
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page

=head1 AUTHOR

Thomas Wieland

