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

my $chipfile   = "";
my $inputfile  = "";
my $outdir     = "";
my $settings   = "";
my $overwrite  = 0;
my $help	   = 0;
my $man		   = 0;
my $logfile    = "pipeline.log";
my $loglevel   = "INFO";


GetOptions(
	"c=s"  => \$chipfile,
	"i=s"  => \$inputfile,
	"o=s"  => \$outdir,
	"se=s" => \$settings,
	"v"    => \$overwrite,
	"h"    => \$help,
	"lf=s" => \$logfile,
	"ll=s" => \$loglevel,
	"man"  => \$man
);

pod2usage( {-exitval => 0  ,-verbose => 1} ) if $help;
pod2usage( {-exitval => 0  ,-verbose => 2} ) if $man;
pod2usage( {-exitval => 1  ,-verbose => 1} ) if $chipfile eq "" || $settings eq "";

Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();

my $params  = Utilities::getParams();
my $bedtools = $params->{programs}->{bedtools}->{path};
my $samtools = $params->{programs}->{samtools}->{path};
my $sissrs   = $params->{programs}->{sissrs}->{path};
my $ref      = $params->{settings}->{$settings}->{reference};

$outdir = dirname($chipfile) if $outdir eq "";

$logger->info("Preparing ChIP-File...");
&prepareFile($chipfile,"$outdir/ChIP.");

my $inputoption = "";
if($inputfile ne ""){
	$logger->info("Preparing Input-File...");
	&prepareFile($inputfile,"$outdir/Input.");
	$inputoption = "-b $outdir/Input.no0.bed";
}


$logger->info("Getting genome size...");
open FAI,"awk 'BEGIN{sum=0}{sum+=\$2}END{print sum}' $ref.fai |" or exit $logger->error("Can't open $ref.fai!");
my $genomesize = <FAI>;
chomp $genomesize;
close FAI;



my $command;

#running sissrs
if(!(-e "$outdir/ChIP.sissrs.out") || $overwrite){
	$command = "perl $sissrs -i $outdir/ChIP.no0.bed -o $outdir/ChIP.sissrs.out -s $genomesize $inputoption";
	$logger->info("Runnings sissrs...");
	$logger->debug("Sissrs command: ".$command);
	system($command)
}


#converting sissrs to BED
if(!(-e "$outdir/ChIP.sissrs.bed") || $overwrite){
	$command = "awk '{if(\$1 ~ /^chr/ ){if(\$2<1){\$2=1};print \$1\"\t\"\$2\"\t\"\$3}}' $outdir/ChIP.sissrs.out > $outdir/ChIP.sissrs.bed";
	$logger->info("Converting sissrs output to BED...");
	$logger->debug("Convert: ".$command);
	system($command)
}

#converting BED to FASTA
if(!(-e "$outdir/ChIP.sissrs.fa") || $overwrite){
	$command = "$prog_path/getFastaFromBed.pl -i $outdir/ChIP.sissrs.bed -r $ref -b $chipfile > $outdir/ChIP.sissrs.fa";
	$logger->info("Converting BED to FASTA...");
	$logger->debug("Convert: ".$command);
	system($command)
}


#calculate statistics
if(!(-e "$outdir/ChIP.sissrs.stats") || $overwrite){
	
	open STATS, ">$outdir/ChIP.sissrs.stats" or exit $logger->error("Can't open $outdir/ChIP.sissrs.stats for writing!");
	print STATS "#Uniquely mapped reads\tPeaks\tNRF(>0.8)\tFRiP(>0.01)\n";
	
	$logger->info("Calculating NRF...");
	open NRF,"cut -f 1,2 $outdir/ChIP.no0.bed | uniq -c | awk '{sum+=\$1}END{print sum\"\t\"(NR/sum)}' |" or exit $logger->error("Can't open $outdir/ChIP.no0.bed!");
	my $tmp = <NRF>;
	chomp $tmp;
	close NRF;
	my ($reads,$nrf) = split("\t",$tmp);

	$logger->info("Calculating Number of Peaks...");
	open PEAKS, "wc -l $outdir/ChIP.sissrs.bed | cut -d \" \" -f 1 |" or exit $logger->error("Can't open $outdir/ChIP.sissrs.bed!");
	my $peaks = <PEAKS>;
	chomp $peaks;
	close PEAKS;
	
	$logger->info("Calculating FRiP...");
	open FRIP, "$bedtools/intersectBed -a $outdir/ChIP.no0.bed -b $outdir/ChIP.sissrs.bed -u -wa | wc -l |" or exit $logger->error("Can't open $outdir/ChIP.sissrs.bed or $outdir/ChIP.no0.bed!");
	my $frip = <FRIP>;
	chomp $frip;
	$frip = $frip/$reads;
	close FRIP;
	
	print STATS "$reads\t$peaks\t$nrf\t$frip\n";
	close STATS;
		
	
}


##############################################
sub prepareFile {
	my $infile    = shift;
	my $outprefix = shift;
	
	
	my $command;
	
	#index input file
	if(!(-e "$infile.bai") || $overwrite){
		$command = "$samtools index $infile";
	 	$logger->debug("Indexing: ".$command);
	 	system($command);
	}
	
	#create BED file with non-0 mapping quality reads
	if(!(-e $outprefix."no0.bed") || $overwrite){
		$command = "$samtools view -bh -q 1 $infile | $bedtools/bamToBed -i stdin > ".$outprefix."no0.bed";
	 	$logger->debug("Converting to BED: ".$command);
	 	system($command);
	}
	
}



=head1 NAME

runSissrs.pl

=head1 SYNOPSIS

runSissrs.pl -c chip.bam -i input.bam 
 
=head1 DESCRIPTION

Wrapper script that calls sissrs for a ChIP-Seq sample and an (optional) Input DNA sample.
Does all the necessary conversions.

=head1 OPTIONS

 -c	<chip.bam>; BAM file of the ChIP experiment. Best choice is a file with duplicate reads (i.e. merged.bam); REQUIRED
 -i	<input.bam>; BaM file of the Input DNA experiment; optional
 -o	</output/dir>; directory where the output files will be generated; default: directory where chip.bam lies
 -se	settings; REQUIRED
 -v	overwrite files that have been already generated
 -lf	<log file>; default: pipeline.log
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page

=head1 AUTHOR

Thomas Wieland

