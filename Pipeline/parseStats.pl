#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
#use DBI;
#use POSIX;
use File::Basename;
use Pod::Usage;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);


my $prog_path = dirname(abs_path($0));
require $prog_path."/Utilities.pm";
require $prog_path."/CheckGender.pm";
require $prog_path."/CollectMetrics.pm";



my $bamfile		= "";
my $elementNum  = 0;
my $outfile	    = "";
my $flagstat    = "";
my $rmdup       = "";
my $calcBaseCov = "";
my $onTarget    = "";
my $verifyBam   = "";
my $exomeDepth  = "";
my $onofftarget = "";
my $checkGender = 0;
my $picardStats = 0;
my $settings    = "default";

my $logfile     = "SCREEN";
my $loglevel    = "INFO";
my $help	    = 0;
my $man         = 0;

GetOptions(
"b=s"  => \$bamfile,
"n=s"  => \$elementNum,
"o=s"  => \$outfile, 
"f=s"  => \$flagstat,
"r=s"  => \$rmdup, 
"c=s"  => \$calcBaseCov, 
"l=s"  => \$onTarget,
"v=s"  => \$verifyBam,
"e=s"  => \$exomeDepth,
"g"    => \$checkGender,
"p"	   => \$picardStats,
"ct=s" => \$onofftarget,
"se=s" => \$settings, 
"lf=s" => \$logfile,
"ll=s" => \$loglevel, 
"man"  => \$man,
"h"    => \$help 
);


pod2usage( {-exitval => 0  ,-verbose => 1} ) if $help;
pod2usage( {-exitval => 0  ,-verbose => 2} ) if $man;
pod2usage( {-exitval => 1  ,-verbose => 1} ) if $bamfile eq "";

Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();
my $outdir = dirname(abs_path($bamfile));

my $params = Utilities::getParams();
my $samtools = $params->{programs}->{samtools}->{path};

#open output file, if required
if($outfile ne ""){
	close STDOUT;
	if($outfile =~ /\.gz$/){
		open(STDOUT,"| bgzip -c > $outfile") or die "Can't open $outfile for writing!\n";
	}else{
		open(STDOUT,">",$outfile) or die "Can't open $outfile for writing!\n";
	}
	
}




my $read1Length = 0;
my $read2Length = 0;
my $itemCount   = 0;
my $tmplen      = 0;
my $sample      = "";

my $duprate     = "NULL";
my $optduprate  = "NULL";

# TODO:
# 
# miSeq experiments have hard/trimmed reads  --> if ALIGNMENT METRIC WAS RUN then can use MEAN_READ_LENGTH
# HISEQ 4k as well
# READ > 2500 reads and use max len 

#open merge file to get read length and sample name
open(MERGE, "$samtools view -h $bamfile |");
while(<MERGE>) {
	chomp;
	if($_ =~ /^\@/){
		if($_ =~ /^\@RG/){
			my @columns = split();
			foreach my $field(@columns){
				if($field =~ /SM:(.+)/){
					$sample = $1;
				}
			}
		}
		next;
	}
	my @columns = split(' ');
	$tmplen = length $columns[9];
	
	if($columns[1] & 64) #if first in pair
	{			
		 $read1Length = ( $read1Length < $tmplen ? $tmplen : $read1Length );  
	} 
	else
	{
		 
		 $read2Length = ( $read2Length < $tmplen ? $tmplen : $read2Length );
	}
	
	$itemCount++;
	
	if($read1Length != 0 && $read2Length != 0 && $itemCount>2500) {
		last;
	}
}
close MERGE;
	



#print header
print "#sample\tduplicates(proportion)\treads\tmapped reads\tmapped(in \%)\tmapped sequence(GB)\tuncovered(in \%)\t>=1X(in \%)\t>=4X(in \%)\t>=8X(in \%)\t>=20X(in \%)\taverage cov\tstd\tmedian coverage\tm-std\tautosomal coverage\ton target(in \%)\tcontamination(in \%)\tRsd\tSRY\tBAM path\tmismatch rate\tlibcomplexity\tavgqual\tavgquallast5\tq30fraction\tonofftarget50percentcoverage\topticalduplicates\tproperlyp\n";
#print "#sample\tduplicates(proportion)\treads\tmapped reads\tmapped(in \%)\tmapped sequence(GB)\tuncovered(in \%)\t>=1X(in \%)\t>=4X(in \%)\t>=8X(in \%)\t>=20X(in \%)\taverage cov\tstd\tmedian coverage\tm-std\tautosomal coverage\ton target(in \%)\tcontamination(in \%)\tRsd\tSRY\tBAM path\tmismatch rate(in \%)\tA(in \%)\tC(in \%)\tG(in \%)\tT(in \%)\tN(in \%)\tcycles with base imbalance diff>average+/-4s\tAT dropout\tGC dropout\tlibrary complexity\tQ average\tcycles with Q<Qavg-4sigma\tQ>10(in \%)\tQ>20(in \%)\tQ>30(in \%)\n";
print "$sample\t";


my $mappedReads = 0;



#TW 14.04.2016: if CollectMetrics::getDuplicates does not return duplicates AND no samtools rmdup log file is given --> extract duplicate rate from flagstat file, if given
if($rmdup ne "" && -e $rmdup &&  CollectMetrics::getDuplicates($bamfile, $outdir, $settings) eq 'NULL' ){
	&parseRmdup($rmdup);
}
elsif ( CollectMetrics::getDuplicates($bamfile, $outdir, $settings) ne 'NULL' ) 
{
	my %dupmetric=CollectMetrics::getDuplicates($bamfile, $outdir, $settings);
		$duprate    = $dupmetric{"DUPLICATION_RATE"};
		$optduprate = $dupmetric{"OPTICAL_DUPLICATION_RATE"};
		print $duprate."\t";
}


if($flagstat ne "" && -e $flagstat){
	&parseFlagstat($flagstat);
}else {
	print "NULL\t" if $rmdup eq "" && $duprate eq "NULL";
	#print "0\t0\t0\t0\t";
	print "NULL\tNULL\tNULL\tNULL\t";
}
if($calcBaseCov ne "" && -e $calcBaseCov) {
	&parseCalcBaseCov($calcBaseCov);
} else {
	#$logger->debug("no calcBaseCov file available -> writing 0 values");
	#print "0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t";
	print "NULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\t";
}

if($onTarget ne "" && -e $onTarget) {
	&parseOnTarget($onTarget);
} else {
	#$logger->debug("no onTarget file available -> writing 0 values");
	#print "0\t";
	print "NULL\t";

}

if($verifyBam ne "" && -e $verifyBam){
	&parseVerifyBam($verifyBam);
}else{
	print "NULL\t";
}

if($exomeDepth ne "" && -e $exomeDepth){
	&parseExomeDepth($exomeDepth);
}else{
	print "NULL\t";
}

if($checkGender){
	 print CheckGender::calcCov($params->{settings}->{$settings}->{sry},$settings,$mappedReads,$bamfile,$params->{settings}->{$settings}->{sryfactor},$logfile,$loglevel)."\t";
}else{
	print "NULL\t";
}

my $bamfilePath = abs_path($bamfile);
if($elementNum > 0){
	$elementNum *= -1;
	my @columns  = split(/\//,$bamfilePath);
	$bamfilePath = join("/",@columns[$elementNum..-1]);
}

print $bamfilePath."\t";

if($picardStats){				#2016-02-18 RB ADDED PICARD STATISTICS
	
	# ALIGNMENT_METRIC
	# 	PF_MISMATCH RATE
	# 	mismatchrate;
	# LIBCOMPLEXITY (#) GET NUMBER
	# 	libcomplexity;
	# meanqualitybycycle	
	# 	avgqual;
	# 	avgquallast5;
	# QDIST	
	# 	q30fraction;
	
	
	#ALN
	my $mismatch_rate 	= CollectMetrics::getAlignmentMetric($bamfile, $outdir, $settings, "PF_MISMATCH_RATE" );
	#BaseD
	my %basedistr     	= CollectMetrics::getBaseDistribution($bamfile, $outdir, $settings);
		my $A=$basedistr{"A"}; my $ASD=$basedistr{"ASD"};
		my $C=$basedistr{"C"}; my $CSD=$basedistr{"CSD"};
		my $G=$basedistr{"G"}; my $GSD=$basedistr{"GSD"};
		my $T=$basedistr{"T"}; my $TSD=$basedistr{"TSD"};
		my $N=$basedistr{"N"}; my $NSD=$basedistr{"NSD"};
		my $badCompReadCount4s = CollectMetrics::getBaseDistributionBadCycles($bamfile, $outdir, $settings, 4);	
	#GCBIAS
	my $gcATDropout		= CollectMetrics::getGCMetric($bamfile, $outdir, $settings, "AT_DROPOUT" );
	my $gcGCDropout		= CollectMetrics::getGCMetric($bamfile, $outdir, $settings, "GC_DROPOUT" );
	#LIBCPLX
	my $libcomplexity 	= CollectMetrics::estimateLibraryComplexity($bamfile, $outdir, $settings);
	my $libcomplexityN	= CollectMetrics::estimateLibraryComplexityNormalised($bamfile, $outdir, $settings);
	#Q
	my %QavgBlock		= CollectMetrics::getMeanQuality($bamfile, $outdir, $settings);
		my $Qavg 		= $QavgBlock{"Q"};
	my $badCyclesQ4s 	= CollectMetrics::getMeanQualityBadCycles($bamfile, $outdir, $settings, 4);
	my %QavgLast5Block  = CollectMetrics::getMeanQualityLastN($bamfile, $outdir, $settings, 5);
		my $QavgLast5 	= $QavgLast5Block{"Q"}; 
	#QDIST
	my @qfract			= CollectMetrics::getQFraction($bamfile, $outdir);
		my $Q10 = $qfract[0];
		my $Q20 = $qfract[1];
		my $Q30 = $qfract[2];
	
	# All items (15)
	#print
	#	$mismatch_rate*100	.
	#	"\t" .	$A  . "\t" .	$C  . "\t" .	$G  . "\t" .	$T  . "\t" .	$N  . "\t" .	$badCompReadCount4s	.
	#	"\t" .	$gcATDropout	. "\t" .	$gcGCDropout	.
	#	"\t" .	$libcomplexity	.
	#	"\t" .	$Qavg	. "\t" .	$badCyclesQ4s	.
	#	"\t" .	$Q10*100	."\t" .	$Q20*100	. "\t" .	$Q30*100	.
	#	"\t";
	
	# Selection to the db:
	print	$mismatch_rate*100	. "\t" .
			$libcomplexityN		. "\t" .
			$Qavg				. "\t" .
			$QavgLast5			. "\t" .
			$Q30				. "\t";
	 
}else{
	print "NULL\tNULL\tNULL\tNULL\tNULL\t";
}

if ($onofftarget ne "" && -e $onofftarget."_ontarget" && -e $onofftarget."_offtarget")
{
	&diffOnOffTarget($onofftarget)
}
else
{
	print "NULL\t";
}

# Optduprate (it's null in case. Separated by its function to maintain compatibility with old versions of the file)
print $optduprate."\t";

# Print properlypaired properlyp
if($flagstat ne "" && -e $flagstat){
	&parseFlagstatPP($flagstat);
}else {
	print "NULL";
}


#############
#subroutines#
#############


############
#parseRmdup#
############
sub parseRmdup {
	my $rmdup = shift;
	
	open(RM, $rmdup) || exit $logger->error("Cannot open $rmdup");
	my $allreads = 0;
	my $rmdups   = 0;
	while(<RM>)
	{
		chomp;
		my $line = $_;
		if($line =~ /library/)
		{ 
			$line =~ /\s(\d+)\s\/\s(\d+)\s/;
			$rmdups   += $1;
			$allreads += $2;
			
		}
	}
	$rmdups /= $allreads;	#27.06.2011 patch TW: samtools rmdup writes a separate entry for each lib --> calculate mean duplicate rate of the libs
	print "$rmdups\t";
	close(RM);
}


###############
#parseFlagstat and parseFlagstatPP
###############
sub parseFlagstat {
	my $flagstat = shift;

	
	open(FL, "$flagstat") || exit $logger->error("Cannot open $flagstat");
	my $read1;
	my $read2;
	my $total      = 1;
	my $duplicates = 0;
	while(<FL>)
	{
		chomp;
		my $line = $_;
		
		#total reads
		if ($line =~ /(\d+)\s*.*\sin\stotal/)
		{
			$total = $1;		#--> don't print now, has to be after duplicates if rmdup file is not given
		}
		elsif ($line =~ /(\d+)\s*.*\sduplicates/){
			$duplicates = $1;
			
					
		}
		elsif ($line =~ /(\d+)\s*.*\smapped\s\(/) #mapped reads and %
		{
			
			if($duplicates == 0){
				$duplicates = "NULL";
			}else{
				$duplicates /= $1;
			}
			print "$duplicates\t" if $rmdup eq "" && $duprate eq "NULL";
			print "$total\t";	
			print "$1\t".sprintf("%.2f",($1/$total)*100)."\t";
			$mappedReads = $1;
			
		}
		elsif ($line =~ /(\d+)\s*.*\sread1/) # number of read 1 !!!!(mapped+unmapped)!!!!
		{
			$read1 = $1;
			
		}
		elsif ($line =~ /(\d+)\s*.*\sread2/) # number of read 2 !!!!(mapped+unmapped)!!!!
		{
			$read2 = $1;
			
		}
	}
	#my $sequenced = (($read1*$read1Length)+($read2*$read2Length)) / 1000000000;
	my $sequenced = ($mappedReads* (($read1Length+$read2Length)/2) ) / 1000000000;
	print "$sequenced\t";
	
	close(FL);
}

sub parseFlagstatPP {
	my $flagstat = shift;

	
	open(FL, "$flagstat") || exit $logger->error("Cannot open $flagstat");
	my $read1;
	my $read2;
	my $total      = 1;
	my $duplicates = 0;
	while(<FL>)
	{
		chomp;
		my $line = $_;
		
		if ($line =~ /(\d+)\s*.*\sin\stotal/)
		{
			$total = $1;		#--> don't print now, has to be after duplicates if rmdup file is not given
		}
		elsif ($line =~ /(\d+)\s*.*\sproperly/) #properly paired
		{
			print sprintf("%.2f",($1/$total)*100)."";			
		}
	}
	close(FL);
}



##################
#parseCalcBaseCov#
##################
sub parseCalcBaseCov {
	my $calcBaseCov = shift;
	open(CB, $calcBaseCov) || exit $logger->error("Cannot open $calcBaseCov");
	while(<CB>)
	{
		chomp;
		my $line = $_;
	
		if ($line =~ /\d+\s\((.+?)\)/) #covered stats
		{  
			my $tmp = $1;
			$tmp =~ s/\%$//;
			print "$tmp\t";
		}
		elsif ($line =~ /(.+?)\saverage\scoverage\s\((.*?)\sstandard\sdeviation\)/)	#average coverage
		{
			my $tmp = $1;
			$tmp =~ s/\%$//;
			print "$tmp\t";
			$tmp = $2;
			$tmp =~ s/\%$//;
			print "$tmp\t";
		}
		elsif ($line =~ /(.+?)\smedian\scoverage\s\((.*?)\)/) 	#median coverage
		{
			my $tmp = $1;
			$tmp =~ s/\%$//;
			print "$tmp\t";
			$tmp = $2;
			$tmp =~ s/\%$//;
			print "$tmp\t";
		}
		elsif ($line =~ /(.+?)\saverage\scoverage\sof\sautosomes/) 	#autosomal coverage
		{
			my $tmp = $1;
			$tmp =~ s/\%$//;
			print "$tmp\t";
		}
	}
	close(CB);
}

#########################
#parseOnTarget          #
#########################
sub parseOnTarget {
	my $onTarget = shift;
	
	open(LO, $onTarget) || exit $logger->error("Cannot open $onTarget");
	while(<LO>)
	{
		chomp;
		my $line = $_;
		if($line =~ /search/)
		{ 
			next;
		}
		else
		{
			if ($line =~ /\d+\s\((.+?)\)\son\starget/) #on target
			{
				my $tmp = $1;
				$tmp =~ s/\%$//;
				print "$tmp\t";
			}
		}
	}
	close(LO);
}


#########################
#parseVerifyBam         #
#########################
sub parseVerifyBam {
	my $verifyBam = shift;
	
	open(IN,$verifyBam) or exit $logger->error("Can't open $verifyBam!");
	my $line = <IN>;
	$line = <IN>;
	my @columns = split("\t",$line);
	#print $columns[-4]."\t";			#old verifyBamID version
	print $columns[6]."\t";			#new verifyBamID version
}


#########################
#parseExomeDepth        #
#########################
sub parseExomeDepth {
	my $exomeDepth = shift;
	open EDSTATS, $exomeDepth or exit $logger->error("Cannot open $exomeDepth!");
	my $line = <EDSTATS>;
	$line    = <EDSTATS>;
	my @columns = split(' ',$line);
	print $columns[3]."\t";	
}

##################
#parseOnOffTarget#
##################
#sub parseOnOffTarget {
#	
#	# Read on and off target coverage histogram and get the top coverage for the 50% of the coding / noncoding regions
#	# 0x 100%, 1X 99.9%, ...etc
#	
#	my $basename = shift;
#		my $ontarget  = $basename."_ontarget";
#		my $offtarget = $basename."_offtarget";
#	my $threshold = 0.5;
#		
#	open(IN_ONTARGET,  $ontarget  ) or exit $logger->error("Can't open $ontarget ! ");
#	open(IN_OFFTARGET, $offtarget ) or exit $logger->error("Can't open $offtarget !");
#	
#	my $ontarget_coverage=0;
#	my $offtarget_coverage=0;
#	
#	while (<IN_ONTARGET>)
#	{
#		my @vals=split(" ", $_);
#		my $frac=$vals[1];
#		
#		if ( $frac < $threshold )
#		{
#			last;
#		}		
#						
#		$ontarget_coverage = $vals[0];		
#	}
#
#	while (<IN_OFFTARGET>)
#	{
#		my @vals=split(" ", $_);
#		my $frac=$vals[1];
#		
#		if ( $frac < $threshold )
#		{
#			last;
#		}		
#						
#		$offtarget_coverage = $vals[0];		
#	}
#	
#	print $ontarget_coverage."\t".$offtarget_coverage."\t";
#}

##################
#diffOnOffTarget#
##################
sub diffOnOffTarget {
	
	# Read on and off target coverage histogram and get the top coverage for the 50% of the coding / noncoding regions
	# 0x 100%, 1X 99.9%, ...etc
	
	my $basename = shift;
		my $ontarget  = $basename."_ontarget";
		my $offtarget = $basename."_offtarget";
	my $threshold = 0.5;
		
	open(IN_ONTARGET,  $ontarget  ) or exit $logger->error("Can't open $ontarget ! ");
	open(IN_OFFTARGET, $offtarget ) or exit $logger->error("Can't open $offtarget !");
	
	my $ontarget_coverage=0;
	my $offtarget_coverage=0;
	
	while (<IN_ONTARGET>)
	{
		my @vals=split(" ", $_);
		my $frac=$vals[1];
		
		if ( $frac < $threshold )
		{
			last;
		}		
						
		$ontarget_coverage = $vals[0];		
	}

	while (<IN_OFFTARGET>)
	{
		my @vals=split(" ", $_);
		my $frac=$vals[1];
		
		if ( $frac < $threshold )
		{
			last;
		}		
						
		$offtarget_coverage = $vals[0];		
	}
	
	print $ontarget_coverage-$offtarget_coverage."\t";
}

=head1 NAME

parseStats.pl

=head1 SYNOPSIS

parseStats.pl -b merged.bam -o summary.stats.tsv -r merged.rmdup.out -f merged.flagstat.out -c calcBaseCov.out -l calcOnTarget.out -g

=head1 DESCRIPTION

This script parses the output files of several scripts:
 -) calcBaseCov_api.pl: custom script to calculate coverage of target regions
 -) calcOnTarget.pl: custom script to calculate number and proportion of reads on target
 -) samtools flagstat: samtools script that parses BAM file and calculates summary metrics
 -) samtools rmdup: the log file of samtools rmdup contains the percentage of duplicate reads that is extracted here
 -) verifyBamID: a script by the University of Michigan to detect putative contamination of a sample
    see: http://genome.sph.umich.edu/wiki/VerifyBamID
 -) ExomeDepth.stats: the ExomeDepth R package detects CNVs from exome data. It also outputs a value "Rsd"
    which is a measure of how similar a given sample is to a set of reference samples.
    see: https://cran.r-project.org/web/packages/ExomeDepth/index.html

The script also opens a given BAM file to get the sample name and the read length in order to calculate the amount of
sequence. It calculates the length for read1 and read2 separately, but only once, so reads with different lengths
(e.g. clipped reads) are not supported. It also checks the coverage of the SRY gene as a measure for gender, if required.

All the parsed information is printed to a single TSV output file which can be imported into the database using
the script statsdb.pl. Any information that is not given is set to NULL or empty string.


=head1 OPTIONS

 -b	<reads.bam> BAM file to extract read length and sample name; REQUIRED
 -n	[number of path components] number of path components of the bam file that should be written to the file.
 	The reads from a BAM file can be accessed using the web interface. The interface supports to show only BAM files relative
 	to a given root directory. For this reason only a relative path of a BAM file can be stored in the database.
 	Which part should be stored can be controlled using this parameter. 
 	Example: -n 4; -b /long/full/qualified/path/to/bam/file.bam 
 	--> only the last 4 elements should be stored: path/to/bam/file.bam
 	Default: take full path
 -o	<output.tsv> default: print to stdout
 -f	<flagstat.out> output file of samtools flagstat to get number of reads 
 -r	<rmdup.out> log file of samtools rmdup. If this file is specified, the percentage of duplicated reads is taken
 	from this file instead of the flagstat file
 -c	<calcBaseCov.out> output file of the calcBaseCov_api.pl script to get the coverage values
 -l	<calcOnTarget.out> output file of the calcOnTarget.pl script to get the percentage of on-target reads
 -v	<verifyBamID.out> output file of the verifyBamID script (see http://genome.sph.umich.edu/wiki/VerifyBamID)
 -e	<exomeDepth.stats> statistics from ExomeDepth (see https://cran.r-project.org/web/packages/ExomeDepth/index.html)
 -g	calculate the normalized coverage of the SRY gene to check the gender of the given sample.
 -p add picard statistics
 -ct <coverageprofile> basename excluded [_on/offtarget] for coverage profile statistics
 -se	name of the settings in the current.config.xml file that holds path to reference genome, 
 	to the annotation file and to possible additional annotation files; use default settings if nothing is given
 -lf	log file; default: print to screen
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page

=head1 AUTHOR

Thomas Wieland

=cut
