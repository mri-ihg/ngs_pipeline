#!/usr/bin/perl

use strict;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);
use Pod::Usage;

#include Utilities.pm
my $prog_path = dirname(abs_path($0));
require $prog_path."/Utilities.pm";


my $outprefix     = "";
my $help          = 0;
my $settings      = "";
my $createIndex	  = 0;

my $infile1       = "";
my $infile2       = "";
my $threads       = -1;
my $bam			  = 0;
my $illumina	  = 0;
my $logfile  	  = "pipeline.log";
my $loglevel 	  = "INFO";

my $readGroupId     = "";
my $readGroupSample = "";
my $readGroupLib    = "LIB";
my $platform        = "Illumina";
my $maxInssize		= 500;
my $man				= 0;
my $mem				= 0;
my $dryrun			= 0;

my $isExternalBam	= 0;

my $removePhiX		= 0;
my $removeAdapters  = 0;

my $helptext      = 
"\n";


GetOptions(
	"o=s"  => \$outprefix, 
	"a=s"  => \$maxInssize,
	"f=s"  => \$infile1, 
	"F=s"  => \$infile2,
	"n"    => \$createIndex,
	"m"    => \$mem,
	"b"    => \$bam, 
	"I"    => \$illumina,
	"i=s"  => \$readGroupId, 
	"s=s"  => \$readGroupSample, 
	"l=s"  => \$readGroupLib, 
	"t=s"  => \$threads,
	"d"    => \$dryrun,
	"man"  => \$man, 
	"h"    => \$help, 
	"lf=s" => \$logfile,
	"ll=s" => \$loglevel, 
	"se=s" => \$settings,
	"externalBAM"    => \$isExternalBam,
	"removeAdapters" => \$removeAdapters,
	"removePhix"     => \$removePhiX );
	
pod2usage( {-exitval => 0  ,-verbose => 1} ) if $help;
pod2usage( {-exitval => 0  ,-verbose => 2} ) if $man;
pod2usage( {-exitval => 1  ,-verbose => 1} ) if $settings eq "" || $infile1 eq "" || $outprefix eq "";


Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();
my $params = Utilities::getParams();

my $inputinfiles  = "inputinfiles";
my $bwa           = $params->{programs}->{bwa}->{path};		
my $bedtools      = $params->{programs}->{bedtools}->{path};			
my $bwareadtrim   = $params->{programs}->{bwa}->{readtrim};		
my $samtools      = $params->{programs}->{samtools}->{path};		
my $sammaxmem     = $params->{programs}->{samtools}->{maxmem};
my $cutadapt	  = $params->{programs}->{cutadapt}->{path};

my $ref           = $params->{settings}->{$settings}->{reference};
my $phiX          = $params->{settings}->{$settings}->{phiX} if defined $params->{settings}->{$settings}->{phiX}; 
my $cutadapt_adapter = $params->{programs}->{cutadapt}->{adaptersequence};

$phiX=$ref if $phiX eq "";

my $bwaref = $ref;
	$bwaref = $phiX if $removePhiX;
	
$bwaref =~ s/.fasta$//;
$bwaref =~ s/.fa$//;

my $removePhixOption = ( $removePhiX ? " grep -v phiX | " : "" );

if($threads == -1){
	$threads = 1;
	if($ENV{NSLOTS}){		#get number of slots given by SGE
		$threads = $ENV{NSLOTS};
		
	}
}

my @tmp         = ();
my @bwafiles    = ();


my $tmp_infile1="";
my $tmp_infile2="";

my $outfile1    = $outprefix."_1";
my $outfile2    = $outprefix."_2";

my $rgId = basename($outprefix);


my $readGroupString = "\'\@RG\tID:$rgId\tSM:$readGroupSample\tLB:$readGroupLib\tPL:$platform\'";

my $command;

# Patch RB 20170606
# To correct convert bam to sam: WE DO NOT KNOW IF WE HAVE AN ALIGNED BAM, IT MUST BE SORTED BY READ NAME, SECONDARY ALIGNMENT MUST BE FILTERED OUT
# samtools view -b -h -F 0x900 test.bam  | samtools sort -n -o - testTemp  | /usr/local/packages/seq/BEDTools/bedtools bamtofastq -i - -fq test1.fq -fq2 test2.fq
# previous version comment: bwa mem can't handle bam files as an input --> convert them to fastq before alignment -> all version of bwa now get converted fastq

# PREPROCESS: bam2fastq
if($bam){			
	
	my $bedCommand = "";
	
	#The filter + sort must be done for external imported files 	
	if ( $isExternalBam )
	{
		my $tmporiginalbamsort=$outprefix."_tmporiginalbamsort";
		# Aligned bam must be filtered for secondary alignment and sorted by read name 
		$bedCommand = "$samtools view -b -h -F 0x900 $infile1 | $samtools sort -n -o - $tmporiginalbamsort | $bedtools/bedtools bamtofastq -i - "; # -fq r1.fastq -fq2 r2.fastq "
	}
	else
	{
		# Our bams are already ok, just convert them to fastq
		$bedCommand = "$bedtools/bedtools bamtofastq -i $infile1 "; # -fq r1.fastq -fq2 r2.fastq "
	}
	
	$infile1 = $outprefix."_R1.fastq";
	$bedCommand .= " -fq $infile1";
	if($infile2 ne ""){
		$infile2 = $outprefix."_R2.fastq";
		$bedCommand .= " -fq2 $infile2";
	}
	
	# Over abundant error messages of no use
	$bedCommand .= " 2> /dev/null ";
		
	$logger->info("Converting BAM file into FASTQ files for BWA");
	$logger->debug($bedCommand);
	system($bedCommand) if !$dryrun;
}
	
# PREPROCESS: Remove adapter
if ($removeAdapters)
{
	my $cutadaptCommand1;
	my $cutadaptCommand2;
		
	#R1
	$tmp_infile1=$infile1;
	$infile1=$outprefix."_CA_R1.fastq";
	$cutadaptCommand1 = "$cutadapt -a $cutadapt_adapter -o $infile1 $tmp_infile1";
			
	#R2 if exists
	if($infile2 ne "")
	{
		$tmp_infile2=$infile2;
		$infile2=$outprefix."_CA_R2.fastq";
		$cutadaptCommand2 = "$cutadapt -a $cutadapt_adapter -o $infile2 $tmp_infile2";
	}
		
	# Build final command:
	# the syntax allows parallel execution of the launched command 
	# system function is blocked by "wait" until both tasks are completed
	my $cutadaptCommand=$cutadaptCommand1;
	$cutadaptCommand.=" & $cutadaptCommand2 & wait;" if $infile2 ne "";
	
	$logger->debug($cutadaptCommand);
	# Launch
	system($cutadaptCommand);
		
}
	

# ALIGNMENT	
if($mem){	
	#bwa mem
	
	$command  = "$bwa mem -R $readGroupString -M -t $threads $bwaref $infile1 ";
	$command .= "$infile2 " if $infile2 ne "";
	$command .= "2> $outprefix.log | ";
}else{		
	#bwa aln
	
	#aligning read 1
	$logger->info("Starting alignment of  $infile1");
	if($illumina){
		$command = "$bwa aln -q $bwareadtrim -t $threads -I $bwaref $infile1 > $outfile1.sai 2> $outfile1.log";
	}else{
		$command = "$bwa aln -q $bwareadtrim -t $threads  $bwaref $infile1 > $outfile1.sai 2> $outfile1.log";
	}
	
	$logger->debug($command);
	system($command) if !$dryrun;
	
	if($infile2 ne ""){
		$logger->info("Starting alignment of  $infile2");
		if($illumina){
			$command = "$bwa aln -q $bwareadtrim -t $threads -I $bwaref $infile2 > $outfile2.sai 2> $outfile2.log";
		}else{
			$command = "$bwa aln -q $bwareadtrim -t $threads  $bwaref $infile2 > $outfile2.sai 2> $outfile2.log";
		}
		
		$logger->debug($command);
		system($command) if !$dryrun;
	
		#running bwa sampe (paired end)
		$command = "$bwa sampe -a $maxInssize $bwaref $outfile1.sai $outfile2.sai $infile1  $infile2 -r $readGroupString 2> $outprefix.log | ";	
		
	}else{
		#running bwa samse (single end)
		$command = "$bwa samse $bwaref $outfile1.sai $infile1 -r $readGroupString 2> $outprefix.log | ";
		
	}

}

$command .= $removePhixOption;
$command .= "$samtools view -Sbh - | ";
$command .= "$samtools sort -\@ $threads -m $sammaxmem  - $outprefix.sort";

$logger->info("Running bwa mem or bwa sams/pe | samtools import | samtools sort: " . $command);
$logger->debug($command);

system($command) if !$dryrun;

if($createIndex){
	$command = "$samtools index $outprefix.sort.bam";
	$logger->info("Creating index for BAM file");
	$logger->debug($command);
	system($command) if !$dryrun;
}

#remove temp fastq files, if needed
#if( ($bam && $mem) || ($removeAdapters && $mem) ){
if( ($bam ) || ($removeAdapters ) ){
	unlink($infile1);
	unlink($infile2) if $infile2 ne "";
}

#if( ($bam && $mem) && ($removeAdapters && $mem) ){
if( ($bam ) && ($removeAdapters ) ){
	unlink($tmp_infile1);
	unlink($tmp_infile2) if $tmp_infile2 ne "";
}


=head1 NAME

 bwa_aln_all_Formats.pl

=head1 SYNOPSIS

 bwa_aln_all_Formats.pl -f read1.fq.gz -F read2.fq.gz -o ./outprefix -se hg19

=head1 DESCRIPTION

This script is a wrapper script for the aligner BWA.

=head1 OPTIONS

 -f	<infile read 1>; REQUIRED
 -F	<infile read 2>; empty if single reads
 -o	<output prefix>; REQUIRED
 -se	settings; REQUIRED
 -n	create index of the BAM file after alignment
 -m	use the "bwa mem" algorithm instead of the old one
 -t	<threads>, default: take number of slots from SGE environment variable NSLOTS or 1 if no SGE job
 -i	<readgroupId>, (should contain flowcell, lane and, if barcoded, bin)
 -s	<sample name>
 -l	<library name>, default LIB
 -a	maximum insert size (for paired-end libraries only); default: 500
 -b	infiles are in .bam format
 -I	infiles are in Illumina sequence.txt format
 -d	dry run; don't run commands just output the commands to logfile
 -lf	log file; default: pipeline.log
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -removePhix		remove PhiX spikes
 -removeAdapters	rempove Illumina adapters
 -man	show man page
 -h	this help

=head1 AUTHOR

Thomas Wieland, Riccardo Berutti

=cut



