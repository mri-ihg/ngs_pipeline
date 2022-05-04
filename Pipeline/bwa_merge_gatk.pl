#!/usr/bin/perl

use strict;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);


my $prog_path = dirname(abs_path($0));
require $prog_path."/Utilities.pm";

my $outdir        = "";
my $infiles       = "";
my @infiles       = ();
my $aligner       = "bwa";
my $help          = 0;
my $settings      = "";
my $removeInfiles = 0;
my $removeSai     = 0;
my $dontRmdp      = 0;
my $rmdupWithPicard = 0;
my $markdupWithPicard = 0;
my $clipMIPPrimer = 0;
my $recal = 0;
my $sampleName = "";

my $logfile  = "SCREEN";
my $loglevel = "INFO";

my $run = 1;

my $helptext      = 
"-i	<bam-files to be merged>; if empty: use all .sort.bam (or .gem.sort.bam) files in the output directory
-o	<output directory>
-r	remove all infiles after merging
-rs	remove all .sai files (BWA index files)
-d	don't remove duplicates, just create link
-picard use picard to remove duplicates   (alternative to -picardmark)
-picardmark use picard to mark duplicates (alternative to -picard)
-c	clip MIP primers; clip the primers from a MIP experiment using the clipMIPPrimer.pl script; only works if -d is set.
-se	settings
-lf	log file; default: pipeline.log
-ll log level: ERROR,INFO,DEBUG; default: INFO
-h this help\n";




GetOptions(
"o=s"  => \$outdir, 
"i=s"  => \$infiles, 
"r"    => \$removeInfiles,
"rs"   => \$removeSai,
"d"    => \$dontRmdp,
"picard"=> \$rmdupWithPicard,
"picardmark" => \$markdupWithPicard,
"c"    => \$clipMIPPrimer,
"recal"=> \$recal,
"lf=s" => \$logfile,
"ll=s" => \$loglevel,
"se=s" => \$settings,
"a=s"  => \$aligner,
"s=s"  => \$sampleName,
"h" => \$help);

Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();

if ($help == 1 || $outdir eq "" || ( $markdupWithPicard && $rmdupWithPicard )) {
	print $helptext;exit(1);
}

my $n = 0;

if ($infiles ne "") { #infiles on command line
	@infiles = split(/\s+/,$infiles);
}
else { # read all sort.bam files from directory
	if ($aligner eq "tophat") {
		@infiles = glob("$outdir/TophatAlignment/accepted_hits.sort.bam");
	}elsif ($aligner eq "gem") {
 		@infiles = glob("$outdir/*.gem.sort.bam");
	} elsif ($aligner eq "fastq2bam") {
		@infiles = glob("$outdir/$sampleName*_L*_001.bam");
	} else {
		@infiles = glob("$outdir/*.sort.bam");
	}
	$n = @infiles;
	foreach (@infiles) {
		$infiles.="$_ "; 
	}
}
$logger->debug("Infiles for merging: $infiles");

if ((scalar (@infiles)) < 1) {
	$logger->error("no files found. exiting ...");
	exit(11);
}

my $params = Utilities::getParams();

my $bedtools = $params->{programs}->{bedtools}->{path};
my $samtools = $params->{programs}->{samtools}->{path};
my $sammem   = $params->{programs}->{samtools}->{maxmem};		#get params
my $java     = $params->{programs}->{java}->{path};
my $gatk     = $params->{programs}->{gatk}->{path};
my $gatktmp  = $params->{programs}->{gatk}->{tmpdir};
my $picard   = $params->{programs}->{picard}->{path};

my $ref      = $params->{settings}->{$settings}->{reference};


my $header = "";
my @headerContent = ();

my $filecount = 1;
foreach my $file (@infiles)
{
	system("$samtools view -H $file > $file.header");
	open(HEAD,"$file.header") || exit $logger->error("Cannot open $file.header");

	if($filecount == 1)
	{
		while(<HEAD>)
		{
			$header .= $_;
		}
	}	
	else
	{
		while(<HEAD>){
			if($_ =~ /\@RG/){
				$header .= $_;
			}
		}
		
	}
	$filecount++;
	close(HEAD);
}

open(OUT, ">$outdir/header.sam")||exit $logger->error("Cannot open $outdir/header.sam");
print OUT "$header";
close OUT;


if($run)
{

	#if(!(-e "$outdir/merged.bam")){			# !!!!!!!!!!! merged.bam only gets created when it doesn't exist already
		system("rm $outdir/merged.bam");
		if ($n > 1) {
			$logger->info("Running samtools merge");
			$logger->debug("CMD: $samtools merge -h $outdir/header.sam $outdir/merged.bam $infiles");
			
			if(system("$samtools merge -h $outdir/header.sam $outdir/merged.bam $infiles") == 0 && $removeInfiles){		#remove infiles if chosen and merge was successful
			#if($removeInfiles){
				$logger->info("Removing input files...");
				foreach(@infiles){
					unlink($_);
				}
			}
		}
		else {
			my $infile = basename($infiles);
			$logger->info("infiles: $infiles : $infile");
			if ($aligner ne "tophat") {
				system("ln -s $infile $outdir/merged.bam");
			} else {
				system("ln -s $infiles $outdir/merged.bam");
			}
		}
	#}
	
	

	
	
	
	
	if($dontRmdp){
		
		if($clipMIPPrimer){
			unlink("$outdir/merged.rmdup.bam");
			my $command = "perl $prog_path/clipMIPPrimer.pl -i $outdir/merged.bam -o $outdir/merged.rmdup.bam -c $outdir/mips.count -se $settings -lf $logfile -ll $loglevel";
			$logger->info("Clipping MIP primers...");
			$logger->debug($command);
			system($command);
		}else{
			$logger->info("Set link for rmdup file");
			system("ln -s merged.bam $outdir/merged.rmdup.bam");
		}
	}else{
		
		# Picard or Samtools
		if($rmdupWithPicard || $markdupWithPicard) {
			# Get option to mark or remove dups
			my $rmdupopt=( $rmdupWithPicard ? "true" : "false" );
			$logger->info("Running Picard MarkDuplicates: $java -Xmx6g -XX:ParallelGCThreads=1 -jar $picard MarkDuplicates INPUT=$outdir/merged.bam OUTPUT=$outdir/merged.rmdup.bam M=$outdir/rmdup.metric REMOVE_DUPLICATES=$rmdupopt ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT");
			system("$java -Xmx6g -XX:ParallelGCThreads=1 -jar $picard MarkDuplicates INPUT=$outdir/merged.bam OUTPUT=$outdir/merged.rmdup.bam M=$outdir/rmdup.metric REMOVE_DUPLICATES=$rmdupopt ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT");  
		} else {
			$logger->info("Running samtools rmdup: $samtools rmdup $outdir/merged.bam $outdir/merged.rmdup.bam 2> $outdir/merged.rmdup.log");
			system("$samtools rmdup $outdir/merged.bam $outdir/merged.rmdup.bam 2> $outdir/merged.rmdup.log");
		}
	}
	
	$logger->info("Running samtools index $outdir/merged.bam");
	system("rm $outdir/merged.bai -rf");
	system("rm $outdir/merged.bam.bai -rf");
	system("$samtools index $outdir/merged.bam");

	$logger->info("Running samtools index $outdir/merged.rmdup.bam");
	system("rm $outdir/merged.rmdup.bai -rf");
	system("rm $outdir/merged.rmdup.bam.bai -rf");
	system("$samtools index $outdir/merged.rmdup.bam");
	
	$logger->info("Running samtools flagstat on merged.rmdup.bam");
	system("$samtools flagstat $outdir/merged.rmdup.bam > $outdir/merged.rmdup.flagstat.out");
	
	$logger->info("Running samtools flagstat on merged.bam");
	system("$samtools flagstat $outdir/merged.bam > $outdir/merged.flagstat.out");


	$logger->info("Remvoing .sai files...") if $removeSai;
	system("rm $outdir/*.sai 2> /dev/null") if $removeSai;


	# Recalibration ( if needed )
	if ( $recal )
	{
		# Recal
		my $recalcommand = "perl $prog_path/recalBam.pl -i $outdir/merged.rmdup.bam -o $outdir/merged.rmdup.bam -m 24g -se $settings -gatk4 -indels_context_size 8 -mismatches_context_size 4 -lf $logfile -ll $loglevel";
		#system( $recalcommand );
		
		if (&Utilities::executeCommand($recalcommand, "Launching recalibration", $logger)) {
			$logger->error("Error executing recalibration");
			exit(100);
		}
		
		
		
		# Index recalibrated file
		system ( "$samtools index $outdir/merged.rmdup.bam" ); # if ( ! -e "$outdir/merged.rmdup.bam.bai" ); 

	}
}
