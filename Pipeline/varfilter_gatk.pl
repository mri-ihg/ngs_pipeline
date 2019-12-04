#!/usr/bin/perl

use strict;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);
use Pod::Usage;
umask(002);


my $prog_path = dirname(abs_path($0));
require $prog_path."/Utilities.pm";

my $outfile    = "";
my $help       = 0;
my $params     = Utilities::getParams();
my $logfile    = "pipeline.log";
my $loglevel   = "INFO";
my $settings   = "";
my $infiles    = "";
my $region     = "";
my $outputAll  = 0;
my $regionFile = "";
my $useUG   	= 0;
my $nt			= 1;
my $nct			= 1;
my $maxRam		= "4g";
my $man			= 0;
my $isArrayJob	= 0;
my $isBEDArrayJob= 0;
my $noVCF		= 0;

# Special Filtering options:
my $dropSoftClipped = 0;
my $isPCRFree 	= 0;
my $hostilePCRCorrection = 0;
my $minBaseQualityScore = 0;
my $minPruning 	= -1;
my $contamination_fraction = 0;

my $Q30 = 0; # Use Q30 bases only!

# GATK4 is from 01.10.2019 default caller
my $usegatk4 = 1;
my $usegatk3 = 0;
my $forceActive = 0;
my $emitRealBam = 0;

# Use GVCF Block compression instead of BP_RESOLUTION
my $useGVCFResolution = 0;

GetOptions(
"i=s"  => \$infiles,
"r=s"  => \$region,
"l=s"  => \$regionFile,
"a"	   => \$outputAll,
"aj"   => \$isArrayJob,
"ajb"  => \$isBEDArrayJob,
"v"    => \$noVCF,
"noVCF"=> \$noVCF,
"fa"   => \$forceActive,
"o=s"  => \$outfile, 
"se=s" => \$settings,
"dropSoftClipped" 	=> \$dropSoftClipped,
"pcrfree" 			=> \$isPCRFree,
"pcrindelcorrection"=> \$hostilePCRCorrection,
"minbq=s"			=> \$minBaseQualityScore,
"contamination_fraction=s" => \$contamination_fraction,
"q30"  => \$Q30,
"ug"   => \$useUG,
"gatk3"=> \$usegatk3,
"gatk4"=> \$usegatk4,
"gvcfblocks" => \$useGVCFResolution,
"emitrealbam" => \$emitRealBam,
"m=s"  => \$maxRam,
"nt=s" => \$nt,
"nct=s"=> \$nct,
"lf=s" => \$logfile,
"ll=s" => \$loglevel,
"man"  => \$man,
"hcminpruning" => \$minPruning,
"h"    => \$help);

# -gatk3 option disables DEFAULT usegatk4
$usegatk4 = 0 if ( $usegatk3 );

pod2usage( {-exitval => 0  ,-verbose => 1} ) if $help;
pod2usage( {-exitval => 0  ,-verbose => 2} ) if $man;
pod2usage( {-exitval => 1  ,-verbose => 1} ) if $outfile eq "" || $settings eq "" || ( $useUG && $usegatk4 );

Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();


# Externally supplied region file
if($region ne "" && $regionFile eq ""){
	$regionFile = $region;
}

my $ref  = $params->{settings}->{$settings}->{reference};
my $gatk = $params->{programs}->{gatk}->{path};
my $java = $params->{programs}->{java}->{path};
my $tmp  = $params->{programs}->{gatk}->{tmpdir};
my $sam  = $params->{programs}->{samtools}->{path};

my $gatk4 = $params->{programs}->{gatk4}->{path};

# Set our version of Java in Path with priority so to overcome any current installation (required by GATK4)
$ENV{'PATH'}=dirname(abs_path($java)).":".$ENV{'PATH'};


my $listPrefix = "";
if($isArrayJob &&  $ENV{SGE_TASK_ID}){
	
	open BED, $params->{settings}->{$settings}->{normalchromosomes} || exit $logger->error("Can't open $params->{settings}->{$settings}->{normalchromosomes}!");
	$. = 0;
	my $line;
	do { $line = <BED> } until $. == $ENV{SGE_TASK_ID} || eof;				#read the line corresponding to this job of the array from the sequence dictionary to extract the chromosome to call
	close BED;
	
	chomp $line;
	my ($chr,$start,$end) = split("\t",$line);
	$regionFile = $chr;
	$listPrefix = $chr.".";
	$outfile = dirname($outfile)."/$chr.".basename($outfile);
	
	
}elsif($isBEDArrayJob &&  $ENV{SGE_TASK_ID} && $params->{settings}->{$settings}->{genome_splits} ){
	$regionFile = $params->{settings}->{$settings}->{genome_splits}."/".$ENV{SGE_TASK_ID}.".bed";
	$listPrefix = $ENV{SGE_TASK_ID}.".";
	$outfile = dirname($outfile)."/".$ENV{SGE_TASK_ID}.".".basename($outfile);
}


if($infiles eq ""){
	$infiles = dirname($outfile)."/merged.rmdup.bam";
}elsif(-d $infiles){
	system("ls $infiles/*.sort.bam > ". dirname($outfile)."/".$listPrefix."sort.bam.list"); 	#create list with all sort.bam files to call variants from
	$infiles = dirname($outfile)."/".$listPrefix."sort.bam.list";
}


# GATK logfile - too much verbose info to save into our log
my $javalog = $outfile;
$javalog    =~ s/vcf$/log/;

my $command="";

my $gvcf =  $outfile;
$gvcf    =~ s/vcf$/gvcf/;
$gvcf   .= ".gz";  # GATK Now allows to produce gzipped gvcf directly

my $isMultisample = 0;


# Variant calling
unless($useUG){

	# GATK HaplotypeCaller
	
	#check if multisample calling should be performed (gVCF mode doesn't support that yet )
	if($infiles =~ /\.list$/){
		open LIST, "$infiles" or exit $logger->error("Can't open $infiles!");
		my %samples;
		while(<LIST>){
			chomp;
			my $line = `$sam view -H $_ | grep \"\@RG\"`;
			chomp $line;
			my @columns = split(" ",$line);
			foreach my $currField(@columns){
				if($currField =~ /^SM:/){
					$samples{$currField} =1;
					last;
				}
			}
		}
		close LIST;
		$isMultisample = 1 if keys %samples > 1;
	}
	
	
	# GATK3 and GATK4 support
	$command = 
		$usegatk4 	?
		"
			$gatk4 --java-options  \"-Xmx$maxRam\" HaplotypeCaller \\
		"
		:
		"
			$java -XX:ParallelGCThreads=1 -Xmx$maxRam -jar $gatk -T HaplotypeCaller   \\
		"; 
		
		
	$command .= "-R $ref -I $infiles \\\n";
	$command .= "-nct $nct \\\n" if (!$usegatk4);
	
	
	$command .= "-A DepthPerAlleleBySample 	\\\n";
	$command .= "-A InbreedingCoeff 		\\\n";
	$command .= "-A HaplotypeScore			\\\n" 	if ( ! $usegatk4 );
	
	$command .= " --logging_level $loglevel -U ALLOW_N_CIGAR_READS \\\n"	if ( ! $usegatk4 );
	$command .= " --verbosity $loglevel \\\n"		if ( $usegatk4 );

	$command .= "--forceActive \\\n"  if ( $forceActive && ! $usegatk4 );
		$command .= "--force-active \\\n"  if ( $forceActive && $usegatk4 );
		
	$command .= "-L $regionFile \\\n" if $regionFile ne "";
		
	
	$command .= " --minPruning $minPruning " if ( $minPruning >0 && !$usegatk4 );
		$command .= " --min-pruning $minPruning " if ( $minPruning >0 && $usegatk4 );
	
	$command .= " --dontUseSoftClippedBases \\\n"	if $dropSoftClipped;
	$command .= " --pcr_indel_model NONE \\\n" 		if ( $isPCRFree && !( $hostilePCRCorrection )) && (! $usegatk4);
	$command .= " --pcr_indel_model HOSTILE \\\n" 	if ( $hostilePCRCorrection )&& (! $usegatk4);
		$command .= " --pcr-indel-model NONE \\\n" 		if ( $isPCRFree && !( $hostilePCRCorrection )) && ($usegatk4);
		$command .= " --pcr-indel-model HOSTILE \\\n" 	if ( $hostilePCRCorrection )&& ($usegatk4);
	
	
	$command .= " --min_base_quality_score $minBaseQualityScore \\\n"	if ( $minBaseQualityScore > 0 && !$usegatk4);
		$command .= " --min-base-quality-score $minBaseQualityScore \\\n"	if ( $minBaseQualityScore > 0  && $usegatk4);
	
	$command .= "--contamination_fraction_to_filter $contamination_fraction " if ( $contamination_fraction > 0 && !$usegatk4 );
		$command .= "--contamination-fraction-to-filter $contamination_fraction " if ( $contamination_fraction > 0 && $usegatk4 );
	
	# To be removed
	#	$command .= "--annotation StrandOddsRatio ";
	#	$command .= "--annotation StrandAlleleCountsBySample ";
	#	$command .= "--annotation StrandBiasBySample ";
		
	$command .= "--annotation OxoGReadCounts " if ( $usegatk4 );	
	
	
	$command .= " --bamOutput ".$gvcf.".bam" if ( (! $usegatk4) && $emitRealBam );
		$command .= " --bam-output ".$gvcf.".bam" if ( $usegatk4 && $emitRealBam );
	
	if($isMultisample){
		$command .= " -o $outfile" if ( ! $usegatk4 );
			$command .= " -O $outfile" if ( $usegatk4 );
		
		$command .= " --output_mode EMIT_ALL_SITES" if ( $outputAll && !$usegatk4 );
			$command .= " --output-mode EMIT_ALL_SITES" if ( $outputAll && $usegatk4 );
	}else{
		#Resolution
		my $resolution = $useGVCFResolution ? "GVCF" : "BP_RESOLUTION";
		
		$command .= " --emitRefConfidence $resolution " if ( ! $usegatk4 );
			$command .= " --emit-ref-confidence $resolution " if ( $usegatk4 );
		
		
			$command .= " --min_base_quality_score 30 " if ( $usegatk4 && $Q30 );
		
		$command .= " --variant_index_type LINEAR --variant_index_parameter 128000 -o $gvcf" if ( ! $usegatk4 );
		
			$command .= " -O $gvcf" if ( $usegatk4 );
	}
	$logger->info("Running GATK HaplotypeCaller...");
	
}else{
	
	# GATK UnifiedGenotyper
	
	# TODO: support GATK4 or remove...we don't need it anymore indeed
	# For compatibility (GATK4 not supported)
	$command = "
			$java -XX:ParallelGCThreads=1 -Xmx$maxRam -jar $gatk \\
				-R $ref -I $infiles -nct $nct  \\
				-A DepthPerAlleleBySample -A InbreedingCoeff -A HaplotypeScore --logging_level $loglevel -U ALLOW_N_CIGAR_READS \\
		";
	
	$command .= "--forceActive \\\n"  if $forceActive;
	$command .= "-L \"$regionFile\" \\\n" if $regionFile ne "";
	
	$command .= "--output_mode EMIT_ALL_SITES \\\n" if $outputAll;
	$command .= "-T UnifiedGenotyper -nt $nt -glm BOTH -o $outfile";
	$logger->info("Running GATK UnifiedGenotyper...");
}

$command .= " >> $javalog 2>&1";

$logger->debug($command);
system($command);



# GenotypeGVCFs Stage 
# 		Produce VCF file
# 		Not needed if UnifiedGenotyper is used

unless($useUG || $isMultisample || $noVCF){
	 my $tabix = $params->{programs}->{tabix}->{path};
	 my $bgzip = $params->{programs}->{bgzip}->{path};

	 #Save time removing external bgzipping and tabix-ing
	 #Avoid tabix wherever possible.

	 ##zip gvcf file
	 #$command = "$bgzip -f $gvcf";
	 #$logger->debug($command);
	 #system($command);
	 
	 #index zipped gvcf file
	 #$gvcf .= ".gz";
	 #$command = "$tabix -fp vcf $gvcf";
	 #$logger->debug($command);
	 #system($command);
	
	 $command = "
		$java -XX:ParallelGCThreads=1 -Xmx$maxRam -jar $gatk \\
		-R $ref -T GenotypeGVCFs   \\
   		--variant $gvcf \\
   		-o $outfile" if ( ! $usegatk4 );

	 $command = "
		$gatk4 --java-options  \"-Xmx$maxRam\" GenotypeGVCFs\\
		-R $ref \\
   		--variant $gvcf \\
   		-O $outfile" if ( $usegatk4 );   		
   
   	 $command .= " >> $javalog 2>&1";
   
    $logger->debug($command);
	system($command);
   
	
}


=head1 NAME

varfilter_gatk.pl

=head1 SYNOPSIS

 varfilter_gatk.pl -i merged.rmdup.bam -l targets.bed -o ontarget.vcf -se hg19_test

=head1 DESCRIPTION

This is a wrapper script for GATK VariantCalling. By default GATK HaplotypeCaller is used and a outfile.gvcf file is
generated as preliminary result which is then translated into vcf. Default version GATK 4.x.
GATK3.x, and the outdated caller UnifiedGenotyper are still supported with appropriate flags.

=head1 OPTIONS

 -i	<list of infiles> or text file that includes bam files (must have .list as filetype); if empty outdir/merged.rmdup.bam will be taken as infile
 	if a directory is given: use all *.sort.bam files in this directory
 -l	<region.bed> BED file with regions for which variants should be called
 -r	<region> where variants should be called; DEPRECATED: region can now be defined via -l
 -o	<outfile>; required
 -se	<settings>; required
 -a	output VCF variant entries for all positions in <region> even if no alternative genotype is called; useful for comparing samples; only works for UnifiedGenotyper
 -ug	run GATK UnifiedGenotyper instead of HaplotypeCaller (default). Requires -gatk3 flag. <SOON TO BE DEPRECATED>
 -gatk3 use GATK3 (DEPRECATED)
 -gatk4 use GATK4 (DEFAULT). Not compatible with UnifiedGenotyper (SOON TO BE DEPRECATED option -ug ). 
 -nt	number of GATK data threads; default: 1 (only for UnifiedGenotyper)
 -nct	number of GATK CPU threads; default: 1; NOT SUPPORTED IF -aj IS CHOSEN; NOT SUPPORTED BY GATK4
 -m	max. memory for each Java GATK job; default 4g
 -fa	sets HaplotypeCaller's --forceActive option. If this option is selected no region will be skipped.
 -aj	is SGE Array job; if this flag is chosen, the script runs as a part of a SGE array job
		this means, that the environmental variable SGE_TASK_ID holds the number of this job
		The script will process only the chromosome in line SGE_TASK_ID of the dict file of the reference genome.
		It will use the chromsome name as prefix for all generated output files.
-ajb	is SGE Array job; if this flag is chosen, the script runs as a part of a SGE array job
		this means, that the environmental variable SGE_TASK_ID holds the number of this job
		The script will only process the file "SGE_TASK_ID.bed" in the folder $params->{settings}->{$settings}->{genome_splits}.
		SGE_TASK_ID will be used as a fileprefix.
 -v	don't generate a VCF file (i.e. only a gvcf file)
 -lf	<log file>; default: pipeline.log
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page

=head1 AUTHOR

Riccardo Berutti, Thomas Wieland

=cut
 
