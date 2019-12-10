#!/usr/bin/perl

# TODO FILTERGATK  GO STATIC

use strict;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;
use Cwd qw(abs_path);
use Pod::Usage;
use Log::Log4perl qw(get_logger :levels);

my $prog_path = dirname(abs_path($0));
require $prog_path."/Utilities.pm";

# Parameters
my $outfile       = "";
my $help          = 0;
my $man			  = 0;
my $params        = Utilities::getParams();
my $logfile       = "pipeline.log";
my $loglevel      = "INFO";
my $settings      = "";
my $infile        = "";
my $maxRam        = "6g";
my $usegnomad	  = 0;  #Todo change to 1 when it will be DEFAULT
my $usegermansnps = 0;  #Use GERMAN Genomes 

my $createModel   = 0;

# Static filter options (for non-human, small targets or legacy samples)
my $snpFilter = 0;
my $indelFilter = 0;

my $isUG = 0;

my $nodatafoundrerun = 0;
my $cmdline = join(" ", @ARGV);

# Dummy compatibility arguments:
my $dummy = 0;


# Getting options
GetOptions(
"i=s"  => \$infile,
"o=s"  => \$outfile, 
"se=s" => \$settings,
"m=s"  => \$maxRam,
"sf"   => \$snpFilter,
"if"   => \$indelFilter,
"lf=s" => \$logfile,
"ll=s" => \$loglevel,
"ug"   => \$isUG,
"nt"   => \$dummy,
"gh"   => \$dummy,
"w"    => \$dummy,
"wg"   => \$dummy,
"cont" => \$dummy,
"hp"   => \$dummy,
"maxGaussians" => \$dummy,
"disableInbreedingCoefficient" => \$dummy,
"customTranches" => \$dummy,
"custom95" => \$dummy,
"nodatafoundrerun" => \$dummy,
"gnomad"		=> \$usegnomad,
"deutschland" 	=> \ $usegermansnps,
"german" 		=> \ $usegermansnps,
"dzhk" 			=> \ $usegermansnps,
"h"    => \$help);

# Help text 
pod2usage( {-exitval => 0  ,-verbose => 1} ) if $help;
pod2usage( {-exitval => 0  ,-verbose => 2} ) if $man;
pod2usage( {-exitval => 1  ,-verbose => 1} ) if ($help == 1 || $outfile eq "" || $settings eq "" || $infile eq ""); 

# Logging init
Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();

# Java logging init
my $javalog = $outfile;
$javalog    =~ s/vcf$/log/;

# Deprecated arguments
$logger->info("Used one or more deprecated arguments. They are now ignored") if ( $dummy );

# Getting parameters
my $ref     = $params->{settings}->{$settings}->{reference};
my $gatk4    = $params->{programs}->{gatk4}->{path};
my $java    = $params->{programs}->{java}->{path};
my $tmp     = $params->{programs}->{gatk}->{tmpdir};

# Define command
my $command = "";

# IsMultisample
my $isMultiSample=0; 

# Get the number of samples from the VCF file
if ( $infile =~ /\.gz$/ )
{
	open IN, "zcat $infile | " or die "Can't open $infile!\n";
    while(<IN>){
    	next if $_ =~ /^##/;
        chomp;
        my @columns = split();
        $isMultiSample = 1 if @columns>17;      #more than 10 samples --> InbreedingCoeff can be used
        	#TODO:
        	#$isMultiSample = 1 if @columns>57;	#more than 50 samples --> Can create recal model
        last;
	}
}
else
{
	open IN, $infile or die "Can't open $infile!\n";
	while(<IN>){
		next if $_ =~ /^##/;
		chomp;
		my @columns = split();
		$isMultiSample = 1 if @columns>17;	#more than 10 samples --> InbreedingCoeff can be used
			#TODO:
			#$isMultiSample = 1 if @columns>57;	#more than 50 samples --> Can create recal model
		last;
	}
}



# Look if at least one SNP and INDEL FILE is defined #TODO: Add GnoMAD SNPS here
my $snpFiles   = ($params->{settings}->{$settings}->{HapMap} || $params->{settings}->{$settings}->{Omni1000G} || $params->{settings}->{$settings}->{gatksnps} || $params->{settings}->{$settings}->{SNPs1000G});
my $indelFiles = ($params->{settings}->{$settings}->{gatksnps} || $params->{settings}->{$settings}->{GoldStandardIndels});

# Set our version of Java in Path with priority so to overcome any current installation (required by GATK4)
$ENV{'PATH'}=dirname(abs_path($java)).":".$ENV{'PATH'};


#Fallback to static filtering for SNPs and INDELs if trackss are missing
if ( $snpFiles eq 0 )
{
	$logger->error("Missing at least one SNP track. Falling back to static SNP filter!");
	$snpFilter=1;
}

if ( $indelFiles eq 0 )
{
	$logger->error("Missing at least one INDEL track. Falling back to static INDEL filter!");
	$indelFilter=1;
}


# Command line
my $command="";


#######################################################
# SNPS:
# Performing recalibration or filtering SNPS
#######################################################

	# Building command for GATK4 - Model Build
	my $maxGaussians=6;

		$command ="$gatk4 --java-options \"-Xmx$maxRam -Xms$maxRam\" 	\\
						VariantRecalibrator 							\\
						--reference 	$ref 							\\
						--variant 		$infile 						\\
						--output 		$outfile.SNP.recal 				\\
						--tranches-file $outfile.SNP.tranches 			\\
						--rscript-file 	$outfile.SNP.plot.R 			\\
					    --verbosity		$loglevel  						\\
					    --trust-all-polymorphic							\\
					    --max-gaussians $maxGaussians \\
						-tranche  100.0 \\
						-tranche  99.95\\
						-tranche  99.9 \\
						-tranche  99.8 \\
						-tranche  99.6 \\
						-tranche  99.5 \\
						-tranche  99.4 \\
						-tranche  99.3 \\
						-tranche  99.0 \\
						-tranche  98.0 \\
						-tranche  97.0 \\
						-tranche  90.0 \\
						-mode SNP 		\\
							-an QD \\
							-an MQRankSum \\
							-an ReadPosRankSum \\
							-an FS \\
							-an MQ \\
							-an SOR \\
							-an DP \\
						";
						
						$command .= "-an InbreedingCoeff \\\n" if $isMultiSample;

			# Resource Tracks (IF DEFINED)
				$command .= "-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $params->{settings}->{$settings}->{HapMap} \\\n" 						if $params->{settings}->{$settings}->{HapMap};
				$command .= "-resource:omni,known=false,training=true,truth=true,prior=12.0 $params->{settings}->{$settings}->{Omni1000G} \\\n" 					if $params->{settings}->{$settings}->{Omni1000G};
				$command .= "-resource:1000G,known=false,training=true,truth=false,prior=10.0 $params->{settings}->{$settings}->{SNPs1000G} \\\n" 					if $params->{settings}->{$settings}->{SNPs1000G}  && ( ! $usegnomad  ) && ( ! $usegermansnps );
				$command .= "-resource:gnomad,known=false,training=true,truth=false,prior=10.0 $params->{settings}->{$settings}->{gnomADsnps} \\\n" 				if $params->{settings}->{$settings}->{gnomADsnps} && (   $usegnomad  ) && ( ! $usegermansnps );
				$command .= "-resource:gnomad,known=false,training=true,truth=false,prior=10.0 $params->{settings}->{$settings}->{germansnps} \\\n" 				if $params->{settings}->{$settings}->{germansnps} && ( ! $usegnomad  ) && (   $usegermansnps );
				$command .= "-resource:dbsnp,known=true,training=false,truth=false,prior=7.0 $params->{settings}->{$settings}->{gatksnps} \\\n" 					if $params->{settings}->{$settings}->{gatksnps};
				## Ready to introduce a CUSTOM dataset from US
				##$command .= "-resource:hmgu,known=true,training=false,truth=false,prior=1.0 /home/berutti/hail2VCF/singeltons_8000_LUEB0080G-2018-01-31.vcf \\\n" if $hmguset;

			# Add java logger
				$command .= " 2>&1 >> $javalog";


	#EXECUTE!

		# Execute (if not SNP static filter)
		if ( ! $snpFilter ){
			$logger->debug($command);
				if (&Utilities::executeCommand($command, "Running VariantRecalibrator for SNPs...", $logger)) {
				$logger->error("Error executing VariantRecalibrator for SNPs. Execute with less gaussians.");
					
					my $maxGaussiansNew = $maxGaussians-2;
					$command =~ s/\-\-max\-gaussians\ $maxGaussians/\-\-max\-gaussians\ $maxGaussiansNew /;
				
					if (&Utilities::executeCommand($command, "Running VariantRecalibrator for SNPs with $maxGaussiansNew gaussians...", $logger)) {
						$logger->error("Error executing VariantRecalibrator for SNPs with less gaussians. That's hopeless. Bye");
						#Exit 
						exit(-1);
					}
				}
			
		
			# Produce Rplots in PDF
			# GATK uses a deprecated argument (legend) which make the plotting fail. Remove it and plot 
			$command = "replace \",legend=FALSE\" \"\" -- $outfile.SNP.plot.R; Rscript $outfile.SNP.plot.R";
			system($command);

		}

#######################################################
# INDELs
# Performing recalibration or filtering INDELs
#######################################################

	# Building command for GATK4 - Model Build
	my $maxGaussiansIndels=4;

		$command ="$gatk4 --java-options \"-Xmx$maxRam -Xms$maxRam\" 	\\
						VariantRecalibrator 							\\
						--reference 	$ref 							\\
						--variant 		$infile 						\\
						--output 		$outfile.INDEL.recal 			\\
						--tranches-file $outfile.INDEL.tranches 		\\
					    --verbosity		$loglevel  						\\
					    --trust-all-polymorphic							\\
						-tranche  100.0 \\
						-tranche  99.95\\
						-tranche  99.9 \\
						-tranche  99.5 \\
						-tranche  99.0 \\
						-tranche  97.0 \\
						-tranche  96.0 \\
						-tranche  95.0 \\
						-tranche  94.0 \\
						-tranche  93.5 \\
						-tranche  93.0 \\
						-tranche  92.0 \\
						-tranche  91.0 \\
						-tranche  90.0 \\
						-mode INDEL 		\\
							-an FS \\
							-an ReadPosRankSum \\
							-an MQRankSum \\
							-an QD \\
							-an SOR \\
							-an DP \\
						--max-gaussians $maxGaussiansIndels \\
						";

			$command .= "-an InbreedingCoeff \\\n" if $isMultiSample;

			$command .= "--resource:mills,known=false,training=true,truth=true,prior=12.0 $params->{settings}->{$settings}->{GoldStandardIndels} \\\n" if $params->{settings}->{$settings}->{GoldStandardIndels};
			$command .= "--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $params->{settings}->{$settings}->{gatksnps} \\\n" if $params->{settings}->{$settings}->{gatksnps};
			
			$command .= " 2>&1 >> $javalog";
			
	#EXECUTE!

		# Execute (if not INDEL static filter)
		if ( ! $indelFilter )
		{
			$logger->debug($command);
			if (&Utilities::executeCommand($command, "Running VariantRecalibrator for INDELs...", $logger)) {
				$logger->error("Error executing VariantRecalibrator for INDELs");
				#Exit 
				exit(-1);
			}
		}	
			
#######################################################
# APPLY RECALIBRATION
# Apply recalibration table for SNPs and INDELs sequentially
#######################################################			

			
		# Apply recalibrations for SNPs (unless SNP Static Filter )
			
		if ( ! $snpFilter)
		{
			$command ="$gatk4 --java-options \"-Xmx$maxRam -Xms$maxRam\" 	\\
					ApplyVQSR \\
					--reference $ref  \\
					--variant $infile \\
					-mode SNP \\
					--truth-sensitivity-filter-level 99.7 \\
					--recal-file $outfile.SNP.recal \\
					--tranches-file $outfile.SNP.tranches \\
					--output $outfile.tmp.vcf \\
					--verbosity $loglevel \\
					";
		
			$command .= " 2>&1 >> $javalog";
		
			$logger->debug($command);
			$logger->info("Apply Recalibration for SNPs...");
			
			if (&Utilities::executeCommand($command, "Running VariantRecalibrator for SNPs...", $logger)) {
					$logger->error("Error executing ApplyRecalibration for SNPs - No data found encountered");
					exit(-1);
			}
			
		}
		else
		{
			# SNP Static	
			
			#run SNP filter  ##### NOTE: logging level is set to ERROR because otherwise warnings will pop up if a field is missing
			$command = "$gatk4 --java-options \"-Xmx$maxRam -Xms$maxRam\" \\
				VariantFiltration \\
				-R $ref \\
				-V $infile \\
				-O $outfile.tmp.vcf  --verbosity ERROR  \\
				--filter-expression \"vc.isSNP() && QD < 2.0\" --filter-name \"QDfilter\" \\
				--filter-expression \"vc.isSNP() && MQ < 40.0\" --filter-name \"MQfilter\" \\
				--filter-expression \"vc.isSNP() && FS > 60.0\" --filter-name \"FSfilter\" \\
				--filter-expression \"vc.isSNP() && MQRankSum < -12.5\" --filter-name \"MQRankSumfilter\" \\
				--filter-expression \"vc.isSNP() && ReadPosRankSum < -8.0\" --filter-name \"ReadPosRankSumfilter\" ";
	
				$command .= "--filter-expression \"vc.isSNP() && HaplotypeScore > 13.0\" --filter-name \"HaplotypeScorefilter\" " if $isUG;

			$command .= " 2>&1 >> $javalog";
			
			$logger->debug($command);
			$logger->info("Apply Static Filter for SNPs...");
			
			if (&Utilities::executeCommand($command, "Running Static Filter for SNPs...", $logger)) {
					$logger->error("Error executing Static Filter for SNPs");
					exit(-1);
			}
			
			
		}
			
		# Apply recalibrations for INDELs (unless INDEL Static Filter)
		
		
		if ( ! $indelFilter )
		{
			$command ="$gatk4 --java-options \"-Xmx$maxRam -Xms$maxRam\" 	\\
					ApplyVQSR \\
					--reference $ref  \\
					--variant $outfile.tmp.vcf \\
					-mode INDEL \\
					--truth-sensitivity-filter-level 99.7 \\
					--recal-file $outfile.INDEL.recal \\
					--tranches-file $outfile.INDEL.tranches \\
					--output $outfile \\
					--verbosity $loglevel \\
					";
					
			$command .= " 2>&1 >> $javalog";
		
			$logger->debug($command);
			$logger->info("Apply Recalibration for INDELs...");
			
			if (&Utilities::executeCommand($command, "Running VariantRecalibrator for INDELs...", $logger)) {
					$logger->error("Error executing ApplyRecalibration for INDELs - No data found encountered");
					exit(-1);
			}
		}
		else
		{
			# INDEL Static
			#run INDEL filter ##### NOTE: logging level is set to ERROR because otherwise warnings will pop up if a field is missing
			$command = "$gatk4 --java-options \"-Xmx$maxRam -Xms$maxRam\" \\
				VariantFiltration \\
				-R $ref \\
				-V $outfile.tmp.vcf \\
				-O $outfile --verbosity ERROR \\
				--filter-expression \"vc.isIndel() && QD < 2.0\" --filter-name \"QDfilter\" \\
				--filter-expression \"vc.isIndel() && FS > 200.0\" --filter-name \"FSfilter\" \\
				--filter-expression \"vc.isIndel() && ReadPosRankSum < -20.0\" --filter-name \"ReadPosRankSumfilter\" ";
		
			$command .= "--filter-expression \"vc.isIndel() && InbreedingCoeff < -0.8\" --filter-name \"InbreedingCoefffilter\" " if $isMultiSample;
			$command .= " >> $javalog 2>&1";
			
			$logger->debug($command);
			$logger->info("Apply Static Filter for INDELs...");
			
			if (&Utilities::executeCommand($command, "Running Static Filter for INDELs...", $logger)) {
					$logger->error("Error executing Static Filter for INDELs");
					exit(-1);
			}
			
			
		}	
			



# Remove tmp files

	unlink("$outfile.tmp.vcf");
	unlink("$outfile.tmp.vcf.idx");
	
# Remove intermediate. For debug keep them, by commenting down here, they are helpful!

	unlink("$outfile.INDEL.vcf");
	unlink("$outfile.INDEL.vcf.idx");
	unlink("$outfile.INDEL.tranches");
	unlink("$outfile.INDEL.recal");
	unlink("$outfile.INDEL.recal.idx");
	#
	unlink("$outfile.SNP.filtered.vcf");
	unlink("$outfile.SNP.filtered.vcf.idx");
	unlink("$outfile.SNP.vcf");
	unlink("$outfile.SNP.vcf.idx");
	unlink("$outfile.SNP.tranches");
	unlink("$outfile.SNP.recal");
	unlink("$outfile.SNP.recal.idx");



=head1 NAME

filterGATK4.pl

=head1 SYNOPSIS

 filterGATK4.pl -i <input.vcf> -o <output.vcf> -se <settings> [ -m <maxram, def 4g> ] 
 
=head1 DESCRIPTION

This script filters a GATK VCF file by running VariantRecalibrator. 
Uses GATK4, does not implement any static filter. Keeping it simple. 
Todo: for static or no files fallback to filterGATK.  

=head1 OPTIONS

 -i	<infile.vcf>; required
 -o	<outfile.vcf>; required
 -se	<settings>; required
 -m	maximum RAM for java virtual machine; default: $maxRam
 -lf	log file; default: pipeline.log
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -ug	legacy UnifiedGenotyper support, deprecated. Give option if UG data
 -sf 	static filter for SNPs
 -if	static filter for Indels
 -gnomad use gnomad SNPs
 -dzhk/-german use internal dzhk test SNP set (only on IHGSEQ12!)
 -h	this help
 -man  man mode

=head1 AUTHOR

Riccardo Berutti 

=cut




