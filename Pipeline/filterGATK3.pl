#!/usr/bin/perl

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
my $maxRam        = "4g";
my $nt		      = -1;
my $isMultiSample = 0;
my $isHC		  = 0;
my $snpFilter     = 0;
my $indelFilter   = 0;
my $isWG		  = 0;
my $lowVariants   = 0;
my $cont          = 0;
my $usegnomad	  = 0;  #Todo change to 1 when it will be DEFAULT
my $usegermansnps = 0;  #Use GERMAN Genomes 
my $maxGaussians  = -1;
my $disableInbreedingCoefficient=0;
my $customTranches = 0;
my $customTranches95 = 0;
my $annotateHomopolymerRun = 0;

my $nodatafoundrerun = 0;
my $cmdline = join(" ", @ARGV);

# Getting options
GetOptions(
"i=s"  => \$infile,
"o=s"  => \$outfile, 
"se=s" => \$settings,
"m=s"  => \$maxRam,
"nt=s" => \$nt,
"gh"   => \$isHC,
"sf"   => \$snpFilter,
"if"   => \$indelFilter,
"w"    => \$isWG,
"wg"   => \$isWG,
"hp"   => \$annotateHomopolymerRun,
"l"    => \$lowVariants,
"lf=s" => \$logfile,
"ll=s" => \$loglevel,
"cont" => \$cont,
"gnomad"=> \$usegnomad,
"deutschland" => \ $usegermansnps,
"german" => \ $usegermansnps,
"hmgu" => \ $usegermansnps,
"maxGaussians=s" => \$maxGaussians,
"disableInbreedingCoefficient" => \$disableInbreedingCoefficient,
"customTranches" => \$customTranches,
"custom95" => \$customTranches95,
"nodatafoundrerun" => \$nodatafoundrerun,
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


# Getting parameters
my $ref     = $params->{settings}->{$settings}->{reference};
my $gatk    = $params->{programs}->{gatk}->{path};
my $java    = $params->{programs}->{java}->{path};
my $tmp     = $params->{programs}->{gatk}->{tmpdir};

# Define command
my $command = "";

# Determining number of threads (adapt to the number of slots given by SGE if applicable)
if($nt == -1){
	$nt = 1;
	if($ENV{NSLOTS}){		#get number of slots given by SGE
		$nt = $ENV{NSLOTS};
	}
}

# Get the number of samples from the VCF file
open IN, $infile or die "Can't open $infile!\n";
while(<IN>){
	next if $_ =~ /^##/;
	chomp;
	my @columns = split();
	$isMultiSample = 1 if @columns>17;	#more than 10 samples --> InbreedingCoeff can be used
	last;
}

# Look if at least one SNP and INDEL FILE is defined #TODO: Add GnoMAD SNPS here
my $snpFiles   = ($params->{settings}->{$settings}->{HapMap} || $params->{settings}->{$settings}->{Omni1000G} || $params->{settings}->{$settings}->{gatksnps} || $params->{settings}->{$settings}->{SNPs1000G});
my $indelFiles = ($params->{settings}->{$settings}->{gatksnps} || $params->{settings}->{$settings}->{GoldStandardIndels});


# If no data found fallback in the settings that most likely will run
$lowVariants  = 1 if $nodatafoundrerun;
$maxGaussians = 2 if $nodatafoundrerun;


#######################################################
# SNPS:
# Performing recalibration or filtering SNPS
#######################################################

if(!$snpFilter && $snpFiles){

	# RECALIBRATION REQUIRED (AND POSSIBLE)
	# Run VariantRecalibrator


		#Building command for GATK<4
		$command ="$java -XX:ParallelGCThreads=1 -Xmx$maxRam -jar $gatk -T VariantRecalibrator \\";
		# Common covariates
			# QD
			# MQRankSum
			# ReadPosRankSum
		    $command .="	-R $ref -input $infile -recalFile $outfile.SNP.recal -tranchesFile $outfile.SNP.tranches -rscriptFile $outfile.SNP.plot.R \\
							--logging_level $loglevel --disable_auto_index_creation_and_locking_when_reading_rods \\
							-nt $nt \\
							-mode SNP \\
							-an QD \\
							-an MQRankSum \\
							-an ReadPosRankSum \\
							-an FS \\
						";

			# DEPRECATED STUFF (Notice to keep here)
						# MQ 
						# $command .= "-an MQ \\\n" if $isWG;
						
						# By default GATK>=3.8 minNumBadVariants is deprecated 
						# MaxGaussians CAN (but not MUST) be omitted for WG as in the GATK new reference information
							#$command .= "--maxGaussians 4 \\\n" if $isWG;
							#$command .= "--maxGaussians 6 --minNumBadVariants 4000 \\\n" if ( $lowVariants || (!$isMultiSample && !$isWG) );
							#$command .= "--maxGaussians 4 --minNumBadVariants 4000 \\\n" if ( $lowVariants || (!$isMultiSample && !$isWG) ); 


			# HaplotypeCaller dependent covariate (HaplotypeScore)			
			$command .= "-an HaplotypeScore \\\n" unless $isHC;


			# Optional Covariates (-switches dependent)

				# SWITCH TO DISABLE
					# InbreedingCoefficient
					$command .= "-an InbreedingCoeff \\\n" if ($isMultiSample && (! $disableInbreedingCoefficient));
			
				# SWITCH TO ENABLE
					# Custom Tranches
					$command .= " -tranche 100 -tranche 99.9 -tranche 99.5 -tranche 99 -tranche 90 " if $customTranches;
					$command .= " -tranche 100 -tranche 99.9 -tranche 99.5 -tranche 99 -tranche 90 -tranche 95.0 " if $customTranches95;
				
					# HomopolymerRun
					$command .= "-an HRun \\\n" if $annotateHomopolymerRun;
					
					# MaxGaussians
					$command .= "--maxGaussians ".($maxGaussians != -1 ? $maxGaussians : 4)." \\\n" if ( $lowVariants || (!$isMultiSample && !$isWG) ); #No more minNumBadVariants required		
						
					$command .= "--maxGaussians $maxGaussians \\\n"			if (($maxGaussians != -1) && !(( $lowVariants || (!$isMultiSample && !$isWG) )));
			
			# Experiment dependent:
				# WG
					$command .= "-an DP \\\n" 	if $isWG;


			# Resource Tracks (IF DEFINED)
				$command .= "-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $params->{settings}->{$settings}->{HapMap} \\\n" 						if $params->{settings}->{$settings}->{HapMap};
				$command .= "-resource:omni,known=false,training=true,truth=true,prior=12.0 $params->{settings}->{$settings}->{Omni1000G} \\\n" 					if $params->{settings}->{$settings}->{Omni1000G};
				$command .= "-resource:1000G,known=false,training=true,truth=false,prior=10.0 $params->{settings}->{$settings}->{SNPs1000G} \\\n" 					if $params->{settings}->{$settings}->{SNPs1000G}  && ( ! $usegnomad  ) && ( ! $usegermansnps );
				$command .= "-resource:gnomad,known=false,training=true,truth=false,prior=10.0 $params->{settings}->{$settings}->{gnomADsnps} \\\n" 				if $params->{settings}->{$settings}->{gnomADsnps} && (   $usegnomad  ) && ( ! $usegermansnps );
				$command .= "-resource:gnomad,known=false,training=true,truth=false,prior=10.0 $params->{settings}->{$settings}->{germansnps} \\\n" 				if $params->{settings}->{$settings}->{germansnps} && ( ! $usegnomad  ) && (   $usegermansnps );
				$command .= "-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $params->{settings}->{$settings}->{gatksnps} \\\n" 					if $params->{settings}->{$settings}->{gatksnps};
				## Ready to introduce a CUSTOM dataset from us
				##$command .= "-resource:hmgu,known=true,training=false,truth=false,prior=1.0 /home/berutti/hail2VCF/singeltons_8000_LUEB0080G-2018-01-31.vcf \\\n" if $hmguset;


			#EXECUTE!

			# Add java logger
			$command .= " 2>&1 >> $javalog";
		

			# Checking if intermediate SNP file was generated already, if $cont continue requested it will keep it
			if ( $cont eq 0 || ( ! -e "$outfile.SNP.recal" ))
			{
				$logger->debug($command);
				if (&Utilities::executeCommand($command, "Running VariantRecalibrator for SNPs...", $logger)) {
					$logger->error("Error executing VariantRecalibrator for SNPs");
					# Here catch the no data found error and restart the whole process BUT ONLY ONCE please!
					
					$logger->info("GATK VariantRecalibrator useless no data found error. Why is it so sloppy? Relaunching with lenient parameters.") if (! $nodatafoundrerun );
					system ( "perl ".$prog_path."/filterGATK.pl -nodatafoundrerun ".$cmdline ) if (! $nodatafoundrerun );
					
					$logger->info("GATK VariantRecalibrator couldn't help himself and failed again. Check manually. ") if ( $nodatafoundrerun );
					# Exit anyway: wether 
					exit(-1);
				}
			}
			else
			{
				$logger->info("$outfile.SNP.recal exists - skipping");
			}
			 
			
			# Produce Rplots in PDF
			# GATK uses a deprecated argument (legend) which make the plotting fail. Remove it and plot 
			$command = "replace \",legend=FALSE\" \"\" -- $outfile.SNP.plot.R; Rscript $outfile.SNP.plot.R";
			system($command);


			# Apply recalibrations for SNPs
			$command = "
							$java -XX:ParallelGCThreads=1 -Xmx$maxRam -jar $gatk \\
							-R $ref -T ApplyRecalibration -input $infile -recalFile $outfile.SNP.recal -tranchesFile $outfile.SNP.tranches \\
							--logging_level $loglevel -o $outfile.tmp.vcf --disable_auto_index_creation_and_locking_when_reading_rods \\
							-mode SNP --ts_filter_level 99.0";
		
			$command .= " 2>&1 >> $javalog";
		
			$logger->debug($command);
			$logger->info("Apply Recalibration for SNPs...");
			if ( $cont eq 0 || ( ! -e "$outfile.SNP.tranches" ))
			{
				if (&Utilities::executeCommand($command, "Running VariantRecalibrator for SNPs...", $logger)) {
					$logger->error("Error executing ApplyRecalibration for SNPs - No data found encountered");
					exit(-1);
				}
				#OLD STX
				#system($command);
			}
			else
			{
				$logger->info("$outfile.SNP.tranches exists - skipping");
			}


}else{

	# NO RECALIBRATION REQUIRED (or POSSIBLE)
	# HARD SNP FILTERING 	
		
		#run SNP filter  ##### NOTE: logging level is set to ERROR because otherwise warnings will pop up if a field is missing
		$command = "$java -XX:ParallelGCThreads=1 -Xmx$maxRam -jar $gatk \\
			-R $ref -T VariantFiltration \\
			--variant $infile \\
			-o $outfile.tmp.vcf  --logging_level ERROR --disable_auto_index_creation_and_locking_when_reading_rods \\
			--filterExpression \"vc.isSNP() && QD < 2.0\" --filterName \"QDfilter\" \\
			--filterExpression \"vc.isSNP() && MQ < 40.0\" --filterName \"MQfilter\" \\
			--filterExpression \"vc.isSNP() && FS > 60.0\" --filterName \"FSfilter\" \\
			--filterExpression \"vc.isSNP() && MQRankSum < -12.5\" --filterName \"MQRankSumfilter\" \\
			--filterExpression \"vc.isSNP() && ReadPosRankSum < -8.0\" --filterName \"ReadPosRankSumfilter\" ";

	$command .= "--filterExpression \"vc.isSNP() && HaplotypeScore > 13.0\" --filterName \"HaplotypeScorefilter\" " unless $isHC;

	$command .= " 2>&1 >> $javalog";

	$logger->debug($command);
	$logger->info("Running SNP filter...");
	system($command);
	
}

#######################################################
# INDELs
# Performing recalibration or filtering INDELs
#######################################################

if(!$indelFilter && $indelFiles){
	# RECALIBRATION REQUIRED (AND POSSIBLE)
	# Run VariantRecalibrator
	
			$command ="
		$java -XX:ParallelGCThreads=1 -Xmx$maxRam -jar $gatk \\
		-R $ref -T VariantRecalibrator -input $outfile.tmp.vcf -recalFile $outfile.INDEL.recal -tranchesFile $outfile.INDEL.tranches \\
		--logging_level $loglevel --disable_auto_index_creation_and_locking_when_reading_rods \\
		-mode INDEL \\
		-an QD -an FS -an ReadPosRankSum -an MQRankSum -nt $nt \\
		";
			#TODO: check SOR
			$command .= "-an SOR \\\n";
			$command .= "-an HRun \\\n" if $annotateHomopolymerRun;
			$command .= "-an DP \\\n" if $isWG;
			$command .= "-an InbreedingCoeff \\\n" if ($isMultiSample && (! $disableInbreedingCoefficient));	
			$command .= "--maxGaussians ".($maxGaussians != -1 ? $maxGaussians : 4)." \\\n" if ( $lowVariants || !$isMultiSample );  # --minNumBadVariants 5000 No more required
			# Custom Tranches
            $command .= " -tranche 100 -tranche 99.9 -tranche 99.5 -tranche 99 -tranche 90 " if $customTranches;

			$command .= "-resource:mills,known=false,training=true,truth=true,prior=12.0 $params->{settings}->{$settings}->{GoldStandardIndels} \\\n" if $params->{settings}->{$settings}->{GoldStandardIndels};
			$command .= "-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $params->{settings}->{$settings}->{gatksnps} \\\n" if $params->{settings}->{$settings}->{gatksnps};
			
			$command .= " 2>&1 >> $javalog";
			
			$logger->debug($command);
			$logger->info("Running VariantRecalibrator for INDELs...");
			system($command);
		
			
			#apply recalibrations for INDELs
			$command ="
		$java -XX:ParallelGCThreads=1 -Xmx$maxRam -jar $gatk \\
		-R $ref -T ApplyRecalibration -input $outfile.tmp.vcf -recalFile $outfile.INDEL.recal -tranchesFile $outfile.INDEL.tranches \\
		--logging_level $loglevel -o $outfile --disable_auto_index_creation_and_locking_when_reading_rods \\
		-mode INDEL --ts_filter_level 99.0";
		
			$command .= " 2>&1 >> $javalog";
		
			$logger->debug($command);
			$logger->info("Apply Recalibration for INDELs...");
			system($command);
			

}else{
	# NO RECALIBRATION REQUIRED (or POSSIBLE)
	# HARD INDEL FILTERING 	
		
			#run INDEL filter ##### NOTE: logging level is set to ERROR because otherwise warnings will pop up if a field is missing
			$command = "$java -XX:ParallelGCThreads=1 -Xmx$maxRam -jar $gatk \\
		-R $ref -T VariantFiltration \\
		--variant $outfile.tmp.vcf \\
		-o $outfile --logging_level ERROR --disable_auto_index_creation_and_locking_when_reading_rods \\
		--filterExpression \"vc.isIndel() && QD < 2.0\" --filterName \"QDfilter\" \\
		--filterExpression \"vc.isIndel() && FS > 200.0\" --filterName \"FSfilter\" \\
		--filterExpression \"vc.isIndel() && ReadPosRankSum < -20.0\" --filterName \"ReadPosRankSumfilter\" ";
		
			$command .= "--filterExpression \"vc.isIndel() && InbreedingCoeff < -0.8\" --filterName \"InbreedingCoefffilter\" " if $isMultiSample;
			
			$command .= " >> $javalog 2>&1";
		
			$logger->debug($command);
			$logger->info("Running INDEL filter...");
			system($command);
			
}



# Remove tmp files

	unlink("$outfile.tmp.vcf");
	unlink("$outfile.tmp.vcf.idx");
	
# Keep the others now, they are helpful!

	#unlink("$outfile.INDEL.vcf");
	#unlink("$outfile.INDEL.vcf.idx");
	#unlink("$outfile.INDEL.tranches");
	#unlink("$outfile.INDEL.recal");
	#unlink("$outfile.INDEL.recal.idx");
	#
	#unlink("$outfile.SNP.filtered.vcf");
	#unlink("$outfile.SNP.filtered.vcf.idx");
	#unlink("$outfile.SNP.vcf");
	#unlink("$outfile.SNP.vcf.idx");
	#unlink("$outfile.SNP.tranches");
	#unlink("$outfile.SNP.recal");
	#unlink("$outfile.SNP.recal.idx");



=head1 NAME

filterGATK.pl

=head1 SYNOPSIS

 filterGATK.pl -i <input.vcf> -o <output.vcf> -se <settings> [ -m <maxram, def 4g> ] [-nt <num_threads>] [-gh] [-w] 
 
=head1 DESCRIPTION

This script filters a GATK VCF file either by running VariantRecalibrator if the files are provided (currently only for humans) or
a static filter in other cases (or if chosen, e.g. because of sample size).

=head1 OPTIONS

 -i	<infile.vcf>; required
 -o	<outfile.vcf>; required
 -se	<settings>; required
 -m	maximum RAM for java virtual machine; default: $maxRam
 -nt	number of GATK data threads; default: take from SGE environment variable NSLOTS (only for VariantRecalibrator)
 -gh	calls have been produced by GATK HaplotypeCaller instead of UnifiedGenotyper (different set of parameters)
 -sf	run static filter for SNPs; default: run VariantRecalibrator
 -if	run static filter for indels; default: run VariantRecalibrator
 -w	the VCF file is from a whole genome experiment --> DP will be used for variant recalibration
 -wg the VCF file is from a whole genome experiment --> DP will be used for variant recalibration (syn)
 -hp annotate homopolymer run
 -l	low number of variants; --> use --maxGaussians & --minNumBadVariants settings although the input is a multisample or whole genome file
	 useful if filtering fails due to low number of variants
 -lf	log file; default: pipeline.log
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -cont continue: do not overwrite SNP temp files if they were already created
 -hmguset use internal HMGU test SNP set (only on IHGSEQ12!)
 -disableInbreedingCoefficient
 -customTranches
 -h	this help
 -man  man mode

=head1 AUTHOR

Riccardo Berutti, Thomas Wieland 

=cut
 




