#!/usr/bin/perl

use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use DBI;
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);
use List::MoreUtils qw/ uniq /;
use Pod::Usage;

my $prog_path = dirname( abs_path($0) );
my $extra_prog_path = $prog_path."/../OtherScripts";
require $prog_path . "/Utilities.pm";


# database
#my $dbh           = "";
#my $sql           = "";
#my $sth           = "";

# Settings
#my $settings  	  = "hg19_wholegenome";

my $vcf           = "";
my $ped           = "";
my $trait         = "";
my $covariates    = "";
my $outdir        = ".";
my $rerun         = 0;

my $allvariants   = 0;	# If set removes the PASS filter (useful to process unfiltered VCFs)

# Cores:
my $cores         = 34;

# Filtering Parameters
my $minDP         = 30;
my $minGQ		  = 30;
my $minMQ		  = 50;
my $minCallRate   = 0.95;
	my $noCallFraction  = 1 - $minCallRate;
	

# Better-to-be-fixed parameters
my $minMAF_Kinship = 0.01;

	#LOF
	my $maxMAF_LOF 		= 0.05;
	
	#LOFrare
	my $maxMAF_LOFrare 	= 0.01;

	#Missense
	my $minMAF_missense	= 0.01;
	my $maxMAF_missense	= 0.05;
	
	#Missenserare
	my $maxMAF_missenserare = 0.01;

	#Nonsyn
	my $maxMAF_nonsyn 	= 0.01;

# What to do one or both
my $burdentest    = 0;
my $singlevarianttest = 0;

my $help          = 0;
my $man           = 0;

# Logger
my $logfile       = "SCREEN";
my $loglevel      = "INFO";

# Util vars:
my $command       = "";

# TODO: 
# Possible evolution
# take multisample VCF from ms calling, and trait as option

GetOptions(
	"msvcf=s"   => \$vcf,
	"vcf=s"     => \$vcf,
	"ped=s"		=> \$ped,
	"trait=s"   => \$trait,
	"cov=s"		=> \$covariates,
	"outdir=s"  => \$outdir,
	"burden"    => \$burdentest,
	"singlevariant"    => \$singlevarianttest,
	"allvariants"=> \$allvariants,
	"mincallrate=s"=> \$minCallRate,
	"minDP=s"		=> \$minDP,
	"minGQ=s"		=> \$minGQ,
	"minMQ=s"		=> \$minMQ,
	"rerun"     => \$rerun,
	"lf=s"      => \$logfile,
	"ll=s"      => \$loglevel,
	"h"         => \$help,
	"man"       => \$man
);


# At least something
pod2usage( { -exitval => 1, -verbose => 1 } ) if ( 
													( ( $burdentest || $singlevarianttest ) == 0 ) || 
													( $vcf eq "" ) ||
													( $ped eq ""   ) ||
													( $trait eq "" )
													);

pod2usage( { -exitval => 0, -verbose => 1 } ) if $help;
pod2usage( { -exitval => 0, -verbose => 2 } ) if $man;


my $params = Utilities::getParams();

# Software
my $EPACTS       = $params->{programs}->{EPACTS}->{path};
my $intersectBed = $params->{programs}->{bedtools}->{intersectbed};
my $filterVCFbyAD= $extra_prog_path."/filterVCFbyAD.pl";
my $rscript      = $params->{programs}->{Rscript}->{path};

# Third class utils:
my $bgzip = $params->{programs}->{bgzip}->{path};
my $tabix = $params->{programs}->{tabix}->{path};

# Required tracks ( to be improved - this is FAR too old)
my $LCR_file = "/data/isilon/seq/analysis/single_analysis/MultiSampleCall/jul_2014/EPACTS/allSamples/LCR-hs37d5.bed";


# Might be useful
#my $exomedb                   = $params->{settings}->{$settings}->{exomedb}->{database};
#my $snvTable                  = $params->{settings}->{$settings}->{exomedb}->{snvtable};
#my $snvsampleTable            = $params->{settings}->{$settings}->{exomedb}->{snvsampletable};

#my $coredb                    = $params->{coredb}->{database};
#my $sampleTable               = $params->{coredb}->{sampletable};

Utilities::initLogger( $logfile, $loglevel );
my $logger = Utilities::getLogger();



# Recalculate quantities needed
$noCallFraction  = 1 - $minCallRate;



# If my outdir is unchanged give a short warning that I'm creating stuff here
if ( $outdir eq "." )
{
    $logger->info("WARNING: no outfolder specified, writing to current working directory. Press CTRL-C to stop");
    sleep(5);
}


# Create EPACTS folder structure:
my $outEPACTS=$outdir."/EPACTS";
	# Subfolders
	my $outEPACTS_burden=$outEPACTS."/burden";
	my $outEPACTS_single=$outEPACTS."/single_variant";

	# Check that it does not exist already
	if ( (-e "$outdir/EPACTS") && (! $rerun) )
	{
		$logger->error("ERROR: $outdir/EPACTS folder exists");
		exit(-9);
	}
	
	# Make
	mkdir $outEPACTS        or die $logger->error("Cannot create folder $outEPACTS. Check.")		if ( ! $rerun || ($rerun && ! -e $outEPACTS) );
	mkdir $outEPACTS_burden or die $logger->error("Cannot create folder $outEPACTS_burden. Check.") if ($burdentest && ( ! $rerun || ($rerun && ! -e $outEPACTS_burden) ) );
	mkdir $outEPACTS_single or die $logger->error("Cannot create folder $outEPACTS_single. Check.") if ($singlevarianttest && ( ! $rerun || ($rerun && ! -e $outEPACTS_single) ) );


#
# Log parameters:
#
my $logger_message_params = "\n".
			"epacts\t".$EPACTS."\n".
			"input_file\t".$vcf."\n".
			"output_folder\t".$outEPACTS."\n".
			"locomplexityregionsfile\t".$LCR_file."\n".
			"selected_variants\t".( $allvariants ? "ALL" : "PASS" )."\n".
			"miDP\t".$minDP."\n".
			"minGQ\t".$minGQ."\n".
			"minMQ\t".$minMQ."\n".
			"minCallRate\t".$minCallRate."\n".
			"minMAFforKinship\t".$minMAF_Kinship."\n".
			"maxMAFLof\t".$maxMAF_LOF."\n".
			"minMAFmissense\t".$minMAF_missense."\n".
			"maxMAFmissense\t".$maxMAF_missense."\n".
			"maxMAFnonsyn\t".$maxMAF_nonsyn."\n".
			"covariates\t".$covariates."\n";

#Print and log:
print "$logger_message_params";	
$logger->info($logger_message_params);


# Common steps:

	# Filter variants
		# Out File
			my $filtered = $outEPACTS."/filtered.vcf.gz";
	
		# Filter for PASS, Remove Low complexity regions, Filter for minimum depth $minDP, GQ $minGQ, 100%ref, callrate 1-$minCallRate
			my $cat = ( $vcf =~ /vcf.gz$/ ) ? "zcat" : "cat";

			# Max No CALL Fraction
			my $maxNoCall = 1 - $minCallRate;
		
			# If allvariants then DO NOT REQUIRE PASS filter set
			if ( $allvariants )
			{
				$command = "$cat $vcf | \
						sed 's/^chr//'| \
						intersectBed -a - -b $LCR_file -v -header  |  \
						/usr/local/packages/seq/bcftools-1.6/bcftools filter -e \"INFO/MQ <= $minMQ\" - | \
						perl /data/isilon/users/scripts/eclipse_workspace_berutti/OtherScripts/filterVCFbyAD.pl -d $minDP -g $minGQ -o $filtered -rn -ncf $noCallFraction \
						";
			}
			else
			{
				$command = "$cat $vcf | \
						awk '{if( \$1 ~ /^#/ || \$7 ~ /PASS/){print}}' | \
						sed 's/^chr//'| \
						intersectBed -a - -b $LCR_file -v -header  |  \
						/usr/local/packages/seq/bcftools-1.6/bcftools filter -e \"INFO/MQ <= $minMQ\" - | \
						perl /data/isilon/users/scripts/eclipse_workspace_berutti/OtherScripts/filterVCFbyAD.pl -d $minDP -g $minGQ -o $filtered -rn -ncf $noCallFraction \
						";
			}


			# If rerun and filtered file not existing then skip creation
			if ( (! $rerun) || ( $rerun && ! -e $filtered ) )
			{ 	
		    	if (&Utilities::executeCommand($command, "Filtering vcf file: $vcf", $logger)) {
		    		$logger->error("Error while filtering vcf file: $vcf");
		    		exit(-1);
		    	}
			}
			else
			{
				$logger->info("RERUN: Filtered file exists. Skipping.")				
			}
		    
		    
		    # If rerun and filtered file index not existing then skip creation 
		    if ( (! $rerun) || ( $rerun && ! -e $filtered.".tbi" ) )
		    {
				$command = "tabix -f $filtered";
				if (&Utilities::executeCommand($command, "Indexing filtered vcf file: $filtered", $logger)) {
		    		$logger->error("Error indexing filtered vcf file: $filtered");
		    		exit(-1);
		    	}
		    }

	# Annotate variants
		# Out File
			my $annotated = $outEPACTS."/annotated.vcf.gz";
	
		# Annotate

			# If rerun and annotated file not existing then skip creation
			if ( (! $rerun) || ( $rerun && ! -e $annotated ) )
			{ 			
				# EPACTS anno engine
				$command = "$EPACTS anno --in $filtered --out $annotated"; 
				if (&Utilities::executeCommand($command, "Annotating filtered vcf file: $filtered to $annotated", $logger)) {
		    		$logger->error("Error while annotating filtered vcf file: $filtered");
		    		exit(-1);
		    	}
			}
			else
			{
				$logger->info("RERUN: Annotated file exists. Skipping.")				
			}		
		
			# If rerun and filtered file index not existing then skip creation 
		    if ( (! $rerun) || ( $rerun && ! -e $annotated.".tbi" ) )
		    {		
				$command = "$tabix -f $annotated";
				if (&Utilities::executeCommand($command, "Indexing annotated vcf file: $annotated", $logger)) {
		    		$logger->error("Error indexing annotated vcf file: $annotated");
		    		exit(-1);
		    	}
		    }
		    

	# Run kinship matrix calculation for EMMAX
		# Out File
			my $kinship = $outEPACTS."/phenotype.kinf";
		
		# Create kinship file	
		
			# If rerun and kinship file not existing then skip creation
			if ( (! $rerun) || ( $rerun && ! -e $kinship ) )
			{ 			
				# EPACTS calculate kinship file
				$command = "$EPACTS make-kin --vcf $annotated --ped  $ped --min-maf $minMAF_Kinship --min-callrate $minCallRate --out $kinship --run $cores";
				if (&Utilities::executeCommand($command, "Calculating kinship matrix: $kinship", $logger)) {
		    		$logger->error("Error while calculating kinship matrix: $kinship");
		    		exit(-1);
		    	}
			}
			else
			{
				$logger->info("RERUN: Kinship file exists. Skipping.")				
			}		


# Burden TEST
if ( $burdentest )
{	
	# Grouping for burden test:

		# Groups:
		
			# Loss of Function
				my $group_lof      = $outEPACTS_burden."/lof.grp";
				
				# EPACTS make-group
					$command = "$EPACTS make-group --vcf $annotated --out $group_lof --format epacts --type Frameshift --type Essential_Splice_Site --type Start_Loss --type Stop_Loss --type Stop_Gain";
					
					if ( (! $rerun) || ( $rerun && ! -e $group_lof ) )
		    		{
						if (&Utilities::executeCommand($command, "Grouping LoF variants", $logger)) {
		    				$logger->error("Error grouping LoF variants");
	    					exit(-1);
			    		}
		    		}			
				
			# Missense				
				my $group_missense = $outEPACTS_burden."/missense.grp";
				
				# EPACTS make-group
					if ( (! $rerun) || ( $rerun && ! -e $group_missense ) )
		    		{
						$command = "$EPACTS make-group --vcf $annotated --out $group_missense --format epacts --type Nonsynonymous";
						if (&Utilities::executeCommand($command, "Grouping missense variants", $logger)) {
		    				$logger->error("Error grouping missense variants");
	    					exit(-1);
						}
			    	}	
				
				
			# Non Synonymous				
				my $group_nonsyn   = $outEPACTS_burden."/nonsyn.grp";
					# EPACTS make-group
					$command = "$EPACTS make-group --vcf $annotated --out $group_nonsyn --format epacts --type Frameshift --type Essential_Splice_Site --type Start_Loss --type Stop_Loss --type Stop_Gain --type Nonsynonymous";
					if ( (! $rerun) || ( $rerun && ! -e $group_nonsyn ) )
					{
						if (&Utilities::executeCommand($command, "Grouping non synonymous variants", $logger)) {
		    				$logger->error("Error grouping non synonymous variants");
	    					exit(-1);
			    		}
					}	

	# Launch tests:
	
		print "Launch EPACTS\n";
		
		# Build covariate strings
		# Transform covariates comma separated argument into the correct syntax for EPACTS
		# --cov SEX --cov PROJECT etc
		$covariates =~ s/,/ --cov /g;
		$covariates = "--cov $covariates" if $covariates ne "";
	
		print "Cov command: ";
		print $covariates."\n";
	
		# Loss of Function
			my $lof_out_base = $outEPACTS_burden."/lof";
		
			# Command. Create makefile and execute it
			$command = "$EPACTS group --vcf $annotated --ped $ped --kin $kinship --test emmaxCMC --pheno $trait --groupf $group_lof --out $lof_out_base  --max-maf $maxMAF_LOF --min-callrate $minCallRate $covariates; \
					make -f $lof_out_base.Makefile -j $cores ";
			# Launch		
			if ( (! $rerun) || ( $rerun && ! -e $lof_out_base.".epacts" ) )
			{		
				if (&Utilities::executeCommand($command, "Burden testing LoF", $logger)) {
		    					$logger->error("Error burden testing LoF");
		    					exit(-1);
				}
			}
			# Create Case Control Matrix
			# get_case_ctl_vcf.pl

		# Loss of Function RARE
			my $rarelof_out_base = $outEPACTS_burden."/rarelof";
		
			# Command. Create makefile and execute it
			$command = "$EPACTS group --vcf $annotated --ped $ped --kin $kinship --test emmaxCMC --pheno $trait --groupf $group_lof --out $rarelof_out_base  --max-maf $maxMAF_LOFrare --min-callrate $minCallRate $covariates; \
					make -f $rarelof_out_base.Makefile -j $cores ";
			# Launch
			if ( (! $rerun) || ( $rerun && ! -e $rarelof_out_base.".epacts" ) )
			{	
				if (&Utilities::executeCommand($command, "Burden testing rare LoF", $logger)) {
		    					$logger->error("Error burden testing rare LoF");
		    					exit(-1);
				}
			}
			
			# Create Case Control Matrix
			# get_case_ctl_vcf.pl
				
		
		# Missense
			my $missense_out_base = $outEPACTS_burden."/missense";
		
			# Command. Create makefile and execute it
			$command = "$EPACTS group --vcf $annotated --ped $ped --kin $kinship --test emmaxCMC --pheno $trait --groupf $group_missense --out $missense_out_base  --min-maf $minMAF_missense --max-maf $maxMAF_missense --min-callrate $minCallRate $covariates; \
					make -f $missense_out_base.Makefile -j $cores ";
			# Launch
			if ( (! $rerun) || ( $rerun && ! -e $missense_out_base.".epacts" ) )
			{						
				if (&Utilities::executeCommand($command, "Burden testing missense", $logger)) {
		    					$logger->error("Error burden testing missense");
		    					exit(-1);
				}
			}
		
		# Missense RARE
			my $raremissense_out_base = $outEPACTS_burden."/raremissense";
		
			# Command. Create makefile and execute it
			$command = "$EPACTS group --vcf $annotated --ped $ped --kin $kinship --test emmaxCMC --pheno $trait --groupf $group_missense --out $raremissense_out_base  --max-maf $maxMAF_missenserare --min-callrate $minCallRate $covariates; \
					make -f $raremissense_out_base.Makefile -j $cores ";
			# Launch
			if ( (! $rerun) || ( $rerun && ! -e $raremissense_out_base.".epacts" ) )
			{		
				if (&Utilities::executeCommand($command, "Burden testing rare missense", $logger)) {
		    					$logger->error("Error burden testing rare missense");
		    					exit(-1);
				}
			}
		
		# Non Synonymous
			my $nonsyn_out_base = $outEPACTS_burden."/nonsyn";
		
			# Command. Create makefile and execute it
			$command = "$EPACTS group --vcf $annotated --ped $ped --kin $kinship --test emmaxCMC --pheno $trait --groupf $group_nonsyn --out $nonsyn_out_base  --max-maf $maxMAF_nonsyn --min-callrate $minCallRate $covariates; \
					make -f $nonsyn_out_base.Makefile -j $cores ";
			# Launch
			if ( (! $rerun) || ( $rerun && ! -e $nonsyn_out_base.".epacts" ) )
			{			
				if (&Utilities::executeCommand($command, "Burden testing non synonymous", $logger)) {
		    					$logger->error("Error burden testing non synonymous");
		    					exit(-1);
				}
			}
		
	

	# Add Case Control counts (lof missense nonsyn)
	#	then apply case control counts to EPACTS files and produce CSVs
	
		print "Annotating files with case control counts\n";
		
		# Launch LoF
		if ( ! -e $lof_out_base.".casectl.csv" )
		{
			$command = "zcat $annotated | perl $extra_prog_path/EPACTS_get_case_ctl_vcf.pl $ped lof $maxMAF_LOF > $lof_out_base.epacts.counts 2> $lof_out_base.epacts.varlist";
			
					if ( ! -e "$lof_out_base.epacts.counts" )
					{
						if (&Utilities::executeCommand($command, "LoF annotation", $logger)) {
	    					$logger->error("Error LoF annotation");
	    					exit(-1);
						}
					}
			#Apply
			$command = "$extra_prog_path/EPACTS_add_case_ctl_vcf.sh $lof_out_base > $lof_out_base.casectl.csv";
			system($command);
		}

		# Launch LoF Rare
		if ( ! -e $rarelof_out_base.".casectl.csv")
		{
			$command = "zcat $annotated | perl $extra_prog_path/EPACTS_get_case_ctl_vcf.pl $ped lof $maxMAF_LOFrare > $rarelof_out_base.epacts.counts 2> $rarelof_out_base.epacts.varlist";
					if ( ! -e "$rarelof_out_base.epacts.counts" )
					{
						if (&Utilities::executeCommand($command, "Rare LoF annotation", $logger)) {
	    					$logger->error("Error Rare LoF annotation");
	    					exit(-1);
						}
					}
			#Apply
			$command = "$extra_prog_path/EPACTS_add_case_ctl_vcf.sh $rarelof_out_base > $rarelof_out_base.casectl.csv";
			system($command); 
		}
		
		# TODO FINISH count outfile protection for reruns
		
		# Launch Missense
		if ( ! -e $missense_out_base.".casectl.csv")
		{
			$command = "zcat $annotated | perl $extra_prog_path/EPACTS_get_case_ctl_vcf.pl $ped missense $maxMAF_missense $minMAF_missense > $missense_out_base.epacts.counts 2> $missense_out_base.epacts.varlist";
			if (&Utilities::executeCommand($command, "Missense annotation", $logger)) {
	    			$logger->error("Error Missense annotation");
	    			exit(-1);
			}
			#Apply
			$command = "$extra_prog_path/EPACTS_add_case_ctl_vcf.sh $missense_out_base > $missense_out_base.casectl.csv";
			system($command);
		} 

		# Launch Missense Rare
		if ( ! -e $raremissense_out_base.".casectl.csv")
		{
			$command = "zcat $annotated | perl $extra_prog_path/EPACTS_get_case_ctl_vcf.pl $ped missense $maxMAF_missenserare > $raremissense_out_base.epacts.counts 2> $raremissense_out_base.epacts.varlist";
			if (&Utilities::executeCommand($command, "Rare Missense annotation", $logger)) {
	    			$logger->error("Error Rare Missense annotation");
	    			exit(-1);
			}
			#Apply
			$command = "$extra_prog_path/EPACTS_add_case_ctl_vcf.sh $raremissense_out_base > $raremissense_out_base.casectl.csv";
			system($command); 
		}


		# Launch Nonsyn
		if ( ! -e $nonsyn_out_base.".casectl.csv")
		{
			$command = "zcat $annotated | perl $extra_prog_path/EPACTS_get_case_ctl_vcf.pl $ped nonsyn $maxMAF_nonsyn > $nonsyn_out_base.epacts.counts 2> $nonsyn_out_base.epacts.varlist";
			if (&Utilities::executeCommand($command, "Nonsyn annotation", $logger)) {
	    			$logger->error("Error Nonsyn annotation");
	    			exit(-1);
			}
			#Apply
			$command = "$extra_prog_path/EPACTS_add_case_ctl_vcf.sh $nonsyn_out_base > $nonsyn_out_base.casectl.csv";
			system($command); 
		}
	
	
	
	# Produce REPORT 	
	#for file in *counts; do echo $file | replace ".counts" ""; cat $file | awk '{VARCA=VARCA+$2;VARCO=VARCO+$3;UNI=UNI+$7}END{print "Variants in cases: "VARCA"\nVariants in controls:"VARCO"\nUnique variants: "UNI}';done

	my $reportfile = $outEPACTS_burden."/reportfile.txt";
	$command = "for type in lof missense nonsyn rarelof raremissense; do echo \${type} >> $reportfile; cat $outEPACTS_burden/\${type}.epacts.counts |  awk \'{VARCA=VARCA+\$2;VARCO=VARCO+\$3;UNI=UNI+\$7}END{print \"Variants in cases: \"VARCA\"\\nVariants in controls:\"VARCO\"\\nUnique variants: \"UNI}\'>> $reportfile; done ";
	
	system($command);


	# Produce QQ-Plots
	$command = "$rscript $prog_path/EPACTS_qqPlots_GeneLabels.R $outEPACTS_burden 2.0";
		system($command);
	
	
	print "\nEPACTS run finished! Check your results at $outEPACTS_burden\n";
	
}

# SINGLE VARIANT TEST
if ( $singlevarianttest )
{
	print "Sorry: single variant test not yet implemented\n";
	print "ZK!!!!! :D \n";
	exit -1;
}


# Done!



=head1 NAME

runEPACTS.pl

=head1 SYNOPSIS

 runEPACTS.pl -vcf <multisample.vcf> -ped <samples.ped> -trait <phenotype_name> -outdir </path/to/out/dir/> -burden [-cov <covariate_name1,covariate_name2,...,covariate_nameN>] [ -allvariants ] [ -mincallrate <min_call_rate> ] [ -minDP <min_DP>] [ -minGQ <min_GQ> ] [ -minMQ <min_MQ> ] [ --rerun ] -lf <logfile.log> -ll <loglevel> 

=head1 DESCRIPTION

This script runs EPACTS for burden test (and in future single variant association test TODO), according to specified filter parameters, and given multisample vcf and ped file

=head1 OPTIONS

Mandatory options

-vcf/-msvcf		multisample VCF 
	-ped			pedigree file
	-trait			name of the disease to run the association test. default: [PHENO]
	-outdir			directory to create the EPACTS data structure
	-burden			execute burden test (mandatory for now)

Filters at the variant level:
	-mincallrate	minimum call rate to call variants, range  [0.00 - 1.00]
	-minDP			minimum depth to call variants					
	-minGQ 			minimum genotype quality to call variants 
	-minMQ			minimum mapping quality to call variants

Optional:
	-cov <cov1,cov2,cov3,...,covN>	covariate names as defined in ped file - comma separated, without spaces (or enclosed into "" )
	-singlevariant					execute single variant test (disabled for now)
	-allvariants					don't take only PASS but all variants
	-lf								log file default: pipeline.log
	-ll								log level: ERROR,INFO,DEBUG; default: INFO

Help:
	-h	this help
	-man show as man page
	
=head1 AUTHOR

Riccardo Berutti

=cut

exit(1);