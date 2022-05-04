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
use File::Copy;

my $prog_path = dirname( abs_path($0) );
require $prog_path . "/Utilities.pm";
#require $prog_path . "/BEDRecord.pm";

# Options
my $input_bam  = "";
my $output_bam = "";
my $ref        = "";
my $maxmemory  = "4g";
my $settings   = "";
my $usebamUtil = 0;
my $fast 	   = 0;
my $downsample = -1;
my $nt         = -1;
my $usegatk3   = 0;
my $usegatk4   = 1;
my $emitIndelQuals = 0;
my $indelcontextsize    = -1;
my $mismatchcontextsize = -1;
my $interval   = "";
my $logfile="SCREEN";
my $loglevel="INFO";
my $help = 0; 
my $man =  0;

# Config
my $default_downsample = 0.10; # 10% downsampling for the fast mode

# Get'em
GetOptions(
	"i=s"	  => \$input_bam,
	"o=s"	  => \$output_bam,
	"ref=s"	  => \$ref,
	"m=s"	  => \$maxmemory,
	"se=s"	  => \$settings,
	"umich"   => \$usebamUtil,
	"bamutil" => \$usebamUtil,
	"fast" => \$fast,
	"downsample=s" => \$downsample,
	"emitIndelQuals" => \$emitIndelQuals,
	"indels_context_size=s"		=> \$indelcontextsize,
	"mismatches_context_size=s" => \$mismatchcontextsize,
	"interval=s" => \$interval,
	"nt=s"	  => \$nt,
	"gatk3"=> \$usegatk3,
	"gatk4"=> \$usegatk4,
	"lf=s" => \$logfile,
	"ll=s" => \$loglevel,
	"h"    => \$help,
	"man"  => \$man
);
	

pod2usage( {-exitval => 0  ,-verbose => 1} ) if $help;
pod2usage( {-exitval => 0  ,-verbose => 2} ) if $man;

# Conditions to crash
pod2usage( {-exitval => 1  ,-verbose => 1} ) if ( 
		( $input_bam  eq ""                       ) ||
		( $settings eq "" && $ref eq ""           ) ||
		( ! -f $input_bam                         )
		);

# Determining number of threads (adapt to the number of slots given by SGE if applicable)
if($nt == -1){
	$nt = 1;
	if($ENV{NSLOTS}){		#get number of slots given by SGE
		$nt = $ENV{NSLOTS};		
	}
}

# Get parameters
my $params     = Utilities::getParams();

# Start logger
Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();

# Get tool paths
$ref  = ( $ref ne "" ? $ref : ( $params->{settings}->{$settings}->{reference} ? $params->{settings}->{$settings}->{reference} : "" ) );
my $gatk = $params->{programs}->{gatk}->{path};
my $java = $params->{programs}->{java}->{path};
my $gatk4 = $params->{programs}->{gatk4}->{path};
my $tabix= $params->{programs}->{tabix}->{path};
my $bamutil=$params->{programs}->{bamutil}->{path};
my $samtools      = $params->{programs}->{samtools}->{path};

# GATK 3 version
$usegatk4 = 0 if ($usegatk3);

# Set our version of Java in Path with priority so to overcome any current installation
$ENV{'PATH'}=dirname(abs_path($java)).":".$ENV{'PATH'};

# Reference?
if ( $ref eq "" )
{
	$logger->error("No default reference for settings $settings and no reference specified. Use -ref");
	exit(-1);
}

# Use only chr[1-22],X,Y,M to recalibrate
my $normalchromosomes="";
	$normalchromosomes = $params->{settings}->{$settings}->{normalchromosomes}  	if ( $settings ne "" &&  $params->{settings}->{$settings}->{normalchromosomes} );

	$normalchromosomes = $interval		if ( $interval ne "" );
	$normalchromosomes = ""				if ( $interval eq "unrestricted");

$logger->debug("Using reference $ref");


# Bam names processing:

# Replace input (saving output)
if ( $input_bam eq $output_bam )
{
	# Move Input to unrecal_$input_bam and index
	my $new_input = dirname( abs_path( $input_bam ) )."/"."unrecal".".".basename($input_bam);
		move ( $input_bam, $new_input );
		
		
		# If m.bam index is m.bam.bai
		move ( $input_bam.".bai", $new_input.".bai" ) if -e $input_bam.".bai";
		
		# If m.bam index is m.bai
		my $input_bai = $input_bam;
			$input_bai =~ s/\.bam$/\.bai/g;
		my $new_input_bai = $new_input;
			$new_input_bai =~ s/\.bam$/\.bai/g;
			#move ( $input_bai, $new_input_bai ) if -e $input_bai; #bad name, bad mess
			system ("rm $input_bai -rf");
			system ("rm $new_input_bai -rf");
			system ("$samtools index $new_input") if ! -e $input_bam.".bai";
		$input_bam=$new_input;
}

# Default output name
if ( $output_bam eq "" )
{
	# Define output file name
	$output_bam = dirname( abs_path( $input_bam ) )."/"."recal".".".basename($input_bam);
}

# Recal tmp:
my $recal_tmp = $output_bam.".recal.grp";

# Downsample settings
$downsample = $default_downsample if (( $downsample == -1 ) && $fast );
my $downsample_opt = (	(($downsample != -1) || $fast) 
						?
						( $usebamUtil 
								?
									" --fast "
								:
									" --downsample_to_fraction $downsample " 
					 	)
					 	:
					 	""
					 );


# The Context options - disabled:
# --indels_context_size			size of the region to be considered around variants for indels (defaults 3)
# --mismatches_context_size		size of the k-mer to be considered around SNPs (defaults 2) - nicer to set 4 to avoid normal error trains
# Increased CTX SIZE (go default)
#--mismatches_context_size 4 \\
#		--indels_context_size 6 \\

my $indelcontext_opt = "";
	$indelcontext_opt = " --indels-context-size $indelcontextsize "  		if ( $indelcontextsize != -1 );
my $snpcontext_opt = "";
	$snpcontext_opt   = " --mismatches-context-size $mismatchcontextsize "	if ( $mismatchcontextsize != -1 );



# Gatk 3.x Syntax
#	
#	java -jar GenomeAnalysisTK.jar \
#	   -T PrintReads \
#	   -R reference.fasta \
#	   -I input.bam \
#	   -BQSR recalibration_report.grp \
#	   -o output.bam
#

# Other opts
# --disable_indel_quals		It creates two extra sets of qualities for every read, increases the bam size too much!

# Build extra options string:
my $extra_options_BaseRec = "";
	$extra_options_BaseRec .= $downsample_opt . $indelcontext_opt . $snpcontext_opt ; 
	$extra_options_BaseRec .= " --intervals $normalchromosomes " if ( $normalchromosomes ne "" );

my $extra_options_PrintReads = "";
	$extra_options_PrintReads .= " --disable_indel_quals " if ( ! $emitIndelQuals  && ! $usegatk4);
	$extra_options_PrintReads .=  " --intervals $normalchromosomes " if ( $normalchromosomes ne "" );
	

###################################################
# GATK command
################################################### 
my $gatk_command = " \\
	$java -XX:ParallelGCThreads=1 -Xmx".$maxmemory." \\
	-jar $gatk \\
		-T BaseRecalibrator \\
		-R $ref \\
		-I $input_bam \\
		-o $recal_tmp \\
		-knownSites $params->{settings}->{$settings}->{HapMap} \\
		-knownSites $params->{settings}->{$settings}->{Omni1000G} \\
		-knownSites $params->{settings}->{$settings}->{SNPs1000G} \\
		-knownSites $params->{settings}->{$settings}->{gatksnps} \\
		-cov ContextCovariate -cov CycleCovariate $extra_options_BaseRec ; \\
	$java -XX:ParallelGCThreads=1 -Xmx".$maxmemory." \\
	-jar $gatk \\
		-nt $nt \\
		-T PrintReads \\
		-R $ref \\
		--validation_strictness LENIENT \\
		-I $input_bam \\
		-BQSR $recal_tmp \\
		-o $output_bam $extra_options_PrintReads \\
";

###################################################
# GATK4 command
###################################################
# --nt $nt \\ only with Spark
my $gatk4_command = " \\
	$gatk4 --java-options  \"-Xmx$maxmemory\" \\
		BaseRecalibrator \\
		--reference $ref \\
		--input $input_bam \\
		--output $recal_tmp \\
		--known-sites $params->{settings}->{$settings}->{HapMap} \\
		--known-sites $params->{settings}->{$settings}->{Omni1000G} \\
		--known-sites $params->{settings}->{$settings}->{SNPs1000G} \\
		--known-sites $params->{settings}->{$settings}->{gatksnps} \\
		$extra_options_BaseRec ; \\
	$gatk4 --java-options  \"-Dsamjdk.compression_level=5 -Xmx$maxmemory\" \\
		ApplyBQSR \\
		--reference $ref \\
		--read-validation-stringency LENIENT \\
		--input $input_bam \\
		--bqsr-recal-file $recal_tmp \\
		--output $output_bam \\
";
 
 #TO ADD OR REMOVE OPT
 #		--static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \\

###################################################
# Umich Abecasis-Lab bamUtils:
# bamUtil command:
###################################################
my $bamUtil_command =" \\
	$bamutil recab 	--in $input_bam --out $output_bam \\
					--log ".$output_bam.".log --noPhoneHome $downsample_opt \\
					--refFile $ref --dbsnp $params->{settings}->{$settings}->{gatksnps} \\
		";

###################################################
# Build command
###################################################
my $command = ( $usebamUtil ? $bamUtil_command : ( $usegatk4 ? $gatk4_command :  $gatk_command ) );


$logger->info("Launching ".( $usebamUtil ? "bamUtil recab" : "GATK PrintReads")." (to recalibrate BQSR) with command: $command");

# If batch already exists do not create it, but issue WARNING
if (&Utilities::executeCommand($command, "Launching ".( $usebamUtil ? "bamUtil recab" : "GATK PrintReads")." (to recalibrate BQSR) on $input_bam to $output_bam", $logger)) {
	$logger->error("Error executing ".( $usebamUtil ? "bamUtil recab" : "GATK PrintReads")." (to recalibrate BQSR) on $input_bam to $output_bam");
	exit(100);
}
else
{
	$logger->info("Recalibration done");
	# TODO:
	# UNLINK ORIGINAL BAM
}

#TODO: quant qual
# quantizeQuals



=head1 NAME

recalBam.pl

=head1 SYNOPSIS

 recalBam.pl -i input.bam [ -o output.bam ] -se settings -r reference.fa [ -m 4g ] 
 
=head1 DESCRIPTION

This script uses GATK (default)/ bamUtil to recalibrate base qualities in a bam file

=head1 OPTIONS

 -i    	<infile> path to the input bam to recalibrate
 -o    	<outfile> path to the output recalibrated bam. If empty <path/to/>recal.<infile.bam>, if equals to <infile.bam> then infile will be moved to </path/to/>unrecal.<infile.bam> 
 -m	   	<memory> max. memory for each Java GATK job; default 4g
 -se	<settings>; required
 -umich   	use bamUtil from Umich gotCloud pipeline (Abecasis Lab). Default GATK 
 -bamutil 	use bamUtil from Umich gotCloud pipeline (Abecasis Lab). Default GATK [Syn to umich]
 -gatk3		use GATK3 [default GATK 4.x]
 -gatk4		use GATK4 [default] for compatibility
 -emitIndelQuals	emit indel qualities into the resulting BAM file. WARNING: It grows >2x in size. Don't do it unless you need them
 -quantizeQuals				quantize qualities as suggested in GATK Best (bad?) Practices 10,20,30 [NOT ENABLED]
 -indels_context_size		size in bp of the region evaluated around indels (defaults to GATK default. Actually 3)
 -mismatches_context_size 	size in bp of the region evaluated around ref mismatches (defaults to GATK default. Actually 2) 
 -interval	specify GATK compatible interval (eg: chr10, chr5:1295102-1951020, intervals.bed), or "unrestricted" for no interval. Defaults to normal chromosomes if defined in the configuration.
 -nt    number of threads
 -lf	<log file>; default: pipeline.log
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h		print this helptext
 -man	show man page
 
=head1 AUTHOR

Riccardo Berutti 

=cut