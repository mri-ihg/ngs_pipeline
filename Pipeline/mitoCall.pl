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
my $output = "/dev/stdout";
my $sample_assigned_name = "";
my $ref        = "";
my $settings   = "mtDNA";
my $interval   = "chrM";
my $calcNoise  = 0;
my $noisePlot  = 0;
my $callGT     = 0;
my $emitDP     = 0; 
my $isSingleReadExperiment = 0;
my $outall     = 0;
my $rawcalc    = 0; # No noise or gaussian filter
my $ignoreIndels = 0; # Do not consider indels
my $dontsortgenos = 0; # Do not sort output genotypes
my $outVCF = 0; # VCF Output

my $logfile="SCREEN";
my $loglevel="INFO";
my $help = 0; 
my $man =  0;


######################################################################################################################################################################################
# TODO:
# This caller works pretty well also for non mito stuff. A global switch for mito / genomic / exomic settings with default settings can be a good idea
######################################################################################################################################################################################

# Default SETTINGS

	# Minimum heteroplasmy accepted
	my $HMIN = 0.0199;

	# Forward reverse imbalance MAX
	my $FWRVIMB = 0.3;

	# How many alleles to call MAX
	my $MAX_ALLELE = 2;
  		$MAX_ALLELE=4;

	# Minimum supporting reads 
	my $MIN_SUPPORTING = 0;	# TODO: set a default for mitoCalling like 10?
		
	# Minimum depth 
	my $MIN_DEPTH = 0;
	
	# Maximum depth
	my $MAX_DEPTH = 20000;	#High enough for mt
	
	# Minimum Mapping quality
	my $MIN_MAPQ = 0;
	
	# Minimum Base Quality (13 samtools mpileup default)
	my $MIN_BQ = 13;
	
	
# Get'em
GetOptions(
	"i=s"	  => \$input_bam,
	"o=s"     => \$output,
	"s=s"     => \$sample_assigned_name,
	"ref=s"	  => \$ref,
	"se=s"	  => \$settings,
	"minhet=s"=> \$HMIN,
	"fwrvimbmax=s" => \$FWRVIMB,
	"max_alleles=s" => \$MAX_ALLELE,
	"minsupporting=s" => \$MIN_SUPPORTING,
	"mindepth=s"	=> \$MIN_DEPTH,
	"maxdepth=s"	=> \$MAX_DEPTH,
	"minMQ=s"=> \$MIN_MAPQ,
	"minBQ=s"=> \$MIN_BQ,
	"noise"=> \$calcNoise,
	"noiseplot" => \$noisePlot,
	"call" => \$callGT,
	"emitDP"=>\$emitDP,
	"single" => \$isSingleReadExperiment,
	"a"    => \$outall,
	"raw"  => \$rawcalc,
	"noindels" => \$ignoreIndels,
	"nogenosort" => \$dontsortgenos,
	"vcf" => \$outVCF,
	"l=s"  => \$interval,
	"lf=s" => \$logfile,
	"ll=s" => \$loglevel,
	"h"    => \$help,
	"man"  => \$man
);

# Post Config

	# Threshold allelecount
	my $ALLELECOUNT_THRESHOLD=4-$MAX_ALLELE;
	
	$calcNoise = 1 if $noisePlot;

pod2usage( {-exitval => 0  ,-verbose => 1} ) if $help;
pod2usage( {-exitval => 0  ,-verbose => 2} ) if $man;

# Conditions to crash
pod2usage( {-exitval => 1  ,-verbose => 1} ) if ( 
		( $input_bam  eq ""  && $sample_assigned_name eq ""  ) ||
		( $settings eq "" && $ref eq ""           )
		);

# Get parameters
my $params     = Utilities::getParams();

# Start logger
Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();

# Get tool paths
$ref  = ( $ref ne "" ? $ref : ( $params->{settings}->{$settings}->{reference} ? $params->{settings}->{$settings}->{reference} : "" ) );
my $samtools = $params->{programs}->{samtools}->{path};

# Analysis folder
my $analysis_folder = $params->{settings}->{$settings}->{analysis}->{folder};


# Reference?
if ( $ref eq "" )
{
	$logger->error("No default reference for settings $settings and no reference specified. Use -ref");
	exit(-1);
}

# Sample name specified? 
if ( $sample_assigned_name ne "" )
{
	 $logger->info("Looking for file $analysis_folder/S0*/$sample_assigned_name/*out/paired-endout/merged.rmdup.bam");
	 my @files=glob($analysis_folder."/S0*/$sample_assigned_name/*out/paired-endout/merged.rmdup.bam");
	 
	 	if ( scalar(@files) != 1 )
		{
			$logger->error("Sample name $sample_assigned_name not found or ambiguous, please check\n");
			exit(-1);
		}
	 
	 $input_bam = $files[0];
}

# Input file?
if ( ! -f $input_bam )
{
	$logger->error("Input file $input_bam does not exist!");
	exit(-1);
}

# Output file:
open (my $fout, ">", $output ) or die $logger->error("Cannot write to output file $output"); 

# IsSingleReadExperiment ? 
# -A enables anomalous paired (and unpaired) reads to be accounted for for the mpileup
my $singleReadExperiment = ($isSingleReadExperiment ? "-A" : "" );


# If not calcnoise, but calling, calculate noise threshold first
my $sample_noise = 0;
my $sample_name = "";

if ( (! $calcNoise) && ($interval eq "chrM") && (! $rawcalc) )
{
	#Debug
	$logger->info("Sample noise calculation: ");
	my $sample_noise_command = "perl $prog_path/mitoCall.pl -i $input_bam -ref $ref -se $settings -noise -l $interval -lf $logfile -ll $loglevel"; 
	   $sample_noise_command = $sample_noise_command . ( $isSingleReadExperiment ? " -single " : "" );
	$logger->debug("Command: ".$sample_noise_command);
	$sample_noise = `$sample_noise_command`;
	$logger->info("Sample noise: ".$sample_noise);
	
	
	#Get Sample Name 
	my $sample_name_command="$samtools view $input_bam -H | grep \"\@RG\" | awk -F\"SM:\" \'{print \$2}\' | awk -F\"\t\" \'{print \$1}\' | tail -n 1 ";
	$sample_name = `$sample_name_command`;
	chomp($sample_name); 
}

$sample_name = "undefined" if ($sample_name eq "" );


# samtools mpileup AEV7F_84737M_FULL.FULL.bam -r chrM -f /data/mirror/goldenpath/hg19_decoy/chromosome/noPAR.hg19_decoy.fa | head
#
#chrM	1	G	67	^],^],^],^],^].^].^].^].^].^].^].^].^].^].^].^9.^],^],^],^],^],^],^],^],^],^],^$,^],^],^],^],^],^],^;,^],^9,^3,^],^B,^],^],^],^!,^],^I,^],^],^],^],^],^],^-,^],^V,^],^],^],^],^],^],^],^>,^9,^],^],^],^!,EE<EECBCCA1CA/CCEEEEEEEEEEEEEEEEEEEEEEEEEEAEEEEEEEE2EEEEEEEECEEEEEE
#chrM	2	A	109	,,c,,,............G...,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,^],^],^],^],^],^],^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].	HH0/H.GBA@C/CA/>GB//CFGHHHHHHHHHGHHHHHHHHHHHHGHGBHHHHHHHHFHHHHHHFAHHHHHGHEEE=EEA?ACCACBC3ACBBBBBCAACA?ABACA>A
#chrM	3	T	195	,,,,,,,.................,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,...............................^],^],^],^],^],^],^],^],^],^],^],^],^],^],^],^],^],^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^E.^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^].^],^E,	H5H05H4HBCCC5BA55>G@?5DFHHHHHHHHHHFHHHHHHHHCHHHGHGFHFHDHHHH2HHFHHHH3HHHHHHHHHGFHHAABEBBCBCAAC3CBBBBCBACBABBBCAAAEEEEEEEEEEEEEEEEECACCCCCCCCBCC?BBBCAEADCCBCEBCBEACBCAADBCBBDDBBAAAACBBBCCA3BCBBBBEE


# In the pileup format (without -u or -g), each line represents a genomic position, consisting of chromosome name, 1-based coordinate, reference base, the number of reads covering the site, 
# read bases, base qualities and alignment mapping qualities. Information on match, mismatch, indel, strand, mapping quality and start and end of a read are all encoded at the read base column. 

# SNPs
# a dot stands for a match to the reference base on the forward strand, 
# a comma for a match on the reverse strand, 
# a '>' or '<' for a reference skip
# `ACGTN' for a mismatch on the forward strand and 
# `acgtn' for a mismatch on the reverse strand.

# INDELs
# A pattern `\\+[0-9]+[ACGTNacgtn]+' indicates there is an insertion between this reference position and the next reference position. The length of the insertion is given by the integer in the pattern, followed by the inserted sequence. 
# Similarly, a pattern `-[0-9]+[ACGTNacgtn]+' represents a deletion from the reference. The deleted bases will be presented as `*' in the following lines. 

# START/END
# Also at the read base column, a symbol `^' marks the start of a read. The ASCII of the character following `^' minus 33 gives the mapping quality. 
# A symbol `$' marks the end of a read segment.

# What are ">" 


# Call Genotype OUT?
my $GTOUT = $callGT ? " GT " : " ";
my $DPOUT = $emitDP ? "DP " : "";

# Print header 
if ( ! $calcNoise ) 
{
	print $fout "#SAMPLE: $sample_name\n";
	print $fout "#MinHet: $HMIN\n";  
	print $fout "#NOISE_LEVEL: $sample_noise\n"		if (! $rawcalc);
	print $fout "#SEQ_NOISE_FILTER: on\n"				if (! $rawcalc);
	print $fout "#GAUSSIAN_NOISE_FILTER: on\n"		if (! $rawcalc);
	print $fout "CHR POS REF".$GTOUT.$DPOUT."A C G T INS DEL FWIMB\n";
	
	$logger->info("Starting genotyping sample $sample_name");
}

my $noise=0; my $noisecount=0; my $noisetot=0; my $noisecounttot=0;


# Interval:

if ( $interval =~ /.*\.bed$/ )
{
	$interval = " -l $interval ";
}
else
{
	$interval = " -r $interval "; 
}

# Opening bam file
# Samtools command
# $singleReadExperiment == -A if single read or it's empty
my $command = "$samtools mpileup $input_bam $interval -f $ref -q $MIN_MAPQ -Q $MIN_BQ  ".$singleReadExperiment." 2>/dev/null ";
open DATA, $command." | ";

# Reading data
while(<DATA>)
{
	# Getting output:
	my @line=split("\t", $_);
	my $isref=0;
	my $currchr=$line[0];
	my $pos=$line[1];
	my $refallele=uc($line[2]);
	my $lcref=lc($refallele);
	#	$read=uc($line[4]);
	my $read=($line[4]);
	
	
		# Extraction of indels first:
		#
		my %insertions; my %deletions;
		
		my $INS=0;	my $INSSEQ=""; 
		my $DEL=0;	my $DELSEQ="";
		
		if ( ! $ignoreIndels )
		{
			# This matches the indel length ($1) 
			my @indelmatches = ( $read =~  m/[+-]([0-9]+)(??{ "[ACGTacgt]{$1}" })/g );			
			
			my $maxCoveredIns = 0;
			my $maxCoveredDel = 0;
			
			# Reparse individually per indel length
			foreach my $match ( uniq(@indelmatches) )
			{
				# Each indel length might contain more indels so split
				foreach my $indel ( $read =~ m/([+-]$match[ACGTacgt]{$match})/g )
        		{
                	#print $indel."\n";
                	my ($sign, $length, $seq)=( $indel =~ m/([+-])([0-9]+)([ACGTacgt]+)/ );
                	
                	
                	#TODO: add FWD/REV balance here!
                	$seq=uc($seq);
                	
                	if ($sign eq "+" )
                	{
                		$insertions{$seq} = ( exists $insertions{$seq} ? $insertions{$seq}+1 : 1  ); 
                		$maxCoveredIns = $insertions{$seq} if $insertions{$seq} > $maxCoveredIns;
                	}
                	elsif ( $sign eq "-" )
                	{
                		$deletions{$seq}  = ( exists $deletions{$seq}  ? $deletions{$seq}+1  : 1  );
                		$maxCoveredDel = $deletions{$seq} if $deletions{$seq} > $maxCoveredDel;
                	}
                	else
                	{
                		$logger->error("Unparsable data from mpileup.");
                		exit(-1);
                	}
        		}
			}
			
			# Clean and purge elements:
			#INS
			foreach my $seq ( keys %insertions )
			{
				delete $insertions{$seq} if ( ($insertions{$seq} < $maxCoveredIns ) || ($insertions{$seq} < $MIN_SUPPORTING ) );
				#Debug
				#print $seq."> seen".$insertions{$seq}." times" if exists( $insertions{$seq} );
				
				#If I didn't purge it save it. I pick the last one if many are represented the same amount of times
				if (  exists( $insertions{$seq} ) ) 
				{
					$INSSEQ=$seq;
					$INS=$insertions{$seq};
				}
			}
			#DEL
			foreach my $seq ( keys %deletions )
			{
				delete $deletions{$seq} if ( ($deletions{$seq} < $maxCoveredDel )   || ($deletions{$seq} < $MIN_SUPPORTING ) );
				#Debug
				#print $seq."> seen".$deletions{$seq}." times" if exists( $deletions{$seq} );
				#If I didn't purge it save it. I pick the last one if many are represented the same amount of times
				if (  exists( $deletions{$seq} ) ) 
				{
					$DELSEQ=$seq;
					$DEL=$deletions{$seq};
				}
				
			}			
								
		}
	
		#Old parsing
		#$read =~ s/[\+\-][0-9]*[ACGTacgt]\.*/./g;
		
		#Strip away indels and continue
		# This syntax matches [+-] one number and as many [ACGTacgt] as the number previously matched
		# This MIGHT NOT WORK with new PERL versions since it uses the Perl $1 variable which is already set.
		$read =~ s/[+-]([0-9]+)(??{ "[ACGTacgt]{$1}" })//g;
		
		#Ignore remaining symbols: $,^Qual, * dels
		$read =~ s/[^ACGTacgt.,]//g;
		$read =~ s/[.]/$refallele/g;
		$read =~ s/[,]/$lcref/g;

		my $A_F = ($read =~ tr/A//);	my $A_R = ($read =~ tr/a//);
		my $C_F = ($read =~ tr/C//);	my $C_R = ($read =~ tr/c//);
		my $G_F = ($read =~ tr/G//);	my $G_R = ($read =~ tr/g//);
		my $T_F = ($read =~ tr/T//);	my $T_R = ($read =~ tr/t//);

		$read=uc($read);
		
		my $A = ($read =~ tr/A//);
		my $C = ($read =~ tr/C//);
		my $G = ($read =~ tr/G//);
		my $T = ($read =~ tr/T//);
		#$INS #Defined and can be zero
		#$DEL #Defined and can be zero
		
		my $TOT=$A+$C+$G+$T+$INS+$DEL;
		
		
		if ( $calcNoise && $TOT>0 ) # TOT>0 PREVENT CRASHES
		{
			$noise=0;
			$noisecount=0;
			
			# Rules for noise: <20% het, and 3rd or 4th allele call
			#if (( $A/$TOT < 0.2 )){	$noise += $A/$TOT;	$noisecount++;	};
			#if (( $C/$TOT < 0.2 )){	$noise += $C/$TOT;	$noisecount++;	};
			#if (( $G/$TOT < 0.2 )){	$noise += $G/$TOT;	$noisecount++;	};
			#if (( $T/$TOT < 0.2 )){	$noise += $T/$TOT;	$noisecount++;	};
			
			# For now do not use INDELs for noise
			
			# Use the two less likely alleles, then divide by two and multiply by three
			if (( $A/$TOT < 0.1 )&&( (($A<$C) + ($A<$G) + ($A<$T))>2) ){	$noise += $A/$TOT;	$noisecount++;	};
			if (( $C/$TOT < 0.1 )&&( (($C<$A) + ($C<$G) + ($C<$T))>2) ){	$noise += $C/$TOT;	$noisecount++;	};
			if (( $G/$TOT < 0.1 )&&( (($G<$A) + ($G<$C) + ($G<$T))>2) ){	$noise += $G/$TOT;	$noisecount++;	};
			if (( $T/$TOT < 0.1 )&&( (($T<$A) + ($T<$C) + ($T<$G))>2) ){	$noise += $T/$TOT;	$noisecount++;	};
			$noise=($noise*3)/2;

			# Add gaussian sqrt(N) noise
			#$noisetot += $noise + sqrt($A+$C+$G+$T)/$TOT;
			$noisetot += $noise;
			
			$noisecounttot += $noisecount;
			$noisecounttot += ( $noisecount > 0 ? 1 : 0 );
			
			#print (( $noisecount == 0 ? 0 : $noise/$noisecount )."\n") if $noisePlot; 
			print (( $noisecount == 0 ? 0 : $noise )."\n") if $noisePlot;
		
			#TODO: FOR NOISE CALC CAN I STOP HERE? 	
		}


		# Min supporting reads:
		$A = 0 if ( $A < $MIN_SUPPORTING || $TOT < $MIN_DEPTH || $TOT > $MAX_DEPTH );
		$C = 0 if ( $C < $MIN_SUPPORTING || $TOT < $MIN_DEPTH || $TOT > $MAX_DEPTH );
		$G = 0 if ( $G < $MIN_SUPPORTING || $TOT < $MIN_DEPTH || $TOT > $MAX_DEPTH );
		$T = 0 if ( $T < $MIN_SUPPORTING || $TOT < $MIN_DEPTH || $TOT > $MAX_DEPTH );
		#$INS = 0 if ( $INS < $MIN_SUPPORTING );	#ALREADY FILTERED OPTIM.
		#$DEL = 0 if ( $DEL < $MIN_SUPPORTING );	#ALREADY FILTERED OPTIM.
		$INS = 0 if ( $TOT < $MIN_DEPTH || $TOT > $MAX_DEPTH );	#ALREADY FILTERED OPTIM.
		$DEL = 0 if ( $TOT < $MIN_DEPTH || $TOT > $MAX_DEPTH );	#ALREADY FILTERED OPTIM.

		$TOT=0 if ( $TOT < $MIN_DEPTH );

		# Remove a gaussian figure of noise (if not raw calculation required)
		if ( ! $rawcalc )
		{
			$A=$A-sqrt($A)-($A*$sample_noise);
			$C=$C-sqrt($C)-($C*$sample_noise);
			$G=$G-sqrt($G)-($G*$sample_noise);
			$T=$T-sqrt($T)-($T*$sample_noise);
			# TODO: A indel noise would be required
		}


		# Minhet should be set to 0 even in case of rawcalc 
		my $HA = 0;
			$HA = ( $A/$TOT > $HMIN ? $A/$TOT : 0) if $TOT>0;
		my $HC = 0;
			$HC = ( $C/$TOT > $HMIN ? $C/$TOT : 0) if $TOT>0;
		my $HG = 0;
			$HG = ( $G/$TOT > $HMIN ? $G/$TOT : 0) if $TOT>0;
		my $HT = 0;
			$HT = ( $T/$TOT > $HMIN ? $T/$TOT : 0) if $TOT>0;
		my $HINS = 0;
			$HINS = ( $INS/$TOT > $HMIN ? $INS/$TOT : 0) if $TOT>0;
		my $HDEL = 0;
			$HDEL = ( $DEL/$TOT > $HMIN ? $DEL/$TOT : 0) if $TOT>0;
				

		my $SB = 0;
		
		my $HA_FIN = ( ($HA>$HC)+($HA>$HG)+($HA>$HT)+($HA>$HINS)+($HA>$HDEL) >=$ALLELECOUNT_THRESHOLD ? $HA : 0);			$SB = ( $HA_FIN == 0 ? 0 : ( abs($A_F-$A_R)/($A_F+$A_R) > $FWRVIMB ? 1 : $SB ));
		my $HC_FIN = ( ($HC>$HA)+($HC>$HG)+($HC>$HT)+($HC>$HINS)+($HC>$HDEL) >=$ALLELECOUNT_THRESHOLD ? $HC : 0);			$SB = ( $HC_FIN == 0 ? 0 : ( abs($C_F-$C_R)/($C_F+$C_R) > $FWRVIMB ? 1 : $SB ));
		my $HG_FIN = ( ($HG>$HC)+($HG>$HA)+($HG>$HT)+($HG>$HINS)+($HG>$HDEL) >=$ALLELECOUNT_THRESHOLD ? $HG : 0);			$SB = ( $HG_FIN == 0 ? 0 : ( abs($G_F-$G_R)/($G_F+$G_R) > $FWRVIMB ? 1 : $SB ));
 		my $HT_FIN = ( ($HT>$HC)+($HT>$HG)+($HT>$HA)+($HT>$HINS)+($HT>$HDEL) >=$ALLELECOUNT_THRESHOLD ? $HT : 0);			$SB = ( $HT_FIN == 0 ? 0 : ( abs($T_F-$T_R)/($T_F+$T_R) > $FWRVIMB ? 1 : $SB ));
		my $HINS_FIN=( ($HINS>$HA)+($HINS>$HC)+($HINS>$HG)+($HINS>$HT)+($HINS>$HDEL) >=$ALLELECOUNT_THRESHOLD ? $HINS : 0);	#TODO # For now no strand balance for indels but: TO DO!
		my $HDEL_FIN=( ($HDEL>$HA)+($HDEL>$HC)+($HDEL>$HG)+($HDEL>$HT)+($HDEL>$HINS) >=$ALLELECOUNT_THRESHOLD ? $HDEL : 0);	#TODO # For now no strand balance for indels but: TO DO!

		# If ONLY one is called (sum of possible alleles with supporting reads is one AND the refallele has reads: then it's ref
		if (
			( ( ($HA_FIN > 0) + ($HC_FIN > 0) + ($HG_FIN > 0) + ($HT_FIN > 0) ) == 1 )
			&&
			(
				( $refallele eq "A" ) && ( $HA_FIN > 0 ) ||
				( $refallele eq "C" ) && ( $HC_FIN > 0 ) ||
				( $refallele eq "G" ) && ( $HG_FIN > 0 ) ||
				( $refallele eq "T" ) && ( $HT_FIN > 0 )
			)
		)
		{
			$isref=1;
		}
		else
		{
			$isref=0;
		}


	if ( (( ! $isref || $outall )) && (! $calcNoise ) )
	{
		$SB = ( $SB == 1 ? "IMB" : "" );
		
		
		#$GTOUT=" ";
		
		#$GTOUT=" ".
		#		($HA_FIN>0 ? "A" : "").
		#		($HC_FIN>0 ? "C" : "").
		#		($HG_FIN>0 ? "G" : "").
		#		($HT_FIN>0 ? "T" : "").
		#		" "
		#		if $callGT;
		   
		#$GTOUT = " N " if ($callGT && $GTOUT eq "  ");
		
		
		#Build genotype string:
		#my @genotypes;
		
		my $extraref = ""; 
			$extraref = $DELSEQ if $HDEL_FIN>0;
		
		my $outref=$refallele.$extraref;
		
		
		#TODO # In case of VCF OUT SUPPRESS REF ALLELE FOR GENOTYPE OUT (TODO!)
		
		#my $AOUT=""; 	$AOUT="A".$outref if ( $HA_FIN>0 );
		#my $COUT=""; 	$COUT="C".$outref if ( $HC_FIN>0 );
		#my $GOUT=""; 	$GOUT="G".$outref if ( $HG_FIN>0 );
		#my $TOUT=""; 	$TOUT="T".$outref if ( $HT_FIN>0 );
		#my $INSOUT="";  $INSOUT=$refallele.$INSSEQ.$outref if ($HINS_FIN>0); 
		#my $DELOUT="";  $DELOUT=$refallele if ($HDEL_FIN>0);
		
		#%genotypes{"A".$outref}=$HA_FIN>0; then sort, push out 
		
		my %genotypes;

		#Set genotype hash
		$genotypes{"A".$extraref}=$HA_FIN if $HA_FIN>0;
		$genotypes{"C".$extraref}=$HC_FIN if $HC_FIN>0;
		$genotypes{"G".$extraref}=$HG_FIN if $HG_FIN>0;
		$genotypes{"T".$extraref}=$HT_FIN if $HT_FIN>0;  
		$genotypes{$refallele.$INSSEQ.$extraref}=$HINS_FIN if $HINS_FIN>0; 
		$genotypes{$refallele}=$HDEL_FIN if $HDEL_FIN>0;
		
		my @genotypesout;
		
		if ($dontsortgenos)
		{
				push @genotypesout, "A".$extraref if $HA_FIN>0;
				push @genotypesout, "C".$extraref if $HC_FIN>0;
				push @genotypesout, "G".$extraref if $HG_FIN>0;
				push @genotypesout, "T".$extraref if $HT_FIN>0;
				push @genotypesout, $refallele.$INSSEQ.$extraref if $HINS_FIN>0;
				push @genotypesout, $refallele	  if $HDEL_FIN>0;
		}
		else
		{
			foreach my $allele ( sort {$genotypes{$b} <=> $genotypes{$a}} keys %genotypes)
			{
				push @genotypesout, $allele;
			}
		}
	 	
	 	my $geno = join ",", @genotypesout;
		
		$GTOUT=" ".$geno." ";
		$GTOUT = " N " if ($callGT && $GTOUT eq "  ");
		$GTOUT = "" if (!$callGT);
		
		$DPOUT = ( $emitDP ? ($callGT ? "":" ").$TOT." " : "" );
		
		print $fout "$currchr $pos $outref".$GTOUT.$DPOUT."$HA_FIN $HC_FIN $HG_FIN $HT_FIN $HINS_FIN $HDEL_FIN $SB\n";
	}

}

if ( $calcNoise && ! $noisePlot )
{	
	$noisecounttot=1 if $noisecounttot == 0;
	print $noisetot/$noisecounttot."\n";	
}



=head1 NAME

mitoCall.pl 

=head1 SYNOPSIS

 mitoCall.pl  

=head1 DESCRIPTION

This script calls mtDNA variants. With tuned options can be used to precisely call any region of the genome.
There are only few statistical tricks that can be disabled so you know why you are getting what! 

NOTE: There will be no genotype 

=head1 OPTIONS


 -i <input.bam>     input bam file. Alternative to -s
 -s <samplename>    sample name. Alternative to -i 
 -o <outfile.vcf>/<outfile.data> output VCF/TSV file optional
 -ref <hg19plus.fa> reference
 -se <mtDNA>        pipeline settings (takes the reference in) 
 -minhet <0.0199>   miminum heteroplasmy to call
 -fwrvimbmax <0.3>  forward reverse imbalance max
 -max_alleles <4>   max number of alleles returned
 -minsupporting <0> minimum number of reads to support a call (disabled by default, enable with option)
 -mindepth <0>      minimum total depth to support a call (disabled by default)
 -maxdepth <20000>  maximum totla depth to accept a call (with 20000 it's almost disabled)
 -noise             calculate noise and nothing else
 -noiseplot         output per position noise calculation
 -a                 output all positions
 -raw               disables noise and gaussian filters ( minhet must be set to 0 additionally)
 -call              callGT (call a Genotype, primitive representation)
 -emitDP            emit depth value per position
 -noindels          forget indels
 -nogenosort        do not sort genotypes for frequency
 -vcf               output VCF
 -l	<chrM>          interval (contig name [chrM] or contig coordinates [chrM:1234-5678])
 -lf <SCREEN>       log file; default: print to screen
 -ll <INFO>         log level: ERROR,INFO,DEBUG; default: INFO
 -h                 print this helptext
 -man               show man page


=head1 AUTHOR

Riccardo Berutti

=cut

