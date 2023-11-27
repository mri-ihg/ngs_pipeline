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
require $prog_path . "/Utilities.pm";
#require $prog_path . "/BEDRecord.pm";

# CONFIG
my $max_gvcf_count = 251;

# Options
my $sample_list = "";
my $gvcf_list = "";
my $ref = "";
my $regionFile  = "";
my $settings = "";
my $analysis_folder = "";
my $batch_in = "";
my $batch_out = "";
my $maxmemory = "24g";
my $isBEDArrayJob = 0;
my $BEDArrayJobNumber = -1;
my $libtypeconstraint = "";
my $libpairconstraint = "";
my $logfile="";
my $loglevel="INFO";
my $help = 0; 
my $man =  0;

# Get'em
GetOptions(
	"samples=s"		=> \$sample_list,
	"gvcfs=s"		=> \$gvcf_list,
	"ref=s"			=> \$ref,
	"l=s"			=> \$regionFile,
	"se=s"			=> \$settings,
	"i=s"			=> \$batch_in,
	"o=s"			=> \$batch_out,
	"m=s"			=> \$maxmemory,
	"ajb"			=> \$isBEDArrayJob,
	"ajboverride=s" => \$BEDArrayJobNumber,
	"p=s"			=> \$analysis_folder,
	"lt=s"			=> \$libtypeconstraint,
	"lp=s"			=> \$libpairconstraint,
	"lf=s" => \$logfile,
	"ll=s" => \$loglevel,
	"h"    => \$help,
	"man"  => \$man
);
	

pod2usage( {-exitval => 0  ,-verbose => 1} ) if $help;
pod2usage( {-exitval => 0  ,-verbose => 2} ) if $man;


#TODO: Conditions to crash :)
pod2usage( {-exitval => 1  ,-verbose => 1} ) if ( 
		( $batch_out eq ""                       ) ||
		( $settings eq "" && $sample_list ne ""  ) ||
		( $sample_list ne "" && $gvcf_list ne "" ) ||
		( $sample_list eq "" && $gvcf_list eq "" ) ||
		( $sample_list ne "" && $gvcf_list ne "" ) ||
		( $settings eq "" && $ref eq ""          ) ||
		( $settings eq "" && $analysis_folder eq "" ) 
		);
		
		# Amended
		#		( $isBEDArrayJob && $gvcf_list ne ""     ) ||

# Get parameters
my $params     = Utilities::getParams();

# Start logger
Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();

# Get tool paths
$ref  = ( $ref ne "" ? $ref : ( $params->{settings}->{$settings}->{reference} ? $params->{settings}->{$settings}->{reference} : "" ) );    
my $gatk = $params->{programs}->{gatk}->{path};
my $gatk4 = $params->{programs}->{gatk4}->{path};

my $java = $params->{programs}->{java}->{path};
my $tabix= $params->{programs}->{tabix}->{path};

# Set our version of Java in Path with priority so to overcome any current installation (required by GATK4)
$ENV{'PATH'}=dirname(abs_path($java)).":".$ENV{'PATH'};

# DB Parameters
	# Hosts (must be one)
	if ( $params->{coredb}->{host} ne $params->{solexadb}->{host} )
	{
		$logger->error("Solexa DB and Core DB must be on the same machine to run this script.");
		exit(-1);
	}
	#Databases
	my $coredb   = $params->{coredb}->{database};
	my $solexadb = $params->{solexadb}->{database};
	#Tables
	my $sampletable 	= $coredb.".".	$params->{coredb}->{sampletable};
	my $projecttable 	= $coredb.".".	$params->{coredb}->{projecttable};
	my $organismtable   = $coredb.".".	$params->{coredb}->{organismtable};
	my $sample2librarytable = $solexadb.".".$params->{solexadb}->{sample2librarytable};
	my $librarytable	= $solexadb.".".$params->{solexadb}->{librarytable};
	my $libtypetable	= $solexadb.".".$params->{solexadb}->{libtypetable};
	my $libpairtable	= $solexadb.".".$params->{solexadb}->{libpairtable};

# Connect to database:
my $dbh = Utilities::connectCoreDB();

# Get input analysis folder:
if ( $analysis_folder eq "" )
{
	if ( $params->{settings}->{$settings}->{analysis}->{folder} )
	{
		$analysis_folder = $params->{settings}->{$settings}->{analysis}->{folder};
	}
	else
	{
		$logger->error("No default analysis folder defined for settings $settings. Specify with -p as input parameter.");
		exit(1); 
	}
}

# Reference?
if ( $ref eq "" )
{
	$logger->error("No default reference for settings $settings and no reference specified. Use -ref");
	exit(-1);
}
$logger->debug("Using reference $ref");

# Batch out should be a filename
if (abs_path($batch_out) eq abs_path(dirname($batch_out)))
{
	$logger->error("Batch out must be a filename, not a folder name");
	exit(-1);
}


# Get SGE BATCH status
my $listPrefix = "";

# SGE JobID
if ( $isBEDArrayJob )
{
	# If parameter passed override SGE_TASK_ID
	$ENV{SGE_TASK_ID} = ( $BEDArrayJobNumber != -1 ? $BEDArrayJobNumber : $ENV{SGE_TASK_ID} );
}


# How to find GVCFs
my $gvcf_subpath=""; # path inside sample folder

# Batch in:
my $batch_in_gvcf = "";

# Defaults for the path inside sample folders:
 	# BedArray
	my $gvcfSubFolder_BEDArray	="HaplotypeCaller";
	my $gvcfSuffix_BEDArray		=".haplotypecaller.gvcf.gz";
	# Normal
	my $gvcfSubFolder_Exome		=".";
	my $gvcfSuffix_Exome		="gatk.ontarget.haplotypecaller.gvcf.gz";

# Get region file and be flexible as a BEDArrayJob or a normal job, create batch_in filename, create batch out filename
if($isBEDArrayJob &&  $ENV{SGE_TASK_ID} && $params->{settings}->{$settings}->{genome_splits} ){
	$regionFile = $params->{settings}->{$settings}->{genome_splits}."/".$ENV{SGE_TASK_ID}.".bed";
	$listPrefix = $ENV{SGE_TASK_ID}.".";
	#$outfile = dirname($outfile)."/".$ENV{SGE_TASK_ID}.".".basename($outfile);
	# Build gvcf subpath
	$gvcf_subpath=$gvcfSubFolder_BEDArray."/".$ENV{SGE_TASK_ID}.$gvcfSuffix_BEDArray;
	# Batch in
	$batch_in_gvcf = ( $batch_in ne "" ? dirname($batch_in)."/".$ENV{SGE_TASK_ID}.".".basename($batch_in) : "" );
	# Batch out
	$batch_out = dirname($batch_out)."/".$ENV{SGE_TASK_ID}.".".basename($batch_out);
}
elsif ( $regionFile ne "" )
{
	# No array job
	$isBEDArrayJob=0;
	# Build gvcf subpath
	$gvcf_subpath=$gvcfSubFolder_Exome."/".$gvcfSuffix_Exome;
	# Batch in
	$batch_in_gvcf = $batch_in;
	# Batch out
	$batch_out = $batch_out;
}
elsif ( $params->{settings}->{$settings}->{normalchromosomes} )
{
	# No array job
	$isBEDArrayJob=0;
	$regionFile = $params->{settings}->{$settings}->{normalchromosomes};
	# Build gvcf subpath
	$gvcf_subpath=$gvcfSubFolder_Exome."/".$gvcfSuffix_Exome;
	# Batch in
	$batch_in_gvcf = $batch_in;
	# Batch out
	$batch_out = $batch_out;
}
else
{
	$logger->error("No target region defined for the sample set.");
	exit(-1);
}

# Make GVCF list
my @sample_name_list;
my @gvcf_list;

# If initial batch present add to the list:
if ( $batch_in_gvcf ne "" )
{
	if (!(-r -f $batch_in_gvcf ))
	{
		$logger->error("Initial batch $batch_in_gvcf does not exist");
		exit (-1);
	}
		
	push(@gvcf_list, $batch_in_gvcf);
	$logger->debug("Included batch $batch_in_gvcf to the merge");
}


# Accumulate libtype and organism to detect errors
my $libtype_hist="";
my $organism_hist="";

# This cycle opens input file
# 	- If it is sample list find sample and get gvcf name
# 	- if it is gvcf list then collect name
# Then check that the gvcfs exist. Exit if one does not exist
my $sample_counter=0;
my $errors=0;
open( IN,  ($sample_list ne "" ? $sample_list : $gvcf_list ))  || die $logger->error("Cannot open input file ".($sample_list ne "" ? $sample_list : $gvcf_list )."\n");
while (my $item = <IN>)
{
	chomp $item;
	# Rename for clarity
	my $sample_name = $item;
	
	# Skip empty lines;
	next if $sample_name eq "";
	
	# Prepare GVCF file name
	my $gvcf="";
	
	# If Sample list
	if ( $sample_list ne "" )
	{
		# Constraint conditions for libtype (genomic, exomic, etc) and for libpair (paired-end, single-end)
		# necessary when a sample has multiple libs with different type or paired states.
		my $libtypecondition= ( $libtypeconstraint ne "") ? "and LT.ltlibtype=\"$libtypeconstraint\"" : "";
		my $libpaircondition= ( $libpairconstraint ne "") ? "and LP.lplibpair=\"$libpairconstraint\"" : "";
		
		# Find sample
		my $sql = qq{
				select 
				    S.name, P.pname, O.orname, group_concat(distinct(LT.ltlibtype)), group_concat(distinct(LP.lplibpair))
				from 
				    $sampletable S 
				    inner join $projecttable P on P.idproject = S.idproject 
				    inner join $organismtable O on O.idorganism=S.idorganism
				    inner join $sample2librarytable S2L on S2L.idsample=S.idsample 
				    inner join $librarytable L on S2L.lid=L.lid  
				    inner join $libtypetable LT on L.libtype=LT.ltid 
				    inner join $libpairtable LP on LP.lpid=L.libpair  
				where 
				    S.name=\"$sample_name\"
				    and S.sbam<>\"\" 
				    $libtypecondition
				    $libpaircondition
				    and S.nottoseq=0 
				    and L.lfailed=0 
				group by 
				    S.name 
		};
		
		my $sth = $dbh->prepare($sql);
		$sth->execute();
		
		#Get libtype; get libpair; get organism
		#  libtype and liborganism serve to check that the batch is made of compatible samples
		#  libtype and libpair can be multiple (group conact distinct, in case they are the script crashes and asks for specification)
		my ($sample, $projectname, $organism, $libtype, $libpair)=$sth->fetchrow_array();

				
		# Sample not found
		if ( !defined($sample) )
		{
			$logger->error("Sample $sample_name not found into database.");
			$errors++;
			next;
		}
				
		# If there are commas in libtype and libpair
		my @libtypes = split(",", $libtype);
		my @libpairs = split(",", $libpair);
		
		if ((@libtypes != 1) || (@libpairs != 1))
		{
			$logger->error("Sample $sample has multiple libtypes ($libtype) or libpairs ($libpair), please specify libtype and libpair as arguments of the analysis");
			$errors++;
			next;
		} 

		# Check organism compatibility
		if ( $organism_hist eq "" )
		{
			$organism_hist=$organism;
		}
		elsif ($organism_hist ne $organism )
		{
			$logger->error("Sample $sample is organism: $organism. First samples were $organism_hist");
			$errors++;
			next;
		}
		
		# Check library type compatibility
		# TODO: Possible exception exome/genome
		if ( $libtype_hist eq "" )
		{
			$libtype_hist=$libtype;
		}
		elsif ($libtype_hist ne $libtype )
		{
			$logger->error("Sample $sample has libtype: $libtype");
			$errors++;
			next;
		}
		
		# Build Gvcf path - Kept into account BEDARRAY JOB ($gvcf_subpath pre calculated)
		$gvcf = $analysis_folder."/".$projectname."/".$sample_name."/".$libtype."out"."/".$libpair."out"."/".$gvcf_subpath;
		
	}
	else
	{
		# GVCF list
		
		# If BEDArray then GVCF Path (try not to miss it)
		if($isBEDArrayJob &&  $ENV{SGE_TASK_ID} && $params->{settings}->{$settings}->{genome_splits} ){
			$gvcf = dirname($item)."/".$ENV{SGE_TASK_ID}.".".basename($item);
			if ( ! -f $gvcf )
			{
				$gvcf = dirname($item)."/".$ENV{SGE_TASK_ID}.".gvcf.gz";
			}
			
			if ( ! -f $gvcf )
			{
				$gvcf = $item."/".$ENV{SGE_TASK_ID}.".gvcf.gz";
			}						
		}
		else # Else GVCF List
		{
			$gvcf=$item;	
		}
		
	}
	
	# Error if expected GVCF not found
	if ( ! ( -r -f $gvcf  ))
	{
		$logger->error("Required gvcf file $gvcf not found.");
		$errors++;
		next;
	}

	push(@gvcf_list, $gvcf);
	$logger->debug("Included sample $sample_name from path $gvcf");
	
	$sample_counter++;
}

# Catch errors. It's not interrupted before to have a list of the problems all in once
if ($errors>0)
{
	$logger->error("Found $errors sample errors. Terminated");
	exit(-1);
}

# Limit to combine max $max_gvcf_count GVCFS
if ( $sample_counter>$max_gvcf_count )
{
	$logger->error("This script can combine maximum $max_gvcf_count samples. $sample_counter present as input. Please subbatch them");
	exit(-1);
}

# Write list
my $batch_out_list = $batch_out.".combinegvcf.list"; 
    $batch_out_list =~ s/\.gvcf\.gz//g; #Remove .gvcf.gz extension if in;

open( GVCF_OUT, ">$batch_out_list" ) || die $logger->error("Cannot open $batch_out_list for writing\n");
foreach my $gvcfsToWrite (@gvcf_list)
{
	print GVCF_OUT "$gvcfsToWrite\n";
}
close(GVCF_OUT);


# Issue CombineGVCFs command! 
my $batch_out_file = $batch_out;
    $batch_out_file .= ".gvcf.gz" if ($batch_out !~ m/\.gvcf\.gz$/);
	
my $batch_out_gatklog = $batch_out.".gatklog";
    $batch_out_gatklog =~ s/\.gvcf\.gz//g; #Remove .gvcf.gz extension if in;


	#Legacy GATK3 command, to remove 
	#my $gatk_command = " \\
	#	$java -XX:ParallelGCThreads=1 -Xmx".$maxmemory." \\
	#	-jar $gatk \\
	#		-R $ref \\
	#		-T CombineGVCFs \\
	#		-L $regionFile \\
	#		--validation_strictness LENIENT \\
	#		--variant $batch_out_list \\
	#		-o $batch_out_file \\
	#";

# GATK4 fully supported
my $gatk_command = " \\
		$gatk4 --java-options  \"-Xmx$maxmemory -Xms$maxmemory\" CombineGVCFs \\
		-R $ref \\
   		-L $regionFile \\
   		--variant $batch_out_list \\
   		--output $batch_out_file";


$logger->info("Launching GATK CombineGVCFs with command: $gatk_command");

# If batch already exists do not create it, but issue WARNING
if ( ! -f $batch_out_file )
{
	# If merging only one file than issue hard link:
	if ( (scalar @gvcf_list) > 1 )
	{
	    if (&Utilities::executeCommand($gatk_command, "Launching GATK CombineGVCFs to merge GVCFs from the list $batch_out_list into $batch_out_file", $logger)) {
	    		$logger->error("Error merging the files from the list $batch_out_list into $batch_out_file");
	    		exit(-1);
	    }
	}
	elsif ( (scalar @gvcf_list) == 1 )
	{
		$logger->warn("You are merging just one file. Linking to the old file");
		my $ln_command = "ln -s ".$gvcf_list[0]." ".$batch_out_file; 
		 if (&Utilities::executeCommand($ln_command, "Hard linking ".$gvcf_list[0]." to $batch_out_file", $logger)) {
	    		$logger->error("Error linking ".$gvcf_list[0]." to $batch_out_file");
	    		exit(-1);
	    }
	}
	else
	{
		$logger->error("Merging no files");
	}
}
else
{
	$logger->warn("The output batch $batch_out_file already exists. It will not be overwritten, delete it if you want to regenerate it.");
}

# Tabix indexing file
my $tabix_command ="$tabix -f -p vcf $batch_out_file";
$logger->info("Launching Tabix Indexing with command: $tabix_command");

my $tabix_out = $batch_out_file.".tbi";

if ( ! -f $tabix_out )
{
	if (&Utilities::executeCommand($tabix_command, "Tabix indexing file $batch_out_file", $logger)) {
			$logger->error("Tabix failed in indexing file $batch_out_file");
			exit(-1);
	}
}
else
{
	$logger->warn("File $batch_out_file already indexed with tabix, skipping");
}


$logger->info("Merged file $batch_out_file ready");


	
#TODO: Sistema roba
#	#$command .= "-L $regionFile \\\n" if $regionFile ne "";



=head1 NAME

combineGVCFs.pl

=head1 SYNOPSIS

 combineGVCFs.pl -samples SAMPLE.list -l region.bed -se hg19_wholegenome -o /path/to/combined_batch (without extension)

=head1 DESCRIPTION

This is a wrapper script for GATK combineGVCFs. It builds a combined GVCF out of the provided sample list. It searches for the samples into the default
folder for the specified settings
 

=head1 OPTIONS

 -samples <sample list>	 text file containing sample names whose GVCFs to merge, one line per sample; required OR -gvcfs
 -gvcfs <gcvf list>  text file containing gvcf paths; required OR -samples. If -ajb GVCF folders
 -l	    <region.bed> BED file with regions for which variants should be merged;  optional
 -se	<settings>; required
 -i     <infile> previous batch to merge data to; if empty no initial file is appended
 -o     <outfile> path to the output file without .gvcf.gz (.<SGE_TASK_ID> will be appended between it and the extension if genome splits are used); required 
 -m	    <memory> max. memory for each Java GATK job; default 4g
 -ajb	is SGE Array job; if this flag is chosen, the script runs as a part of a SGE array job
		this means, that the environmental variable SGE_TASK_ID holds the number of this job
		The script will only process the file "SGE_TASK_ID.bed" in the folder $params->{settings}->{$settings}->{genome_splits}.
		SGE_TASK_ID will be used as a fileprefix.
 -p		<projectfolder> Folder in which S0xyz/Sample/...out/...out/*gvcf files are located
 -lt	<libtype> (exomic, genomic, mtDNA, etc), necessary if some of the samples have multiple libraries of different type
 -lp	<libpair> (paired-end, single-end, etc), necessary if some of the samples have multiple libraries with different paired status
 -lf	<log file>; default: pipeline.log
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page

=head1 AUTHOR

Riccardo Berutti 

=cut
 
