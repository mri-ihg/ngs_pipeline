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
use DateTime;

my $prog_path = dirname( abs_path($0) );
require $prog_path . "/Utilities.pm";
#require $prog_path . "/BEDRecord.pm";


# database
my $dbh           = "";
my $sql           = "";
my $sth           = "";
my $logfile       = "SCREEN";
my $loglevel      = "INFO";

# Options
my $inputlist     = "";
my $gvcfout		  = "";
my $isBEDArray    = 0;
my $batch_name    = "";
my $timestamp_mysql = "";
my $multisettings = "";
my $libtypeconstraint = "";
my $isbatch		  = 0;
my $isdatafreeze  = 0;
my $maxRam        = "UNDEF";
my $help          = 0;
my $man           = 0;



GetOptions(
	"samples=s" 		=> \$inputlist,
	"gvcfs=s"			=> \$inputlist,
	"o=s"		        => \$gvcfout,
	"ajb"				=> \$isBEDArray,
	"date=s"			=> \$timestamp_mysql,
	"name=s"			=> \$batch_name,
	"multisettings=s"	=> \$multisettings,
	"batch"				=> \$isbatch,
	"datafreeze" 		=> \$isdatafreeze,
	"m=s"				=> \$maxRam,
	"lt=s"			    => \$libtypeconstraint,
	"lf=s" 				=> \$logfile,
	"ll=s" 				=> \$loglevel,
	"h"    => \$help,
	"man"  => \$man
);

pod2usage( { -exitval => 0, -verbose => 1 } ) if $help;
pod2usage( { -exitval => 0, -verbose => 2 } ) if $man;
pod2usage( { -exitval => 0, -verbose => 1 } ) if (( $isbatch eq 0 && $isdatafreeze eq 0 ) || ( $isbatch eq 1 && $isdatafreeze eq 1 ) );

# Start Logging
Utilities::initLogger( $logfile, $loglevel );
my $logger = Utilities::getLogger();

# MaxRam:
# If supplied use specified value, otherwise differentiate for batch (48g) and datafreeze (128g) 
$maxRam = 
			$maxRam ne "UNDEF" ?	$maxRam :
									$isbatch ? "64g" : "128g" ; 

# Getting parameters:
my $params = Utilities::getParams();

	# Tools    
	my $gatk  = $params->{programs}->{gatk}->{path};
	my $gatk4 = $params->{programs}->{gatk4}->{path};
	my $java  = $params->{programs}->{java}->{path};
	my $tabix = $params->{programs}->{tabix}->{path}; 
	
	# JAVA override system install
	# Set our version of Java in Path with priority so to override any current installation (required by GATK4)
	$ENV{'PATH'}=dirname(abs_path($java)).":".$ENV{'PATH'};

	# Get the CoreDB
	my $coredb                    = $params->{coredb}->{database};
	my $sampletable               = $params->{coredb}->{sampletable};

	# Multisample db and tables
	my $multisampledb			  = $params->{multisampledb}->{database};
	my $batchtable                = $params->{multisampledb}->{batchtable};
    my $datafreezetable           = $params->{multisampledb}->{datafreezetable};
    my $sample2batchtable         = $params->{multisampledb}->{sample2batchtable};
    my $settingstable             = $params->{multisampledb}->{settingstable};

	# Get the multisample settings
	my $coredbh = Utilities::connectCoreDB();
	
	my $settings_sql = qq{select idmssettings, selectquery, settings, bedarray from $multisampledb.$settingstable where name=?};
	my $settings_sth = $coredbh->prepare($settings_sql);
	$settings_sth->execute($multisettings) || die $logger->error("Error in getting settings: $settings_sql");
	
	my ($idmultisettings, $sampleselect_sql, $settings, $perBEDArray)=$settings_sth->fetchrow_array();
	
	# No multisample settings found
	if ( ! defined $idmultisettings )
	{
		print "No multisample-calling settings named $multisettings. Create settings first\n";
		exit 1;
	}
	
	# Refgenome
	my $ref  = $params->{settings}->{$settings}->{reference};
	 
	# Genome splits:
	my $genome_splits			  = $params->{settings}->{$settings}->{genome_splits} if ($perBEDArray);
	my @bedArrayItems			  = glob($genome_splits."/*.bed") if ($perBEDArray);
	
	# Get multisample path:
	my $multisample_basepath      = $params->{settings}->{$settings}->{analysis}->{multisample} if defined $params->{settings}->{$settings}->{analysis}->{multisample} or die ("No multisample path defined for pipeline settings $settings");

	
	# Add suffix (solve mess)
	$gvcfout.=".gz" if ($gvcfout =~ m/\.gvcf$/);	
	$gvcfout.=".gvcf.gz" if ($gvcfout !~ m/\.gvcf\.gz$/);


	# Errorfile
	my $errorfile = dirname($gvcfout)."/ERROR";
	
	# If errorfile exists already save time  and do not run
	if ( -f $errorfile )
	{
		$logger->error("Errorfile $errorfile found. Not executing to save time.");
		exit -1;
	}
	
	# Test that inputfile exists:
	if ( ! -f $inputlist )
	{
		$logger->error("Sample/GVCF list $inputlist not found.");
		exit -1;		
	}

	# Define $ENV{SGE_TASK_ID} if not.
	$ENV{SGE_TASK_ID}="" if (! defined $ENV{SGE_TASK_ID} );
	
# COMBINEGVCF
#GetOptions(
#	"samples=s"		=> \$sample_list,
#	"gvcfs=s"		=> \$gvcf_list,
#	"ref=s"			=> \$ref,
#	"l=s"			=> \$regionFile,
#	"se=s"			=> \$settings,
#	"i=s"			=> \$batch_in,
#	"o=s"			=> \$batch_out,
#	"m=s"			=> \$maxmemory,
#	"ajb"			=> \$isBEDArrayJob,
#	"p=s"			=> \$analysis_folder,
#	"lt=s"			=> \$libtypeconstraint,
#	"lp=s"			=> \$libpairconstraint,
#	"lf=s" => \$logfile,
#	"ll=s" => \$loglevel,
#	"h"    => \$help,
#	"man"  => \$man

	
# Build CombineGVCFCommand:

my 	$command  = "perl ".$prog_path."/combineGVCFs.pl ";
	$command .= " -samples $inputlist " if $isbatch;			# Samples if batch, gvcfs if datafreeze
	$command .= " -gvcfs   $inputlist " if $isdatafreeze;		#
	$command .= " -o $gvcfout ";
	$command .= " -se $settings "; 
	$command .= " -ajb " if $isBEDArray;
	$command .= " -m $maxRam";
	$command .= " -lt $libtypeconstraint" if $libtypeconstraint ne "";
	$command .= " -lf $logfile -ll $loglevel ";
	
	
# Track error
my $error = 0;

# Execute command and grab output
# TABIX is included by CombineGVCFs,pl
if (&Utilities::executeCommand($command, "CombineGVCFs.pl for inputlist: $inputlist, BEDArray (status: $isBEDArray - TaskID(ifapp): ".$ENV{SGE_TASK_ID}.") - Combining for ". ($isbatch ? "BATCH " : "DATAFREEZE"), $logger)) 
{
		$logger->error("Error in CombineGVCFs.pl Command line : $command");
		$error=1;
		
		#ON ERROR output ERROR file in batch folder to mark it for deletion after the whole thing is over.
		CreateErrorFile();	
}


# DATAFREEZE
# Genotype GVCF for Datafreeze 
if ( $isdatafreeze && ($error eq 0) )
{
	# Build name if bedarraytast
	my $gvcfout_bedarrayaware;
		
	if($isBEDArray &&  $ENV{SGE_TASK_ID} && $params->{settings}->{$settings}->{genome_splits} )
	{
		$gvcfout_bedarrayaware = dirname($gvcfout)."/".$ENV{SGE_TASK_ID}.".".basename($gvcfout) ;
	}
	else
	{
		$gvcfout_bedarrayaware = $gvcfout;
	}
		
	# GenotypeGVCF
	my $vcfout = $gvcfout_bedarrayaware;
		$vcfout =~ s/.gvcf.gz/.vcf.gz/g;
	
	# GATK3 legacy command - to be removed
	#my $command = "
	#				$java -XX:ParallelGCThreads=1 -Xmx$maxRam -jar $gatk \\
	#				-R $ref -T GenotypeGVCFs \\
   	#				--variant $gvcfout_bedarrayaware \\
   	#				-o $vcfout";
   					
   	# GATK4
   	# TODO: Use genotypeGVCF.pl from pipeline
	my $command = "
					$gatk4 --java-options  \"-Xmx$maxRam\" GenotypeGVCFs\\
					-R $ref \\
					--variant $gvcfout_bedarrayaware \\
					-o $vcfout";
					# --max-alternate-alleles 3 \\ # if necessary can be added here.
   					
    $logger->debug($command);
	if (&Utilities::executeCommand($command, "GenotypeGVCF for datafreeze $gvcfout_bedarrayaware", $logger)) 
	{
   		$logger->error("Error genotyping GVCF $gvcfout_bedarrayaware");
   		$error=1;
   		CreateErrorFile();
	}
	
	
	# Tabix indexing file
	my $tabix_command ="$tabix -f -p vcf $vcfout";

	$logger->info("Launching Tabix Indexing with command: $tabix_command");
	if (&Utilities::executeCommand($tabix_command, "Tabix indexing file $vcfout", $logger)) {
		$logger->error("Tabix failed in indexing file $vcfout");
		$error=1;
		CreateErrorFile();
	}
	
}


# Reconnect to mySQL:
# it's harmless. It will likely have timed out after the long execution.
$coredbh = Utilities::connectCoreDB();


# Log into BATCH or DATAFREEZE 
# (on duplicate key update -> BEDARRAY -> On error log ERROR for job)
# Depending on what I am running
my $table = $isbatch ? 
			"$multisampledb.$batchtable" :
			"$multisampledb.$datafreezetable";
			
my $idname = $isbatch ? 
				"idmsbatch" :
				"idmsdatafreeze" ;
			
$sql = "INSERT INTO $table 
				( 
					name, date, idmssettings, gvcf_path, errors, list 
				)
				VALUES
				( 
					\"$batch_name\", \"$timestamp_mysql\", $idmultisettings, \"$gvcfout\", $error, \"$inputlist\"  
				)
			ON DUPLICATE KEY UPDATE errors=(errors+$error), date=\"$timestamp_mysql\", $idname=LAST_INSERT_ID($idname)
			"; 

# Insert or update
$logger->debug($sql);
$sth = $coredbh->prepare($sql)
		  || $logger->error("Can't prepare statement: $DBI::errstr");
	$sth->execute() || $logger->error($DBI::errstr);
	# Get Inserted ID
    my $batchid = $coredbh->last_insert_id(undef, undef, qw($table $idname)) or exit $logger->error("Can't retrieve last inserted id on table $table! - SQL: $sql");


# If it's a batch, then insert mssample2batch 
if ( $isbatch )
{
	# Open sample file:
	open SAMPLELIST, "< $inputlist" or die $logger->error("Sample/GVCF list $inputlist cannot be opened.");
	
	while(my $sample_name = <SAMPLELIST>)
	{
		# Insert IGNORE (when necessary to FIX old batches)
		chomp $sample_name;	
		$sql = "INSERT IGNORE INTO $multisampledb.$sample2batchtable ( idsample, idmsbatch ) values ( (SELECT idsample from $coredb.$sampletable where name=\"$sample_name\") , $batchid )";
		$logger->debug($sql);					
		$sth = $coredbh->prepare($sql)
		  || $logger->error("Can't prepare statement: $DBI::errstr");
		$sth->execute() || $logger->error("Cannot insert sample $sample_name into batch $batch_name. Sample does not exist. DBIErr: ".$DBI::errstr );
	}
}


# Return 0 if no error, -1 if error
exit(-$error); 
 
 
sub CreateErrorFile
{
		if ( ! -f $errorfile )
		{
			$command="touch $errorfile";
			(&Utilities::executeCommand($command, "Creating ERROR file $errorfile to mark directory for deletion", $logger));
		} 	
}


=head1 NAME

multisample_batch.pl

=head1 SYNOPSIS


multisample_batch.pl 

This is a service script to launch datafreeze or batches. No user launch


=head1 DESCRIPTION

This is a service script to launch datafreeze or batches. No user launch

=head1 OPTIONS

 -h	          print this helptext
 -man	      show man page

=head1 AUTHOR

Riccardo Berutti

=cut

