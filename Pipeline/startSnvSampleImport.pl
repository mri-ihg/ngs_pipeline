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

my $run = 1;

# database
my $dbh           = "";
my $sql           = "";
my $sth           = "";
my $logfile       = "SCREEN";
my $loglevel      = "INFO";

my $settings  	  = "hg19_wholegenome";
my $doJob     	  = 0;
my $snvOnly       = 0;
my $processLimit  = 9999;
my $help          = 0;
my $man           = 0;
my $counter       = 0;
my $numberOfSamples = 0;

# DZHKomics 1000 German Genomes samples are never imported/shown unless required
my $DZHKomics 		  = 0;

# Lock import to one sample (when needed to import one fast)
my $oneSampleImport   = "";


GetOptions(
	"se=s"      => \$settings,
	"do"        => \$doJob,
	"snvonly"   => \$snvOnly,
	"limit=s"   => \$processLimit,
	"lf=s"      => \$logfile,
	"ll=s"      => \$loglevel,
	"DZHKomics" => \$DZHKomics,
	"sample=s"	=> \$oneSampleImport,
	"h"         => \$help,
	"man"       => \$man
);

pod2usage( { -exitval => 0, -verbose => 1 } ) if $help;
pod2usage( { -exitval => 0, -verbose => 2 } ) if $man;

my $params = Utilities::getParams();


my $exomedb                   = $params->{settings}->{$settings}->{exomedb}->{database};
my $snvTable                  = $params->{settings}->{$settings}->{exomedb}->{snvtable};
my $snvsampleTable            = $params->{settings}->{$settings}->{exomedb}->{snvsampletable};

my $coredb                    = $params->{coredb}->{database};
my $sampleTable               = $params->{coredb}->{sampletable};


Utilities::initLogger( $logfile, $loglevel );
my $logger = Utilities::getLogger();


# Check if another instance is running
my $amIRunning=`ps aux | grep startSnvSampleImport.pl | grep -v grep | grep -v edit | wc -l`;
if ( $amIRunning > 1 )
{
	$logger->error("Another startSnvSampleImport.pl instance is running.");
	exit(1);	
}

# Connect to Exome DB
$dbh = Utilities::connectExomeDB($settings);

# Preparing query: if one sample specified then lock to one sample:
$oneSampleImport = "S.name=\"".$oneSampleImport."\" and " if ($oneSampleImport ne "");

# Select for each sample only the most recent pipeline (2017-05-23 RB)
# and if DZHKomics option selected then show only the ones from the project DZHKomics, otherwise show the others
my $DZHKomicsSQL = ($DZHKomics ? "=" : "<>");
my $getImportFiles_sql = qq{
		select 
			S.idsample, S.name, P.currentsettings, P.snvsampleImportFile, P.idpipeline 
		from 
			$coredb.pipeline P 
		inner join 
			$coredb.sample S on S.idsample=P.idsample  
		where 
			$oneSampleImport
			P.ready2importsnvsample=? and 
			P.snvsampleimportFinished=? and 
			currentsettings=? and 
			endtime<>0 and
			S.idproject $DZHKomicsSQL (  select idproject from $coredb.project where pdescription="DZHKomics" ) and
			S.idproject $DZHKomicsSQL (  select idproject from $coredb.project where pdescription="DZHK_Omics_Plus" ) and  
			S.idproject NOT IN (  select idproject from $coredb.project where pdescription="Controls" ) and
			P.idpipeline in ( select max(idpipeline) from $coredb.pipeline group by idsample,currentsettings ) 
		order by P.idpipeline desc 
		limit $processLimit 
		};
		
my $getImportFiles_sth = $dbh->prepare($getImportFiles_sql) || $logger->error("Can't prepare statement: $DBI::errstr");

# Update pipeline					
my $updatePipeline_sql = qq{update $coredb.pipeline set ready2importsnvsample=?,snvsampleimportFinished=? where idsample=? and currentsettings=?};
my $updatePipeline_sth = $dbh->prepare($updatePipeline_sql) || $logger->error("Can't prepare statement: $DBI::errstr");

# Delete SNVs for the sample
my $deleteSNVs_sql = qq{delete from $exomedb.$snvsampleTable where idsample=?};
my $deleteSVNs_sth = $dbh->prepare($deleteSNVs_sql) || $logger->error("Can't prepare statement: $DBI::errstr");


# Check if sample id is valid
$getImportFiles_sth->execute(1,0,$settings) || $logger->error($DBI::errstr);
$numberOfSamples = $getImportFiles_sth->rows;

# Do import:
if ($doJob)
{
	#Prepare DB to the import
	$dbh->do("SET \@\@session.unique_checks = 0");
	$dbh->do("SET \@\@session.foreign_key_checks = 0");
	$dbh->do("SET \@\@session.sql_log_bin=0");
}
else
{
	printf "\n";
	printf "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n";
	printf "\| %-15s \| %-15s \| %-25s \| %-110s \| %-110s \|\n", "SampleName", "idPipeline", "Settings", "SNPs","SVs";
	printf "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n";
}

# Samples whose import failed
my @failedsamples;

# Process each pipeline
while (my ($sampleID, $sampleName, $pipelineSettings, $importFile, $idpipeline) =  $getImportFiles_sth->fetchrow_array()) {
	
	$counter++;
	
	
	# my $importFile; 
	my $mergedSV = dirname( abs_path($importFile) )."/merged.sv.vcf";
	my $mergedSVText = ( -e $mergedSV ? $mergedSV : "MISSING SVs"); 

	if($doJob)
	{
		# Check that files exist or don't waste time:
		
		my $missingfile = 0;
		
		if ( ! -f $importFile )
		{
			$logger->error("SNVs not found for sample $sampleName at $importFile");
			$missingfile++;
		}
		
		if ( ! -f $mergedSV )
		{	
			$logger->error("SVs not found for sample $sampleName at $mergedSV");
			$missingfile++;
		}
		
		if ( $missingfile  ){
			push( @failedsamples, $sampleName);
		}
		else
		{	
		
			$logger->info("Processing sample $sampleName ... \t($counter / $numberOfSamples)");
	
			#Not necessary - checks that the software is not running multiple times
			#$updatePipeline_sth->execute(0,0,$idpipeline) || $logger->error($DBI::errstr);
			
			#delete
			$deleteSVNs_sth->execute($sampleID) || $logger->error($DBI::errstr);
			
			# Importing tsv variant file	
			my $failed=0;
			$logger->info("Importing file $importFile into $exomedb.$snvsampleTable");
			my $sql="LOAD DATA LOCAL INFILE '".$importFile."' INTO TABLE $exomedb.$snvsampleTable";
			$logger->info("Executing statement: $sql");
			$dbh->do($sql) or $failed=1;
			
			if ($failed) {
	
				# FAILURE: 
				$logger->error("Import of $importFile FAILED! Continuing with the next file, if present...");
				push( @failedsamples, $sampleName); 
	
			} else {
	
				# SUCCESS: 
				$logger->info("Executing statement: $sql FINISHED");
				$logger->info("Import of $importFile OK");
				
				#Calculate variantstat
				my $command = $prog_path . "/updateVariantStat.pl -s $sampleName -se $pipelineSettings -pipe $idpipeline -lf $logfile -ll $loglevel";
				if (&Utilities::executeCommand($command, "Updating variantstat for $sampleName", $logger)) {
					$logger->error("Updating variantstat for $sampleName FAILED! Continuing with the next sample...");
					push( @failedsamples, $sampleName);
				} else {
					$logger->info("Updating variantstat for $sampleName OK");
				
							
					# Import SV
					if ( -f $mergedSV )
					{
					     if ( ! $snvOnly )
					     {
						my $command = $prog_path . "/insertSV.pl -i $mergedSV -s $sampleName -se $pipelineSettings -d -lf $logfile -ll $loglevel";
						if (&Utilities::executeCommand($command, "Inserting SVs for $sampleName", $logger)) {
							$logger->error("SVs import for sample $sampleName at $mergedSV FAILED! Continuing with the next sample...");
							push( @failedsamples, $sampleName);
						} else {
						
							$logger->info("SVs import for sample $sampleName at $mergedSV OK");	
							
							# Only if SNV and SV import successful:
							# Reset current and past pipelines for the sample 
							$updatePipeline_sth->execute(1,1,$sampleID, $pipelineSettings) || $logger->error($DBI::errstr);
							
						}
					     }
					     else
					     {
					     	$logger->info("Skipping SVs");
						$updatePipeline_sth->execute(1,1,$sampleID, $pipelineSettings) || $logger->error($DBI::errstr);
					     }
					}
					else
					{
						$logger->error("SVs not found for sample $sampleName at $mergedSV");
						push( @failedsamples, $sampleName);
					}
				}
			}
		}
	}
	else
	{
		# Or show just the sample(s) to be processed
		printf "\| %-15s \| %15s \| %-25s \| %-110s \| %-110s \|\n", $sampleName, $idpipeline, $pipelineSettings, $importFile, $mergedSVText;
	}
}

#Import/display finished
if($doJob)
{
	$logger->info("Import finished");
	
	
	if ( ( scalar @failedsamples ) > 0 )
	{
		print "\nThere are ".(scalar @failedsamples)." failed samples:\n";	
		foreach my $failedsample (@failedsamples)
		{
			printf $failedsample."\n"; 
		}
	}
}
else
{
	printf "\| %-287s \|\n", "No samples to be inserted!" if $counter == 0;
	printf "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n";
	printf "\n"
}


=head1 NAME

startSnvSampleImport.pl

=head1 SYNOPSIS

startSnvSampleImport.pl [ -se pipeline_settings ] [ -do ]

=head1 DESCRIPTION

This script checks the pipeline-table if new snvsample-import files are ready to import and issues the mysqlimport command for each of them.
SVs are also imported.

=head1 OPTIONS

 -se	settings; default: hg19_wholegenome
 -sample sample name to import
 -do	start insert; (default: no); if -do is not set, only samples that would be inserted are listed
 -limit max samples to process: default NO LIMIT
 -lf	log file; default: print to screen
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -DZHKomics  uses only DZHKomics project samples (without the switch they are excluded)
 -h	print this helptext
 -man	show man page

=head1 AUTHORS

Riccardo Berutti

=cut
