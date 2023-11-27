#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use DBI;
use File::Basename;
use File::Path qw(make_path);
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);
use List::MoreUtils qw/ uniq /;
use Pod::Usage;
use DateTime;

my $prog_path = dirname( abs_path($0) );
require $prog_path . "/Utilities.pm";
require $prog_path . "/pipeconfig.pm";
#require $prog_path . "/BEDRecord.pm";

#Config
my $max_batch = 200; # Grouping, max 100 samples per multisample calling batch

# database
my $dbh           = "";
my $sql           = "";
my $sth           = "";
my $logfile       = "";
my $loglevel      = "INFO";

# Options
my $multisettings = "";
my $checkOnly     = 0;
my $noStart       = 0;
my $dataFreeze    = 0;
my $libtypeforce  = "";
my $help          = 0;
my $man           = 0;

my $isBedArray    = 1; # Defaults for genome

GetOptions(
	"mse=s" 		=> \$multisettings,
	"checkonly" 	=> \$checkOnly,
	"nostart"		=> \$noStart,
	"datafreeze"	=> \$dataFreeze,
	"forcelibtype=s"=> \$libtypeforce,
	"lf=s" => \$logfile,
	"ll=s" => \$loglevel,
	"h"    => \$help,
	"man"  => \$man
);

pod2usage( { -exitval => 0, -verbose => 1 } ) if $help;
pod2usage( { -exitval => 0, -verbose => 1 } ) if $multisettings eq "";
pod2usage( { -exitval => 0, -verbose => 2 } ) if $man;


# Logger for checkOnly must be on screen, but for pipeline must be in pipeline folder, so once i get multisettings i should start logging where proper!
if ( $checkOnly )
{
	$logfile  = "SCREEN";
	$loglevel = "ERROR";
}

# Getting run parameters:
my $params = Utilities::getParams();

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
	$settings_sth->execute($multisettings);
	
	my ($idmultisettings, $sampleselect_sql, $settings, $perBEDArray)=$settings_sth->fetchrow_array();
		 
	# No multisample settings found
	if ( ! defined $idmultisettings )
	{
		print "No multisample-calling settings named $multisettings. Create settings first";
		exit 1;
	}
	 
	# Get the settings:
	my $exomedb                   = $params->{settings}->{$settings}->{exomedb}->{database} if defined $params->{settings}->{$settings}->{exomedb}->{database} or die ("No pipeline settings named $settings.");
	my $snvtable                  = $params->{settings}->{$settings}->{exomedb}->{snvtable};
	my $snvsampletable            = $params->{settings}->{$settings}->{exomedb}->{snvsampletable};

	# Genome splits:
	my $genome_splits			  = $params->{settings}->{$settings}->{genome_splits} if ($perBEDArray);
	my @bedArrayItems			  = glob($genome_splits."/*.bed") if ($perBEDArray);
	
	# Get multisample path:
	my $multisample_basepath      = $params->{settings}->{$settings}->{analysis}->{multisample} if defined $params->{settings}->{$settings}->{analysis}->{multisample} or die ("No multisample path defined for pipeline settings $settings");

# Get a timestamp for the pipeline submission - it's important to date batches and datafreeze
# Formats: 
# - batch_timestamp 20170701_120004
# - mysql_timestamp	2017-07-01 12:00:04
my $datetime_now  = DateTime->now;   # Stores current date and time as datetime object
	my $datetime_date = $datetime_now->ymd;   # Retrieves date as a string in 'yyyy-mm-dd' format
	my $datetime_time = $datetime_now->hms;

	my $timestamp_file = $datetime_date."_".$datetime_time;
		$timestamp_file =~ s/\-//g;
		$timestamp_file =~ s/\://g;
	my $timestamp_mysql= $datetime_date." ".$datetime_time;

#TODO CLEANUP
#print "Timestamp:\nFile:\t$timestamp_file\nSql:\t$timestamp_mysql\n";
#exit;
	

# BUILD STRUCTURE:

# Make path for multisample calling if not present
#
#				make_path("$outdir/$project/$sample/$folder/$subfolder/HaplotypeCaller", {				#create directory for single files
#					mode => 0775
#				});

# Create folder structure:
#
	my $multisettings_folder 	= $multisample_basepath."/".$multisettings;		# Base for mssettings
	my $batches_folder 			= $multisettings_folder."/"."batches"; 			# where batches gvcfs are saved
	my $datafreeze_folder		= $multisettings_folder."/"."datafreeze";		# where datafreeze gvcfs are saved
	my $pipelines_folder		= $multisettings_folder."/"."pipelines";		# where pipeline logs are saved always

	# Main multisample folder
	if ( ! -d $multisample_basepath )
	{
		make_path($multisample_basepath, { mode => 0775 }) or die "Cannot create base path folder $multisample_basepath";	
	}
	
	# Make multisettings folder
	if ( ! -d $multisettings_folder )
	{
		make_path( $multisettings_folder, { mode => 0775 }) or die "Cannot create folder $multisettings_folder";	
	}

	# Make batches folder
	if ( ! -d $batches_folder )
	{
		make_path( $batches_folder, { mode => 0775 }) or die "Cannot create folder $batches_folder";	
	}	
	
	# Make datafreeze folder 
	if ( ! -d $datafreeze_folder )
	{
		make_path( $datafreeze_folder, { mode => 0775 }) or die "Cannot create folder $datafreeze_folder";	
	}

	# Make pipelines folder 
	if ( ! -d $pipelines_folder )
	{
		make_path( $pipelines_folder, { mode => 0775 }) or die "Cannot create folder $pipelines_folder";	
	}


# Logfile if not specified (and not checkonly - defaults to screen)
if ( $logfile eq "" )
{
	$logfile = $pipelines_folder."/".$timestamp_file."_pipeline.log";
}

# Start Logging
Utilities::initLogger( $logfile, $loglevel );
my $logger = Utilities::getLogger();


# Build some stats:
	# How many samples in this multisample ? 
	my $totalsamplelist_sth=$coredbh->prepare($sampleselect_sql);
	$totalsamplelist_sth->execute();
	my $stats_total_samples_in_msquery=$totalsamplelist_sth->rows();

	# How many samples/batches/datafreezes already processed? 
	my $totalsamples_in_batches_sql=qq{ 
		select 
			count(distinct(MSS2B.idsample)),
			count(distinct(MSS2B.idmsbatch)),
			count(distinct(MSD.idmsdatafreeze))
		from 
			$multisampledb.$sample2batchtable MSS2B 
			inner join $multisampledb.$batchtable MSB on MSB.idmsbatch=MSS2B.idmsbatch 
			left  join $multisampledb.$datafreezetable MSD on MSD.date>MSB.date 
		where
			MSB.idmssettings=$idmultisettings
			and MSB.errors = 0
	};
	my $totalsamples_in_batches_sth=$coredbh->prepare($totalsamples_in_batches_sql);
	$totalsamples_in_batches_sth->execute();	
	my ($total_samples_processed, $total_batches_processed, $total_datafreezes_processed)=$totalsamples_in_batches_sth->fetchrow_array();


# New samples to be processed
#
# 	Samples which are:
# 	- selected according to mssettings query
# 	- but not in batches for mssettings
# 	- have been run at least once for pipeline settings specified in mssettings
#
#TODO: pipeline table name change to some $paramz->
#

my $newsampleselect_sql=qq{
	select 
		S.idsample, S.name
	from 
		( $sampleselect_sql ) SL
		inner join $coredb.$sampletable S on S.idsample=SL.idsample  
	where 
		S.idsample
			not in 
			(
				select 
					idsample 
				from 
					$multisampledb.$sample2batchtable MSS2B
					inner join $multisampledb.$batchtable MSB on MSB.idmsbatch =MSS2B.idmsbatch
					inner join $multisampledb.$settingstable MSS on MSS.idmssettings=MSB.idmssettings
				where
					MSS.name="$multisettings"
					and MSB.errors = 0
			)
			and
		S.idsample 
			in
			(
				select
					distinct( idsample )
				from
					$coredb.pipeline P 
				where 
					currentsettings="$settings"
		  	)
};


my $newsamplelist_sth=$coredbh->prepare($newsampleselect_sql);
$newsamplelist_sth->execute();

# Ready for cycling through samples
my ($sampleid, $sample_name);
my $counter = 0;
my @batches; #Batches structure goes here
my @samples; 
my @samples_inbatch;
	my $batch_count = 1;
	my $batch_sample_count 		= 0;
	my $total_sample_count 		= 0;
	my $discarded_sample_count 	= 0;
	my @discarded_samples;
	
# Loop for samples
while ( ($sampleid, $sample_name) = $newsamplelist_sth->fetchrow_array()  )
{
	print " > ".$sample_name."\n";
	
	# Reset gvcf name
	my $gvcf_name = "";
	
	# Check that all the GVCF files exist, otherwise problems are  
	# Process according to perBEDArray or not 0 or 1;
	if ( $perBEDArray )
	{
		# Track errors while checking for the gvcfs to exist:
		my $bedarrayerror=0;
		
		# If perBEDArray check that all GVCFs exist:		
		for (my $bedArrayItem=1; $bedArrayItem<=(scalar @bedArrayItems); $bedArrayItem++ )
		{
			# Gvcf name
			if ( $libtypeforce ne "" )
			{
				$gvcf_name = Utilities::getGVCFPath("$sample_name", "$settings", "$bedArrayItem", "$libtypeforce");
			}
			else
			{
				$gvcf_name = Utilities::getGVCFPath("$sample_name", "$settings", "$bedArrayItem");
			}

			# Find sample if exists GVCF otherwise log info and go on. Do not CRASH! Think how many samples must be incomplete around...!
			if ( (! -f $gvcf_name)  || ($gvcf_name eq "") )
			{
				# Discarded stats and list	
				$discarded_samples[$discarded_sample_count] = $sample_name;
				$discarded_sample_count++;
				#Error
				$logger->error("Sample excluded: $sample_name. GVCF file $gvcf_name not found where expected ");
				#Track it to exit the main array
				$bedarrayerror=1;
				#End cycling in the sub-bed-array
				last;
			}
		}
		
		# Previous last was in the sub-bed array, in case of error go to the next sample
		next if ($bedarrayerror == 1);
				
	}
	else
	{
		
		# Gvcf name
		if ( $libtypeforce ne "" )
		{
			$gvcf_name = Utilities::getGVCFPath("$sample_name", "$settings", "0", "$libtypeforce");
		}
		else
		{		
			$gvcf_name = Utilities::getGVCFPath("$sample_name", "$settings", "0");
		}
		
		# Find sample if exists GVCF otherwise log info and go on. Do not CRASH! Think how many samples must be incomplete around...!
		if ( (! -f $gvcf_name)  || ($gvcf_name eq "") )
		{
			# Discarded stats and list
			$discarded_samples[$discarded_sample_count] = $sample_name;
			$discarded_sample_count++;
			# Error
			$logger->error("Sample excluded: $sample_name. GVCF file $gvcf_name not found where expected ");
			# Next sample
			next;
		}
	}


	# Increment total count of samples and count of samples in current batch
	# I want only existing samples here, that's why I checked $gvcf_name for existence before!
	$total_sample_count++;
	$batch_sample_count++;	

	# If batch full, new batch
	if ( $batch_sample_count > $max_batch )
	{
		$batch_sample_count = 1;
		$batch_count++;
	}
	
	# Now I can push back
	# 	Do I want to use gvcf name? For bedArrayJob it would be a mess. CombineGVCFs makes good use of sample names <3 
	#	$batches[$batch_count][$batch_sample_count]=$gvcf_name if $gvcf_name ne "";
	
	$batches[$batch_count][$batch_sample_count]=$sample_name if $gvcf_name ne "";
	$samples[$total_sample_count]="$sample_name";
	$samples_inbatch[$total_sample_count]="$batch_count";
	
	# TODO remove
	#print "\$samples[\$total_sample_count]=\$samples[$total_sample_count]=".$samples[$total_sample_count]."\n";
}

# Batch starts from 1 
$batch_count = 0 if ( $total_sample_count == 0 );

# If check only print info and exit 
if ($checkOnly)
{
	# Print how many samples, and how many batches I'm adding
	print "\n";
	print "====================================================================\n";
	print "Settings $multisettings\n";
	print "====================================================================\n";
	print "Query:\t$sampleselect_sql\n";
	print "Pipeline Settings:\t$settings\n";
	print "Total samples:\t$stats_total_samples_in_msquery\n";
	print "====================================================================\n";
	print "Samples already present to multisample for settings $multisettings\n";
	print "====================================================================\n";
	print "Sample count:\t$total_samples_processed\n";
	print "Batch count:\t$total_batches_processed\n";
	print "Datafr count:\t$total_datafreezes_processed\n";
	print "====================================================================\n";
	print "New samples to be added to multisample for settings $multisettings\n";
	print "====================================================================\n";
	print "New sample count:\t$total_sample_count\n";
	print "Batch count:\t$batch_count\n";
	print "Discarded count:\t$discarded_sample_count\n"; 
	print "====================================================================\n";
	print "Sample list to be added and in which batch\n";
	print "====================================================================\n";
	# Print sample list;
	for (my $loop=1; $loop<=$total_sample_count; $loop++)
	{
		print $samples[$loop]."\t".$samples_inbatch[$loop]."\n";
	}
	print "====================================================================\n";
	print "Discarded samples with missing GVCF files\n";
	print "====================================================================\n";
	foreach my $discarded_sample_name ( @discarded_samples )
	{
		print $discarded_sample_name."\n";
	}
	print "====================================================================\n\n";
	exit 0;
}

# End here if no samples to be added
if ( $total_sample_count == 0 )
{
	print "No new samples to be added to multisample calling $multisettings."; 
	
	if ( !  $dataFreeze )
	{
	 	print " Bye!\n";
		exit 0;
	}
	else
	{
		print " Checking if old batches can be merged for datafreeze.\n";
	}
}


# Create pipeline cfg file
my $pipe_name        = "MultiSample_".$timestamp_file;
my $pipe_config_file = $pipelines_folder."/".$timestamp_file."_pipeline.cfg";

# Use pipeconfig utility
my $pipecfg = pipeconfig->new($pipe_name, $pipe_config_file, "Multisample Calling Pipeline: $pipe_name ");

	# Header
	$pipecfg->header(
		{ "start_time" 			=> $timestamp_mysql 		},
		{ "timestamp"  			=> $timestamp_file 			},
		{ "total_batches" 		=> $batch_count 			},
		{ "total_samples"   	=> $total_sample_count 		},
		{ "discarded_samples"	=> $discarded_sample_count 	},
		{ "multisettings"		=> $multisettings			},
		{ "settings"			=> $settings				},{}
	);


# Newly generated Batches (names and gvcfs)
my @new_batches_names;
my @new_batches_gvcfs; 

# For each batch schedule the batch processing:
# Batches
for ( my $batch_id=1; $batch_id<=$batch_count; $batch_id++ )
{
	# Define batch ident:
	#  Date, time, settings, batchid  20170623_091201_genomems_0001 (sprintf pad 4 dig year, 2 dig month, day hour min sec, 3 dig batchid )
	# OR GET BATCH COUNT FOR SETTINGS FROM DB and ADD to BATCH_ID into making name!
	
	# indexes for the array $batches
		# $batch_id, batch number
	my $batch_name = $timestamp_file."_".$multisettings."_".( sprintf("%04d", $batch_id) );
	
	# Make batch dir
	my $current_batch_folder = $batches_folder."/".$batch_name;
	if ( ! -d $current_batch_folder )
	{
		make_path( $current_batch_folder, { mode => 0775 }) or die $logger->error("Cannot create folder $current_batch_folder");
	}
	else
	{
		die $logger->error("Batch folder $current_batch_folder already exists. Something is wrong with that");	
	}	
	
	# Open batch file:
	my $batch_list = $current_batch_folder."/".$batch_name.".list"; 
	open(BATCH_LIST, ">$batch_list");
	
	# Samples each batch has max_batch samples except the last one 
	for ( my $sample_id=1; $sample_id<=( ($batch_id == $batch_count) ? ($total_sample_count - ($max_batch*($batch_count-1))) : $max_batch ); $sample_id++ )
	{
		# indexes for the array $batches
		# 	$sample_id, internal sample count into batch
		
		# Give sample name not GVCF name for combineGVCFs.pl
		my $current_sample = $samples[($batch_id-1)*$max_batch + $sample_id];
	
		# Write entry
		print BATCH_LIST $current_sample."\n";	

	}
	
	# Close batch file
	close(BATCH_LIST);
	
	# Create batch submission entry
	# Pass batch ident, batch settings
	
	# Batch gvcf:
	my $batch_gvcf = $current_batch_folder."/".$batch_name.".gvcf.gz"; 
	
	# New batch names and gvcfs
	push(@new_batches_names, $batch_name);
	push(@new_batches_gvcfs, $batch_gvcf);
	
	# Add program
	$pipecfg->program(
				"multisample_batch.pl",
				"runs combineGVCFs.pl on the batch and writes down the batch ",
				1, 
				"",
				( $perBEDArray ? "perBEDArray" : 8 ),
				{ "samples" 		=> "$batch_list" },
				{ "o" 				=> "$batch_gvcf" },
				( $perBEDArray ? { "ajb" => ""} : {"#" => "#" } ), 
				{ "name" 	=> "$batch_name" },
				{ "multisettings" 	=> "$multisettings" },  
				{ "date" 			=> "\"$timestamp_mysql\"" },
				{ "batch" => ""},{}
			);
}


# If datafreeze
if ( $dataFreeze )
{
	
	# Get last datafreeze
	my $df_sql = qq{select idmsdatafreeze, name, date, gvcf_path from $multisampledb.$datafreezetable where idmssettings=? and errors=0 order by date desc limit 1};
	my $df_sth = $coredbh->prepare($settings_sql);
	$df_sth->execute($idmultisettings);
	
	# Get id, date, gvcf path
	my ($idmsdatafreeze, $old_datafreeze_name, $old_datafreeze_date, $old_datafreeze_gvcf)=$settings_sth->fetchrow_array();
	
	# If old datafreeze batch must be selected for merging only if newer than last datafreeze
	my $datafreeze_sql = "";
	# If older datafreeze exist then
	if ( defined $idmsdatafreeze )
	{
		$datafreeze_sql = " and date > \"$old_datafreeze_date\" ";
	}

	# Store here all batches > date + current
	my $batches_to_merge = 0;
	my @merge_batches_names;
	my @merge_batches_gvcfs;
	
	# Get all batches > date(datafreeze) and paths but not current
	# batches with errors are not merged, their samples are also excluded on the beginning
	my @old_batches_names;
	my @old_batches_gvcfs;																					#datafreeze_sql > date
	my $batches_sql = qq{select idmsbatch, name, date, gvcf_path, list from $multisampledb.$batchtable where idmssettings=? $datafreeze_sql and errors=0 order by date asc};
	my $batches_sth = $coredbh->prepare($batches_sql) or die $logger->error("SQL Error");
	$batches_sth->execute($idmultisettings) or die $logger->error("SQL Error");

	# Foreach batch
	while ( my ($old_idmsbatch, $old_batch_name, $old_batch_date, $old_batch_gvcf, $old_batch_list) = $batches_sth->fetchrow_array() )
	{
		print "Batch: $old_batch_name\n";
		
		# Check if GVCF exists if date <> today's batch $timestamp_mysql
		# If it does not exist - rerun batch and log error
		# Old batches must be in place already
		if ( ! &gvcfexists($old_batch_gvcf) )
		{
			$pipecfg->program(
				"multisample_batch.pl",
				"runs combineGVCFs.pl on the batch and writes down the batch ",
				1, 
				"",
				( $perBEDArray ? "perBEDArray" : 20 ),
				{ "samples" 				=> "$old_batch_list" },
				{ "o" 					 	=> "$old_batch_gvcf" },
				( $perBEDArray ? 	{ "ajb" => ""} 
								: 	{"#" 	=> "#" } ), 
				{ "name" 					=> "$old_batch_name" },
				{ "multisettings" 			=> "$multisettings" },  
				{ "date" 					=> "\"$old_batch_date\"" },
				{ "batch" 					=> ""},{}
			);
		}
		
		# Push
		push(@old_batches_names, $old_batch_name);
		push(@old_batches_gvcfs, $old_batch_gvcf);
		
	}
	
	# Add LATEST datafreeze if present
	if ( defined $idmsdatafreeze ) 
	{
		push(@merge_batches_names, "$old_datafreeze_name");
		push(@merge_batches_names, "$old_datafreeze_gvcf");
	}
	
	# Merge names (old and new)
	foreach my $bname ( @new_batches_names, @old_batches_names )
	{
		push(@merge_batches_names, $bname);
	}
	
	# Merge gvcfs (old and new)
	foreach my $bgvcf ( @new_batches_gvcfs, @old_batches_gvcfs )
	{
		push(@merge_batches_gvcfs, $bgvcf);
		$batches_to_merge++ ; 
	}
	
	# If nothing to do exit
	if ( $batches_to_merge ==  0 )
	{
		$logger->info("No batches to merge - no need to run me. Next time run me with -checkonly option to see if there's something to do");
		exit;
	}
	
	
	# Datafreeze creation
	
	# Datafreeze name
	my $new_datafreeze_name = $timestamp_file."_".$multisettings;
	
	# Define datafreeze folder
	my $new_datafreeze_folder = $datafreeze_folder."/".$new_datafreeze_name;
	
	# Datafreeze GVCF
	my $new_datafreeze_gvcf = $new_datafreeze_folder."/".$new_datafreeze_name.".gvcf.gz";
	
	# Datafreeze VCF
	my $new_datafreeze_vcf = $new_datafreeze_folder."/".$new_datafreeze_name.".vcf.gz";
	
	# Make datafreeze Folder 
	if ( ! -d $new_datafreeze_folder )
	{
		make_path( $new_datafreeze_folder, { mode => 0775 }) or die $logger->error("Cannot create folder $new_datafreeze_folder");
	}
	
	# Make datafreeze gvcf merge file	
	my $new_datafreeze_gvcf_list = $new_datafreeze_folder."/".$new_datafreeze_name.".list";
		
	# Create merge file
	open(FREEZE_LIST, ">$new_datafreeze_gvcf_list");
	foreach my $gvcf ( @merge_batches_gvcfs )
	{
		print FREEZE_LIST $gvcf."\n";
	}
	close(FREEZE_LIST);
	
	# DATAFREEZE!
	# Add program multisample_datafreeze.pl 
	# -> it is a link to multisample_batch.pl 
	# it has to wait multisample_batch.pl and the different name
	# is meant for when we have BEDArrayJob otherwise every job waits for each other
	$pipecfg->program(
			"multisample_datafreeze.pl",
			"runs combineGVCFs.pl between last datafreeze and recent batches ",
			1, 
			"multisample_batch.pl",
			( $perBEDArray ? "perBEDArray" : 1 ),
			{ "gvcfs" => "$new_datafreeze_gvcf_list" },
			{ "o" => "$new_datafreeze_gvcf" },
			( $perBEDArray ? { "ajb" => ""} : {"#" => "#" } ), 
			{ "name" => "$new_datafreeze_name" },
			{ "multisettings" => "$multisettings" },  
			{ "date" => "\"$timestamp_mysql\"" },
			{ "datafreeze" => ""},{}
		);
	
	
	# Concatenate GVCF/VCF CHUNKS if necessary
	# If BEDARRAY I must combine the VCFs [ and opt the GVCFs ]
	# GIVE THEM THE FINAL NAME
	if ( $perBEDArray )
	{
		
		#concatVCF.pl -i . -o concat.vcf.gz -e gvcf.gz -se hg19_wholegenome -n
		
		# Optional Combine GVCF (why doing that)
		#TODO: ADD RUNS AS OPTION
		$pipecfg->program(
			"concatVCF.pl",
			"concats VCF for BEDArray",
			1, 
			"multisample_datafreeze.pl",
			1,
			{ "i" => "$new_datafreeze_folder" },
			{ "o" => "$new_datafreeze_gvcf" }, 
			{ "e" => "gvcf.gz" },
			{ "se"=> "$settings" },
			{ "n" => "" },{}
		);
		
		# Combine called VCFs
		$pipecfg->program(
			"concatVCF.pl",
			"concats VCF for BEDArray",
			1, 
			"multisample_datafreeze.pl",
			1,
			{ "i" => "$new_datafreeze_folder" },
			{ "o" => "$new_datafreeze_vcf" }, 
			{ "e" => "vcf.gz" },
			{ "se"=> "$settings" },
			{ "n" => "" },{}
		);
	}

	# Filter VCF
	
	
	# Genotype phasing

}


# Write pipeline configuration
$pipecfg->write() or die $logger->error("Cannot create pipeline config file $pipe_config_file");
    

# Start parallelpipeline from settings (in pipeline folder)

my $pipe_command="$prog_path/parallelpipeline.pl -i $pipe_config_file -sge custom.q  -pri \"-500\" -pe pipeline -otherpipe";
  
if ( $noStart )
{
	print "********************WARNING********************\n";
	print "You chose not to launch the pipeline after data\n";
	print "collection and batching.  To launch it manually\n";
	print "use the command below.\n";
	print "***********************************************\n";
	print "$pipe_command\n";
	print "**********************END**********************\n";
	
}
else
{
	#START
	system($pipe_command);
	print "[JOBS LAUNCHED] Ya know that they'll fail. ;) Good Luck!";
}

 
# Update SNV PASS status from the datafreeze

 
 
 

 
exit(0);

# Test if a GVCF exists
# 
# if bedArray tests using the chromosome_split prefixes
#
# if normal tests the file 
#
sub gvcfexists()
{
	my $gvcf = shift;
	my $exists=1;
	
	if ( $isBedArray )
	{
		# Defined up
		for (my $bedArrayItem=1; $bedArrayItem<=(scalar @bedArrayItems); $bedArrayItem++ )
		{
			my $item = dirname($gvcf)."/".$bedArrayItem.".".basename($gvcf);
			$exists = $exists && -f $item;
		}
	}
	else
	{
		$exists = -f $gvcf;
	}
	
	return $exists;
}

 

=head1 NAME

startMultiSampleCalling.pl

=head1 SYNOPSIS

startMultiSampleCalling.pl -mse hg19WGS [ -checkOnly ] [ -dodatafreeze ] 

=head1 DESCRIPTION

This scripts creates batches for multisample calling and starts a data freeze according to the multisample settings specified as input

=head1 OPTIONS

 -mse	       multisample sample selection settings 
 -datafreeze   run the multisample calling on the batches
 -checkonly    check status  
 -forcelibtype force a library type ( for old samples that have genome and exome with the same sample ID)
 -lf	       log file; default: print to screen
 -ll	       log level: ERROR,INFO,DEBUG; default: INFO
 -h	           print this helptext
 -man	       show man page

=head1 AUTHOR

Riccardo Berutti

=cut



# TABLES:
# create table `mssettings` ( `idmssettings` int(11) unsigned NOT NULL AUTO_INCREMENT, `name` varchar(30) not null, `selectquery` blob not null, `settings` varchar(30) not null, `bedarray` tinyint(1) unsigned not null, `description` blob, primary key (`idmssettings`), unique key (`name`) ) ENGINE=InnoDB DEFAULT CHARSET=latin1;
# create table `mssample2batch` ( `idmssample2batch` int(11) unsigned NOT NULL AUTO_INCREMENT, `idsample` int(11) unsigned not null, `idmsbatch` int(11) unsigned not null, primary key (`idmssample2batch`), unique key `sample2batch` (`idsample`,`idmsbatch`) ) ENGINE=InnoDB DEFAULT CHARSET=latin1;

# create table `msbatch` ( `idmsbatch` int(11) unsigned NOT NULL AUTO_INCREMENT, `name` varchar(45) not null, `date` timestamp not null, `idmssettings` int(11) not null, `gvcf_path` varchar(255) not null, `errors` int(11) not null, `list` varchar(255) not null, primary key (`idmsbatch`), unique key `namemssettings` (`name`, `idmssettings`) ) ENGINE=InnoDB DEFAULT CHARSET=latin1;
# create table `msdatafreeze` ( `idmsdatafreeze` int(11) unsigned NOT NULL AUTO_INCREMENT, `name` varchar(45) not null, `date` timestamp not null, `idmssettings` int(11) not null, `gvcf_path` varchar(255) not null, `errors` int(11) not null, `list` varchar(255) not null, primary key (`idmsdatafreeze`), unique key `namemssettings` (`name`, `idmssettings`) ) ENGINE=InnoDB DEFAULT CHARSET=latin1;
