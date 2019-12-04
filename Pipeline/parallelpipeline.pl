#!/usr/bin/perl

use strict;
#use warnings;
#use diagnostics;
use Getopt::Long;
use Log::Log4perl qw(get_logger :levels);
umask(022);
use DateTime;
use DBI;
use Text::ParseWords;

use File::Basename;
use Cwd qw(abs_path);

use File::HomeDir;
my $homedir = File::HomeDir->my_home;

my $pgrpath = dirname(abs_path($0));
require $pgrpath."/Utilities.pm";


# Tim M Strom May 2010
# Script to automate next generation sequencing analysis
# Sebastian Eck Jun 2010
# added:
# -log file and error handling
# -optional program list
# Thoams Wieland July 2011
# submitting every script as a single SGE job including dependencies


my $line     = "";
#my $pgrpath = "/usr/local/packages/tools/solexa";
my $type     = "";
my $item1    = "";
my $item2    = "";
my $pgr     = "";
my $param    = "";
my @lane     = ();
my $help     = 0;
my $cfgfile  = "pipeline.cfg";
my $logfile  = "";
my $loglevel = "INFO";
my $programs = "";
#my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $exec = 1; #execute flag
my $sgequeue    = "";
my $debug	    = 0;
my $sleep	    = 0;
my $skipSleep   = 0;
my $firstTime   = "";
my $otherpipe   = 0;
my $clearOutput = 0;
my $finishedFile = "";
my $parallelEnv  = "pipeline";
my $priority     = "";
my $dependsOnSample = "";

GetOptions(
"i=s" => \$cfgfile, 
"p=s" => \$programs,
"pri=s" => \$priority,
"sge=s" => \$sgequeue,
"pe=s"  => \$parallelEnv,
"f=s" => \$finishedFile,
"d"     => \$debug,
"sl=s" => \$sleep,
"n=s"  => \$skipSleep,
"st=s" => \$firstTime,
"co"   => \$clearOutput,
"otherpipe" => \$otherpipe,
"do=s" => \$dependsOnSample,
"lf=s" => \$logfile,
"ll=s" => \$loglevel,
"h" => \$help
);
if ($help == 1) {
print 
"-i	<cfgfile>
-p	<programs> executes only specified programs
-f	<finishedFile.txt> file to write information about finished pipeline; DEPRECATED: information on starting & stopping is now inserted into the database
-sge	<sge_queue> submit pipeline job to specified SGE queue
-pri	priority of the new jobs. Overwrites standard values (-10,-100). Only negative values are possible. NOTE: negative values must be in apostrophes.
-pe		<parallel_environment> use specified SGE parallel environment to submit multithreaded jobs
-st	time at which the job(s) should start; format: DD.MM.YYYY,hh:mm; default: now
-sl	<minutes_to_sleep> time the script waits until it submits the next SGE job, default: 0; note that the sleep and start times only apply to jobs without dependencies
-n	number of jobs that should be started at the beginning without sleeping (i.e. there is no sleeping time between them)
-co	remove all files from OUT folder; useful if pipeline was run before and files should be cleaned
-d	debug; don't start jobs, just print the qsub call
-otherpipe this is not a sample alignment pipeline but another generic pipeline which uses the parallelpipeline launching script (no pipeline on db, no lastjob)
-lf	log file; default: pipeline.log
-ll	log level: ERROR,INFO,DEBUG; default: INFO
-h	this help\n";
exit(1);
}

if($logfile eq ""){
	$logfile = dirname($cfgfile)."/pipeline.log";
}

Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();

my $params     = Utilities::getParams();


my @pgrList = ();
my %pgrHash = ();
if($programs ne "") {
	$exec = 0;
	@pgrList = split(" ",$programs);
	foreach (@pgrList) {
		$pgrHash{$_} = 1;
		$logger->debug("$_");
	}
}

my @removes;
my %dependson;
my %runningJobs;
my $sample;
my $jobCount       = 0;
my $outpath;
my $addedSleepTime = 0;
my $libtype = "";
my $libpair = "";
my $allJobNames = "";
my $project     = "";
my $settings    = "";
my $firstRun    = 1;
my $runscripts  = "";
my $slots       = "1";

my $pipename    = "";



open(IN, "$cfgfile") || exit $logger->error("Cannot open configfile $cfgfile");
#my $finishedfile = "";

$logger->info("########################################################");
$logger->info("pipeline job submission started!");
$logger->info("########################################################");

while (<IN>) {
	$_ =~ s/^\s+//;
	if (/^#/) {next;}
	chomp;
	s/\s+//g;
	$line = $_;
	#($type,$item1,$item2)=split(/\:/,$line);
	($type,$item1,$item2)=quotewords(":",0,$line);
	
	# Parsing cfg file items:
	if ($type eq "outdir"){
		$outpath = $item1;
	}
	elsif ($type eq "libtype"){
		$libtype = $item1;
		$logger->info("LIBTYPE: $libtype");
	}
	elsif ($type eq "libpair"){
		$libpair = $item1;
		#$finishedfile =dirname($cfgfile)."/$sample.$libtype.$libpair.finished";		#delete old finished file
		#unlink($finishedfile);
		$logger->info("LIBPAIR: $libpair");
	}
	elsif ($type eq "sample"){
		$sample = $item1;
		$logger->info("SAMPLE: $sample");
	}elsif ($type eq "settings"){
		$settings = $item1;
		$logger->info("SETTINGS: $settings");
	}
	elsif ($type eq "project"){
		$project = $item1;
	}
	elsif ($type eq "pgr") {
		# Reset the param string
		$pgr  = "";
		$param = "-lf $logfile -ll $loglevel ";
		@lane  = ();
		$pgr  = $item1;
		if ($programs ne "" and exists($pgrHash{$pgr})) {
			$exec = 1;
		}
		elsif ($programs ne "" and !exists($pgrHash{$pgr})) {
			$exec = 0;
		}#dirname($cfgfile)
	}
	elsif ($type eq "param") {
		$param .= " -" . $item1 . " " . $item2;
	}
	elsif ($type eq "lane") {
		push(@lane,$item1);
	}
	elsif ($type eq "dependson"){
		$dependson{$pgr} = $item1;
	}
	elsif($type eq "slots"){
		$slots = $item1;
	}
	elsif($type eq "pipename"){
		# Then it's not a sample analysis but a generic pipeline
		$pipename = $item1;
		$otherpipe=1;
	}
	#programs without lanes
	elsif ($type eq "run" and $exec) {
		#$logger->info("START: pipeline $pgr $param");
		#system("$pgrpath/$pgr $param");
		
		if($firstRun){
			&getRunningJobs();
			$firstRun = 0;
		}
		
		my $jobName = "";
		
		if($runningJobs{$pgr}){		#if there is already a job with this name running --> add this job into the list
			my $currentRunning = $runningJobs{$pgr};
			my @columns = split(",",$currentRunning);
			$jobName = $columns[-1];
			$jobName =~ /(\d+)$/;
			my $newCount = $1+1;
			$jobName =~ s/(\d+)$/$newCount/;	#substitute last number for higher number
			$currentRunning .= ",".$jobName;
			$runningJobs{$pgr} = $currentRunning;
		
		}else{						#construct jobname from scratch
			my $subPgr = $pgr;
			$subPgr =~ s/\./_/g;
			#$subPgr =~ s/_//g;
			$jobName = namePrefix()."_".$subPgr."1";
			$runningJobs{$pgr} = $jobName;
		}
		
		my $dependString = "";					#construct the string of jobs this job depends on
		if($dependson{$pgr}){
			my @dependencies = split(",",$dependson{$pgr});
			foreach(@dependencies){
				if($runningJobs{$_}){
					if($dependString eq ""){
						$dependString = $runningJobs{$_};
					}else{
						$dependString .= ",".$runningJobs{$_};
					}
				}
			}
		}
		#if sample depends on other sample:
		if ($dependsOnSample) {
			if($dependString ne ""){
				$dependString .= ",";
			}
			$dependString .= "'*_" . $dependsOnSample . "_*'";
		}
		
		#clear output dir if requested
		if($clearOutput){
			if($debug){
				#print "outpath: $outpath\n";
			}else{
				system("find $outpath -type f -delete");
			}
		}
		
		#generate qsub string
		my $qsubstring = "qsub -q $sgequeue -N $jobName ";
		$allJobNames .= $jobName.",";
		
		
		
		if($slots eq "perChrArray"&& $params->{settings}->{$settings}->{normalchromosomes}){				#start as an array job with one job per chromosome
			my $normalChroms = $params->{settings}->{$settings}->{normalchromosomes}; #get number of chromosomes, i.e. number of jobs in array, from dict file
			my $chromosomes = `cat $normalChroms | wc -l`;
			chomp $chromosomes;
			$qsubstring .= "-t 1-$chromosomes ";
			
			
		}elsif($slots eq "perBEDArray" && $params->{settings}->{$settings}->{genome_splits}){ #start as an array job with one job per bed file in a folder
			my $dir = $params->{settings}->{$settings}->{genome_splits};
			my @beds = glob("$dir/*.bed");
			$qsubstring .= "-t 1-".@beds." ";
			
		}elsif($parallelEnv ne "" && (length($slots)>1 || $slots>1)){ #add parallel environment
			$qsubstring .= "-pe $parallelEnv $slots ";
		}
		
		#set error/out path
		#$qsubstring .= "-e ".dirname($cfgfile)."/GE.error -o ".dirname($cfgfile)."/GE.out ";
		
		if($dependString eq ""){				#if a job has no dependencies the start delay applies
			
			
			
				if ( $jobCount >= $skipSleep ) { #start the first "skipSleep" jobs right now
					$addedSleepTime +=
					  $sleep;    #add "sleep" minutes to the starttime of all other jobs
				}
			
				my $starttime = &getStartTime($addedSleepTime);
				
				
				$qsubstring .= "-a $starttime ";
				
				if($priority eq ""){
					$qsubstring .= "-p -100 ";
				}else{
					$qsubstring .= "-p $priority ";
				}
				
				
			
				$jobCount++;
		}else{
			$qsubstring .= "-hold_jid $dependString ";
			if($priority eq ""){
				$qsubstring .= "-p -10 ";					# --> priority of jobs with dependency is higher than for other jobs --> a pipeline should finish first before another pipeline starts
			}else{
				$qsubstring .= "-p $priority ";
			}
			
		}
		
		#create temp script to submit
		open (OUT, ">$homedir/$jobName.sgeSubmTmp.sh") || exit $logger->error("Cannot open $homedir/$jobName.sgeSubmTmp.sh");
		print OUT "#\$  -S /bin/bash\n";
		print OUT "perl $pgrpath/$pgr $param";
		close OUT;
		
		$qsubstring .= "$homedir/$jobName.sgeSubmTmp.sh";
		
		if($debug){
			print $qsubstring."\n";
		}else{
			$logger->info("SUBMITTING: $qsubstring: $pgrpath/$pgr $param");
			
			system($qsubstring);
			#add name of program to list of run scripts
			$runscripts .= "$pgr:";
		}
		system("rm $homedir/$jobName.sgeSubmTmp.sh"); 	#remove tmp script

		
		#$logger->info("END: pipeline $pgr $param");
	}
	elsif ($type eq "rm"){
		push(@removes,$item1);
	}
}

########################
#post into finished file that a new job has been submitted (to track jobs that are started but never finished)
#open  OUT, ">>$finishedFile" || exit $logger->error("Cannot open $finishedFile for appending!");
#my @timeData = localtime(time);
#print OUT "START\t$sample\t$project\t$libtype\t$libpair\tNA\t$logfile\t".($timeData[5]+1900)."/".($timeData[4]+1)."/$timeData[3] $timeData[2]:$timeData[1]:$timeData[0]\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
#close OUT;


if($allJobNames eq ""){
	$logger->info("No jobs to submit for ". ( $otherpipe == 0 ? "sample" : "pipeline") ." $sample! - exiting");
	exit(0);
}



# Last steps will run ONLY if SAMPLE PIPELINE 
if ( $otherpipe == 0  )
{
	
	#######################
	# Write Pipeline information into the database IF Sample Pipeline	
	#UPDATE: TW, 03.03.2014 --> information on start and stop of pipeline is now also stored in the database
	my $dbh         = Utilities::connectCoreDB();
	my $params      = Utilities::getParams();
	my $solexadb    = $params->{solexadb}->{database};
	my $sampleTable = $params->{coredb}->{sampletable};
	
	$runscripts =~ s/:$//;
	
	
	
	#get GIT version of the pipeline
	chdir($pgrpath);
	open GIT,"git log | head -n 1 | sed 's/commit //'|";
	my $gitversion = <GIT>;
	chomp $gitversion;
	close GIT;
	
	
	#insert new entry into pipeline file
	my $sql = "insert into pipeline (idsample,idlibtype,idlibpair,logfile,pipelinegitversion,currentsettings,executedscripts, outputdir)
	values ((select idsample from $sampleTable where name='$sample'),(select ltid from $solexadb.libtype where ltlibtype='$libtype'),(select lpid from $solexadb.libpair where lplibpair='$libpair'),'$logfile','$gitversion','$settings','$runscripts', '$outpath');";
	
	$dbh->do($sql);
	my $idpipeline = $dbh->last_insert_id(undef, undef, qw(pipeline idpipeline)) or die $logger->error($DBI::errstr);


	########################
	# LastJob is submitted for sample processing pipeline
	#submit lastJob.pl for cleanup purpose

	$allJobNames =~ s/,$//g;
	my $jobName = namePrefix()."_lastJob_pl";
	my $qsubstring = "qsub -q $sgequeue -N $jobName -hold_jid $allJobNames -p -1 ";
	#create temp script to submit
	open (OUT, ">$homedir/$jobName.sgeSubmTmp.sh") || exit $logger->error("Cannot open $homedir/$jobName.sgeSubmTmp.sh");
	print OUT "#\$  -S /bin/bash\n";
	print OUT "perl $pgrpath/lastJob.pl -ip ".$idpipeline." -s $sample -se $settings -lt $libtype -lp $libpair -p $project -lf $logfile -ll $loglevel";
	close OUT;
	
	$qsubstring .= "$homedir/$jobName.sgeSubmTmp.sh";
	
	if($debug){
		print $qsubstring."\n";
	}else{
		$logger->info("SUBMITTING LAST JOB: $qsubstring");
		#print "\n\nSUBMITTING: $qsubstring\n";
		
		system($qsubstring);
		#system("sh $pgrpath/sgeSubmTmp.sh");
	}
	system("rm $homedir/$jobName.sgeSubmTmp.sh"); 	#remove tmp script
}


#################################################################################################
# Functions
#################################################################################################
sub getRunningJobs {
	#initialise the list of running jobs to (possibly) wait for with jobs for this sample/library type already running on the SGE
	
	# Il Libtype and LibPair are defined build prefix with them otherwise use pipename
	my $prefix=namePrefix();
	
	open QSTAT, "qstat -r | grep \"Full jobname\" | grep $prefix |";
	while(<QSTAT>){
		chomp;
		$_ =~ s/^\s+//;
		my ($dummy,$dummy2,$jobName) = split(' ',$_);
		
		#build script name from jobname
		if($jobName =~/$prefix\_(\w*)\d+$/){
			my @columns = split("_",$1);
			my $progName = join("_",@columns[0..(@columns-2)]).".".$columns[-1];
			
			if($runningJobs{$progName}){
				$runningJobs{$progName} = $runningJobs{$progName}.",".$jobName;
			}else{
				$runningJobs{$progName} = $jobName;
			}
			
		}
	}
	
}



#################################################################################################
sub getStartTime {
	my $minOffset = shift;

	my $dt;
	my $tz = DateTime::TimeZone->new( name => 'local' );
	if ( $firstTime eq "" ) {
		$dt = DateTime->now;
		$dt->set_time_zone($tz);
	}
	else {
		my ( $date, $time ) = split( ",", $firstTime );
		my ( $day, $month, $year ) = split( /\./, $date );
		my ( $hour, $minute ) = split( ":", $time );
		$dt = DateTime->new(
			year      => $year,
			month     => $month,
			day       => $day,
			hour      => $hour,
			minute    => $minute,
			time_zone => $tz
		);
	}

	$dt->add( minutes => $minOffset );

	my $hour   = sprintf( "%02d", $dt->hour );
	my $minute = sprintf( "%02d", $dt->minute );
	my $date   = $dt->ymd("");

	return "$date$hour$minute";

}

#################################################################################################
# Build prefix for pipeline name (was repeated in the code before)
sub namePrefix {
	
	my $libtype_in  = shift;
	my $libpair_in  = shift;
	my $sample_in   = shift;
	my $pipename_in = shift;
	my $outPrefix   = "";
	
	$libtype_in  = $libtype  if $libtype_in  eq "";
	$libpair_in  = $libpair  if $libpair_in  eq "";
	$sample_in   = $sample   if $sample_in   eq "";
	$pipename_in = $pipename if $pipename_in eq "";
	
	# Il Libtype and LibPair are defined build prefix with them otherwise use pipename
	if( $libtype_in ne "" && $libpair_in ne "" )
	{
		$outPrefix = (uc substr($libtype_in,0,1)).(uc substr($libpair_in,0,1))."_".$sample_in;
	}
	else
	{
		$outPrefix = "PL_".$pipename_in;
	}
	
	
	return $outPrefix;
}



#remove files
#foreach(@removes){
#	$logger->info("Removing file(s) $_");
#	system("rm $_");
#}

#$logger->info("########################################################");
#$logger->info("pipeline finished!");
#$logger->info("########################################################");


#open(OUT, ">$finishedfile") || exit $logger->error("Cannot open $finishedfile");
#my @timeData = localtime(time);
#my $month = $timeData[4]+1;
#my $year  = $timeData[5]+1900;
#print OUT "Pipeline finished at $timeData[2]:$timeData[1]:$timeData[0] on $timeData[3].$month.$year\n";
#my $errors = qx/grep ERROR $logfile | wc -l/;
#chomp $errors;
#print OUT "Logfile contains $errors errors! \n";
#close OUT;

