#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use POSIX;

use File::Basename;
use File::chmod::Recursive;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);
use DBI;

my $prog_path = dirname(abs_path($0));
require $prog_path."/Utilities.pm";


my $logfile  = "pipeline.log";
my $loglevel = "INFO";
my $idpipeline = -1;
my $sample = "";
my $project = "";
my $libtype = "";
my $libpair = "";
my $settings = "";
my $outDir = "";

#program flags
my $help       = 0;

GetOptions(
"ip=s" => \$idpipeline,
"s=s"  => \$sample,
"p=s"  => \$project,
"se=s" => \$settings,
"lt=s" => \$libtype,
"lp=s" => \$libpair,
"lf=s" => \$logfile,
"ll=s" => \$loglevel, 
"h" => \$help);

Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();


if ($help == 1) {
print 
"
This script is meant to run as last job in a batch of SGE jobs. Its main purpose is handling of a batchs status.
It writes some basic information about time, errors, ... into the database table \"pipeline\"

-ip	idpipeline; --> entry in database table pipeline that should be updated
-p	<project name>
-s	<sample name>
-lt	<library type>
-lp <library pair>
-se <settings>
-lf	<log file>; default: pipeline.log
-ll	<log level>: ERROR,INFO,DEBUG; default: INFO
-h	this help\n\n";
exit(1);
}


#parse log file for times and errors
open LOG, $logfile or exit $logger->error("Can't open $logfile!");
my $numberOfErrors = 0;
my $submTime       = "0";
my $firstTime      = "0";
my $lastTime	   = "0";
while(<LOG>){	#read in whole log file but throw away everything of the older runs
	chomp;
	my ($date,$time,$type,$script,@rest) = split(' ');
	$lastTime = $date." ".$time if $date && $time;
	if($_ && $_ =~ /pipeline job submission started/){
		$numberOfErrors = 0;
		$firstTime      = "";
		$submTime       = $lastTime;
	}
	if($firstTime eq "" && !($script =~ /parallelpipeline/)){		#firstTime is the time where the first script actually started running, not the submission time
		$firstTime = $lastTime;
	}
	if($type && $type =~ /ERROR/){
		$numberOfErrors++;
	}
}

#Query database
my $params = Utilities::getParams();
my $database       = $params->{coredb}->{database};
my $statTable      = $params->{coredb}->{stattable};
my $sampleTable    = $params->{coredb}->{sampletable};
my $solexadb       = $params->{solexadb}->{database};

my $exomedb        = $params->{settings}->{$settings}->{exomedb}->{database};
my $snvsampleTable = $params->{settings}->{$settings}->{exomedb}->{snvsampletable};

my $dbh = Utilities::connectCoreDB();


#udpate pipeline information in database
#first: get versions of programs
if($idpipeline != -1){
	my $versionstring = "";
	foreach my $program(keys %{$params->{programs}}){
		if($params->{programs}->{$program}->{getversion}){
			open IN,"$params->{programs}->{$program}->{getversion}|" or next;
			my $version = <IN>;
			chomp $version;
			$versionstring .= "$program: $version\n";
			
		}
	}
	
	my $varquery = "";
	my $importquery = "";
	my $sql;
	my $out;
	
	#getoutputpath
	$sql = qq{select outputdir, DATE_FORMAT(starttime, '%Y%m%d%H%i%S') from pipeline where idpipeline=$idpipeline;};
	$logger->debug($sql);
	$out = $dbh->prepare($sql) || exit $logger->error($DBI::errstr);
	$out->execute() || exit $logger->error($DBI::errstr);
	my ($outputDirectory, $starttime) = $out->fetchrow_array();
	$outDir = $outputDirectory;
	
	
	# TODO: IMPORTANT: DO NOT HARDCODE HERE -> CODE IT IN SETTINGS -> ASYNCHRONOUS IMPORT -> SOMETHING LIKE:
	# if ( $params->{settings}->{$settings}->{exomedb}->{asyncimport} == "true"  )
	if ($settings eq "hg19_wholegenome" || $settings eq "hg19_genomedeepvariant") {
				
		#Merge snvsample.$sample.* files
		opendir(DIR, $outputDirectory);
		my $counter=0;
		my $snvsampleImportFile = "snvsample.merged.tsv";
		my $command = "";
		#remove $snvsampleImportFile if existing before creating it
		if (-e "$outputDirectory/$snvsampleImportFile") {
			if (&Utilities::executeCommand("rm $outputDirectory/$snvsampleImportFile", "Removing merged snvsample file $outputDirectory/$snvsampleImportFile before (re)creating it", $logger)) {
				$logger->error("Error removing file $outputDirectory/$snvsampleImportFile")
			}
		}
		while(my $f = readdir(DIR)) {
			if ($f =~ m/^snvsample\.$sample/) {	#The $sample is here in the regex to avoid merging snvsample.merged.tsv to itself
				$counter++;
				
				# TODO: HARDCODING AGAIN - REMOVE
				if ( $settings eq "hg19_genomedeepvariant" && $f eq "snvsample.$sample.gatk.ontarget.haplotypecaller.filtered.dbSNP.plus.checked.vcf.tsv"){
					$counter--;
					next;
				}
				# TODO: HARDCODING AGAIN - REMOVE
				if ( $settings eq "hg19_wholegenome" && $f eq "snvsample.$sample.deepvariant.dbSNP.plus.checked.vcf.tsv"){
					$counter--;
					next;
				}
				
				if ($counter<=1) {
					$command = "cat $outputDirectory/$f > $outputDirectory/$snvsampleImportFile";
				} else {
					$command = "cat $outputDirectory/$f >> $outputDirectory/$snvsampleImportFile";
				}
				if (&Utilities::executeCommand($command, "Adding file $outputDirectory/$f to final SNV-import-file", $logger)) {
					$logger->error("Error adding snv file $outputDirectory/$f to the final import file $outputDirectory/$snvsampleImportFile, i.e. $outputDirectory/$f will not be imported to DB!!!")
				}
			}
		}
		closedir(DIR);
		
		#TODO: check if any of the snsample.$sample files have been changed by this pipeline
		
		#if snvsample import file is found
		if (-e "$outputDirectory/$snvsampleImportFile") {
			#first of all, check date
			my $lastFileChangeDate = POSIX::strftime("%Y%m%d%H%M%S", localtime((stat("$outputDirectory/$snvsampleImportFile"))[9]));
			#check if file is empty
			my $lineCounts = `wc -l < $outputDirectory/$snvsampleImportFile`;
			if ($lastFileChangeDate > $starttime && $lineCounts > 0) { 		#check if this pipeline run changed the snvsample import file and if the file has entries; if yes, it should be importet; if not, no import
				$importquery = ",snvsampleImportFile='$outputDirectory/$snvsampleImportFile', ready2importsnvsample=1, snvsampleimportFinished=0";
			}
			
			#TODO: INVALIDATE PREVIOUS IMPORT INSTANCES FOR THE SAME SAMPLE WITH THE SAME SETTINGS
		}
	} elsif ($exomedb) {
		$varquery = ",numofvars=(SELECT count(DISTINCT x.idsnvsample) FROM $exomedb.$snvsampleTable x WHERE p.idsample=x.idsample and (caller='samtools' or caller='gatk' or caller='deepvariant' ) )
		,numofsvs=(SELECT count(DISTINCT x.idsnvsample) FROM $exomedb.$snvsampleTable x WHERE p.idsample=x.idsample and (caller='pindel') )
		,numofcnvs=(SELECT count(DISTINCT x.idsnvsample) FROM $exomedb.$snvsampleTable x WHERE p.idsample=x.idsample and (caller='exomedepth') )
	 	";
	}
	
	$sql = "update pipeline p set status='finished',endtime=NOW()$varquery $importquery,errors=$numberOfErrors,programversions='$versionstring',
	seq=(select seq from $statTable where idsample=p.idsample and idlibtype=p.idlibtype and idlibpair=p.idlibpair),
	duplicates=(select duplicates from $statTable where idsample=p.idsample and idlibtype=p.idlibtype and idlibpair=p.idlibpair),
	cov20x=(select cov20x from $statTable where idsample=p.idsample and idlibtype=p.idlibtype and idlibpair=p.idlibpair),
	exomedepthrsd=(select exomedepthrsd from $statTable where idsample=p.idsample and idlibtype=p.idlibtype and idlibpair=p.idlibpair),
	sry=(select sry from $statTable where idsample=p.idsample and idlibtype=p.idlibtype and idlibpair=p.idlibpair),
	mix=(select mix from $statTable where idsample=p.idsample and idlibtype=p.idlibtype and idlibpair=p.idlibpair)
	where p.idpipeline=$idpipeline;";
	
	$logger->debug($sql);
	$out = $dbh->prepare($sql) || exit $logger->error($DBI::errstr);
	$out->execute || exit $logger->error($DBI::errstr);
}

if ($libtype eq "RNA") {			#TODO: execute for every lib-type
	chmod_recursive(
	    {
	        dirs  => 0775,       # Mode for directories
	        files => 0664,
	    },
	    "$outDir"
	);
}

# Report pipeline finished into logfile
$logger->info("=========================");
$logger->info("Pipeline Finished!");
$logger->info("SAMPLE : $sample");
$logger->info("LIBTYPE: $libtype");
$logger->info("LIBPAIR: $libpair");
$logger->info("=========================");


exit; # Disable custom email

########################################
### Send Mail for important projects ###
########################################
use Mail::Send;

my $sql = "SELECT proj.pname, s.name 
FROM exomehg19.pipeline p
INNER JOIN exomehg19.sample s ON s.idsample=p.idsample
INNER JOIN exomehg19.project proj ON proj.idproject=p.idproject
WHERE p.idpipeline='$idpipeline'";
my $out = $dbh->prepare($sql) || die print "$DBI::errstr\n";
$out->execute || die print "$DBI::errstr\n";
my $projectName = "";
my $sampleName = "";
while ( ($projectName, $sampleName) = $out->fetchrow_array ) {
	chomp $projectName;
	chomp $sampleName;
	if($projectName eq "S0248") {		#Schunkert
		my $msg = Mail::Send->new(Subject => "Schunkert Sample Analyzed");
		$msg->to($params->{misc}->{email_ts});
		my $fh = $msg->open;
		print $fh "Sample $sampleName analyzed";
		$fh->close or die "couldn't send whole message: $!\n";
	}
}
########################################
########################################



