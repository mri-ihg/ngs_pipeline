#!/usr/bin/perl

use strict;
use Getopt::Long;
use File::Find;
use File::stat;
use DBI;
use XML::Simple;
use DateTime;
use Data::Dumper;
use File::Basename;
use Term::ANSIColor;
use Term::ANSIColor qw(:constants);
use Mail::Send;
use Cwd qw(abs_path);

# CONFIG:
my $dbhost	  = "DBHOST";
	my $dbuser= "DBUSERNAME";
	my $dbpw  = "DBPWD";

my $rundir    = "RUNDIR";
my $bamdir    = "RUNDIR/SequenceBams";
my $startFrom = &getStartTime();
my $help      = 0;
my $exclude   = "PATH/exclude_flowcells.txt";
my $forInitAnalysis = 0;
my $perSample = 0;
my $mailto    = "";
my $noOutput  = 0;
my $mailFile  = "";
my $onlySample= 0;
my $prog_path = dirname(abs_path($0));

my %seqMap;						# to map sequencer serial numbers to "lab names"
$seqMap{1}    = "Seq1";
$seqMap{2} = "Seq2";
$seqMap{3}   = "Seq3";

$seqMap{4}   = "Seq4";
$seqMap{5}   = "Seq5";

$seqMap{6}   = "Seq6";

GetOptions(
	"r=s"  => \$rundir,
	"b=s"  => \$bamdir,
	"s=s"  => \$startFrom,
	"e=s"  => \$exclude,
	"f=s"  => \$mailFile,
	"m=s"  => \$mailto,
	"n"    => \$noOutput,
	"i"	   => \$forInitAnalysis,
	"s"	   => \$perSample,
	"os"   => \$onlySample,
	"h"    => \$help
);

if($help){
	print "
This script looks for lowcells that have been run or are currently running. It only considers
new flowcells (CASAVA version >= 1.8). It shows flowcells in four categories:
1)	Sequencintg in progress (NO RTAComplete.txt)
2)	Sequencing finished, no demultiplexing (RTAComplete.txt, NO DemultiplexedBustardSummary.xml)
3)	Demultiplexing finished, no analysis (DemultiplexedBustardSummary.xml, NO Bamfiles in -b $bamdir)
4)	Analysis running. Only looks for the flowcellname in the jobs currently in a SGE queue. Therefore
	only working on hosts that have access to the Sun Grid Engine (currently only ihgseq). 
	UPDATE: Also includes information per sample, if specified

UPDATE 2013/06/19: The script now also sends e-mails if the state of a flowcell changes. To achieve this
	you have to specify a list of receiver adresses (-m) and a file in which information on already sent
	mails is stored (-f). If the script is run as a cron job you can specify the option -n to prevent 
	output on the commandline.

-r	<rundir>; default: $rundir
-b	<bamdir>; default: $bamdir
-s	<startFrom> month from which to start looking for new flowcells; format: YYMM; default: $startFrom
-e	<exclude_file> file that includes flowcellnames (one per line) that should be excluded from search;
	if \"STDIN\": take flowcellnames from STDIN; default: $exclude
-s	include run information per sample
-os	ONLY print sample information without header (to parse in checkPipelineRuns.pl)
-i	output only flowcell names to analyse --> can be piped directly into initAnalysis.pl
-m	<mailto1\@server.com,mailto2\@server.com> comma seperated list of mail addresses that should receive changes
-f	<got_mail.txt> file that stores information on mails that have been sent, such that in case of a cron job a mail gets only sent once
-n	no output on commandline
-h	print this help
";
	exit(0);
}

my @dirs = glob($rundir.'/*XX');
push(@dirs, glob($rundir.'/*XY'));
my @mailTo = split(",",$mailto);

#read in $mailFile
my %mailSent;
if($mailFile ne ""){
	open IN,$mailFile or die "Can't open $mailFile!\n";
	while(<IN>){
		chomp;
		my ($fc,$state) = split;
		$mailSent{$fc}->{$state} = 1;
	}
}

my $runningPrint = "";
my $multiplPrint = "";
my $analysePrint = "";
my $curranaPrint = "";

my %excludedFCs;
if($exclude ne ""){
	if($exclude ne "STDIN"){	
		open(IN,$exclude) or die "Can't open $exclude!\n";
		while(<IN>){
			chomp;
			$_ =~ s/\s*#.*//;
			$excludedFCs{$_} = 1;
		}
	}else{
		while(<STDIN>){
			chomp;
			$_ =~ s/\s*#.*//;
			$excludedFCs{$_} = 1;
		}
	}
}




#print output
if($onlySample){
	print &getSampleInfo();
}else{
	
	foreach my $dir (@dirs){
		my @columns = split(/\//,$dir);
		my ($date,$sequencer,$run,$flowcell) = split("_",$columns[-1]);
		
		
		my $aorb  = substr($flowcell,0,1);
		$flowcell = substr( $flowcell, ( ( length $flowcell ) - 9 ), length $flowcell );
		
		next if substr($date,0,4) < $startFrom;	#skip old flowcells	
		next if $excludedFCs{$flowcell};		#skip flowcells from list
		#print "date,sequencer,run,flowcell: $date,$sequencer,$run,$flowcell\n";
		
		#reformat date 130319,
		$date = substr($date,4,2).".".substr($date,2,2).".20".substr($date,0,2);
		
		my $sql = "SELECT distinct s.name
					FROM exomehg19.sample s
					INNER JOIN sample2library sl on sl.idsample=s.idsample
					INNER JOIN library l on sl.lid=l.lid
					INNER JOIN library2pool lp on lp.lid=l.lid
					INNER JOIN pool p on p.idpool=lp.idpool
					INNER JOIN lane la on la.idpool=p.idpool
					INNER JOIN run r on r.rid=la.rid
					WHERE r.rname='$flowcell'";
		my $dbh = &connectDB();
		my $out = $dbh->prepare($sql) || die print "$DBI::errstr\n";
		$out->execute || die print "$DBI::errstr\n";
		my $grepString = "";
		my $sampleName = "";
		while ( $sampleName = $out->fetchrow_array ) {
			chomp $sampleName;
			$grepString .= " -e $sampleName"."_lastJob_pl";
		}
		#open(IN,"qstat -r 2> /dev/null | grep $flowcell |");  #old
		open(IN,"qstat -u \"*\" -r 2> /dev/null | grep $grepString |");	#new
		
		if(<IN>){ # flowcell is currently analysed
			$curranaPrint .= "\n".$flowcell." - ".&getFlowcellInfo($flowcell);	
		}
		else{
			if(-e $dir."/RTAComplete.txt"){
				
				
				if(-e $dir."/Demultiplexed/Basecall_Stats_".$flowcell."/Demultiplex_Stats.htm"  || -e $dir."/Demultiplexed/Stats/DemultiplexingStats.xml"){		#TW 04.02.2016: new demultiplexing software version puts Stats file to different direction
				#if(-e $dir."/Demultiplexed/DemultiplexedBustardSummary.xml"){
					if(-d $bamdir."/$flowcell"){
						my $fastqfiles = 0;
						my $isPairedEnd = 0;

						foreach(glob($dir."/Demultiplexed/Project*")){
							find(sub{if(/\.fastq\.gz$/){$fastqfiles++} if(/_R2_/){$isPairedEnd=1}}, $_);	# count fastq files
						}
						my $bamfiles   = 0;
						find(sub{if(/\.bam$/){$bamfiles++}}, $bamdir."/$flowcell");				# count bam files
						
						if( ($isPairedEnd && $fastqfiles > ($bamfiles*2)) || (!$isPairedEnd && $fastqfiles > $bamfiles) ){					# if not all fastq files where converted --> analysis not finished
							if($forInitAnalysis){
								print $flowcell."\n" unless $noOutput;
							}else{
								$analysePrint .= "\n".$flowcell." - readfiles: $fastqfiles:$bamfiles"." - ".&getFlowcellInfo($flowcell);
							}
							&sendMail($flowcell,3, "", $dir);
						}
						
					}else{									# --> not yet analaysed
						if($forInitAnalysis){
							print $flowcell."\n" unless $noOutput;
						}else{
							$analysePrint .= "\n".$flowcell." - ".&getFlowcellInfo($flowcell);
						}
						&sendMail($flowcell,3, "", $dir);
					}
					
				}elsif(!$forInitAnalysis){
					$multiplPrint .= "\n".$flowcell." - ".&getFlowcellInfo($flowcell);		# --> not yet demultiplexed
					&sendMail($flowcell,2);
				}
			} elsif(!$forInitAnalysis){											# --> running
				
				#get some additional information from Status.xml file
				if(-e $dir."/Data/reports/Status.xml"){
					my $xml    = new XML::Simple;
					my $params = $xml->XMLin($dir."/Data/reports/Status.xml");
				
					$runningPrint .= "\n$flowcell running on $sequencer($seqMap{$sequencer}) slot $aorb, Run started: $date, Read started: $params->{RunStarted}, Called Cycles: $params->{CallCycle}/$params->{NumCycles} - ".&getFlowcellInfo($flowcell);
					&sendMail($flowcell,1,"$sequencer($seqMap{$sequencer})" );
				}else{
					my $xml    = new XML::Simple;
					next unless -e $dir."/RunInfo.xml";
					my $params = $xml->XMLin($dir."/RunInfo.xml");
					#print Dumper($params);
					my $runtype = "High Output";
					$runtype = "Rapid" if $params->{Run}->{FlowcellLayout}->{LaneCount} == 2;
					my $totCycles = 0;
					if (ref($params->{Run}->{Reads}->{Read}) =~ m/ARRAY/) {
						foreach my $read(@{$params->{Run}->{Reads}->{Read}}){
							$totCycles += $read->{NumCycles};
						}
					} else {
						$totCycles = $params->{Run}->{Reads}->{Read}->{NumCycles};
					}
					#get current cycle
					my @cycles = glob($dir."/Data/Intensities/BaseCalls/L001/C*");
					my $maxCycle = 0;
					foreach my $cycle(@cycles){
						$cycle = basename($cycle);
						$cycle =~ s/^C//;
						$cycle =~ s/\.1$//;
						$maxCycle = $cycle if $cycle>$maxCycle;
					}
					
					#print "test: ".$params->{Run}->{Reads}->{Read}[0]->{NumCycles}."\n";
					$runningPrint .= "\n$flowcell running on $sequencer($seqMap{$sequencer}) slot $aorb, Run type: $runtype, Run started: $date, RTA Cycles: $maxCycle/$totCycles - ".&getFlowcellInfo($flowcell);
					&sendMail($flowcell,1,"$sequencer($seqMap{$sequencer})");
				}
			}
		}
		close IN;
	}
	
	
	unless($forInitAnalysis || $noOutput){
		print BOLD WHITE ON_BLACK "\n===========================================================================================\n";
		
		if($runningPrint ne ""){
			#print color 'bold red';
			print BOLD RED ON_BLACK "\nFlowcells on sequencers:";
			print CLEAR RED ON_BLACK $runningPrint;
		}
		if($multiplPrint ne ""){
			print BOLD BLUE ON_BLACK "\n\nFlowcells not yet demultiplexed:";
			print CLEAR BLUE ON_BLACK $multiplPrint;
		}
		if($analysePrint ne ""){
			print BOLD GREEN ON_BLACK "\n\nFlowcells not yet analysed:";
			print CLEAR GREEN ON_BLACK $analysePrint;
		}
		if($curranaPrint ne ""){
			print BOLD YELLOW ON_BLACK "\n\nFlowcells currently analysed (jobs in SGE):";
			print CLEAR YELLOW ON_BLACK $curranaPrint;
		}
		if($perSample){
			#check for samples that are currently analysed
			print BOLD WHITE ON_BLACK "\n\nCurrently analysed samples:\n";
			print CLEAR WHITE ON_BLACK &getSampleInfo();
		}
	
		print BOLD WHITE ON_BLACK "\n\n===========================================================================================\n\n";
	
	
	}
}

###############################################################
sub getSampleInfo(){
	open(IN,"qstat -u \"*\" -r 2> /dev/null | ");
	
	my %samples;       
					   
	my $samplePrint = "==================================================================\nType\tSample\tProject\t\tCurrent Stage\tRunning?\n------------------------------------------------------------------\n" if $onlySample == 0;
	while(<IN>){ # flowcell is currently analysed
		chomp;
		next if !($_ =~ /Full jobname/);
		my ($dummy,$dummy2,$name) = split(' ');
		#next if !($name =~ /^\w{2}_.{5,7}_/);
		my @columns;
		my @columns_tmp = split('_',$name);
				
		my $script_items=0;
		my $script = $columns_tmp[-2].".pl";

		# RB: 20160713 patch
		#
		# Job names can conflict in parsing with sample names:
		# Underscores allowed both in sample name and script name
		# JOBTYPE_SAMPLE_NAME_SCRIPT_NAME_extension[0-9]
		# 
		# script name is parsed first and backwards until a script with such a name exists
		while(!(-e  "$prog_path/$script" ))
		{
			$script_items++;
			if ( $script_items > (scalar(@columns_tmp)-2))
			{
				last;
			}
			
			$script = $columns_tmp[-2-$script_items]."_".$script;
		}
		
		# Fill the fields of the columns as TYPE, SAMPLE, SCRIPT, EXTENSION
		$columns[0]=$columns_tmp[0];
		$columns[1]="";
		$columns[2]=$script;
		$columns[3]=$columns_tmp[-1];
		
		# Get sample name using all the fields minus the ones building the name of the script
		for(my $spname=1; $spname<(scalar(@columns_tmp)-$script_items-2); $spname++)
		{
			if($spname>1){ $columns[1].="_"; }
			$columns[1] .= $columns_tmp[$spname];
		}
				
		my $type = "";
		if($_ =~ /:\s*(.+)_$columns[1]/){
			$type = $1;
		}

		unless($samples{$type."_".$columns[1]}){
			my $next = <IN>;
			my $running = "No";
			#print "test: $next\n";
			if($next =~ /Master Queue/){
				$running = "Yes";
			}
			my $currentJob = join("_",@columns[2..($#columns -1)]);
			
			# Build a key with # character instead of "_"
			$samples{$type."#".$columns[1]} = $currentJob.":".$running;
			
		}
		
	}
	
	# Print sample information
	foreach my $key(keys %samples){
		# Split the key just built before with the hash #  
		my ($type,$sample) = split("#",$key);
		my $sql = "SELECT p.pdescription FROM exomehg19.sample s INNER JOIN exomehg19.project p ON s.idproject=p.idproject WHERE name='$sample';";
		my $dbh = &connectDB();
		my $out = $dbh->prepare($sql) || die print "$DBI::errstr\n";
		$out->execute || die print "$DBI::errstr\n";
		my $project = $out->fetchrow_array;
		if($project){
			my ($currentJob,$running) = split(":",$samples{$key});
			$samplePrint .= "$type\t$sample\t$project\t$currentJob\t$running\n";
		}
	}
	return $samplePrint;
}


###############################################################
sub getFlowcellInfo(){
	my $flowcell = shift;
	
	my $sql = "SELECT run.rcomment, group_concat(DISTINCT pool.odescription),group_concat(lane.alane), run.rfailed,lane.aread1failed,lane.aread2failed
FROM run
INNER JOIN lane on run.rid=lane.rid
INNER JOIN pool on lane.idpool=pool.idpool
WHERE run.rname='$flowcell'
GROUP BY run.rcomment, run.rfailed,lane.aread1failed,lane.aread2failed
ORDER BY lane.aread2failed";
	
	my $dbh = &connectDB();
	my $out = $dbh->prepare($sql) || die print "$DBI::errstr\n";
	$out->execute || die print "$DBI::errstr\n";
	my ( $runComment, $pools, $lanes, $runFailed, $read1Failed, $read2Failed );
	my $failedLanes;
	my $runInfo;
	
	while ( ( $runComment, $pools, $lanes, $runFailed, $read1Failed, $read2Failed ) = $out->fetchrow_array ) {
		
		$runInfo  = "$runComment; Pools: $pools; $runFailed";
		if($read1Failed eq "T" || $read2Failed eq "T"){
			$runInfo .= "; failed lanes: $lanes";
		}
	}

	return $runInfo;
}

###############################################################
sub connectDB(){
	my $dbh = DBI->connect("DBI:mysql:database=solexa;host=$dbhost;port=3306", "$dbuser", "$dbpw",{PrintError => 0, RaiseError => 0	});
	unless($dbh) {
		$dbh = DBI->connect("DBI:mysql:database=solexa;host=localhost;port=3306", "$dbuser", "$dbpw",{PrintError => 0, RaiseError => 0	}) 		#connection over IP address doesn't work locally --> try localhost connection
		|| 
		die print "$DBI::errstr\n"};
	
	return $dbh;
}

#################################################################################################
sub getStartTime {

	my $tz = DateTime::TimeZone->new( name => 'local' );
	my $dt = DateTime->now;
	$dt->set_time_zone($tz);
	
	$dt->subtract( months => 2 );


	return substr($dt->year(),2,2)."".sprintf( "%02d", $dt->month() );

}

###############################################################
sub sendMail{
	my $flowcell = shift;
	my $state    = shift;
	my $sequencer= shift;	#sequencer the flowcell is on, if $state=1
	my $flowcelldir = shift;
	
	return if %mailSent && $mailSent{$flowcell}->{$state};
	
	####
	# states:
	# 1 --> Sequencing in progress (NO RTAComplete.txt)
	# 2 --> Sequencing finished, no demultiplexing (RTAComplete.txt, NO DemultiplexedBustardSummary.xml)
	# 3 --> Demultiplexing finished, no analysis (DemultiplexedBustardSummary.xml, NO Bamfiles in -b $bamdir)
	# 4 --> Analysis running --> Don't send an e-mail
	
	
		
	if($mailto ne ""){
				
		my $subject = "$flowcell - ";
		my $body	= "";
		if($state == 1){
			$subject .= "Run started on $sequencer";
		}elsif($state == 2){
			$subject .= "Run finished. You can start demultiplexing";
		}elsif($state == 3){
			$subject .= "Demultiplexing finished. You can start the pipeline.";
			$body = "Standard call:\nperl /data/isilon2/users/scripts/eclipse_workspace_wieland/Pipeline/initAnalysis.pl -fs $flowcell -sge custom.q -b /data/runs/Runs/SequenceBams/\nperl /data/isilon2/users/scripts/eclipse_workspace_wieland/Pipeline/initAnalysis.pl -fs $flowcell -sge custom.q -v vcf -var";
		
			# check if demultiplexing worked properly, i.e. if specified indices are correct
			my $demuxdir  = $flowcelldir . "/Demultiplexed";
			my @undeterminedFiles = glob($demuxdir.'/Undetermined*');
			my $demuxError = 0;
			my $novaSeqOffset = 1;			#Undetermined files > 2GB
			if ($flowcelldir =~ m/DMXX/) {
				$novaSeqOffset = 2.5;		#Undetermined files > 5GB
			}
			foreach my $undeterminedFile (@undeterminedFiles) {
				my $stat = stat($undeterminedFile);
				if ($stat->size > (2*1024*1024*1024 * $novaSeqOffset)) {
					$demuxError = 1;
				}
			}
			if ($demuxError) {
				my $msg = Mail::Send->new(Subject => "Demultiplexing Error for Flowcell $flowcell");
				$msg->to('user@url.de');
				my $fh = $msg->open;
				print $fh "Problem in directory $demuxdir";
				$fh->close or die "couldn't send whole message: $!\n";
			}
			
		}
		
		#send e-mail
		my $msg = Mail::Send->new(Subject => $subject);
		$msg->to(@mailTo);
		
		my $fh = $msg->open;
		$fh->close or die "couldn't send whole message: $!\n";
		
		if($mailFile ne ""){
			$mailSent{$flowcell}->{$state} = 1;
			open OUT, ">>$mailFile" or die "Can't open $mailFile to append!\n";
			print OUT "$flowcell\t$state\n";
			close OUT;
		}
	}
}

