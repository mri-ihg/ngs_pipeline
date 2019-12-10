#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Basename;
use DateTime;
use DBI;
use Log::Log4perl qw(get_logger :levels);

my $prog_path = dirname( abs_path($0) );
require $prog_path . "/Utilities.pm";

my $params     = Utilities::getParams();
# Config
my $dbhost = $params->{settings}->{hg19_test}->{exomedb}->{host};
my $dbuser = $params->{settings}->{hg19_test}->{exomedb}->{user};
my $dbpw   = $params->{settings}->{hg19_test}->{exomedb}->{password};


my $flowcell       = "";
my $runfolder      = $params->{settings}->{hg19_test}->{analysis}->{runfolder};
my $outputfolder   = "";
my $help           = 0;
my $convertBam     = 0;
my $tmpdir         = "";
my $checkOnly      = 0;
my $qualitFormat   = "Illumina";
my $checkQual      = 0;
my $sgequeue       = "";
my $sleep          = 0;
my $skipSleep      = 0;
my $firstTime      = "";
my $jobCount       = 0;
my $addedSleepTime = 0;
my $logfile        = "";
my $loglevel       = "INFO";
my $convertSingleSample = 0;
my $singleSampleName = "";
my $libtype = "";

GetOptions(
	"f=s"   => \$flowcell,
	"r=s"   => \$runfolder,
	"o=s"   => \$outputfolder,
	"h"     => \$help,
	"b"     => \$convertBam,
	"t=s"   => \$tmpdir,
	"c"     => \$checkOnly,
	"qf=s"  => \$qualitFormat,
	"sge=s" => \$sgequeue,
	"sl=s"  => \$sleep,
	"n=s"   => \$skipSleep,
	"st=s"  => \$firstTime,
	"lf=s"  => \$logfile,
	"ll=s"  => \$loglevel,
	"css"   => \$convertSingleSample,
	"ssn=s" => \$singleSampleName,
	"lt=s"  => \$libtype
);

if ( $qualitFormat eq "Check" ) {
	$checkQual = 1;
}

if ( $help == 1 ) {
	print "
-f	<flowcellname> OR empty: do it for all flowcells in <runfolder> OR \"gerald\": get a list of GERALD folders from <STDIN>
-r	<runfolder>
-o	<outputfolder>
-b	convert file into BAM Format; if not chosen only a hardlink is created
-t	<temdir> used for bam conversion; default $outputfolder/tmp
-c	just look for GERALD folders; don't do anything
-qf	<quality format>; possible formats \"Illumina\", \"Solexa\", \"Standard\" or \"Check\", default: \"$qualitFormat\"; if \"Check\": the script
	looks for the quality format using FastQC for each GERALD folder
-sge	<sge_queue> submit pipeline job to specified SGE queue
-st	time at which the job(s) should start; format: DD.MM.YYYY,hh:mm; default: now
-sl	<minutes_to_sleep> time the script waits until it submits the next SGE job, default: 0
-n	number of jobs that should be started at the beginning without sleeping (i.e. there is no sleeping time between them)
-css convert only a single sample
-ssn name of this single sample
-lf	log file; default: $outputfolder/pipeline.log
-ll	log level: ERROR,INFO,DEBUG; default: INFO
-h	print this helptext
\n";
	exit(0);
}

if ( $logfile eq "" ) {
	$logfile = "$outputfolder/pipeline.log";
}

#get params
my $params = Utilities::getParams();
my $java   = $params->{programs}->{java}->{path};
Utilities::initLogger( $logfile, $loglevel );
my $logger = Utilities::getLogger();

if ( $tmpdir eq "" ) {
	system("mkdir $outputfolder\/tmp 2>\/dev\/null");
	$tmpdir = $outputfolder . "\/tmp";
}

if ( $flowcell eq "gerald" ) {    #get list of gerald folders from STDIN
	my $fc = "";
	while (<STDIN>) {
		chomp;
		if ( $_ =~
/_((a|b|c|d|e|f|g|h|i|j|k|l|m|n|o|p|q|r|s|t|u|v|w|x|y|z|A|B|C|D|E|F|G|H|I|J|K|L|M|N|O|P|Q|R|S|T|U|V|W|X|Y|Z|0|1|2|3|4|5|6|7|8|9)+(XX|XY))/
		  )
		{
			$fc = $1;
			&linkFiles( $fc, $_ );
		}
	}
} elsif ($convertSingleSample && $singleSampleName ne "") {
	if ($flowcell eq "") {    #Safety feature: do not convert entire /data/runs/Runs/ folder if $flowcell is empty but $convertSingleSample is set
		print STDERR "Tried to convert single sample into BAM but no samplename was specified\n";
	} else {
		&linkFiles($flowcell, "");
	}
} elsif ( $flowcell ne "") {    #just for a single flowcell
	&linkFiles( $flowcell, "" );
}
else {                         #get flowcells from file list

	open( DIRS, "ls -d $runfolder/* | grep \"/.*/[0123456789]*_.*\"|" );
	while (<DIRS>) {
		chomp;
		my @columns = split("_");
		$flowcell = $columns[-1];
		if ( $flowcell eq "PE" ) {
			$flowcell = $columns[-2];
		}
		if ( ( length $flowcell ) > 9 )
		{ #some of the first HiSeq flowcell names within the dirname contain the A.../B... that marks the flowcellposition within the HiSeq --> no part of our name
			$flowcell = substr( $flowcell, ( ( length $flowcell ) - 9 ), 9 );
		}

		&linkFiles( $flowcell, "" );
	}

}

##############################################################
sub linkFiles {
	my $fc     = shift;
	my $gerald = shift;
	$logger->info("Start processing $fc...");

	my $fcpath = "";
	if ( $gerald eq "" )
	{    #if no gerald folder is given try to find it in the Run folder

		$fcpath = qx/ls -d $runfolder\/*$fc\/ 2> \/dev\/null/;

		#print "test: $fc : $fcpath";

		chomp $fcpath;

		if ( $fcpath eq "" ) {
			$logger->error("No folder found for Flowcell $fc!");
			return;
		}

		my @columns = split( /\//, $fcpath );
		my $rgname  = $columns[-1];

		my $demultiplexed = 0;
		
		if(glob("$fcpath/Data/Intensities/BaseCalls/".$singleSampleName."*.fastq.gz")){		#TW 22.09.2014: in the case of a MiSeq Run --> fastq.gz files are simply in BaseCalls folder
			$gerald = "$fcpath/Data/Intensities/BaseCalls/";
			if ($convertSingleSample) {
				&createLink( $fc, $gerald, "$outputfolder/", $rgname, "", 1 );
			} else {
				&createLink( $fc, $gerald, "$outputfolder/$fc/", $rgname, "", 1 );
			}
			
		}elsif ( -e "$fcpath/Demultiplexed/Project_all" )
		{ #CASAVA 1.8+ --> demultiplexed files are in these subdirectories grouped by samples
			
			open( IN,
				"find $fcpath/Demultiplexed/Project_all -name \"Sample_" . $singleSampleName . "*\" |" );

			my @samples;
			while (<IN>) {

				#print "test: $_";
				chomp;
				push( @samples, $_ );
			}
			close IN;

			foreach (@samples) {
				my @path = split( /\//, $_ );
				#my ( $dummy, $sample ) = split( "_", $path[-1] );
				my $sample = $path[-1];
				$sample =~ s/Sample_//;
				
				if ($sample =~ m/(.*)_/g) {			#for 10x libraries
					$sample = $1;
				}
				
				if ($convertSingleSample) {
					&createLink( $fc, $_, "$outputfolder/",	$rgname, $sample, 1 );
				} else {
					unless ($checkOnly) {
						system("mkdir $outputfolder/$fc 2> /dev/null");
					}
					&createLink( $fc, $_, "$outputfolder/$fc/$path[-1]/",
						$rgname, $sample, 1 );
				}
			}
		}
		elsif (
			-e "$fcpath/Data/Intensities/BaseCalls/Demultiplexed/SamplesDirectories.csv"
		  )
		{    #Case 1: if this file exists the flowcell was indexed.
			open( IN,
"$fcpath/Data/Intensities/BaseCalls/Demultiplexed/SamplesDirectories.csv"
			);
			my %samples;
			while (<IN>) {
				chomp;
				my @columns = split(",");
				next if $columns[0] eq "FCID";

				#next
				#  if -d "$outputfolder/$fc/$columns[-1]"
				#  ;    #if links have already been generated

				$gerald =
qx/ls -dt $fcpath\/Data\/Intensities\/BaseCalls\/Demultiplexed\/\/$columns[-1]\/GERALD*  2> \/dev\/null/;
				chomp $gerald;
				unless ($checkOnly) {
					system("mkdir $outputfolder/$fc 2> /dev/null");
				}
				$rgname .= "_$columns[-1]";
				unless ( $samples{ $columns[2] } ) {
					&createLink( $fc, $gerald,
						"$outputfolder/$fc/$columns[-1]/",
						$rgname, $columns[2], 0 );
					$samples{ $columns[2] } = 1;
				}

			}
			unless ($checkOnly) {
				system(
"ln $fcpath/Data/Intensities/BaseCalls/Demultiplexed/SamplesDirectories.csv $outputfolder/$fc/ 2> /dev/null"
				);
			}
		}
		else {
			$gerald =
qx/ls -dt $fcpath\/Data\/Intensities\/Bustard*\/GERALD*  2> \/dev\/null/
			  ;    #Case 2: RTA went wrong

			if ( $gerald eq "" ) {
				$gerald =
qx/ls -dt $fcpath\/Data\/*Firecrest*\/Bustard*\/GERALD*  2> \/dev\/null/
				  ;    #Case 3: RTA went wrong

			}
			if ( $gerald eq "" ) {
				$gerald =
qx/ls -dt $fcpath\/Data\/*IPAR*\/Bustard*\/GERALD*  2> \/dev\/null/
				  ;    #Case 3: RTA went wrong

			}

			if ( $gerald eq "" ) {    #Case 4: no demultiplexing, RTA ok
				$gerald =
qx/ls -dt $fcpath\/Data\/Intensities\/BaseCalls\/GERALD* 2> \/dev\/null/;

				if ( $gerald eq "" )
				{ # --> no GERALD directory found --> Illumina pipeline has not yet run or sequencing is in progress
					$logger->error(
						"No GERALD dir found for $fc in path $fcpath!");
					return;
				}

			}
			chomp $gerald;
			&createLink( $fc, $gerald, "$outputfolder/$fc/", $rgname, "", 0 );
		}

	}
	else {

		#print "$fc, $gerald, $outputfolder/$fc/, , , 0\n";
		my $rgname = $fc;
		if($gerald =~ /\/(\d*\w*XX\w*)\//){
			$rgname = $1;
		}
		&createLink( $fc, $gerald, "$outputfolder/$fc/", $rgname, "", 0 );
	}

}

######################################################################################
sub createLink {
	my $fc       = shift;
	my $gerald   = shift;
	my $outdir   = shift;
	my $rgname   = shift;
	my $sample   = shift;
	my $casava18 = shift;

	my $platform         = "Illumina";
	my $sequencingCenter =
	  "'Helmholtz Center Munich, Institute of Human Genetics'";

	my @columns = split( ' ', $gerald );
	if ( @columns > 1 ) {
		$logger->error(
			"WARNING! More than one GERALD dir found! Only newest one linked!");
		$gerald = $columns[0];
	}
	
	unless ($checkOnly) {
		system("mkdir $outdir 2> /dev/null");

		my $query;

		if ($convertBam) {
			if ( $sample ne "" )
			{    #no sample id given --> no multiplexed flowcell
				$query =
"select exomehg19.sample.name, library.lname, lane.alane, ldescription
	from run 
	left join lane on lane.rid=run.rid 
	left join pool on lane.idpool=pool.idpool 
	left join library2pool on library2pool.idpool=pool.idpool
	left join library on library2pool.lid=library.lid
	left join sample2library on sample2library.lid=library.lid
	left join exomehg19.sample on exomehg19.sample.idsample=sample2library.idsample
	where run.rname='$fc' and exomehg19.sample.name='$sample'";
			}
			else {

#no multiplexed flowcell --> group by lane, because it could be a pool (more than one libs per lane) without multiplexing
				$query =
"select group_concat(exomehg19.sample.name), group_concat(library.lname), lane.alane, group_concat(ldescription)
	from run 
	left join lane on lane.rid=run.rid 
	left join pool on lane.idpool=pool.idpool 
	left join library2pool on library2pool.idpool=pool.idpool
	left join library on library2pool.lid=library.lid
	left join sample2library on sample2library.lid=library.lid
	left join exomehg19.sample on exomehg19.sample.idsample=sample2library.idsample
	where run.rname='$fc'
	group by lane.alane;";
			}
			
			$logger->debug($query);
			my $dbh = &connectDB();

			my $out = $dbh->prepare($query)
			  || exit $logger->error("$DBI::errstr");
			$out->execute || exit $logger->error("$DBI::errstr");
			my $foundInDB = 0;
			while ( my ( $sn, $ln, $lane, $descr ) = $out->fetchrow_array ) {
				$foundInDB = 1;
				unless ($descr) {
					$descr = "";
				}
				unless ($sn)
				{    #No connection to sample in the DB for older flowcells
					my $dummy;
					( $sn, $dummy ) = split( /(\s|_)/, $ln );
				}

				if ($libtype =~ m/scATAC-Seq/i && $casava18) {
					my $sampleName10x = basename($gerald);
					$sampleName10x =~ s/Sample_//;
					my $command = "zip -j $outputfolder/$sampleName10x.zip";
					my $formatedLane = sprintf( "%03d", $lane );
					open(IN, "find $gerald -maxdepth 1 -name \"*.fastq.gz\" |");
					while (<IN>) {
						chomp;
						$command .= " $_";
					}
					&submitSGEJob($command, "ZIP" . $sample, $outputfolder, $fc);
					close IN;
				} elsif ($casava18)	{    #in Casava 1.8 there can be more than one file per lane
					my $formatedLane = sprintf( "%03d", $lane );
					open( IN, "find -L $gerald -maxdepth 1 -name \"*\_L$formatedLane\_R1\_*.fastq.gz\" |");
					while (<IN>) {
						chomp;
						my $command = "$java -Xmx4g -XX:ParallelGCThreads=1 -jar /usr/local/packages/seq/picard/picard.jar FastqToSam ";
						$command .= "FASTQ=$_ ";
						my $fastq2 = $_;
						$fastq2 =~ s/_R1_/_R2_/;
						if(-e $fastq2){
							$command .= "FASTQ2=$fastq2 ";		#if it is a paired end run
						}
						

						$fastq2 =~ s/_R2//;
						$fastq2 =~ s/fastq.gz/bam/;

						$fastq2 =~ /_(\d+).bam/;
						my $rgIndex = $1;
						$fastq2 = basename($fastq2);

						$command .=
						  "QUALITY\_FORMAT=Standard OUTPUT=$outdir/$fastq2 ";
						
						$command .=
						    "READ\_GROUP\_NAME=$rgname\_s\_" . $lane
						  . "\_$rgIndex SAMPLE\_NAME='$sn' LIBRARY\_NAME='$ln' ";
						$command .=
"PLATFORM=$platform SEQUENCING\_CENTER=$sequencingCenter DESCRIPTION='"
						  . $descr . "'";

						if ( $tmpdir ne "" ) {
							$command .= " TMP_DIR=" . $tmpdir;
						}

						my $test = qx/ls $outdir\/$fastq2 2> \/dev\/null/;
						

						#print "$command\n";

						unless ($test)
						{    #tests if file exists; -e doesn't work ?!

							#print "Didn't find $fastq2\n";
							#print $test;
							#print "$command\n";
							$logger->debug($command);
							if ( $sgequeue eq "" ) {
								system($command);
							}
							else {
								&submitSGEJob( $command, $fastq2, $outdir,
									$fc );
							}

						}

					}

				}
				else {
					&buildCommand($outdir,$gerald,$lane,$qualitFormat,$rgname,$sn,$ln,$sequencingCenter,$platform,$descr,$fc);
				}

			}
			
			unless($foundInDB){			#if flowcell was not found in DB --> old flowcell; just convert it without information
				
				my @miseq = glob("$gerald/*_L001_R1_*.fastq.gz");
				
				if(@miseq>1){		#apparently miseq run with samples not in database
					foreach(@miseq){
						my $command =
"$java -Xmx4g -XX:ParallelGCThreads=1 -jar /usr/local/packages/seq/picard/picard.jar FastqToSam ";
						$command .= "FASTQ=$_ ";
						my $fastq2 = $_;
						$fastq2 =~ s/_R1_/_R2_/;
						if(-e $fastq2){
							$command .= "FASTQ2=$fastq2 ";		#if it is a paired end run
						}
						

						$fastq2 =~ s/_R2//;
						$fastq2 =~ s/fastq.gz/bam/;

						$fastq2 =~ /_(\d+).bam/;
						my $rgIndex = $1;
						$fastq2 = basename($fastq2);
						
						my ($sn,@tmp) = split("_",$fastq2);

						$command .=
						  "QUALITY\_FORMAT=Standard OUTPUT=$outdir/$fastq2 ";
						
						$command .=
						    "READ\_GROUP\_NAME=$rgname\_s1\_$rgIndex SAMPLE\_NAME='$sn' LIBRARY\_NAME='$sn\_LIB1' ";
						$command .=
"PLATFORM=$platform SEQUENCING\_CENTER=$sequencingCenter DESCRIPTION=''";

						if ( $tmpdir ne "" ) {
							$command .= " TMP_DIR=" . $tmpdir;
						}

						my $test = qx/ls $outdir\/$fastq2 2> \/dev\/null/;
						

						unless ($test)
						{    #tests if file exists; -e doesn't work ?!

							#print "Didn't find $fastq2\n";
							#print $test;
							#print "$command\n";
							$logger->debug($command);
							if ( $sgequeue eq "" ) {
								system($command);
							}
							else {
								&submitSGEJob( $command, $fastq2, $outdir,
									$fc );
							}

						}
					}
				}else{
					for(my $i=1; $i <= 8; $i++){	#try it for lanes 1 - 8
						&buildCommand($outdir,$gerald,$i,$qualitFormat,$rgname,"UNDEF","UNDEF",$sequencingCenter,$platform,"UNDEF",$fc);
					}
				}
				
			}
			
			
			if ($checkQual) {
				$qualitFormat = "Check";
			}

		}
		else {
			system("ln $gerald/* $outdir 2> /dev/null");
		}
	}

}
#################################################################################################
sub buildCommand {
	my $outdir           = shift;
	my $gerald           = shift;
	my $lane             = shift;
	my $qualitFormat     = shift;
	my $rgname           = shift;
	my $sn               = shift; 
	my $ln	             = shift;
	my $sequencingCenter = shift;
	my $platform         = shift;
	my $descr            = shift;
	my $fc               = shift;
	
	 
	my $command =
	  "$java -Xmx2g -XX:ParallelGCThreads=1 -jar /usr/local/packages/seq/picard/picard.jar FastqToSam ";

	my $fastq;

	if ( -e ( "$gerald/s\_" . $lane . "\_1\_sequence.txt" ) ) {  #paired end run
		    #print "paired end: s\_" . $lane . "\_1\_sequence.txt\n";
		$command .=
		    "FASTQ=$gerald/s\_" . $lane
		  . "\_1\_sequence.txt FASTQ2=$gerald/s\_"
		  . $lane
		  . "\_2\_sequence.txt ";
		$fastq = "$gerald/s\_" . $lane . "\_1\_sequence.txt";
	}
	elsif ( -e ( "$gerald/s\_" . $lane . "\_1\_sequence.fastq.gz" ) )
	{   #SPECIAL CASE (only for 1-2 flowcells) where we had to recreate the data
		$command .=
		    "FASTQ=$gerald/s\_" . $lane
		  . "\_1\_sequence.fastq.gz FASTQ2=$gerald/s\_"
		  . $lane
		  . "\_2\_sequence.fastq.gz ";
		$fastq = "$gerald/s\_" . $lane . "\_1\_sequence.fastq.gz";
	}
	else {
		$command .=
		  "FASTQ=$gerald/s\_" . $lane . "\_sequence.txt ";    #single end run
		$fastq = "$gerald/s\_" . $lane . "\_sequence.txt";

	}

	if ( $qualitFormat eq "Check" ) {    #if the quality format isn't defined
		my $dummy =
qx/ cat $outdir\/s\_*\_fastqc\/fastqc\_data.txt 2> \/dev\/null |grep Encoding/;

		unless ($dummy) {                #fastqc hasnt been run
			system(
				"$prog_path/FastQC/fastqc -o $outdir --extract -t 2 $fastq" );
			$dummy =
qx/ cat $outdir\/s\_*\_fastqc\/fastqc\_data.txt 2> \/dev\/null |grep Encoding/;
		}
		chomp $dummy;
		$logger->info("Detected Quality Format: $dummy");

		#print "dummy: $dummy\n";
		( $dummy, $dummy, $qualitFormat ) =
		  split( ' ', $dummy );

		#print "qualit: $qualitFormat\n";
		if (   ( $qualitFormat eq "1.5" )
			|| ( $qualitFormat eq "1.3" ) )
		{
			$qualitFormat = "Illumina";
		}
		else {
			$qualitFormat = "Solexa";
		}
	}

	$command .=
	    "QUALITY\_FORMAT=$qualitFormat OUTPUT=$outdir/s\_" . $lane
	  . "\_sequence.bam ";
	$command .=
	    "READ\_GROUP\_NAME=$rgname\_" . $lane
	  . " SAMPLE\_NAME='$sn' LIBRARY\_NAME='$ln' ";
	$command .=
	  "PLATFORM=$platform SEQUENCING\_CENTER=$sequencingCenter DESCRIPTION='"
	  . $descr . "'";
	if ( $tmpdir ne "" ) {
		$command .= " TMP_DIR=" . $tmpdir;
	}

	my $outfile = "$outdir/s\_" . $lane . "\_sequence.bam ";
	my $test    = qx/ls $outfile 2> \/dev\/null/;

	unless ($test) {    #tests if file exists; -e doesn't work ?!

		#print "Didn't find $outfile\n";
		#print $test;
		#print "$command\n";
		$logger->debug($command);
		if ( $sgequeue eq "" ) {
			system($command);
		}
		else {
			my $bin = "";
			if ( $rgname =~ /_(0\d\d)/ ) {
				$bin = $1;
			}
			&submitSGEJob( $command, "s\_" . $lane . "\_sequence.bam ",
				$outdir, $fc . $bin );
		}

	}
}
#################################################################################################
sub connectDB() {
	my $dbh =
	  DBI->connect( "DBI:mysql:database=solexa;host=$dbhost;port=3306",
		"$dbuser", "$dbpw", { PrintError => 0, RaiseError => 0 } );
	unless ($dbh) {
		$dbh =
		  DBI->connect( "DBI:mysql:database=solexa;host=localhost;port=3306",
			"$dbuser", "$dbpw", { PrintError => 0, RaiseError => 0 } )
		  or #connection over IP address doesn't work locally --> try localhost connection

		  exit $logger->error("$DBI::errstr");
	}
	return $dbh;
}

#################################################################################################
sub submitSGEJob {
	my $command     = shift;
	my $outfilename = shift;
	my $outdir      = shift;
	my $fc          = shift;

	$outfilename =~ s/.bam//;
	my $jobName = "BAM" . $outfilename . $fc;
	$jobName =~ s/\s//g;

	#generate qsub string
	my $qsubstring = "qsub -q $sgequeue -N $jobName ";

	if ( $jobCount >= $skipSleep ) { #start the first "skipSleep" jobs right now
		$addedSleepTime +=
		  $sleep;    #add "sleep" minutes to the starttime of all other jobs
	}

	my $starttime = &getStartTime($addedSleepTime);

	$qsubstring .= "-a $starttime ";

	$jobCount++;

	#create temp script to submit
	open( OUT, ">$outdir/sgeSubmTmp.sh" )
	  || exit $logger->error("Cannot open $outdir/sgeSubmTmp.sh");
	print OUT "#\$  -S /bin/bash\n";
	print OUT $command;
	print OUT "\nif [ \$? -ne 0 ]
then
  $command
else";
	close OUT;

	$qsubstring .= "$outdir/sgeSubmTmp.sh";

	if ($checkOnly) {
		$logger->debug($qsubstring);
	}
	else {

		#$logger->info("SUBMITTING: $qsubstring: $pgrpath/$pgr $param");
		#print "\n\nSUBMITTING: $qsubstring\n";
		$logger->info("SUBMITTING: $qsubstring");
		$logger->debug($command);
		system($qsubstring);

		#system("sh $pgrpath/sgeSubmTmp.sh");
	}
	system("rm $outdir/sgeSubmTmp.sh");    #remove tmp script
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
