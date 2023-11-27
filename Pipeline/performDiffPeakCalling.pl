#!/usr/bin/perl

use strict;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);
use DateTime;
umask(002);

my $prog_path = dirname(abs_path($0));
require $prog_path."/Utilities.pm";

my $help     = 0;
my $dbh = "";
my $query = "";
my $params = Utilities::getParams();
my $logfile  = "pipeline.log";
my $loglevel = "INFO";
my @caseSampleIds;
my @controlSampleIds;
my $outputFolder = "";
my $experiment = "";
my $Rexecutable = "R";
my $out;
my $sampleSheet = "";
my $diffBind = "";
my $tissueName = "";
my $ipName = "";
my $peakCaller = "MACS2";
my $peakFormat = "bed";
my $broad      = 0;
my $species    = "";
my $loc        = "";
my $helptext      = 
"
This scripts starts the Differential-Peak-Calling Pipeline for ChIP-seq samples

-ca	<case sample id's>
-co	<control sample id's>
-o	<output directory>
-e	<experimentname>
-s	<species> either 'human' or 'mouse'
-b if broad peaks should be used
-lf	<log file>	default: pipeline.log
-ll	<log level>	ERROR,INFO,DEBUG; default: INFO
-h	this help text\n";

GetOptions(
"o=s" => \$outputFolder,
"e=s" => \$experiment,
"ca=s" => \@caseSampleIds,
"s=s" => \$species,
"b"   => \$broad,
"co=s" => \@controlSampleIds,
"lf=s" => \$logfile,
"ll=s" => \$loglevel,
"h" => \$help);

if ($help) {
	print $helptext;
	exit 11;
}

#create outputfolder if it doesn't exist
unless (-d $outputFolder) {
	my $ret = system("mkdir $outputFolder");
	if ($ret) {
		print STDERR "CMD died with return value: $ret\n";
		exit(11);
	}
}

if (! $species) {
	print "No species (-s <SPECIES>) specified! Possible values are 'human' or 'mouse'. Exiting ...\n";
	exit 11;
}

Utilities::initLogger($outputFolder."/".$logfile,$loglevel);
my $logger = Utilities::getLogger();

#connect to database
my $db     = $params->{coredb}->{database};
my $host   = $params->{coredb}->{host};
my $port   = $params->{coredb}->{port};
my $user   = $params->{coredb}->{user};
my $passwd = $params->{coredb}->{password};

unless ( $dbh = DBI->connect( "DBI:mysql:database=$db;host=$host;port=$port", $user, $passwd ) ) {
	DBI->connect( "DBI:mysql:database=$db;host=$host;port=$port", $user, $passwd ) || die print "$DBI::errstr\n";
}

if ($experiment eq "") {
	$experiment = "NONAME";
}

$sampleSheet = "$outputFolder/$experiment" . "_samplesheet.csv";
open OUT, ">$sampleSheet" || die print "Error opening samplesheet $sampleSheet for writing\n";
#prepare samle sheet
print OUT "SampleID,Tissue,Factor,Condition,Treatment,Replicate,bamReads,ControlID,bamControl,Peaks,PeakCaller,PeakFormat\n"; 		#print header line
#insert lines for case samples
&createSampleSheet(\@caseSampleIds,"case", $broad);
&createSampleSheet(\@controlSampleIds,"control", $broad);
close OUT;

##DiffBind analysis
$diffBind = $params->{programs}->{chipseq}->{DiffBind};
my $command = "$Rexecutable --no-restore --no-save --args $outputFolder $sampleSheet $species < $prog_path/$diffBind > $outputFolder/DiffBind.Rlog 2>&1";
if (&Utilities::executeCommand($command, "Running: DESeq2 R script for differential expression analysis", $logger)) {
	exit(11);
}

#######
## subs:
########
sub createSampleSheet {
	my $sampleIdsRef = shift;
	my $group = shift;
	my $br = shift;
	
	for my $caseId (@$sampleIdsRef) {
		$query = "select t.name, s.chipseqcontrol
from exomehg19.sample s 
left join exomehg19.tissue t on t.idtissue=s.idtissue
inner join solexa.sample2library sl on sl.idsample=s.idsample
inner join solexa.library l on l.lid=sl.lid
inner join solexa.libtype lt on lt.ltid=l.libtype
where lt.ltlibtype like 'ChIP-Seq%'
and s.name like '$caseId'";			#QUESTION: in DB, should we store name or idsample in chipseqcontrol??
		$out = $dbh->prepare($query) || exit $logger->error("$DBI::errstr");
		$out->execute || exit $logger->error("$DBI::errstr");
		if ($out->rows == 1) {
			($tissueName, $ipName) = $out->fetchrow_array;
			print OUT "$caseId,$tissueName,transcriptionFactor,$group,normal,1,";
			my $c = "find -L /data/isilon/seq/analysis/exomehg19/S*/$caseId/ChIP-Seqout/paired-endout/ -type f -name 'merged.bam'";
			$logger->info("CMD: $c");
			open FIND, "$c | " or exit $logger->error("Error opening $c |");
			my $loc = <FIND>;
			chomp($loc);
			close FIND;
			if ($loc eq "") {
				print("Couldn't find bam file on system for sample $caseId\n");
				exit(11);
			} else {
				print OUT "$loc,";
			}
			print OUT "$ipName,";
			$c = "find -L /data/isilon/seq/analysis/exomehg19/S*/$ipName/ChIP-Seqout/paired-endout/ -type f -name 'merged.bam'";
			$logger->info("CMD: $c");
			open FIND, "$c | " or exit $logger->error("Error opening $c |");
			$loc = <FIND>;
			chomp($loc);
			close FIND;
			chomp($loc);
			print OUT "$loc,";
			if ($br) {
				$c = "find -L /data/isilon/seq/analysis/exomehg19/S*/$caseId/ChIP-Seqout/paired-endout/ -type f -name 'ChIP.macs_peaks.broadPeak'";
			} else {
				$c = "find -L /data/isilon/seq/analysis/exomehg19/S*/$caseId/ChIP-Seqout/paired-endout/ -type f -name 'ChIP.macs_peaks.narrowPeak'";
			}
			$logger->info("CMD: $c");
			open FIND, "$c | " or exit $logger->error("Error opening $c |");
			$loc = <FIND>;
			chomp($loc);
			close FIND;
			chomp($loc);
			if ($loc eq "") {
				print("Couldn't find peak file on system for sample $caseId");
				exit(11);
			} else {
				print OUT "$loc,";
			}
		}
		print OUT "$peakCaller,$peakFormat\n";
		
		$loc        = "";
		$tissueName = "";
		$ipName     = "";
	}
}