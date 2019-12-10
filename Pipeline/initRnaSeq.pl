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
my $settings = "hg19";
my $query = "";
my $params = Utilities::getParams();
my $logfile  = "rnaseq.pipeline.log";
my $loglevel = "INFO";
my $infiles  = "";
my @caseSampleIds = {};
my @controlSampleIds = {};
my $caseSampleIdsFile = "";
my $controlSampleIdsFile = "";
my $bamFileLocation = "";
my $caseCountFile = "";
my $controlCountFile = "";
my $outputFolder = "";
my $experiment = "";
my $countFile = "";
my $fpkmFile = "";
my $conditionFile = "";
my $sgequeue = "";
my $annotationSource = "ucsc";
my $Rexecutable = "R";
my $fileFolder = "";
my $noCreateHeatMap = 0;
my $deSeq = "";
my $doRuv = 0;
my $analysis = "all";
my $fileExtension = "htseqcounts";
my $bamFileName = "merged.rmdup.bam";
my $noFPKMMerge = "0";
my $helptext      = 
"
This scripts starts the RNA-seq pipeline

-as	<case sample id's>
-us	<control sample id's>
-o	<output directory>
-e	<experiment>
-a	<analysis>	which analysis should be performed; possible values are deseq, deseq2, edger, dexseq and all (default: $analysis)
-ff 	<file folder>	folder where to look for the bam files
-nohm		do not create heatmap with all genes included (likely to take long) (default: no)
-nofm		no fpkm file merge
-ruv			remove unwanted variation with RUVr (for now for edgeR only)
-se	<settings>	required for genesymbol translation of the resulting table (default: $settings)
-fe <fileextension>	fileextionsion of count files (default: $fileExtension)
-bfn <bamfilename>	default: $bamFileName
-lf	<log file>	default: pipeline.log
-ll	<log level>	ERROR,INFO,DEBUG; default: INFO
-casef	<casefiles.tmp>	
-ctrlf	<controlfiles.tmp>
Alternatively use:
tmp file format:	SAMPLEID1\tBAMFILELOCATION1
			SAMPLEID2\tBAMFILELOCATION2
					...
					
-h	this help text\n";


GetOptions(
"o=s" => \$outputFolder,
"e=s" => \$experiment,
"as=s" => \$caseSampleIdsFile,
"us=s" => \$controlSampleIdsFile,
"lf=s" => \$logfile,
"se=s" => \$settings,
"a=s" => \$analysis,
"ll=s" => \$loglevel,
"ff=s" => \$fileFolder,
"nohm"	=> \$noCreateHeatMap,
"bfn=s" => \$bamFileName,
"fe=s" => \$fileExtension,
"ruv" => \$doRuv,
"casef=s" => \$caseCountFile,
"nofm" => \$noFPKMMerge,
"controlf=s" => \$controlCountFile,
"h" => \$help);

if ($help == 1 || $experiment eq "" || (($caseSampleIdsFile eq "" || $controlSampleIdsFile eq "") && ($caseCountFile eq "" || $controlCountFile eq ""))) {
	print $helptext;
	exit(1);
}

#create outputfolder if it doesn't exist
unless (-d $outputFolder) {
	my $ret = system("mkdir $outputFolder");
	if ($ret) {
		print STDERR "CMD died with return value: $ret\n";
		exit(11);
	}
}

Utilities::initLogger($outputFolder."/".$logfile,$loglevel);
my $logger = Utilities::getLogger();

#connect to database
my $db     = $params->{coredb}->{database};
my $host   = $params->{coredb}->{host};
my $port   = $params->{coredb}->{port};
my $user   = $params->{coredb}->{user};
my $passwd = $params->{coredb}->{password};

unless ( $dbh = DBI->connect( "DBI:mysql:database=$db", $user, $passwd ) ) {
	DBI->connect( "DBI:mysql:database=$db;host=$host;port=$port", $user, $passwd ) || die print "$DBI::errstr\n";
}

my $fileMerger = $params->{programs}->{rnaseq}->{htseqfilemerger};
my $translateId = $params->{programs}->{rnaseq}->{translateidcolumn};

#fill @caseSampleIds and @controlSampleIds
@caseSampleIds = &parseSampleIds($caseSampleIdsFile);
@controlSampleIds = &parseSampleIds($controlSampleIdsFile);

if ($experiment eq "") {
	$query = "select sbam from exomehg19.sample where name like '$caseSampleIds[0]'";
	my $out = $dbh->prepare($query) || exit $logger->error("$DBI::errstr");
	$out->execute || exit $logger->error("$DBI::errstr");
	if ($out->rows == 1) {
		$experiment = $out->fetchrow_array;
		$experiment =~ s/[\s-\/\\\(\)'"\^°§$%&=´`~*#<>|@]//;       #substitute \s, -, /, \, (, ), ... with emptystring
	} else {
		$experiment = "NONAME";
	}
}

$countFile = $experiment . ".all.counts";
$fpkmFile = $experiment . ".all.fpkms";
$conditionFile = $countFile . ".conditions";

if (!($caseCountFile && $controlCountFile)) {
	$caseCountFile = &getFileLocations(\@caseSampleIds, "case", $fileFolder);
	$controlCountFile = &getFileLocations(\@controlSampleIds, "control", $fileFolder);
}

#merge the count files with htseqOutputFileMerger.pl; kind of a safety-feature and easier for the later R script
my $command = "perl $prog_path/$fileMerger -as $caseCountFile -us $controlCountFile -cfn $countFile -condfn $conditionFile -o $outputFolder -fe $fileExtension";
if (&Utilities::executeCommand($command, "Running: $fileMerger to merge count files in $caseCountFile and $controlCountFile", $logger)) {
	exit(11);
}
#merge the fpkm files with htseqOutputFileMerger.pl like above
unless($noFPKMMerge) {
	my $command = "perl $prog_path/$fileMerger -as $caseCountFile -us $controlCountFile -cfn $fpkmFile -condfn $conditionFile -o $outputFolder -fe fpkm";
	if (&Utilities::executeCommand($command, "Running: $fileMerger to merge count files in $caseCountFile and $controlCountFile", $logger)) {
		exit(11);
	}
}

if ($analysis eq "all" || $analysis eq "deseq2") {
	##DESeq2 analysis
	my $deseq2OutputFolder = "$outputFolder/DESeq2";
	$command = "mkdir -p $deseq2OutputFolder";
	if (&Utilities::executeCommand($command, "making output directory for DESeq2", $logger)) {
		exit(11);
	}
	$deSeq = $params->{programs}->{rnaseq}->{RScripts}->{DESeq2_local};
	my $command = "$Rexecutable --no-restore --no-save --args $outputFolder/$countFile $outputFolder/$conditionFile $deseq2OutputFolder $experiment $outputFolder/$fpkmFile $noCreateHeatMap < $prog_path/$deSeq > $deseq2OutputFolder/DESeq2.Rlog 2>&1";
	if (&Utilities::executeCommand($command, "Running: DESeq2 R script for differential expression analysis", $logger)) {
		exit(11);
	}
	#Gene Ontology Analysis - goseq
	if (-f "$deseq2OutputFolder/session.RData" ) {
		my $goseqOutputFolder = "$deseq2OutputFolder/go_enrichment";
		$command = "mkdir -p $goseqOutputFolder";
		if (&Utilities::executeCommand($command, "making output directory for goseq", $logger)) {
			exit(11);
		}
		my $goseq = $params->{programs}->{rnaseq}->{RScripts}->{goseq};
		my $command = "$Rexecutable --no-restore --no-save --args $goseqOutputFolder/ $deseq2OutputFolder/session.RData $settings < $prog_path/$goseq > $goseqOutputFolder/goseq.Rlog 2>&1";
		
		if (&Utilities::executeCommand($command, "Running: GOseq script for Gene Ontology enrichment analysis", $logger)) {
			exit(11);
		}
	} else {
		print "Neccessary file 'session.RData' missing in DESeq2 output folder\n";
		exit(11);
	}
	#Pathway Analysis
	if(-f "$deseq2OutputFolder/session.RData" ) {
		my $pathwayOutputFolder = "$deseq2OutputFolder/pathway_enrichment";
		$command = "mkdir -p $pathwayOutputFolder";
		if (&Utilities::executeCommand($command, "making output directory for pathway enrichment analysis", $logger)) {
			exit(11);
		}
		my $gagePathwayCommand = $params->{programs}->{rnaseq}->{RScripts}->{pathway};
		my $command = "$Rexecutable --no-restore --no-save --args $settings $deseq2OutputFolder/".$experiment.".csv $pathwayOutputFolder $experiment < $prog_path/$gagePathwayCommand > $pathwayOutputFolder/pathway.Rlog 2>&1";
		if (&Utilities::executeCommand($command, "Running: GAGEPathway script for pathway enrichment analysis", $logger)) {
			exit(11);
		}
	} else {
		print "Neccessary file 'session.RData' missing in DESeq2 output folder\n";
		exit(11);
	}
}



sub parseSampleIds {
	my $file = shift;
	my @tmp = ();
	
	if ((-e $file) && (-f $file)) {
		open IN, $file || exit $logger->error("Couldn't open file $file");
		while(<IN>) {
			chomp($_);
			s/\s//;			#substitute possible space character with empty string
			unless (($_ =~ m/^#/) || ($_ =~ m/^$/)) {
				push(@tmp, $_);
			}
		}
	} else {		#sampleid passed directly
		if ($file =~ /[A-Z]*[0-9]*/) {
			push(@tmp, $file);
		}
	}
	if (@tmp < 1) {
		$logger->error("Couldn't find file $file!");
		exit(11);
	}
	return @tmp;
}

sub getFileLocations {
	my $sampleIdsRef = shift;
	my $name = shift;
	my $ff = shift;
	my $bamFileLocation = "";
	my $resultFileLocation = "$outputFolder/$name.tmp";
	
	open OUT, ">$resultFileLocation";
	if ($ff eq "") {
		foreach my $sId (@$sampleIdsRef) {
			my $c = "find -L " . $params->{settings}->{hg19_test}->{analysis}->{folder} . "/S*/ -maxdepth 1 -type d -name $sId 2> /dev/null";			#maxdepth 2, otherwise it'd be possible that backupfolders like /old/sId/ could be found instead of the actual one
			$logger->info("CMD: $c");
			open FIND, "$c | " or exit $logger->error("Error opening $c |");
			my $loc = <FIND>;
			chomp($loc);
			close FIND;
			if ($loc eq "") {
				$logger->error("Couldn't find bam file on system: $sId");
				exit(11);
			} else {
				$loc .= "/RNAout/paired-endout/$bamFileName";
				print OUT "$sId\t$loc\n";
			}
		}
	}
	close OUT;
	return $resultFileLocation;
}
