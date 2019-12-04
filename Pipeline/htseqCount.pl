#!/usr/bin/perl

use strict;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);
use File::stat;
umask(002);

#include Utilities.pm
my $prog_path = dirname(abs_path($0));
require $prog_path."/Utilities.pm";

my $settings      = "hg19_rna";
my $logfile  	  = "";
my $loglevel 	  = "";
my $inputBAM	  = "";
my $help		  = "";
my $command		  = "";
my $keepTmp		  = "";
my $htseqMode	  = "intersection-nonempty";
my $isStrandedRNA = 0;
my $aligner = "";
my $outputFileName = "";
my $sampleName     = "";
my $sortOption = "pos";

my $helptext      = 
"
-i	<infile>	bam file required
-m  <htseq mode>	possible values: \"union\", \"intersection-strict\" or \"intersection-nonempty\" (default: $htseqMode)
-a	<aligner>	mandatory if -sr is specified, star, tophat, gem (needed to choose the correct stranded option for htseq-count) 
--lf	<log file>		the log file for the pipeline logging (default: pipeline.log)
--ll	<log level>		the log level for the pipeline logging; available options: ERROR, INFO, DEBUG; (default: $loglevel)
--se	<settings>		(default: $settings)
-sr				if stranded RNA-seq
-k				keep temporary files
-h				show help text\n
";

GetOptions(
"h" => \$help,
"a=s" => \$aligner,
"i=s" => \$inputBAM,
"s=s" => \$sampleName,
"lf=s" => \$logfile,
"ll=s" => \$loglevel,
"m=s" => \$htseqMode,
"se=s" => \$settings,
"sr" => \$isStrandedRNA,
"k" => \$keepTmp
);

if (($help == 1) || ($isStrandedRNA && $aligner eq "")) {
	print $helptext;
	exit(1);
}

my $params = Utilities::getParams();
my $samtools      = $params->{programs}->{samtools}->{path};
my $sammaxmem     = $params->{programs}->{samtools}->{maxmem};

my $htseqCount	  = $params->{programs}->{htseq}->{count};
my $featureFile	  = $params->{settings}->{$settings}->{gemannotation};	#use same gtf file as for the alignment!

#init the logger
Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();

$outputFileName = "$inputBAM.htseqcounts";

#check file size (htseq-count has problems with large BAM files that are not name sorted -> buffer overflow)
my $stat = stat($inputBAM);
if ($stat->size > 2147483648 ) {		#rough estimate, everything larger than 2GB is too much
	$command = "samtools sort -n -@ 5 -m 2G $inputBAM $inputBAM.namesorted";
	if (&Utilities::executeCommand($command, "Namesorting input BAM file", $logger, "htseq-count-sort")) {
		exit(11);
	}
	$inputBAM = "$inputBAM.namesorted.bam";
	$sortOption = "name";
}

$command = "$htseqCount -f bam -t exon -i gene_id -m $htseqMode -r $sortOption ";
#$command .= "-a 0 ";

if (!$isStrandedRNA) {
	#connect to database
	my $db     = $params->{coredb}->{database};
	my $host   = $params->{coredb}->{host};
	my $port   = $params->{coredb}->{port};
	my $user   = $params->{coredb}->{user};
	my $passwd = $params->{coredb}->{password};
	
	my $dbh;
	unless ( $dbh = DBI->connect( "DBI:mysql:database=$db;host=$host;port=$port", $user, $passwd ) ) {
		DBI->connect( "DBI:mysql:database=$db;host=$host;port=$port", $user, $passwd ) || die print "$DBI::errstr\n";
	}
	my $query = "select a.name from exomehg19.sample s inner join solexa.sample2library sl on sl.idsample=s.idsample inner join solexa.library l on l.lid=sl.lid inner join solexa.assay a on a.idassay=l.idassay where s.name like '$sampleName';";
	my $out = $dbh->prepare($query) || exit $logger->error("$DBI::errstr");
	$out->execute || exit $logger->error("$DBI::errstr");
	if ($out->rows == 1) {
		my $assayName = $out->fetchrow_array;
		if (($assayName =~ m/Stranded Total RNA$/) || ($assayName =~ m/Stranded mRNA$/)) {
			$isStrandedRNA=1;
		}
	}
}
if ($isStrandedRNA) {
	if ($aligner eq "tophat") {
		$command .= "--stranded=yes ";
	} else {
		$command .= "--stranded=reverse ";
	}
} else {
	$command .= "--stranded=no ";
}




$command .= "$inputBAM $featureFile > $outputFileName 2> $outputFileName.log";
if (&Utilities::executeCommand($command, "Running htseq count of $inputBAM.namesorted", $logger, "htseq-count")) {
	exit(11);
}
