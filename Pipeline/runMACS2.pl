#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;
use Getopt::Long qw(:config no_ignore_case);
use Cwd qw(abs_path);
use File::Basename;
use Pod::Usage;
use Scalar::Util;
umask(002);

my $start_time = time();

my $prog_path = dirname( abs_path($0) );
require $prog_path."/Utilities.pm";

my $chipfile   = "";
my $inputfile  = "";
my $outdir     = "";
my $settings   = "";
my $help	   = 0;
my $man		   = 0;
my $broad      = 0;
my $logfile    = "pipeline.log";
my $Rexecutable = "Rscript";
my $loglevel   = "INFO";
my $organism   = "";
my $dbh 	   = "";
my $query      = "";
my $out        = "";
my $sampleName = "";
my $inputSample = "";
my $broadPeak  = "";

GetOptions(
	"c=s"  => \$chipfile,
	"i=s"  => \$inputfile,
	"o=s"  => \$outdir,
	"se=s" => \$settings,
	"r=s"  => \$organism,
	"h"    => \$help,
	"b"    => \$broad,
	"lf=s" => \$logfile,
	"ll=s" => \$loglevel,
	"sn=s" => \$sampleName,
	"man"  => \$man
);

pod2usage( {-exitval => 0  ,-verbose => 1} ) if $help;
pod2usage( {-exitval => 0  ,-verbose => 2} ) if $man;
pod2usage( {-exitval => 1  ,-verbose => 1} ) if $chipfile eq "" || $settings eq "";

Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();
$logger->info("Started runMACS2");

my $params    = Utilities::getParams();
my $ref       = $params->{settings}->{$settings}->{reference};
my $sppScript = $params->{programs}->{chipseq}->{spp};

#connect to database
my $db     = $params->{coredb}->{database};
my $host   = $params->{coredb}->{host};
my $port   = $params->{coredb}->{port};
my $user   = $params->{coredb}->{user};
my $passwd = $params->{coredb}->{password};

unless ( $dbh = DBI->connect("DBI:mysql:database=$db;host=$host;port=$port", $user, $passwd ) ) {
	$dbh = DBI->connect( "DBI:mysql:database=$db;host=$host;port=$port",$user, $passwd )  || die print "$DBI::errstr\n";
}

$outdir = dirname($chipfile) if $outdir eq "";

my $inputoption = "";
if($inputfile ne ""){
	$inputoption = "-c $inputfile";
} else {
	$query = "select s.chipseqcontrol
from exomehg19.sample s 
inner join solexa.sample2library sl on sl.idsample=s.idsample
inner join solexa.library l on l.lid=sl.lid
inner join solexa.libtype lt on lt.ltid=l.libtype
where lt.ltlibtype like 'ChIP-Seq%'
and s.name like '$sampleName'";			#QUESTION: in DB, should we store name or idsample in chipseqcontrol??
	$out = $dbh->prepare($query) || exit $logger->error("$DBI::errstr");
	$out->execute || exit $logger->error("$DBI::errstr");
	if ($out->rows == 1) {
		$inputSample = $out->fetchrow_array;
		my $c = "find -L /PATHTO/S*/$inputSample/ChIP-Seqout/paired-endout/ -type f -name 'merged.bam'";			#maxdepth 2, otherwise it'd be possible that backupfolders like /old/sId/ could be found instead of the actual one
		$logger->info("CMD: $c");
		open FIND, "$c | " or exit $logger->error("Error opening $c |");
		my $loc = <FIND>;
		chomp($loc);
		close FIND;
		if ($loc eq "") {
			$logger->error("Couldn't find bam file on system for sample $sampleName");
		} else {
			$inputfile = $loc;
			unless (-e $inputfile . ".filtered.bam") {
				my $filterCommand = "samtools view -bh -q 5 $inputfile > " . $inputfile . ".filtered.bam";
				&Utilities::executeCommand($filterCommand, "Filtering low quality mapping before peak calling ...", $logger);
			}
			$inputoption = "-c " . $inputfile . ".filtered.bam";
		}
	} else {
		$logger->info("No ChIP-seq control is specified. Calling peaks without control ...");
	}
}

unless(Scalar::Util::looks_like_number($organism)){
	if($organism eq "human"){
		$organism = "hs";
	}elsif($organism eq "mouse"){
		$organism = "mm";
	}elsif($organism eq "C. elegans"){
		$organism = "ce";
	}elsif($organism eq "fruitfly"){
		$organism = "dm";
	}else{
		$logger->info("Getting genome size...");
		open FAI,"awk 'BEGIN{sum=0}{sum+=\$2}END{print sum}' $ref.fai |" or exit $logger->error("Can't open $ref.fai!");
		$organism = <FAI>;
		chomp $organism;
		close FAI;
	}
}

## M. Heinig (ICB) suggested to filter read with bad mapping quality  
my $filterCommand = "samtools view -bh -q 5  $chipfile > " . $chipfile . ".filtered.bam";
&Utilities::executeCommand($filterCommand, "Filtering low quality mapping before peak calling ...", $logger);

##run spp script -p=<NUMTHREADS> (uses snow library!)
my $sppCommand="$Rexecutable $prog_path/$sppScript -rf -c=$chipfile.filtered.bam -tmpdir=$outdir -odir=$outdir -savp=$outdir/spp_phantom_peak.pdf -out=$outdir/spp_phantom_qc.txt &> $outdir/spp_phantom.log ";
&Utilities::executeCommand($sppCommand, "Running: $sppScript", $logger);

my $command;
#running macs2
$command  = "macs2 callpeak --outdir $outdir -f BAMPE -n ChIP.macs -g $organism --buffer-size 1000000 -t " . $chipfile . ".filtered.bam $inputoption";
&Utilities::executeCommand($command, "Running MACS2 norrow peak calling ...", $logger);

$command .= " --broad";
&Utilities::executeCommand($command, "Running MACS2 broad peak calling ...", $logger);


#create links
if($broad){
	unlink("$outdir/ChIP.macs.bed");
	system("ln -s ChIP.macs_peaks.broadPeak $outdir/ChIP.macs.bed");
}else{
	unlink("$outdir/ChIP.macs.bed");
	system("ln -s ChIP.macs_peaks.narrowPeak $outdir/ChIP.macs.bed");
}

my $end_time = time();
$logger->info("runMACS2 finished in " . &Utilities::seconds_to_ddhhmmss($end_time - $start_time) . " (ddd:hh:mm:ss)");


=head1 NAME

runMACS2.pl

=head1 SYNOPSIS

runMACS2.pl -c chip.bam -i input.bam 
 
=head1 DESCRIPTION

Wrapper script that calls MACS2 for a ChIP-Seq sample and an (optional) Input DNA sample.
Does all the necessary conversions. macs2 must be in path!

=head1 OPTIONS

 -c	<chip.bam>; BAM file of the ChIP experiment. Best choice is a file with duplicate reads (i.e. merged.bam); REQUIRED
 -i	<input.bam>; BaM file of the Input DNA experiment; optional
 -o	</output/dir>; directory where the output files will be generated; default: directory where chip.bam lies
 -r	organism; used for calculation of effective genome size; currently supported by MACS2: human, mouse, C. elegans and fruitfly
 	The effective genome size can also be given directly as number.
 	if empty: use the fasta index file of the reference genome specified in settings to calculate genome size
 -b	use "broad" peak calling mode of MACS2
 -se	settings; REQUIRED
 -sn    sample name (DB field: exomehg19.sample.name)
 -lf	<log file>; default: pipeline.log
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page

=head1 AUTHOR

Thomas Wieland

