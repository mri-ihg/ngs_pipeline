#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use File::Basename;
use Cwd qw(abs_path);
umask(002);

#include Utilities.pm
my $prog_path = dirname(abs_path($0));
require $prog_path . "/Utilities.pm";

my $help		= "";
my $logfile  = "rnaseq.pipeline.log";
my $loglevel = "INFO";
my $outputFile = "";
my $caseFiles = "";
my $controlFiles = "";
my $countFileName = "";
my $conditionFileName = "";
my $outputFolder = "";
my %countFiles = ();
my $fileExtension = "htseqcounts";
my %fp;

my $helptext = "
-as <case file>			a file where all case files are located; files must be output from HTSeq, and have the same structure
-us <control file>		equivalent to case file for the controls
-cfn <count file location>	name the resulting count file
-condfn <condition file name>	name the resulting condition file
-o <output folder>
-fe <file extension>		the file extionsion of the files to merge (default: $fileExtension)
-h				show help
";

GetOptions(
	"as=s" => \$caseFiles,
	"us=s" => \$controlFiles,
	"cfn=s" => \$countFileName,
	"condfn=s" => \$conditionFileName,
	"o=s" => \$outputFolder,
	"fe=s" => \$fileExtension,
	"h" => \$help
);

if ($help) {
	print $helptext;
	exit(1);
}

unless ($caseFiles && $controlFiles && $countFileName && $conditionFileName) {
	print STDERR "not all required parameters specified\n";
	exit(11);
}

&Utilities::initLogger($outputFolder."/".$logfile, $loglevel);
my $logger = &Utilities::getLogger();

&fillCountFileHash($caseFiles, \%countFiles, 1);	#1 for cases and 0 for controls
&fillCountFileHash($controlFiles, \%countFiles, 0);

if ((keys %countFiles) < 1) {
	$logger->error("No files found from the input file $caseFiles and/or $controlFiles");
	exit(11);
} else {
	$logger->info("Found " . (keys %countFiles) . " files to merge");
}

#check if all files have the same line count
$logger->info("Calculating file line counts");
my $lastLineCount = "-1";
my $actLineCount = "-1";
my $i=0;
foreach my $key ( keys %countFiles )
{
	open (IN, "wc -l $countFiles{$key}{'file'} |") || exit $logger->error("Error executing wc -l $countFiles{$key}{'file'}");
	$actLineCount = <IN>;
	$actLineCount =~ s/^([0-9]*)(\s)(.*)/$1/g;
	if (($i != 0) && ($actLineCount ne $lastLineCount)) {
		$logger->error("All files must have same line count in file $countFiles{$key}{'file'}! Use of different annotation files (.gtf)?");
		exit (1);
	}
	$lastLineCount = $actLineCount;
	close(IN);
	$i++;
}
$logger->info("Lenght check successfull! All files have the same line count: $lastLineCount");

#open the files and fill the filehandle hash 
$logger->info("Opening ".(keys %countFiles)." files");
my $countTableHeaderLine = "transcriptid";
foreach my $key ( sort keys %countFiles ) {
	open($fp{$key}, "$countFiles{$key}{'file'}")|| exit $logger->error("Can't open file $countFiles{$key}{'file'}");
	$countTableHeaderLine .= "\t$key";
}

open CF, ">$outputFolder/$countFileName" || exit $logger->error("Couldn't open $outputFolder/$countFileName");
print CF "$countTableHeaderLine\n";

#use first file the synchronize all others
my $guideKey = (sort keys %fp)[0];
my $guideFP = $fp{$guideKey};
delete $fp{$guideKey};

my $line = "";
my $countTableLine = "";
my $lineCount = 0;
my $tmpLine = "";

$logger->info("Start merging files");
while (<$guideFP>) {
	$lineCount++;
	my @columns = split("\t", $_);
	#$columns[1] =~ s/\n$//;
	chomp(@columns);
	#check if second column is a number!!
	if ($columns[1] !~ m/^[0-9\.]+$/) {
		$logger->info("Input files must be in the format: TRANSCRIPTID		212! Found: \"@columns\"\n");
		exit(1);
	}
	$countTableLine .= "\"$columns[0]\"\t$columns[1]";
	
	foreach my $key (sort keys %fp) {
		$tmpLine = readline($fp{$key});
		$line .= "\t" . (split("\t", $tmpLine))[1];
		$line =~ s/\n$//;
		$countTableLine .= "\t" . (split("\t", $tmpLine))[1];
		$countTableLine =~ s/\n$//;
		#check if all Transcript id's are equal
		if ($columns[0] ne (split("\t", $tmpLine))[0]) {
			$logger->error("Transcript id of sample $guideKey and sample $key in line $lineCount are not equal ($columns[0] != " . (split("\t", $tmpLine))[0]);
			exit (1);
		}								
	}
	$countTableLine .= "\n";
	print CF "$countTableLine";
	$countTableLine = "";
}

#close the files
foreach my $key (keys %fp) {
	close($fp{$key});
}
close(CF);

if ($fileExtension =~ /htseqcounts/) {
	#write the condition file
	open COND, ">$outputFolder/$conditionFileName" || exit $logger->error("Couldn't open $outputFolder/$conditionFileName");
	foreach my $key (sort keys %countFiles) {
		print COND "$key\t" . $countFiles{$key}{"condition"} . "\n";
	}
	close COND;
}

$logger->info("Finished successful");


sub fillCountFileHash {
	my $f = shift;
	my $s = shift;
	my $class = shift;
	
	open FILE, "$f";
	while (<FILE>) {
		chomp($_);
		my ($sId, $fLocation) = split("\t");
		unless (-e $fLocation) {
			$logger->error("File \"$fLocation\" doesn't exit");
			exit(11);
		}
		$fLocation .= ".". $fileExtension;		#htseq file extension which will be assigned in htseqCount.pl
		my %h = ("condition" => $class, "file" => $fLocation);
		$s->{$sId} = {%h};
	}
	close(FILE);
}