#!/usr/bin/perl

use strict;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);
use Pod::Usage;
umask(002);


my $prog_path = dirname(abs_path($0));
require $prog_path."/Utilities.pm";

my $outfile    = "";
my $help       = 0;
my $man		   = 0;
my $params     = Utilities::getParams();
my $logfile    = "pipeline.log";
my $loglevel   = "INFO";
my $settings   = "";
my $indir      = ".";
my $ending	   = "vcf";
my $numbered   = 0;


GetOptions(
"i=s"  => \$indir,
"e=s"  => \$ending,
"o=s"  => \$outfile, 
"se=s" => \$settings,
"n"    => \$numbered,
"lf=s" => \$logfile,
"ll=s" => \$loglevel,
"man"  => \$man,
"h"    => \$help);


pod2usage( {-exitval => 0  ,-verbose => 1} ) if $help;
pod2usage( {-exitval => 0  ,-verbose => 2} ) if $man;
pod2usage( {-exitval => 1  ,-verbose => 1} ) if $outfile eq "" || ($numbered == 0 && $settings eq "");

Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();

my $bgzip = $params->{programs}->{bgzip}->{path};
my $tabix = $params->{programs}->{tabix}->{path};
my $cat = "cat";
$cat    = "$bgzip -cd" if $ending =~ /\.gz$/;

my $isZipped = 0;
if($outfile =~ /\.gz$/){
	$isZipped = 1;
	$outfile =~ s/\.gz$//;
	
}

#first get header file from arbitrary input file
my $command = "ls $indir/*.$ending | head -n 1 | xargs -i $cat {} | awk '{if(\$1 ~ /^#/){print}else{exit 0} }' > $outfile ";
#print $command."\n";

$logger->debug($command);
system($command);

if($numbered){
	#get number of files
	my $fileNumber = `ls $indir/*.$ending | wc -l`;
	chomp $fileNumber;
	for (my $i = 1; $i<= $fileNumber; $i++){
		$command = "find $indir -name \"$i.*.$ending\" | xargs -i $cat {} | awk '{if(\$1 !~ /^#/){print}}' >> $outfile";
		$logger->debug($command);
		system($command);
	}
	
}else{
	my $dict = $params->{settings}->{$settings}->{reference};
	$dict    =~ s/fa$/dict/;
	$dict    =~ s/fasta$/dict/;
	
	open DICT, "grep \"\@SQ\" $dict |" || exit $logger->error("Can't open $dict!");
	while(<DICT>){
		if($_ =~ /SN:(.+?)\s+/) {
			$command = "find $indir -name \"$1.*.$ending\" | xargs -i $cat {} | awk '{if(\$1 !~ /^#/){print}}' >> $outfile";
			
			#print $command."\n";
			$logger->debug($command);
			system($command);
		}
	}
	close DICT;

}

if($isZipped){
	$command = "$bgzip -f $outfile";
	$logger->debug($command);
	system($command);

	$command = "$tabix -fp vcf $outfile.gz";
	$logger->debug($command);
	system($command);
}




=head1 NAME

concatVCF.pl

=head1 SYNOPSIS

 concatVCF.pl -i /path/to/folder/ -e vcf -se hg19_test -o out.vcf

=head1 DESCRIPTION

This is script concats VCF/gVCF files according to the sorting order of a sequence dictionary file from a given ref genome.
Files have to start with chromosomename. or they are numbered 1-n. if "-n" is specified. 
They have to end with file ending (-e).

=head1 OPTIONS

 -i	/input/folder; folder to look for files; default: .
 -e	file ending; default: vcf
 -o	<outfile>; required
 -se	<settings>; required
 -n	files are numbered instead of named by chromosome
 -lf	<log file>; default: pipeline.log
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page

=head1 AUTHOR

Thomas Wieland

=cut
 