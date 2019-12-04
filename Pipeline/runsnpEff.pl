#!/usr/bin/perl

use Getopt::Long qw(:config no_ignore_case);
use Pod::Usage;
use File::Basename;
use Cwd qw(abs_path);

my $prog_path = dirname( abs_path($0) );
require $prog_path."/Utilities.pm";

my $help        = 0;
my $man         = 0;
my $outfile     = "";
my $infile      = "";
my $settings    = "";
my $logfile  	= "pipeline.log";
my $loglevel 	= "INFO";

GetOptions(
"i=s" => \$infile,
"o=s" => \$outfile,
"se=s"=> \$settings,
"lf=s" => \$logfile,
"ll=s" => \$loglevel,
"h"   => \$help,
"man" => \$man
);

pod2usage( {-exitval => 1  ,-verbose => 1} ) if $infile eq "" || $settings eq "";
pod2usage( {-exitval => 0  ,-verbose => 1} ) if $help;
pod2usage( {-exitval => 0  ,-verbose => 2} ) if $man;

if($outfile eq ""){
	$outfile = $infile;
	$outfile =~ s/\.vcf$/\.snpEff\.vcf/;
}

chdir(dirname($outfile));
$outfile = basename($outfile);

Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();

my $params  = Utilities::getParams();
my $snpEff  = $params->{programs}->{snpeff}->{path};
my $java    = $params->{programs}->{java}->{path};
my $datadir = $params->{settings}->{$settings}->{snpeff}->{datadir};
my $dbname  = $params->{settings}->{$settings}->{snpeff}->{dbname};

my $command = "$java -XX:ParallelGCThreads=1 -jar $snpEff -dataDir $datadir -hgvs $dbname $infile > $outfile";
$logger->info("Running snpEff...");
$logger->debug($command);
system($command);



=head1 NAME

 runsnpEff.pl

=head1 SYNOPSIS

 runsnpEff.pl -i ontarget.vcf -se hg19_test

=head1 DESCRIPTION

This script is a wrapper script for the snpEff annotation tool.

=head1 OPTIONS

 -i	input file; REQUIRED
 -o	output file; default: infile.snpEff.vcf
 -se	settings; REQUIRED
 -lf	<log file>; default: pipeline.log
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	show help text
 -man	show man page

=head1 AUTHOR

Thomas Wieland

=cut