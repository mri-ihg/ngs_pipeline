#!/usr/bin/perl

use strict;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;
use File::Copy;
use Cwd qw(abs_path);
use Data::Dumper;
use IO::File;
use Pod::Usage;
umask(002);

my $prog_path = dirname( abs_path($0) );
require $prog_path."/Utilities.pm";


my $bamfile = "";
my $help	= 0;
my $man     = 0;
my $logfile  	  = "pipeline.log";
my $loglevel 	  = "INFO";
my $settings  = "";


GetOptions(
"b=s"  => \$bamfile,
"se=s" => \$settings,
"lf=s" => \$logfile,
"ll=s" => \$loglevel,
"man"  => \$man,
"h"	   => \$help
);

pod2usage( {-exitval => 0  ,-verbose => 1} ) if $help;
pod2usage( {-exitval => 0  ,-verbose => 2} ) if $man;
pod2usage( {-exitval => 1  ,-verbose => 1} ) if $bamfile eq "" || $settings eq "" ;

Utilities::initLogger($logfile,$loglevel);

my $logger = Utilities::getLogger();

my $outfile = $bamfile . ".tdf";

my $params 		  = Utilities::getParams();
my $IGVtools      = $params->{programs}->{IGVtools}->{path};
my $java          = $params->{programs}->{java}->{path};
my $ref      	  = $params->{settings}->{$settings}->{reference};



# Extract the discordant paired-end alignments.
my $command = "$java -Xmx1500m -XX:ParallelGCThreads=1 -jar $IGVtools count -z 5 -w 25 $bamfile $outfile $ref";
$logger->info("Generate coverage track for IGV...");
$logger->debug("Command: ".$command);
system($command);








=head1 NAME

runIGVtools.pl

=head1 SYNOPSIS

 runIGVtools.pl -b merged.rmdup.bam -se hg19_plus

=head1 DESCRIPTION

This script is a wrapper script to run IGVtools in order to generate coverage track in IGV.

=head1 OPTIONS

 -b	<bamfile> to run IGVtools on; required
 -se settings; required
 -lf	<log file>; default: pipeline.log
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page

=head1 AUTHOR

Thomas Wieland

=cut