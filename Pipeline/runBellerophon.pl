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
my $logfile  	  = "SCREEN";
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


my $params 		  = Utilities::getParams();
my $bellerophon	  = $params->{programs}->{bellerophon}->{path};
my $dust_bed_file = $params->{settings}->{$settings}->{bellerophon}->{dust};

#generate blat.conf file
open BLAT,">".dirname($bamfile)."/blat.conf" or exit $logger->error("Can't open ".dirname($bamfile)."/blat.conf");
print BLAT "blat_port=".$params->{programs}->{bellerophon}->{blatconf}->{$settings}->{port}."\n";
print BLAT "blat_server=".$params->{programs}->{bellerophon}->{blatconf}->{$settings}->{server}."\n";
print BLAT "blat_path=".$params->{programs}->{bellerophon}->{blatconf}->{$settings}->{path}."\n";
print BLAT "minScore=".$params->{programs}->{bellerophon}->{blatconf}->{$settings}->{minscore}."\n";
print BLAT "minIdentity=".$params->{programs}->{bellerophon}->{blatconf}->{$settings}->{minidentity}."\n";
close BLAT;


#run bellerophon
my $command = "
export PERL5LIB=\$PERL5LIB:".dirname($bellerophon)." 

$bellerophon --input_bam $bamfile --blat_config ".dirname($bamfile)."/blat.conf --dust $dust_bed_file";

$logger->info("Run Bellerophon...");
$logger->debug("Command: ".$command);
system($command);






=head1 NAME

runBellerophon.pl

=head1 SYNOPSIS

 runBellerophon.pl -b merged.rmdup.bam -se hg19_plus

=head1 DESCRIPTION

This script is a wrapper script to run Bellerophon. It requires a Blat server running on ihgseq3,
as specified in the current.config.xml file.

=head1 OPTIONS

 -b	<bamfile> to run bellerophon on; required
 -se settings; required
 -lf	<log file>; default: pipeline.log
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page

=head1 AUTHOR

Thomas Wieland

=cut
