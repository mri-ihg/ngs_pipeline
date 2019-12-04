#!/usr/bin/perl

use strict;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);
use Sys::Hostname;
use Pod::Usage;




my $prog_path = dirname(abs_path($0));
require $prog_path."/Utilities.pm";


my $help       = 0;
my $man		   = 0;
my $params     = Utilities::getParams();
my $logfile    = "SCREEN";
my $loglevel   = "INFO";
my $settings   = "default";
my $sample     = "";
my $outdir     = ".";
my $R		   = "R";





GetOptions(
"se=s" => \$settings,
"s=s"  => \$sample,
"o=s"  => \$outdir,
"lf=s" => \$logfile,
"ll=s" => \$loglevel,
"man"  => \$man,
"h"    => \$help);

pod2usage( { -exitval => 0, -verbose => 1 } ) if $help;
pod2usage( { -exitval => 0, -verbose => 2 } ) if $man;
pod2usage( { -exitval => 1, -verbose => 1 } ) if $sample eq "";

Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();



my $ref      = $params->{settings}->{$settings}->{reference};
my $script   = $params->{programs}->{homozygosity}->{script};
my $host     = $params->{settings}->{$settings}->{exomedb}->{host};

if($host eq hostname ){		#set hostname to localhost if this script is run on the host where the database is located; otherwise the R script might have problems to connect
	$host = "localhost";
}
									

my $database = $params->{settings}->{$settings}->{exomedb}->{database};
my $user     = $params->{settings}->{$settings}->{exomedb}->{user};
my $password = $params->{settings}->{$settings}->{exomedb}->{password};
my $coredb   = $params->{coredb}->{database};


my $command  = "$R --no-restore --no-save --args $user $password $database $host $sample $coredb <$prog_path/$script > $outdir/homozygosity.log 2>&1";
$logger->debug($command);
$logger->info("Calculating stretches of homozygosity...");
system($command);



=head1 NAME

homozygosity.pl

=head1 SYNOPSIS

homozygosity.pl -s SAMPLENAME

=head1 DESCRIPTION

This script is a wrapper script for homozygosity.R which calculates stretches of homozygosity and inserts
them directly into the database.

=head1 OPTIONS

 -s	<samplename>; REQUIRED
 -o	<outdir>; directory where log files are saved; default: .
 -se	name of the settings in the current.config.xml file that holds path to reference genome, 
 	to the annotation file and to possible additional annotation files; use default settings if nothing is given
 -lf	log file; default: print to screen
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page
 

=head1 AUTHOR

Thomas Wieland

=cut
