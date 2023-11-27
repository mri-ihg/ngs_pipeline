#!/usr/bin/perl
package CheckGender;

use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use DBI;
use POSIX;
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);
use Bio::DB::Sam;
use List::Util qw(sum);

my $prog_path = dirname( abs_path($0) );
require $prog_path . "/Utilities.pm";
require $prog_path . "/CheckGender.pm";

#use lib '/home/eck/work/programming/perlScripts/modules';
#use frequentSubroutines;

#
#calculate average coverage of targeted bases of enrichment experiment
#

my $infile  = "";

#my $refFa = "";
my $settings = "";
my $bam      = "";
my $logfile  = "pipeline.log";
my $loglevel = "INFO";
my $stats	 = "";
my $factor	 = -1;

#program flags
my $help = 0;

GetOptions(
	"t=s"  => \$infile,
	"se=s" => \$settings,
	"s=s"  => \$stats,
	"f=s"  => \$factor,
	"b=s"  => \$bam,
	"lf=s" => \$logfile,
	"ll=s" => \$loglevel,
	"h"    => \$help,
);


if($stats eq ""){
	$stats = "summary.stats.csv";
	if($bam ne ""){
		my $outdir = dirname($bam);
		$stats = $outdir."/".$stats;
	}
}

if ( $help == 1 || $bam eq "" || $settings eq "") {
	print
"\n
Calculates the coverage of the target region and normalises it with the number of total mapped reads.

Required:
-b	<bam file> to check
-se	settings

Optional:
-t	<targetfile> targeted regions (in UCSC or BED format); if empty: check default region on chrY (SRY gene): SRY
-s	<summary.stats.csv> file to take reads on target from; default: $stats
-f	<factor> factor to multiply normalised value; default: FACTOR
-lf	log file; default: pipeline.log
-ll	log level: ERROR,INFO,DEBUG; default: INFO
-h	this help\n\n";
	exit(1);
}

my $avg = CheckGender::calcCov($infile,$settings,$stats,$bam,$factor,$logfile,$loglevel);

print "Normalised Average Coverage: $avg \n";


