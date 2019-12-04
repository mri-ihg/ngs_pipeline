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
my $outdir  = "";
my $help	= 0;
my $man     = 0;
my $logfile  	  = "pipeline.log";
my $loglevel 	  = "INFO";
my $bins    = 500;
my $settings  = "";


GetOptions(
"b=s"  => \$bamfile,
"o=s"  => \$outdir,
"i=s"  => \$bins,
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

$outdir = dirname($bamfile) if $outdir eq "";

my $params 		  = Utilities::getParams();
my $lumpy  		  = $params->{programs}->{lumpy}->{path};
my $sam    		  = $params->{programs}->{samtools}->{path};
my $sammaxmem     = $params->{programs}->{samtools}->{maxmem};
my $normalChrs 	  = $params->{settings}->{$settings}->{normalchromosomes};

my $discordantBAM = $bamfile;
$discordantBAM    =~ s/bam$/discordant/;

my $splitBAM      = $bamfile;
$splitBAM         =~ s/bam$/splitters/;

# Extract the discordant paired-end alignments.
my $command = "$sam view -h -L $normalChrs -b -F 1294 $bamfile | $sam sort -m $sammaxmem - $discordantBAM";
$logger->info("Extracting discordant readpairs...");
$logger->debug("Command: ".$command);
system($command);
$discordantBAM .= ".bam";

# Extract the split-read alignments
$command = "$sam view -h -L $normalChrs -h $bamfile \\
| $lumpy/scripts/extractSplitReads_BwaMem -i stdin \\
| $sam view -Sb - \\
| $sam sort -m $sammaxmem - $splitBAM";
$logger->info("Extracting split readpairs...");
$logger->debug("Command: ".$command);
system($command);
$splitBAM .= ".bam";


#run lumpy
$command = "$lumpy/bin/lumpyexpress \\
-B $bamfile \\
-S $splitBAM \\
-D $discordantBAM \\
-o $outdir/lumpy.vcf";
$logger->info("Run lumpy-sv...");
$logger->debug("Command: ".$command);
system($command);








=head1 NAME

runLumpy.pl

=head1 SYNOPSIS

 runLumpy.pl -b merged.rmdup.bam -se hg19_plus

=head1 DESCRIPTION

This script is a wrapper script to run LUMPY-SV (https://github.com/arq5x/lumpy-sv).
The input BAM file (-b) must be aligned using BWA MEM, because LUMPY-SV requires split reads.
Only tested for whole genome data. It calls variants only on chromosomes 1-22,X,Y and M.

=head1 OPTIONS

 -b	<bamfile> to run cnvnator on; required
 -se settings; required
 -o	outdir; default: directory where bamfile lies in
 -lf	<log file>; default: pipeline.log
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page

=head1 AUTHOR

Thomas Wieland

=cut