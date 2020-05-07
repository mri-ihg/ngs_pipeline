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
my $settings  = "";
my $threads       = -1;



GetOptions(
"b=s"  => \$bamfile,
"o=s"  => \$outdir,
"t=s"  => \$threads,
"se=s" => \$settings,
"lf=s" => \$logfile,
"ll=s" => \$loglevel,
"man"  => \$man,
"h"	   => \$help
);

pod2usage( {-exitval => 0  ,-verbose => 1} ) if $help;
pod2usage( {-exitval => 0  ,-verbose => 2} ) if $man;
pod2usage( {-exitval => 1  ,-verbose => 1} ) if $bamfile eq "" || $outdir eq "" || $settings eq "" ;

Utilities::initLogger($logfile,$loglevel);

my $logger = Utilities::getLogger();

$outdir = dirname($bamfile) if $outdir eq "";

my $params 		  = Utilities::getParams();
my $whamg  		  = $params->{programs}->{whamg}->{path};
my $whamg_filterprogram   = $params->{programs}->{whamg}->{filterprogram};
my $whamg_annotationprog  = $params->{programs}->{whamg}->{annotationprogram};
my $normalChrs  	  = $params->{settings}->{$settings}->{normalchromosomes};
my $ref		 	  = $params->{settings}->{$settings}->{reference};


if($threads == -1){
	$threads = 1;
	if($ENV{NSLOTS}){		#get number of slots given by SGE
		$threads = $ENV{NSLOTS};
		
	}
}

my $command = "$whamg \\
-x $threads \\
-c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY \\
-a $ref \\
-f $bamfile \\
> $outdir/whamg.vcf \\
2> $outdir/whamg.err";


&Utilities::executeCommand($command,"Running whamg", $logger);


$command = "cat $outdir/whamg.vcf | \\
/usr/bin/perl $whamg_filterprogram \\
> $outdir/whamg.filtered.vcf";

&Utilities::executeCommand($command,"Filter whamg results", $logger);


$command = "python $whamg_annotationprog $outdir/whamg.filtered.vcf > $outdir/whamg.filtered.annotated.vcf";

&Utilities::executeCommand($command,"Annotate whamg results", $logger);




=head1 NAME

runWhamg.pl

=head1 SYNOPSIS

 runWhamg.pl -b merged.rmdup.bam -o whamgdir/ -se hg19_plus

=head1 DESCRIPTION

This script is a wrapper script to run Whamg
Only tested for whole genome data. It calls variants only on chromosomes 1-22,X,Y and M.

=head1 OPTIONS

 -b	<bamfile> to run cnvnator on; required
 -o	outdir; required
 -se settings; required
 -t	<threads>, default: take number of slots from SGE environment variable NSLOTS or 1 if no SGE job
 -lf	<log file>; default: pipeline.log
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page

=head1 AUTHOR

Thomas Wieland

=cut
