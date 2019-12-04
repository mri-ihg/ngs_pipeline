#!/usr/bin/perl
use strict;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;
use File::Copy;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);


######################################################
# vcfsorter.pl
#
# Copyright (C) 2011 German Gaston Leparc
#
# sorts VCF by reference genome
#
# Adapted by Thomas Wieland, 24.07.2014
#
######################################################


my $prog_path = dirname(abs_path($0));
require $prog_path."/Utilities.pm";

my $outfile       = "";
my $help          = 0;
my $params        = Utilities::getParams();
my $logfile       = "pipeline.log";
my $loglevel      = "INFO";
my $settings      = "";
my $vcf_file      = "";
my $replace       = 0;

my $helptext      = 
"
This script sorts a VCF file according to the sequence dictionary of a reference genome.

-i	<infile.vcf>; default: take from stdin
-o	<outfile.vcf>; required if -r is not chosen
-se	<settings>; required
-r	replace infile.vcf by sorted file; a temporary file is created next to infile.vcf
-lf	log file; default: pipeline.log
-ll	log level: ERROR,INFO,DEBUG; default: INFO
-h	this help\n";


GetOptions(
"i=s"  => \$vcf_file,
"o=s"  => \$outfile, 
"se=s" => \$settings,
"r"    => \$replace,
"lf=s" => \$logfile,
"ll=s" => \$loglevel,
"h"    => \$help);



$outfile = dirname($vcf_file)."/tmp".time.".vcf" if $replace && $vcf_file ne "";





if ($help == 1  || $settings eq "") {
	print $helptext;exit(1);
}

Utilities::initLogger( $logfile, $loglevel );
my $logger = Utilities::getLogger();

if($vcf_file ne ""){
	close STDIN;
	$vcf_file =~ s/(.*\.gz)\s*$/gzip -dc < $1|/;
	
	open(STDIN,$vcf_file) or exit $logger->error("Can't open $vcf_file!");
	
	
}

if($outfile ne ""){
	close STDOUT;
	if($outfile =~ /\.gz$/){
		open(STDOUT,"| gzip -c > $outfile") or die "Can't open $outfile for writing!\n";
	}else{
		open(STDOUT,">",$outfile) or die "Can't open $outfile for writing!\n";
	}
	
}




my $dict_file = $params->{settings}->{$settings}->{reference};
$dict_file    =~ s/fa$/dict/;


#---------------------------------------- LOAD IN FASTA DICT INTO MEMORY
open(DICT,$dict_file) or exit $logger->error("Can't open $dict_file!");
my @contig_order;
my $c=0;
while(<DICT>)
{
if($_=~ /\@SQ/)
	{
	my ($contig) = $_ =~ /SN:(\S+)/;
	$contig_order[$c]=$contig;
	++$c; 
	#print $contig,"\n";
	}
}
close(DICT);

#---------------------------------------- PARSE VCF FILE & OUTPUT SORTED VCF

#open(VCF,$vcf_file) or exit $logger->error("Can't open $vcf_file!");

my %vcf_hash;
my $header;

while(<STDIN>)
{
if($_=~/^#/){ $header .= $_; } # store header and comment fields
chomp($_);

my @data = split(/\t/,$_);
my $contig = $data[0];
my $start = $data[1];
my $variant = $data[4]."to".$data[5];
my $line = $_;

#print $contig,":",$start,"\n";

$vcf_hash{$contig}{$start}{$variant}=$line;

}
close(STDIN);

#------------------ print out the VCF in the order of the reference genome

#print standard VCF header
print $header;

foreach my $contig (@contig_order) # sort by contig order
	{
	#print $contig,"\n";
	foreach my $start (sort {$a <=> $b} keys %{$vcf_hash{$contig}}) # sort numerically by coordinates
		{
		#print $start,"\n";
		foreach my $variant (keys %{$vcf_hash{$contig}{$start}}) # if overlapping mutation, print each variant
			{
			print $vcf_hash{$contig}{$start}{$variant},"\n";
			}	
		}
		
	}
	


move($outfile,$vcf_file) if $replace;