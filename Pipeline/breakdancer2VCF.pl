#!/usr/bin/perl

use strict;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;
use Cwd qw(abs_path);
#use Log::Log4perl qw(get_logger :levels);
use IO::File;
use Pod::Usage;
umask(002);

my $prog_path = dirname( abs_path($0) );
require $prog_path."/Utilities.pm";


my $infile  = "";
my $outfile = "";
my $help	= 0;
my $man     = 0;
my $logfile  	  = "SCREEN";
my $loglevel 	  = "INFO";



GetOptions(
"i=s"  => \$infile,
"o=s"  => \$outfile,
"lf=s" => \$logfile,
"ll=s" => \$loglevel,
"man"  => \$man,
"h"	   => \$help
);

pod2usage( {-exitval => 0  ,-verbose => 1} ) if $help;
pod2usage( {-exitval => 0  ,-verbose => 2} ) if $man;

Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();


if($infile ne ""){
	close STDIN;
	$infile =~ s/(.*\.gz)\s*$/gzip -dc < $1|/;
	
	open(STDIN,$infile) or die "Can't open $infile!\n";
	
	
}

if($outfile ne ""){
	close STDOUT;
	if($outfile =~ /\.gz$/){
		open(STDOUT,"| gzip -c > $outfile") or die "Can't open $outfile for writing!\n";
	}else{
		open(STDOUT,">",$outfile) or die "Can't open $outfile for writing!\n";
	}
	
}



#print VCF header
print "##fileformat=VCFv4.2
##source=breakdancer
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total number of discordant read pairs.\">
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">
##INFO=<ID=OR1,Number=1,Type=String,Description=\"Read orientations at breakpoint 1\">
##INFO=<ID=OR2,Number=1,Type=String,Description=\"Read orientations at breakpoint 2\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Number of discordant read pairs.\">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT";

while(<STDIN>){
	chomp;
	if($_ =~ /^#/){
		if($_ eq "#Library Statistics:"){
			my $line = <STDIN>;					#get sample name out of library name
			chomp $line;
			my @columns = split("\t",$line);
			foreach my $currCol(@columns){
				my ($key,$value) = split(":",$currCol);
				if($key eq "library"){
					my ($sample,$lib) = split("_",$value);
					print "\t$sample\n";
					last;
				}
			}
		}
	}else{
		my @columns = split("\t");
		next unless ($columns[6] eq "DEL" || $columns[6] eq "INS" || $columns[6] eq "INV");
		
		$columns[7] *= -1 if ($columns[6] eq "DEL" && $columns[7] > 0)  || ($columns[6] ne "DEL" && $columns[7] < 0);
		
		print "$columns[0]\t$columns[1]\t.\tN\t<$columns[6]>\t$columns[8]\t.\tEND=$columns[4];SVTYPE=$columns[6];SVLEN=$columns[7];DP=$columns[9];OR1=$columns[2];OR2=$columns[5]\tGT:GQ:DP\t./.:$columns[8]:$columns[9]\n";
	}
}


=head1 NAME

breakdancer2VCF.pl

=head1 SYNOPSIS

 breakdancer2VCF.pl -i breakdancer.out -o breakdancer.vcf

=head1 DESCRIPTION

This script converts breakdancer output files to VCF file format. Currently it only converts deletions, insertions
and inversions (so no translocations).

=head1 OPTIONS

 -i	<breakdancer.out> breakdancer output file; default: STDIN
 -o	<breakdancer.vcf> VCF output file; default: STDIN
 -lf	<log file>; default: print to screen
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page

=head1 AUTHOR

Thomas Wieland

=cut