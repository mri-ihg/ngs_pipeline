#!/usr/bin/perl

use Getopt::Long qw(:config no_ignore_case);
use Pod::Usage;
use Cwd qw(abs_path);
use File::Basename;

my $prog_path = dirname( abs_path($0) );
require $prog_path . "/Utilities.pm";

my $help     = 0;
my $man      = 0;
my $input    = "";
my $output   = "";

GetOptions(
"i=s" => \$input,
"o=s" => \$output,
"h"   => \$help,
"man" => \$man
);

pod2usage( {-exitval => 0  ,-verbose => 1} ) if $help;
pod2usage( {-exitval => 0  ,-verbose => 2} ) if $man;


#open input/output files, if specified
if($input ne ""){
	close STDIN;
	$input =~ s/(.*\.gz)\s*$/gzip -dc < $1|/;
	
	open(STDIN,$input) or die "Can't open $input!\n";
	
	
}

if($output ne ""){
	close STDOUT;
	if($output =~ /\.gz$/){
		open(STDOUT,"| gzip -c > $output") or die "Can't open $output for writing!\n";
	}else{
		open(STDOUT,">",$output) or die "Can't open $output for writing!\n";
	}
	
}

while(<STDIN>){
	if($_ =~ /^#/){		#print comments
		print $_;
		next;
	}
	
	chomp;
	my @columns = split("\t");
	
	@columns[3..4] = Utilities::getMinimalIndel($columns[3],$columns[4]) if length($columns[3])>1 && length($columns[4])>1;	#calc minimal indel
	
	print join("\t",@columns)."\n";
	
}

close STDIN;
close STDOUT;




=head1 NAME

 getMinimalIndels.pl

=head1 SYNOPSIS

 getMinimalIndels.pl -i input.vcf -o output.vcf

=head1 DESCRIPTION

This script takes a VCF file and calculates the minimial indels (i.e. without surrounding repeats).
This indel format is used for instance by samtools. The script only uses the 4th and 5th column of 
the VCF file.

=head1 OPTIONS

 -i	<input.vcf>; default: take from STDIN
 -o	<output.vcf>; default: print to STDOUT
 -h	show help text
 -man	show man page

=head1 AUTHOR

Thomas Wieland

=cut