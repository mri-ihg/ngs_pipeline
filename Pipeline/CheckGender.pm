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

use constant FACTOR => 1000000000;
use constant SRY    => "chrY:2654896-2655792";


my $prog_path = dirname( abs_path($0) );
require $prog_path . "/Utilities.pm";




################################################################################
sub calcCov {
	my $infile   = shift;	#if file: take regions from file; if region: take region; if empty: take SRY constant
	my $settings = shift;
	my $stats    = shift;	#if file: get reads from stats file; else: reads are given directly
	my $bam      = shift;
	my $factor   = shift;	#if empty: take constant FACTOR
	my $logfile  = shift;
	my $loglevel = shift;
	#my $logger	 = shift;
	
	if(!defined $factor || $factor == -1){
		$factor = FACTOR;
	}
	
	Utilities::initLogger( $logfile, $loglevel );
	my $logger = Utilities::getLogger();
	
	my @coverage;
	my $params   = Utilities::getParams();
	my $ref      = $params->{settings}->{$settings}->{reference};
	
	my $sam = Bio::DB::Sam->new(
		-fasta => $ref,
		-bam   => $bam
	);

	if(-e $infile){
		open(IN,$infile) or exit $logger->error("Can't open $infile!");	#infile given
		while(<IN>){
			chomp;
			push(@coverage,&getCoverage($_,$sam));
		}
	}elsif($infile ne ""){
		push(@coverage,&getCoverage($infile,$sam));		#region given directly
		
	}else{	# no infile given --> standard region
		push(@coverage,&getCoverage(SRY,$sam));
	}
	
	my $mean       = sum(@coverage)/@coverage;
	my $reads	   = $stats;
	$reads 		   = &getReads($stats,$logger) if -e $stats;	#if stats file given read it
	my $normalised = ($mean/$reads)*$factor; 
	
	return $normalised;

}

################################################################################
sub getCoverage{
	my $region = shift;
	my $sam    = shift;
	
	
	my $chr;
	my $start;
	my $end;
	
	if($region =~ /:/){ #UCSC format
		my $tmp;
		($chr,$tmp)   = split(":",$region);
		($start,$end) = split("-",$tmp);
	}else{
		($chr,$start,$end) = split("\t",$region);	#BED format
	}
	
	my ($coverage) = $sam->features(
		-type   => 'coverage',
		-seq_id => $chr,
		-start  => $start,
		-end    => $end
	);
		
	my @cov      = $coverage->coverage;
	my $range    = @cov;

	if ( $range == 0 ) {    #if no reads overlap this region @cov would be empty
		$range = $end - $start ;
		@cov[ 0 .. $range  ] = 0;
	}

	return @cov;	
	
}

################################################################################
sub getReads{
	my $stats  = shift;
	my $logger = shift;
	
	open(STATS,$stats) or exit $logger->error("Can't open $stats!");
	my $header  = <STATS>;
	my $data    = <STATS>;
	my @columns = split("\t",$data);
	return $columns[3];
}