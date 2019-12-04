#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use DBI;
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);
use List::MoreUtils qw/ uniq /;
use Pod::Usage;


my $prog_path = dirname( abs_path($0) );
require $prog_path . "/Utilities.pm";

my $run = 1;

# database
my $dbh           = "";
my $sql           = "";
my $sth           = "";
my $logfile       = "SCREEN";
my $loglevel      = "INFO";


my $infile    = "";
my $targetfile= "";
my $outfile   = "";
my $help      = 0;
my $man		  = 0;


GetOptions(
	"i=s"  => \$infile,
	"t=s"  => \$targetfile,
	"o=s"  => \$outfile,
	"lf=s" => \$logfile,
	"ll=s" => \$loglevel,
	"h"    => \$help,
	"man"  => \$man
);

pod2usage( { -exitval => 0, -verbose => 1 } ) if $help;
pod2usage( { -exitval => 0, -verbose => 2 } ) if $man;
pod2usage( { -exitval => 1, -verbose => 1 } ) if $infile eq "" || $targetfile eq "";


my $params   = Utilities::getParams();
my $samtools = $params->{programs}->{samtools}->{path};


Utilities::initLogger( $logfile, $loglevel );
my $logger = Utilities::getLogger();

#open output file, if required
if($outfile ne ""){
	close STDOUT;
	if($outfile =~ /\.gz$/){
		open(STDOUT,"| bgzip -c > $outfile") or die "Can't open $outfile for writing!\n";
	}else{
		open(STDOUT,">",$outfile) or die "Can't open $outfile for writing!\n";
	}
	
}

my ($totalReads,$ontarget) = &calcOnTarget($infile,$targetfile);

my $targetRegion = `cat $targetfile | wc -l`;
chomp $targetRegion;
my $targetBases  = `awk '{sum+=(\$3-\$2)}END{print sum}' $targetfile`;
chomp $targetBases;

print "$targetRegion distinct regions targeted, total $targetBases bases\n";
print "Found $totalReads total reads\n";
print "$ontarget (".sprintf("%.2f",(($ontarget/$totalReads)*100) )."\%) on target\n";

###### sub #######
sub calcOnTarget {
	my $infile     = shift;
	my $targetfile = shift;

	#calculate total reads on chromosomes that contain at least a single target region
	my $totalReads = 0;
	my $ontarget   = 0;
	
	open CHRS, "cut -f1 $targetfile | sort | uniq|" or exit $logger->error("Can't execute command: cut -f1 $targetfile | sort | uniq|!");
	$logger->info("Getting total number of reads...");
	while(<CHRS>){
		chomp;
		$logger->debug("Counting reads on chromosome $_");
		my $reads = `$samtools view -F 1024 -c $infile $_`;
		chomp $reads;
		$totalReads += $reads;
		
		$reads = `$samtools view -F 1024 -L $targetfile -c $infile $_`;
		chomp $reads;
		$ontarget   += $reads;
		
		
	}
	close CHRS;


	return ($totalReads,$ontarget);
	
}


=head1 NAME

calcOnTarget.pl

=head1 SYNOPSIS

calcOnTarget.pl -i merged.rmdup.bam -t targets.bed -o ontarget.out   

=head1 DESCRIPTION

This script calculates the number and proportion of reads that are on target. To calculate the 
proportion it calculates the total number of reads on all chromosomes that harbor at least one
target region.
NOTE: If the BAM file contains reads marked as duplicates, these reads are discarded.

=head1 OPTIONS

 -i	<input.bam> BAM file to calculate on target proportion; REQUIRED
 -t	<targets.bed> BED file containing the target regions; REQUIRED
 -o	<ontarget.out> file to output reads and proportion; default: stdout
 -lf	log file; default: print to screen
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page

=head1 AUTHOR

Thomas Wieland

=cut