#!/usr/bin/perl

use strict;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);
use Pod::Usage;
use DBI;
use Scalar::Util qw(looks_like_number);

#include Utilities.pm
my $prog_path = dirname(abs_path($0));
require $prog_path."/Utilities.pm";


my $outprefix     = "";
my $help          = 0;
my $settings      = "";

my $folder        = "";

my $logfile  	  = "pipeline.log";
my $loglevel 	  = "INFO";

my $readGroupLib  = "";
my $man			  = 0;

my $helptext      = 
"\n";


GetOptions(
	
	"folder=s"=>\$folder,
	"man"  => \$man, 
	"h"    => \$help, 
	"lf=s" => \$logfile,
	"ll=s" => \$loglevel, 
	"se=s" => \$settings
);
	
pod2usage( {-exitval => 0  ,-verbose => 1} ) if $help;
pod2usage( {-exitval => 0  ,-verbose => 2} ) if $man;
pod2usage( {-exitval => 1  ,-verbose => 1} ) if $settings eq "";


Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();
my $params = Utilities::getParams();

my $inputinfiles  = "inputinfiles";
my $bwa           = $params->{programs}->{bwa}->{path};		
my $bedtools      = $params->{programs}->{bedtools}->{path};			
my $bwareadtrim   = $params->{programs}->{bwa}->{readtrim};		
my $samtools      = $params->{programs}->{samtools}->{path};		
my $sammaxmem     = $params->{programs}->{samtools}->{maxmem};
my $cutadapt	  = $params->{programs}->{cutadapt}->{path};
my $java          = $params->{programs}->{java}->{path};
my $picard        = $params->{programs}->{picard}->{path};

my $ref           = $params->{settings}->{$settings}->{reference};
my $phiX          = $params->{settings}->{$settings}->{phiX} if defined $params->{settings}->{$settings}->{phiX}; 
my $cutadapt_adapter = $params->{programs}->{cutadapt}->{adaptersequence};

my $coredb                    = $params->{coredb}->{database};
my $solexadb                  = $params->{solexadb}->{database};
my $optduptable               = $params->{solexadb}->{opticalduplicatestable};


# Generate per lane metrics

my @bams=glob($folder."/[0-9][0-9][0-9][0-9][0-9][0-9]\_*_[0-9][0-9][0-9][0-9]_*_s_[0-9]_001*bam");

foreach my $bamfile (@bams)
{
	print $bamfile."\n";
	
	if( $bamfile =~ /[0-9][0-9][0-9][0-9][0-9][0-9]\_.*_[0-9][0-9][0-9][0-9]_(.*)_s_([0-9])_001/ )
	{
		my ($flowcell,$lane)=($1,$2);
		$flowcell =~ s/^[AB]//g;
		
		# In case that procedure is executed out of that script get LIB from bam
		$readGroupLib=`$samtools view $bamfile -H | grep "\@RG" | awk -F"LB:" '{print \$2}' | awk -F"\t" '{print \$1}'`;
		chomp($readGroupLib);		
		
		print $readGroupLib."\n";
		
		my $dbh = Utilities::connectExomeDB($settings);
		
		# Run MarkDup, discard output
		$logger->info("Running Picard MarkDuplicates: $java -Xmx6g -XX:ParallelGCThreads=1 -jar $picard MarkDuplicates INPUT=$bamfile OUTPUT=/dev/null M=/dev/stdout REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT 2>/dev/null");
		my ($dups,$optdups) = split(",", `$java -Xmx6g -XX:ParallelGCThreads=1 -jar $picard MarkDuplicates INPUT=$bamfile OUTPUT=/dev/null M=/dev/stdout REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT 2>/dev/null | grep -v "#" | grep -v UNPAIRED_READS_EXAMINED | grep $readGroupLib | awk '{print (2*\$6)/((2*\$3)+\$2)","(2*\$7)/((2*\$3)+\$2)}'`);
	
		# If is numeric
		if ( looks_like_number($optdups) )
		{
			# Delete Old entry 
			my $sql = "DELETE from solexa.opticalduplicates where rname=\"$flowcell\" and lane=\"$lane\" and lname=\"$readGroupLib\";";
			my $sth = $dbh->prepare($sql) || $logger->error("Can't prepare statement: $DBI::errstr");
			$sth->execute() || $logger->error($DBI::errstr);
		
			my $sql = "INSERT INTO solexa.opticalduplicates (rname, lane, lname, duplicates, opticalduplicates) values ( \"$flowcell\", $lane, \"$readGroupLib\", $dups, $optdups);";
			$sth = $dbh->prepare($sql) || $logger->error("Can't prepare statement: $DBI::errstr");
			$sth->execute() || $logger->error($DBI::errstr);
		
			print "\n\nFC: $flowcell LANE: $lane LIB: $readGroupLib\n\n"; 
			print "$optdups\n";
			
		}
		
	}
	
}