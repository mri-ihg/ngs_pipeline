#!/usr/bin/perl

use strict;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);
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
my $idsample   = 0;
my $vcf		   = "";
my $outfile    = "";
my $insertdb   = 0;
my $recombrate = "1e-8"; # Can be a string

GetOptions(
"se=s" => \$settings,
"s=s"  => \$sample,
"i=s"  => \$vcf,
"o=s"  => \$outfile,
"insert" => \$insertdb,
"recombrate=s" => \$recombrate,
"lf=s" => \$logfile,
"ll=s" => \$loglevel,
"man"  => \$man,
"h"    => \$help);

pod2usage( { -exitval => 0, -verbose => 1 } ) if $help;
pod2usage( { -exitval => 0, -verbose => 2 } ) if $man;
pod2usage( { -exitval => 1, -verbose => 1 } ) if ( $sample eq "" && $insertdb ) || ( $vcf eq "" ) || ( $outfile eq "" );

Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();


# Tools
my $bcftools= $params->{programs}->{bcftools}->{path};
my $bgzip=    $params->{programs}->{bgzip}->{path};
my $tabix=    $params->{programs}->{tabix}->{path};

	$tabix = Utilities::getProgram("tabix");
	

# Config
my $ref      = $params->{settings}->{$settings}->{reference};
my $gnomADfreqsROH = $params->{settings}->{$settings}->{gnomADfreqsROH};


# Get DBs and connect:					
my $exomedb  = $params->{settings}->{$settings}->{exomedb}->{database};
my $rohtable = $params->{settings}->{$settings}->{exomedb}->{homozygosity};

my $coredb   = $params->{coredb}->{database};
my $sampletable = $params->{coredb}->{sampletable};

my $dbh = Utilities::connectCoreDB();
my $sql; my $sth;

# Check if sample exists and collect ID (if sample name provided)
if ( $sample ne "" )
{
	$sql = "SELECT idsample from $coredb.$sampletable where name=\"$sample\"";
	$sth = $dbh->prepare($sql) || die "Error in preparing $sql";
	$sth->execute() || die "Error in executing $sql";
	
	my ($tmpidsample ) = $sth->fetchrow_array();
	
	if ( ! defined $tmpidsample || $tmpidsample eq ""  )
	{
		$logger->error("Sample $sample does not exist");
		exit(-1);
	}
	
	$idsample = $tmpidsample;
}

# Run the ROH tool:

my $command = "";

# Indexed VCF file is required. Checking

my $tmpvcf = $vcf;

if ( $vcf =~ /^*.vcf$/ )
{
	$tmpvcf = $vcf;
	$tmpvcf =~ s/.vcf$/.tmproh.vcf.gz/g;
	$command = "cat $vcf | $bgzip -c > $tmpvcf";
	if (&Utilities::executeCommand($command, "Creating tmp compressed vcf for roh", $logger)) {
		$logger->error("Creating tmp compressed vcf for roh failed");
		exit(-1);
	}
}
elsif ( $vcf =~ /^*.vcf.gz$/ )
{
	$tmpvcf = $vcf;
}
else
{
	$logger->error("Input file $vcf not a VCF");
	exit(-1);
}

# Check that the vcf is indexed
if ( ! -e $tmpvcf.".tbi" )
{
	# Run Tabix and index the file (vcf.gz)
	$command = "$tabix -p vcf $tmpvcf";
	if (&Utilities::executeCommand($command, "Indexing compressed vcf $tmpvcf for roh", $logger)) {
		$logger->error("Indexing compressed vcf for roh failed");
		exit(-1);
	}
}


## 
##  # bcftools runs of homozygosity syntax
##  bcftools roh --AF-file (ALLELEFREQ FILE) -M (RECOMB RATE /BP) -o (OUTFILE)  (INFILE).vcf.gz 
##

my $AFoption = ""; 
   $AFoption = "--AF-file $gnomADfreqsROH " if $gnomADfreqsROH ne "";		# ROH can run without the AF option if no freq table is available 

$command = "
	$bcftools	roh 				\\
		-M $recombrate				\\
		-o $outfile					\\
		$AFoption $tmpvcf
";	

$logger->debug($command);

#if ( ! -e $outfile )
#{
	if (&Utilities::executeCommand($command, "Calculating stretches of homozygosity for $sample", $logger)) {
			$logger->error("Calculating stretches of homozygosity for $sample failed");
			exit(100);
	}else{
			$logger->debug("Calculating stretches of homozygosity for $sample OK");
	}
#}

# Delete temporary file
if ( $tmpvcf ne $vcf )
{
	# Remove $tmpvcf and index
	unlink($tmpvcf);
	unlink($tmpvcf.".tbi");
}


# Insert results
if ( $insertdb )
{

	# Open roh
	open(IN, "$outfile") || exit $logger->error("Cannot open $outfile!");
	
	$sql = "DELETE from $exomedb.$rohtable where idsample = $idsample";
	$sth = $dbh->prepare($sql) || die "Error in preparing $sql";
	$sth->execute() || die "Error in executing $sql";
	
	while (<IN>)
	{
		chomp;
		# Skip header and sliding window reports
		next if ( $_ =~ /^#.*$/ || $_ =~ /^ST.*$/ );
		
		# Get data to insert
		#print ">>  ".$_."\n";
		my ( $type, $sampleinvcf, $chrom, $start, $end, $length, $supporting, $score ) = split ("\t", $_ );

		# Skip X/Y
                next if $chrom =~ /X/;
                next if $chrom =~ /Y/;

		print $sampleinvcf."\n";
		if ( $sampleinvcf ne $sample )
		{
			$logger->error("Sample name in VCF does not match provided sample name. Can't insert data");
			exit(-1);
		} 

		# Record in roh table: 
		# idsample chrom start end count score
		$sql = "INSERT into $exomedb.$rohtable ( idsample, chrom, start, end, count, score ) values ( $idsample, \"$chrom\", $start, $end, $supporting, $score  )";
		$sth = $dbh->prepare($sql) || die "Error in preparing $sql";
		$sth->execute() || die "Error in executing $sql";			
	}
	
	
}



=head1 NAME

homozygosity.pl

=head1 SYNOPSIS

homozygosity.pl -s SAMPLENAME -i INPUT.vcf (or vcf.gz) -o OUTPUT_FILE [ -insert ] [ -recombrate ] [ -lf ] [ -ll ]

=head1 DESCRIPTION

This script executes bcftools roh to detect stretches of homozygosity and optionally inserts data into the database

=head1 OPTIONS

 -s	<samplename> 	required if insert results into database -insert enabled 
 -i	<vcf>			input vcf file
 -o <outfile>		output file for the bcftools roh tool
 -se				name of the settings in the current.config.xml file that holds path to reference genome, 
 					to the annotation file and to possible additional annotation files; use default settings if nothing is given
 -insert 			insert results into database
 -recombrate		recombination rate /bp (default 1e-4)
 -lf	log file; default: print to screen
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page
 

=head1 AUTHOR

Riccardo Berutti

=cut
