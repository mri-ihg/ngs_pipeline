#!/usr/bin/perl 

# Ignore vcftools warnings--------------------------------------------------------------------------------
# Note: in the original template insertSV, tabix was used and silenced
BEGIN { $SIG{__WARN__} = sub {warn $_[0] unless( $_[0] =~ m/Vcf.pm/) }; };

# Load general perl modules/libraries --------------------------------------------------------------------

use strict;
use Getopt::Long;
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);
use Pod::Usage;
use Data::Dumper;
use List::Util qw(min max);

# Load pipeline modules/libraries ------------------------------------------------------------------------

# load assisting pipeline perl functions
my $prog_path = dirname( abs_path($0) );
require $prog_path . "/Utilities.pm"; #e.g. for parsing the current.config.xml

# library for parsing a vcf file
my $params = Utilities::getParams();
require $params->{programs}->{vcftools}->{pm};

# Setup prerequisites ---------------------------------------------------------------------------

# define further variables
my $settings = "";

my $infile = "";
my $sample = "";

my $deleteOnly = "";

my $help = 0;
my $man = 0;
my $logfile = "SCREEN";
my $loglevel = "INFO";


# get command line arguments
GetOptions(
	"i=s"  => \$infile,
	"se=s" => \$settings,
	"s=s"  => \$sample,
	"d"    => \$deleteOnly,
	"lf=s" => \$logfile,
	"ll=s" => \$loglevel,
	"h"    => \$help,
	"man"  => \$man
);
# setup the logging system
Utilities::initLogger( $logfile, $loglevel );
my $logger = Utilities::getLogger();


pod2usage( {-exitval => 0  ,-verbose => 1} ) if ( $infile eq "" || $settings eq "" || $sample eq "" || ! -e $infile  );
pod2usage( {-exitval => 0  ,-verbose => 1} ) if $help;
pod2usage( {-exitval => 0  ,-verbose => 2} ) if $man;



# Setup connection to the database (wrapper function creates a DBI object) and get names of the database and tables where to store the results and look up sample information
my $dbh = Utilities::connectExomeDB($settings);

my $exomedb =			$params->{settings}->{$settings}->{exomedb}->{database}; # database where to store results (e.g. 'genomegatk' for setting "hg19_wholegenome"
my $expansionTable =		$params->{settings}->{$settings}->{exomedb}->{expansiontable}; # typical value "expansion"
my $expansionSampleTable =	$params->{settings}->{$settings}->{exomedb}->{expansionsampletable}; # typical value "expansionsample"
my $coredb      = $params->{coredb}->{database}; # typical value "exomehg19": needed to link sample name in $sample with its ID in the database
my $sampleTable = $params->{coredb}->{sampletable}; # typical value "sample": contains the aforementioned link


# Check id sample and if exists
my $query =  qq{select idsample from $coredb.$sampleTable where name = ? };
my @parameters = ( $sample );
my $idsample = runDbQuery($dbh, $query, $logger, \@parameters, "fetch");
if ( !defined( $idsample ) ) {
	$logger->error("sample name $sample not found in $coredb.$sampleTable - exiting");
	exit(100);
}


# Start transaction: all the insertion can only succeed together
my @parameters = ();
runDbQuery($dbh, "START TRANSACTION", $logger, \@parameters );


# Delete previous entries
$logger->info("Removing previous expansion alleles from database for sample...");
$query = "delete from $expansionSampleTable where idsample= ? ";
@parameters = ( $idsample );
runDbQuery($dbh, $query, $logger, \@parameters);

if ( $deleteOnly ) # If delete only then stop here
{
	@parameters = ();
	runDbQuery($dbh, "COMMIT", $logger, () );
	exit(0);
}


# Check: delete orphan expansions
# $query = "delete from $expansionTable where (SELECT count(distinct EST.idexpansion) from $expansionSampleTable EST where EST.idexpansion=$expansionTable.idexpansion)=0";
# runDbQuery($dbh, $query, $logger, "fetch");


# Read in the vcf results file from ExpansionHunter, which contains repeat expansion data only for a single sample $sample
# create new object for the vcf file and remove the header
my $vcf = Vcf->new( file => $infile );
$vcf->parse_header();

# process each "line" of the vcf object, check whether the very expansion was already registered (or create that) and (if not done already) add the expansion for the given $sample.
$logger->info("Inserting Expansions...");
while ( my $vcfline = $vcf->next_data_hash() ) {
	# define variables
	my ($chrom, $start, $end, $refallele, $repeatunit, $refcount, $length,$count1,$count2,$cilo1, $cihi1, $cilo2, $cihi2,$coverage, $filter) = ("NULL") x 15;
	
	# get the values from the hash into new variables for the table $expansionTable: note that the entry for the frequency of the reference expansion has to be calculated for the database (see further below)
	$chrom = $vcfline->{CHROM} if defined $vcfline->{CHROM};
	$start = $vcfline->{POS} if defined $vcfline->{POS};
	$end = $vcfline->{INFO}->{END} if defined $vcfline->{INFO}->{END};
	$refallele = $vcfline->{REF} if defined $vcfline->{REF};
	
	$repeatunit = $vcfline->{INFO}->{RU} if defined $vcfline->{INFO}->{RU};
	$refcount = $vcfline->{INFO}->{REF} if defined $vcfline->{INFO}->{REF};
	$length = $vcfline->{INFO}->{RL} if defined $vcfline->{INFO}->{RL};
	# get the values from the hash into new variables for the table $expansionSampleTable
	
	# TODO: merged.rmdup hardcoded -> might have to be just the filename without .bam 	
	($count1,$count2) = split(/[\/|]/,$vcfline->{gtypes}->{"merged.rmdup"}->{REPCN}) if defined $vcfline->{gtypes}->{"merged.rmdup"}->{REPCN}; # Repeat count for each allele.
	next if $count1 == "."; # is a no call due to low depth
	($cilo1, $cihi1, $cilo2, $cihi2) = split(/[\-\/\|]/,$vcfline->{gtypes}->{"merged.rmdup"}->{REPCI}) if defined $vcfline->{gtypes}->{"merged.rmdup"}->{REPCI}; # bounds of the confidence intervals
	$coverage = $vcfline->{gtypes}->{"merged.rmdup"}->{LC} if defined $vcfline->{gtypes}->{"merged.rmdup"}->{LC}; # Average coverage
	
	$filter = join( ',', @{ $vcfline->{FILTER} } ); #Filter
	$filter = "PASS" if ( $filter eq "." );


	# TODO: Collapse into INSERT ... ON DUPLICATE IGNORE and GET LAST ID.
	# Check whether the variant exists already in the database
	my $query = qq{select idexpansion from $exomedb.$expansionTable ET where ET.chrom = ? and ET.start = ? and ET.end = ? and ET.refallele = ? and ET.repeatunit = ? and ET.refcount = ? and ET.length = ?};
	@parameters = ($chrom, $start, $end, $refallele, $repeatunit, $refcount, $length );
	my $idExpansion = runDbQuery($dbh, $query, $logger, \@parameters, "fetch");
	if ( !defined( $idExpansion ) ) {
		#New expansion insert	
		$query = qq{insert into $exomedb.$expansionTable (chrom, start, end, refallele, repeatunit, refcount, length) values (?,?,?,?,?,?,?)};
		@parameters = ($chrom, $start, $end, $refallele, $repeatunit, $refcount, $length );
		runDbQuery($dbh, $query, $logger, \@parameters);
		
		#Now again look up the current expansion to know the value for the key "idexpansion"
		my $query = qq{select idexpansion from $exomedb.$expansionTable ET where ET.chrom = ? and ET.start = ? and ET.end = ? and ET.refallele = ? and ET.repeatunit = ? and ET.refcount = ? and ET.length = ?};
		@parameters = ($chrom, $start, $end, $refallele, $repeatunit, $refcount, $length );
		$idExpansion = runDbQuery($dbh, $query, $logger, \@parameters, "fetch");	
	}
	# End of TODO
	
	# Insert ExpansionSample
	$query = qq{insert into $exomedb.$expansionSampleTable (idexpansion,idsample,count1,count2,cilo1,cihi1,cilo2,cihi2,coverage,filter) values (?,?,?,?,?,?,?,?,?,?)};
	@parameters = ( $idExpansion, $idsample, $count1, $count2, $cilo1, $cihi1, $cilo2, $cihi2, $coverage, $filter );
	runDbQuery($dbh, $query, $logger, \@parameters);
	
	#update the registered database frequency of the reference allele expansion
	$query = qq{select count(*) from $exomedb.$expansionSampleTable EST where EST.idexpansion = ? and EST.count1 = ? and EST.count2 = ? };
	@parameters = ( $idExpansion, $refcount, $refcount ); 
	my $countHomoRefSample = runDbQuery($dbh, $query, $logger, \@parameters, "fetch");
	
	$query = qq{select count(*) from $exomedb.$expansionSampleTable EST where EST.idexpansion = ? and (EST.count1 = ? or EST.count2 = ? ) };
	@parameters = ( $idExpansion, $refcount, $refcount ); 
	my $countHomoHetRefSample = runDbQuery($dbh, $query, $logger, \@parameters, "fetch");	
	
	$query = qq{select count(*) from $exomedb.$expansionSampleTable EST where EST.idexpansion = ? };
	@parameters = ( $idExpansion );
	my $countTotalSample = runDbQuery($dbh, $query, $logger, \@parameters, "fetch");
	
	# Reference Frequency
	my $refFreq = (2*$countHomoRefSample+($countHomoHetRefSample-$countHomoRefSample)); #/$countTotalSample;
	$query = qq{update $exomedb.$expansionTable ET set ET.freq = ? where ET.idexpansion = ?};
	@parameters = ( $refFreq, $idExpansion); 
	runDbQuery($dbh, $query, $logger, \@parameters);
		
};

# Finished - commit
@parameters=();
runDbQuery($dbh, "COMMIT", $logger, \@parameters );



$logger->info("Finished inserting Expansions.");

#ERIK: added to close db connection
$dbh->disconnect;

sub runDbQuery {
	my $dbObj = shift; #The database object.
	my $dbQuery = shift; #The database query character string.
	my $loggerObj = shift; #The logger object.
	my $parameters = shift;
	my $fetch = shift; #Indicator for fetching a first row. 
	$loggerObj->debug($dbQuery);
	my $dbQueryObj = $dbObj->prepare($dbQuery) || $loggerObj->error("Can't prepare statement: $DBI::errstr");

	#$dbQueryObj->execute(@$parameters) || $loggerObj->error($DBI::errstr);
	
	$dbQueryObj->execute(@$parameters) || die $DBI::errstr;

	if ($fetch eq "fetch") { #Note: Perl interpretes undefined variables as FALSE.
		my $outVal = $dbQueryObj->fetchrow_array(); 
		return $outVal;		
				
	};
	
}

=head1 NAME

insertExpansionHunter.pl

=head1 SYNOPSIS

 insertExpansionHunter.pl -i [EXPANSION-HUNTER-VCF-FILE] -s [SAMPLE] -se settings 

=head1 DESCRIPTION

This script inserts expansion calls from ExpansionHunter output into the database.

=head1 OPTIONS

 -s	sample name; REQUIRED
 -i	VCF file containing ExpansionHunter calls, REQUIRED
 -se	name of the settings in the current.config.xml file that holds path to reference genome, 
 	to the annotation file and to possible additional annotation files; REQUIRED
 -lf	log file; default: print to screen
 -d	exit after deleting previous sample's expansions (no insertion of new data)
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page

=head1 AUTHOR

Erik Tilch, Riccardo Berutti

=cut
