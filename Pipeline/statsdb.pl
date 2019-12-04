#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use DBI;
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);
use Pod::Usage;

my $prog_path = dirname(abs_path($0));
require $prog_path."/Utilities.pm";

my $sql	       = "";
my $sth        = "";

my $infile   = "";
my $assay    = "";                                                
my $help     = 0;
my $man		 = 0;
my $nocommit = 0;

my $logfile  = "SCREEN";
my $loglevel = "INFO";
my $settings = "";
my $libtype	 = "exomic";
my $libpair  = "paired-end";
my $covXFile = "";

my ($sample,$dup,$reads,$mapped,$mappedPer,$seq,$uncov,$x1,$x4,$x8,$x20,$avgCov,$std,$medianCov,$mstd,$autosomalCov,$onBait);
my $idsample = 0;

GetOptions(
"i=s"  => \$infile, 
"a=s"  => \$assay,
"d=s"  => \$covXFile,
"lt=s" => \$libtype,
"lp=s" => \$libpair,
"lf=s" => \$logfile,
"se=s" => \$settings,
"ll=s" => \$loglevel, 
"h"    => \$help,
"man"  => \$man,
"nocommit" => \$nocommit);

pod2usage( {-exitval => 0  ,-verbose => 1} ) if $help;
pod2usage( {-exitval => 0  ,-verbose => 2} ) if $man;
pod2usage( {-exitval => 1  ,-verbose => 1} ) if $infile eq "";


Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();

my $params = Utilities::getParams();

my $database       = $params->{coredb}->{database};
my $statTable      = $params->{coredb}->{stattable};
my $sampleTable    = $params->{coredb}->{sampletable};
my $percentcovTable= $params->{coredb}->{percentcoveragetable};
my $solexadb	   = $params->{solexadb}->{database};
my $libtypetable   = $params->{solexadb}->{libtypetable};
my $libpairtable   = $params->{solexadb}->{libpairtable};
my $assaytable     = $params->{solexadb}->{assaytable};

	# Add stats	from PICARDtools		RB 2016.02.18 
	# 	mismatchrate
	# 	libcomplexity	
	# 	avgqual
	# 	avgquallast5
	# 	q30fraction
	#   avg difference of coverage btw coding and nondcoding regions
	#   support for picard rmdup 
	
	# Added transactions
	# Add SIMULATION MODE (ONLY TEST QUERIES BUT NO COMMIT)

# Open input stats file
open(IN, "$infile") || exit $logger->error("Cannot open $infile!");

# Connect to DB
my $dbh = Utilities::connectCoreDB();

# Everything is a transaction. On failure no commit will be issued
$sql="START TRANSACTION;";
$dbh->do($sql) || die print "$DBI::errstr";


while (<IN>)
{
	chomp;
	$_ =~ s/\%//g;
	my ($sample,$dup,$reads,$mapped,$mappedPer,$seq,$uncov,$x1,$x4,$x8,$x20,$avgCov,$std,$medianCov,$mstd,$autosomalCov,$onBait,$mix,$exomedepthrsd,$gender,$bampath,$mismatchrate,$libcomplexity,$avgqual,$avgquallast5,$q30fraction,$avgdiffdepth,$opticalduplicates,$properlyp) = split("\t",$_);
	
	next if $sample =~ /^#/;	#skip header

	$dup=sprintf("%.2f",$dup) if $dup ne "NULL";
	$seq = sprintf("%.2f",$seq) if $seq ne "NULL";

	# Set last fields to null if accidentally run on old stats file
	#Might be undefined or empty
	
	$mismatchrate = 'NULL' if ! defined $mismatchrate;
	$libcomplexity= 'NULL' if ! defined $libcomplexity;
	$avgqual = 'NULL' if ! defined $avgqual;
	$avgquallast5 = 'NULL' if ! defined $avgquallast5;
	$q30fraction = 'NULL' if ! defined $q30fraction;
	$avgdiffdepth = 'NULL' if ! defined $avgdiffdepth;
	$opticalduplicates = 'NULL' if ! defined $opticalduplicates;
	$properlyp = 'NULL' if ! defined $properlyp;
	
	$mismatchrate 	= (		$mismatchrate 		eq "" ? 'NULL'	:	$mismatchrate  	);
	$libcomplexity	= (		$libcomplexity		eq "" ? 'NULL'	:	$libcomplexity  );
	$avgqual	    = (		$avgqual			eq "" ? 'NULL'	:	$avgqual  		);
	$avgquallast5 	= (		$avgquallast5 		eq "" ? 'NULL'	:	$avgquallast5  	);
	$q30fraction 	= (		$q30fraction 		eq "" ? 'NULL'	:	$q30fraction  	);
	$avgdiffdepth 	= (		$avgdiffdepth 		eq "" ? 'NULL'	:	$avgdiffdepth  	);
	$opticalduplicates= (	$opticalduplicates  eq "" ? 'NULL'	:	$opticalduplicates	);
	$properlyp		= (     $properlyp          eq "" ? 'NULL'  :   $properlyp);		

	$idsample = &getIdsample($sample);
	$logger->debug("$idsample,$sample,$dup,$reads,$mapped,$mappedPer,$seq,$uncov,$x1,$x4,$x8,$x20,$avgCov,$std,$medianCov,$mstd,$onBait,$mix,$exomedepthrsd,$gender,$bampath,$mismatchrate,$libcomplexity,$avgqual,$avgquallast5,$q30fraction,$avgdiffdepth,$opticalduplicates,$properlyp");
	&insert($idsample,$dup,$reads,$mapped,$mappedPer,$seq,$uncov,$x1,$x4,$x8,$x20,$avgCov,$onBait,$libtype,$libpair,$std,$medianCov,$mstd,$mix,$exomedepthrsd,$gender,$bampath,$mismatchrate,$libcomplexity,$avgqual,$avgquallast5,$q30fraction,$avgdiffdepth,$opticalduplicates,$properlyp);
}

# If here with no error then commit
if ( ! $nocommit )
{
	$sql="COMMIT;";
	$dbh->do($sql) || die print "$DBI::errstr";
}
else
{
	$logger->info("Not committed. Remove option -nocommit");
}
my $elapsed = (time - $^T)/60;
$logger->info("finished in $elapsed minutes");

#############
#subroutines#
#############

sub getIdsample {

	my $patId = shift;
	my $sql = "";
	my $sth = "";
	my $id = "";
	
	
	#check if patient id is valid
	$sql = qq{select idsample from $database.$sampleTable where name = '$patId' };
	$sth = $dbh->prepare( $sql ) || $logger->error("Can't prepare statement: $DBI::errstr");
	$sth->execute();
	($id)=$sth->fetchrow_array();
	if(!defined($id))
	{
		print "patient name $patId not found in $database.$sampleTable - exiting\n"; 
		exit(1);
	}
	
	return $id;

}

sub insert {


	my ($idsample,$dup,$reads,$mapped,$mappedPer,$seq,$uncov,$x1,$x4,$x8,$x20,$avgCov,$onBait,$libtype,$libpair,$std,$medianCov,$mstd,$mix,$exomedepthrsd,$gender,$bampath,$mismatchrate,$libcomplexity,$avgqual,$avgquallast5,$q30fraction,$avgdiffdepth,$opticalduplicates,$properlyp) = @_;
	
	my $sql = "";
	my $sth = "";
	
	#get idlibtype
	$sql = qq{select ltid from $solexadb.$libtypetable where replace(ltlibtype,' ','')='$libtype'};
	$logger->debug($sql);
	$sth = $dbh->prepare( $sql ) || $logger->error("Can't prepare statement: $DBI::errstr");
	$sth->execute();
	my $idlibtype = $sth->fetchrow_array();
	
	#get idlibpair
	$sql = qq{select lpid from $solexadb.$libpairtable where replace(lplibpair,' ','')='$libpair'};
	$logger->debug($sql);
	$sth = $dbh->prepare( $sql ) || $logger->error("Can't prepare statement: $DBI::errstr");
	$sth->execute();
	my $idlibpair = $sth->fetchrow_array();
	
	#get idassay
	my $idassay = "NULL";
	if($assay ne ""){
		$sql = qq{select idassay from $solexadb.$assaytable where name='$assay'};
		$logger->debug($sql);
		$sth = $dbh->prepare( $sql ) || $logger->error("Can't prepare statement: $DBI::errstr");
		$sth->execute() || $logger->error("Can't execute statement: $DBI::errstr");
		if ($sth->rows) {
			$idassay = $sth->fetchrow_array();
		} else {
		$idassay = "NULL";
		}
	}

	
	
	# SQL INSERT OR UPDATE
	$sql = qq{
		insert into $database.$statTable 
		( 
			idsample, idassay, duplicates, `reads`, mapped, percentm, seq, onbait, avgcov, uncovered, cov1x, cov4x, cov8x, cov20x, idlibtype,sry,idlibpair,mix,avgcovstd,mediancov,mediancovstd,exomedepthrsd, 
			mismatchrate, libcomplexity, avgqual, avgquallast5, q30fraction, avgdiffdepth, opticalduplicates, properlyp
		) 
		values
		(
			'$idsample', $idassay, '$dup', '$reads', '$mapped', '$mappedPer', '$seq', '$onBait', '$avgCov', '$uncov', '$x1', '$x4', '$x8', '$x20',$idlibtype,$gender,$idlibpair,$mix,$std,$medianCov,$mstd,$exomedepthrsd,
		      $mismatchrate,$libcomplexity,$avgqual,$avgquallast5,$q30fraction,$avgdiffdepth,$opticalduplicates,$properlyp
		) 
	          ON DUPLICATE KEY UPDATE idassay=$idassay, duplicates='$dup', `reads`='$reads', mapped='$mapped', percentm='$mappedPer', seq='$seq', onbait='$onBait', avgcov='$avgCov', uncovered='$uncov', 
	          cov1x='$x1', cov4x='$x4', cov8x='$x8', cov20x='$x20',sry=$gender,mix=$mix,avgcovstd=$std,mediancov=$medianCov,mediancovstd=$mstd,exomedepthrsd=$exomedepthrsd,
			  mismatchrate=$mismatchrate, libcomplexity=$libcomplexity, avgqual=$avgqual, avgquallast5=$avgquallast5, q30fraction=$q30fraction,avgdiffdepth=$avgdiffdepth,opticalduplicates=$opticalduplicates,properlyp=$properlyp
	};
		
	$logger->debug($sql);
	$sth = $dbh->prepare( $sql ) || $logger->error("Can't prepare statement: $DBI::errstr");
	$sth->execute() || $logger->error("Can't execute statement: $DBI::errstr");
	
	
	if($bampath ne "")
	{
		$sql = qq{update $database.$sampleTable set sbam='$bampath' where idsample=$idsample;};
		$logger->debug($sql);
		$sth = $dbh->prepare( $sql ) || $logger->error("Can't prepare statement: $DBI::errstr");
		$sth->execute() || $logger->error("Can't execute statement: $DBI::errstr");
	}
	
	if(-e $covXFile)
	{
		open COVX,$covXFile || exit $logger->error("Cannot open $covXFile!");
		
		#get exomestat id
		$sql = qq{select idexomestat from $database.$statTable where idsample=$idsample and idlibtype=$idlibtype and idlibpair=$idlibpair};
		$logger->debug($sql);
		$sth = $dbh->prepare( $sql ) || $logger->error("Can't prepare statement: $DBI::errstr");
		$sth->execute();
		
		if(my $idexomestat = $sth->fetchrow_array())
		{
			while(<COVX>)
			{
				chomp;
				next if $_ =~ /^#/;
				my ($covx,$percent) = split();
				$sql = qq{insert into $database.$percentcovTable (idexomestat,covx,percent) values($idexomestat,$covx,$percent) };
				$logger->debug($sql);
				$sth = $dbh->prepare( $sql ) || $logger->error("Can't prepare statement: $DBI::errstr");
				$sth->execute();
			}
		}
	}

}


=head1 NAME

statsdb.pl

=head1 SYNOPSIS

statsdb.pl -i summary.stats.tsv

=head1 DESCRIPTION

This script inserts the tab separated statistics file produced by parseStats.pl into the database. 


=head1 OPTIONS

 -i	<summary.stats.tsv> output file of the parseStats.pl script; REQUIRED
 -a	[assay] type of the assay (e.g. exome kit) used to prepare the library; default: empty
 -lt	[libtype] type of the library (e.g. exomic, genomic,...); default: exomic
 -lp	[libpair] type of the run (e.g. paired-end, single-end); default: paired-end
 -d	<coverageDistribution.out> coverage distribution file produced by calcCoverage.pl
 -se	name of the settings in the current.config.xml file that holds path to reference genome, 
 	to the annotation file and to possible additional annotation files; use default settings if nothing is given
 -lf	log file; default: print to screen
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page
 -nocommit tests the script without modification to the db 

=head1 AUTHOR

Thomas Wieland

=cut
