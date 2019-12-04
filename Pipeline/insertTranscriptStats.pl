#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use DBI;
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);


my $prog_path = dirname(abs_path($0));
require $prog_path."/Utilities.pm";




my $infile     = "";
my $assay      = "";                                                
my $help       = 0;
my $logfile    = "pipeline.log";
my $loglevel   = "INFO";
my $sample     = "";
my $libtype	   = "exomic";
my $libpair    = "paired-end";
my $delete     = 0;


GetOptions(
"i=s"  => \$infile, 
"a=s"  => \$assay,
"s=s"  => \$sample,
"d"    => \$delete,
"lt=s" => \$libtype,
"lp=s" => \$libpair,
"lf=s" => \$logfile,
"ll=s" => \$loglevel, 
"h" => \$help);


if ($help == 1 || $sample eq "") {
print 
"
This script inserts stats per transcript into the database (depth, mapping quality,...)


-i	<infile>; required
-s	<samplename>; required
-d	delete (possible) old entries before inserting new ones
-a	assay
-lt	libtype; default: $libtype
-lp	libpair; default: $libpair
-lf	log file; default: pipeline.log
-ll	log level: ERROR,INFO,DEBUG; default: INFO
-h	this help\n";
exit(1);
}

Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();


my $params = Utilities::getParams();

my $transcriptstattable = $params->{coredb}->{transcriptstattable};
my $transcripttable     = $params->{coredb}->{transcripttable};
my $genetable           = $params->{coredb}->{genetable};
my $sampletable         = $params->{coredb}->{sampletable};
my $organismtable       = $params->{coredb}->{organismtable};
my $solexadb			= $params->{solexadb}->{database};
my $libtypetable        = $params->{solexadb}->{libtypetable};
my $libpairtable        = $params->{solexadb}->{libpairtable};
my $assaytable          = $params->{solexadb}->{assaytable};
my $dbh                 = Utilities::connectCoreDB();





#get idsample
my $sql;
my $sth;
my ($idsample,$organism);

if($organismtable){
	$sql = qq{select s.idsample,o.orname from $sampletable s inner join $organismtable o on s.idorganism=o.idorganism where s.name = '$sample' };
	$sth = $dbh->prepare( $sql ) || $logger->error("Can't prepare statement: $DBI::errstr");
	$sth->execute() || $logger->error("Can't execute statement: $DBI::errstr");
	($idsample,$organism)=$sth->fetchrow_array();
	if($organism ne "human"){
		$logger->info("no transcripts inserted because organism \"$organism\" is not \"human\"");
		exit(1);
	}
}else{
	$sql = qq{select s.idsample from $sampletable s where s.name = '$sample' };
	$sth = $dbh->prepare( $sql ) || $logger->error("Can't prepare statement: $DBI::errstr");
	$sth->execute() || $logger->error("Can't execute statement: $DBI::errstr");
	$idsample=$sth->fetchrow_array();
}
if(!defined($idsample))
{
	$logger->error("sample name $sample not found in $sampletable - exiting"); 
	exit(1);
}




#get idlibtype
$sql = qq{select ltid from $solexadb.$libtypetable where replace(ltlibtype,' ','')='$libtype'};
$logger->debug($sql);
$sth = $dbh->prepare( $sql ) || $logger->error("Can't prepare statement: $DBI::errstr");
$sth->execute() || $logger->error("Can't execute statement: $DBI::errstr");
my $idlibtype = $sth->fetchrow_array();

#get idlibpair
$sql = qq{select lpid from $solexadb.$libpairtable where replace(lplibpair,' ','')='$libpair'};
$logger->debug($sql);
$sth = $dbh->prepare( $sql ) || $logger->error("Can't prepare statement: $DBI::errstr");
$sth->execute() || $logger->error("Can't execute statement: $DBI::errstr");
my $idlibpair = $sth->fetchrow_array();

#get idassay
my $idassay = "NULL";
$sql = qq{select idassay from $solexadb.$assaytable where name='$assay'};

$logger->debug($sql);
$sth = $dbh->prepare( $sql ) || $logger->error("Can't prepare statement: $DBI::errstr");
$sth->execute() || $logger->error("Can't execute statement: $DBI::errstr");
if ($sth->rows) {
	$idassay = $sth->fetchrow_array();
} else {
	$idassay = "NULL";
}

#delete old entries from db
if($delete){
	$sql = qq{delete from $transcriptstattable where idsample=$idsample and idlibtype=$idlibtype and idlibpair=$idlibpair};
	$sth = $dbh->prepare( $sql ) || $logger->error("Can't prepare statement: $DBI::errstr");
	$sth->execute() || $logger->error("Can't execute statement: $DBI::errstr");
}

open(IN, "$infile") || exit $logger->error("Cannot open $infile!");

my $notInsertedTranscripts = 0;
my $insertedTranscripts    = 0;

while (<IN>)
{
	chomp;
	next if $_ =~ /^#/;
	my ($region,$regionchr,$exonStarts,$exonEnds,$cov20x,$covAvg,$mqAvg,$scov20x,$scovAvg,$smqAvg) = split();
	my ($transcript,$idtranscript) = split(',',$region);
	
	$scov20x = "" unless $scov20x;
	$scovAvg = "" unless $scovAvg;
	$smqAvg  = "" unless $smqAvg;
	$cov20x  = 0  if(!defined($cov20x) || $cov20x eq "");
	$covAvg  = 0  if(!defined($covAvg) || $covAvg eq "");
	$mqAvg   = 0  if(!defined($mqAvg) || $mqAvg eq "");
	



	#insert transcriptstats into database
	$sql = qq{insert into $transcriptstattable (idtranscript,idsample,idlibtype,idlibpair,avgdepth,qdepth,avgmapqual,avgdepthtotal,qdepthtotal,avgmapqualtotal,idassay) values ($idtranscript,$idsample,$idlibtype,$idlibpair,'$scovAvg','$scov20x','$smqAvg',$covAvg,$cov20x,$mqAvg,$idassay)};
	#print $sql."\n";
	$logger->debug($sql);
	$sth = $dbh->prepare( $sql ) || $logger->error("Can't prepare statement: $DBI::errstr");
	$sth->execute() || $logger->error("Can't execute statement: $DBI::errstr");
	
	$insertedTranscripts++;
}

close IN;

$logger->info("Succesfully inserted $insertedTranscripts entries");
$logger->info("Couldn't insert $notInsertedTranscripts entries, because the referenced transcript is not in the database");