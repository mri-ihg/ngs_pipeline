#!/usr/bin/perl 


use strict;
use Getopt::Long;
use DBI;
use warnings;
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);

my $prog_path = dirname(abs_path($0));
require $prog_path."/Utilities.pm";
require $prog_path."/CollectMetrics.pm";




my $sample   = "";
my $logfile  = "pipeline.log";
my $loglevel = "INFO";
my $help	 = 0;
my $bamfile  = "";
my $libtype  = ""; 

GetOptions(
	"s=s"   => \$sample,
	"b=s"	=> \$bamfile,
	"l=s"	=> \$libtype,
	"lf=s" => \$logfile,
	"ll=s" => \$loglevel,
	"h"	   => \$help
);

if ( $help == 1 ||  $bamfile eq "") {
	print "insertSizeIntoDB.pl
-s	sample name; required if library isn't specified correctly in \@RG header
-l	libtype; required if library isn't specified correctly in \@RG header
-b	</path/to/file.bam> - bam file for which the insert size is calculated; the directory of this file is used as output directory; required
-lf	log file; default: pipeline.log
-ll	log level: ERROR,INFO,DEBUG; default: INFO
-h	this help\n";
	exit(1);
}


my $params         = Utilities::getParams();
my $coredb         = $params->{coredb}->{database};

my @columns;

Utilities::initLogger( $logfile, $loglevel );
my $logger = Utilities::getLogger();

my $dbh = Utilities::connectSolexaDB();

my $outdir = dirname($bamfile);

my ($mean,$sd,$pointer) = CollectMetrics::getInsertSize($bamfile,$outdir);

my %libs = %$pointer;

foreach my $key(keys %libs){		# for each library --> insert 

	my $sql = "UPDATE library
SET linsertsize=$libs{$key}->[0], linsertsizesd=$libs{$key}->[1]
";
	
	my $checksql="SELECT count(lid) FROM library WHERE lname='$key';";
	my $out = $dbh->prepare($checksql) || $logger->error($DBI::errstr);			
	$out->execute || $logger->error($DBI::errstr);
	
	if($out->fetchrow_array() == 1){		# if the correct library ID was in the @RG header
		$sql .= "WHERE lname='$key';";
	}else{
		$sql .= "WHERE libtype=(SELECT ltid FROM libtype WHERE ltlibtype='$libtype')
AND lid IN (SELECT sl.lid FROM sample2library sl INNER JOIN ".$coredb.".sample s ON sl.idsample=s.idsample WHERE s.name='$sample');";
	}
	
	$logger->debug($sql);
	
	$out = $dbh->prepare($sql) || $logger->error($DBI::errstr);			
	$out->execute || $logger->error($DBI::errstr);


	
}