#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use DBI;
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);

my $prog_path = dirname( abs_path($0) );
require $prog_path . "/Utilities.pm";


#deletes all entries in the snvsample table for a given sample



my $logfile  = "pipeline.log";
my $loglevel = "INFO";
my $patId	 = "";
my $help     = 0;
my $settings = "";
my $caller   = "";

GetOptions(
	"p=s"  => \$patId,
	"c=s"  => \$caller,
	"lf=s" => \$logfile,
	"ll=s" => \$loglevel,
	"se=s" => \$settings,
	"h"    => \$help,
);

my $params         = Utilities::getParams();
my $coredb         = $params->{coredb}->{database};
my $sampleTable    = $params->{coredb}->{sampletable};

my $exomedb        = $params->{settings}->{$settings}->{exomedb}->{database};
my $snvsampleTable = $params->{settings}->{$settings}->{exomedb}->{snvsampletable};


Utilities::initLogger( $logfile, $loglevel );
my $logger = Utilities::getLogger();

if ( $help == 1 ) {
	print "snvdbExomeDelete.pl
-p <patientID>
-se	settings
-c	<caller> if defined only variants from the specified caller will be deleted
-lf log file; default: pipeline.log
-ll log level: ERROR,INFO,DEBUG; default: INFO
-h this help\n";
	exit(1);
}

my $dbh = Utilities::connectExomeDB($settings);


my $sql =
qq{delete from $exomedb.$snvsampleTable where $exomedb.$snvsampleTable.idsample=(select $coredb.$sampleTable.idsample from $coredb.$sampleTable where $coredb.$sampleTable.name='$patId')};

$sql .= " AND find_in_set('$caller',caller)>0" if $caller ne "";

$logger->debug("query: $sql");

my $sth = $dbh->prepare($sql)
  || $logger->error("Can't prepare statement: $DBI::errstr");
$sth->execute() || $logger->error($DBI::errstr);