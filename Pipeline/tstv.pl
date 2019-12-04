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


my $sample    = "";
my $settings  = "default";
my $help      = 0;
my $man		  = 0;


GetOptions(
	"s=s"  => \$sample,
	"se=s" => \$settings,
	"lf=s" => \$logfile,
	"ll=s" => \$loglevel,
	"h"    => \$help,
	"man"  => \$man
);

pod2usage( { -exitval => 0, -verbose => 1 } ) if $help;
pod2usage( { -exitval => 0, -verbose => 2 } ) if $man;
pod2usage( { -exitval => 1, -verbose => 1 } ) if $sample eq "";


my $params         = Utilities::getParams();


my $exomedb                   = $params->{settings}->{$settings}->{exomedb}->{database};
my $snvTable                  = $params->{settings}->{$settings}->{exomedb}->{snvtable};
my $snvsampleTable            = $params->{settings}->{$settings}->{exomedb}->{snvsampletable};

my $coredb      = $params->{coredb}->{database};
my $sampleTable = $params->{coredb}->{sampletable};
my $statTable   = $params->{coredb}->{stattable};

Utilities::initLogger( $logfile, $loglevel );
my $logger = Utilities::getLogger();

$dbh = Utilities::connectExomeDB($settings);

$sql = "update $coredb.$statTable e
inner join $coredb.$sampleTable s on s.idsample=e.idsample
set e.tstv=
(select ROUND((select count(v.idsnv) 
               from $snvTable v inner join $snvsampleTable ss on ss.idsnv=v.idsnv 
               where ss.idsample=s.idsample and v.class='snp' and ss.filter='PASS' 
               and ( (refallele='A' AND allele='G') OR  (refallele='G' AND allele='A') OR (refallele='C' AND allele='T') OR (refallele='T' AND allele='C'))
               and (find_in_set('missense',v.func) or find_in_set('syn',v.func) or find_in_set('nonsense',v.func) or find_in_set('stoploss',v.func))
              )/
             (select count(v.idsnv) 
             from $snvTable v inner join $snvsampleTable ss on ss.idsnv=v.idsnv 
             where ss.idsample=s.idsample and v.class='snp' and ss.filter='PASS' 
             and NOT ( (refallele='A' AND allele='G') OR  (refallele='G' AND allele='A') OR (refallele='C' AND allele='T') OR (refallele='T' AND allele='C'))
             and (find_in_set('missense',v.func) or find_in_set('syn',v.func) or find_in_set('nonsense',v.func) or find_in_set('stoploss',v.func))),
           2))
where s.name='$sample';
";

$logger->debug($sql);
$sth = $dbh->prepare($sql) || $logger->error("Can't prepare statement: $DBI::errstr");
$sth->execute() || $logger->error($DBI::errstr);

=head1 NAME

tstv.pl

=head1 SYNOPSIS

tstv.pl -s SAMPLE  

=head1 DESCRIPTION

This script calculates the Ts/Tv ratio of the coding variants of a given sample and updates the respective
value in the exomestat table.

=head1 OPTIONS

 -s	[samplename] REQUIRED
 -se	name of the settings in the current.config.xml file that holds path to reference genome, 
 	to the annotation file and to possible additional annotation files; use default settings if nothing is given
 -lf	log file; default: print to screen
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page

=head1 AUTHOR

Thomas Wieland

=cut
