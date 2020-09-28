#!/usr/bin/perl 

use strict;
use Getopt::Long;
use DBI;
use warnings;
use File::Basename;
use Cwd qw(abs_path);
use Pod::Usage;


my $prog_path = dirname(abs_path($0));
require $prog_path."/Utilities.pm";

my ($name,$chrom,$strand,$txStart,$txEnd,$cdsStart,$cdsEnd,$exonCount,$exonStarts,$exonEnds);
my ($query,$out);
my $seq = "";
my @exonStarts = ();
my @exonEnds = ();
my $settings  = "default";
my $logfile   = "SCREEN";
my $loglevel  = "INFO";
my $genetable = "";
my $database  = "";
my $cdstable  = "";
my $help      = 0;
my $man		  = 0;
my $transcript= "";
my $delete    = 0;

GetOptions(
"se=s" => \$settings,
"g=s"  => \$genetable,
"d=s"  => \$database,
"c=s"  => \$cdstable,
"t=s"  => \$transcript,
"e"    => \$delete,
"lf=s" => \$logfile,
"ll=s" => \$loglevel, 
"man"  => \$man,
"h"    => \$help

);

pod2usage( { -exitval => 0, -verbose => 1 } ) if $help;
pod2usage( { -exitval => 0, -verbose => 2 } ) if $man;

Utilities::initLogger( $logfile, $loglevel );
my $logger = Utilities::getLogger();


#connect
my $dbh = Utilities::connectVariationDB($settings);

if($genetable eq "" || $database eq "" || $cdstable eq ""){
	my $params = Utilities::getParams();
	$database       = $params->{settings}->{$settings}->{variationdb}->{database};
	$genetable      = $params->{settings}->{$settings}->{variationdb}->{genetable};
	$cdstable       = $params->{settings}->{$settings}->{variationdb}->{cds};
}


#create table (if required)
$query ="CREATE TABLE $cdstable (
  `name` varchar(25) NOT NULL,
  `cds` longblob NOT NULL,
  PRIMARY KEY (`name`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1;";
$logger->debug($query);

$out = $dbh->prepare($query) || $logger->error("Can't prepare statement: $DBI::errstr");
$out->execute || $logger->error($DBI::errstr);


if($delete){
	$query  = "delete from $database.$cdstable";
	$query .= " where name='$transcript'" if $transcript ne "";
	$logger->debug($query);
	
	$out = $dbh->prepare($query) || $logger->error("Can't prepare statement: $DBI::errstr");
	$out->execute || $logger->error($DBI::errstr);
}

$query  = "Select name,chrom,strand,txStart,txEnd,cdsStart,cdsEnd,exonCount,exonStarts,exonEnds 
	FROM $database.$genetable where cdsStart < cdsEnd";
$query .= " and name='$transcript'" if $transcript ne ""; 
$logger->debug($query);

$out = $dbh->prepare($query) || $logger->error("Can't prepare statement: $DBI::errstr");
$out->execute || $logger->error($DBI::errstr);
while (($name,$chrom,$strand,$txStart,$txEnd,$cdsStart,$cdsEnd,$exonCount,$exonStarts,$exonEnds) = $out->fetchrow_array)
{
	(@exonStarts) = split(/\,/,$exonStarts);
	(@exonEnds)   = split(/\,/,$exonEnds);
	$seq = getCds($chrom,$strand,$cdsStart,$cdsEnd,$exonCount,\@exonStarts,\@exonEnds);

	&insert($name,$seq);	
	
} 

$query = "select count(*) from $database.$cdstable where length(cds)%3<>0";
$query .= " and name='$transcript'" if $transcript ne "";
$logger->debug($query);
$out = $dbh->prepare($query) || $logger->error("Can't prepare statement: $DBI::errstr");
$out->execute || $logger->error($DBI::errstr);
if(my $count = $out->fetchrow_array){
	$logger->info("$count transcripts with a length NOT mod 3");
}
my $elapsed = (time - $^T)/60;
$logger->info("finished in $elapsed minutes");


#############
#subroutines#
#############
#get coding sequence

sub getCds {

my $chrom        = shift;
my $strand       = shift;
my $cdsStart     = shift;
my $cdsEnd       = shift;
my $exonCount    = shift;
my $exonStarts   = shift;
my $exonEnds     = shift;

my $cds = "";

# UCSC tracks have 'chr' prefix - but reference might not:
$chrom =~ s/^chr//g if (! Utilities::referenceHasChr($settings) );
 

if ($strand eq "+") 
{
	
	$cdsStart++;
	for (my $i=0;$i<$exonCount;$i++) 
	{  
		$$exonStarts[$i]++;
		#single exon gene
		if (($cdsStart >= $$exonStarts[$i]) and ($cdsEnd <= $$exonEnds[$i])) 
		{ 
			$cds .= Utilities::twoBit($chrom,$cdsStart,$cdsEnd,$strand,$settings);
		}
		#exon with cdsStart
		elsif (($cdsStart >= $$exonStarts[$i]) and ($cdsStart <= $$exonEnds[$i])) 
		{ 
			$cds .= Utilities::twoBit($chrom,$cdsStart,$$exonEnds[$i],$strand,$settings);
		}
		#other coding exons
		elsif (($cdsStart < $$exonStarts[$i]) and ($cdsEnd > $$exonEnds[$i]))  
		{ 
			$cds .= Utilities::twoBit($chrom,$$exonStarts[$i],$$exonEnds[$i],$strand,$settings);
		}
		#exon with cdsEnd
		elsif (($cdsEnd >= $$exonStarts[$i]) and ($cdsEnd <= $$exonEnds[$i]))  
		{ 
			$cds .= Utilities::twoBit($chrom,$$exonStarts[$i],$cdsEnd,$strand,$settings);
		}
	}
}
if ($strand eq "-") 
{
	$cdsStart++;
	
	for (my $i=$exonCount-1;$i>=0;$i--) 
	{ 
		$$exonStarts[$i]++;
		#single exon gene 
		if (($cdsStart >= $$exonStarts[$i]) and ($cdsEnd <= $$exonEnds[$i])) 
		{ 
			$cds .= Utilities::twoBit($chrom,$cdsStart,$cdsEnd,$strand,$settings);
		}
		#exon with cdsStart
		elsif (($cdsStart >= $$exonStarts[$i]) and ($cdsStart <= $$exonEnds[$i])) 
		{ 
			$cds .= Utilities::twoBit($chrom,$cdsStart,$$exonEnds[$i],$strand,$settings);
		}
		#internal coding exons
		elsif (($cdsStart < $$exonStarts[$i]) and ($cdsEnd > $$exonEnds[$i]))  
		{ 
			$cds .= Utilities::twoBit($chrom,$$exonStarts[$i],$$exonEnds[$i],$strand,$settings);
		}
		#exon with cdsEnd
		elsif (($cdsEnd >= $$exonStarts[$i]) and ($cdsEnd <= $$exonEnds[$i]))  
		{ 
			
			$cds .= Utilities::twoBit($chrom,$$exonStarts[$i],$cdsEnd,$strand,$settings);
		}
	}
}

return $cds;

#end sub
}



########
#insert#
########
sub insert {

my $name = shift;
my $cds = shift;

#print "$cds\n";
my $sql = qq{insert into $database.$cdstable values ('$name','$cds') };
$logger->debug($sql);

my $sth = $dbh->prepare( $sql ) || $logger->error("Can't prepare statement: $DBI::errstr");
$sth->execute();


}




=head1 NAME

fillCodingSequenceTable.pl

=head1 SYNOPSIS

fillCodingSequenceTable.pl 

=head1 DESCRIPTION

This script extracts the coding sequence of each coding transcript and stores it in the database for
the use in the annotation script annotateVCF.pl.

=head1 OPTIONS

 -g	genetable; default: take from settings
 -d	database; default: take from settings
 -c	cdstable; default: take from settings
 -t	transcript name; do insert (and delete if -e is defined) only for the given transcript
 -e	delete from cds table before inserting
 -se	name of the settings in the current.config.xml file that holds path to reference genome, 
 	to the annotation file and to possible additional annotation files; use default settings if nothing is given
 -lf	log file; default: print to screen
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page
 
=head1 AUTHOR

Thomas Wieland

=cut


