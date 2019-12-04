#!/usr/bin/perl

use strict;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;
use File::Copy;
use Cwd qw(abs_path);
use Data::Dumper;
use IO::File;
use Pod::Usage;
umask(002);

my $prog_path = dirname( abs_path($0) );
require $prog_path."/Utilities.pm";


my $bamfile = "";
my $help	= 0;
my $man     = 0;
my $logfile  	  = "SCREEN";
my $loglevel 	  = "INFO";
my $settings  = "";
my $sample    = "";
my $delete    = 0;

GetOptions(
"p=s"  => \$bamfile,
"se=s" => \$settings,
"s=s"  => \$sample,
"d"    => \$delete,
"lf=s" => \$logfile,
"ll=s" => \$loglevel,
"man"  => \$man,
"h"	   => \$help
);

pod2usage( {-exitval => 0  ,-verbose => 1} ) if $help;
pod2usage( {-exitval => 0  ,-verbose => 2} ) if $man;
pod2usage( {-exitval => 1  ,-verbose => 1} ) if $bamfile eq "" || $settings eq "" || $sample eq "";

Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();

my $params 		  = Utilities::getParams();

my $dbh          = Utilities::connectExomeDB($settings);
my $vardb        = $params->{settings}->{$settings}->{variationdb}->{database};
my $varGeneTable = $params->{settings}->{$settings}->{variationdb}->{genetable};

my $coredb       = $params->{coredb}->{database};
my $sampleTable  = $params->{coredb}->{sampletable};

unless($params->{settings}->{$settings}->{exomedb}->{translocationtable}){
	$logger->error("No translocation table defined for settings $settings - Exiting!");
	exit(1);
}

my $translocationTable 	= $params->{settings}->{$settings}->{exomedb}->{translocationtable};
my $geneTable 			= $params->{settings}->{$settings}->{exomedb}->{genetable};

#get idsample
my $sql = qq{select s.idsample from $coredb.$sampleTable s where s.name = '$sample' };
$logger->debug("Query: $sql");
my $sth = $dbh->prepare($sql)
	  || $logger->error("Can't prepare statement: $DBI::errstr");
$sth->execute() || $logger->error($DBI::errstr);
my $sampleId;
unless (  $sampleId = $sth->fetchrow_array() ) {
	$logger->error(
		"patient name $sample not found in $coredb.$sampleTable - exiting");
	exit(1);
}

if($delete){
	$sql = qq{delete from $translocationTable where idsample=$sampleId };
	$logger->debug("Query: $sql");
	$sth = $dbh->prepare($sql) || $logger->error("Can't prepare statement: $DBI::errstr");
	$sth->execute() || $logger->error($DBI::errstr);
}

my %infiles;
$infiles{$bamfile.".unbalanced.txt"} = "unbalanced";
$infiles{$bamfile.".balanced.txt"} = "balanced";
$infiles{$bamfile.".unclassified.txt"} = "unclassified";
$infiles{$bamfile.".inter_insertions.txt"} = "inter_insertions";

foreach(keys %infiles){
	open IN, $_ or exit $logger->error("Can't open $_");
	my $tmp = <IN>;
	$tmp    = <IN>;
	while(my $line = <IN>){
		next if $line =~ /^-/;
		chomp $line;
		next if $line =~ /^$/;
		
		my @columns = split("\t",$line);
		my ($gene1,$inGene1) = &findGene($columns[0],$columns[1]);
		my ($gene2,$inGene2) = &findGene($columns[2],$columns[3]);
		if($gene1 ne "NULL" && $gene2 ne "NULL"){
			
			#insert the translocation
			$sql = "insert into $translocationTable (idsample,chrom1,pos1,chrom2,pos2,idgene1,isingene1,idgene2,isingene2,varianttype,num_discordant,sc1,sc2) 
			values ($sampleId,'$columns[0]',$columns[1],'$columns[2]',$columns[3],(select idgene from $geneTable where geneSymbol='$gene1'),$inGene1,(select idgene from $geneTable where geneSymbol='$gene2'),$inGene2,'".$infiles{$_}."',$columns[6],$columns[7],$columns[8] );";
			$logger->debug("Query: $sql");
			#print $sql."\n";
			$sth = $dbh->prepare($sql) 	  || $logger->error("Can't prepare statement: $DBI::errstr");
			$sth->execute() || $logger->error($DBI::errstr);
my $sampleId;
		}
	}
	
	close IN;
}






###################################################################################
sub findGene {
	my $chrom = shift;
	my $pos   = shift;
	
	my $query = "Select geneSymbol FROM $vardb.$varGeneTable					
	WHERE chrom='$chrom' 
	AND $pos >= txStart  
	AND $pos <= txEnd 
	";
	$logger->debug( "Query: " . $query );
	my $out = $dbh->prepare($query) || exit $logger->error($DBI::errstr); #get overlapping genes
	$out->execute || exit $logger->error($DBI::errstr);
	
	if( my $gene = $out->fetchrow_array){				#there is an overlapping gene -> return it
		return ($gene,1);						
	}
	
	#search nearest gene if no gene overlaps
	
	#get upstream
	$query = "Select geneSymbol, $pos-txEnd FROM $vardb.$varGeneTable					
	WHERE chrom='$chrom' 
	AND $pos > txEnd
	ORDER BY txEnd DESC
	LIMIT 1;
	";
	$logger->debug( "Query: " . $query );
	$out = $dbh->prepare($query) || exit $logger->error($DBI::errstr); #get overlapping genes
	$out->execute || exit $logger->error($DBI::errstr);
	my ($usGene,$usDist) = $out->fetchrow_array;
	
	
	#get downstream
	$query = "Select geneSymbol, txStart-$pos FROM $vardb.$varGeneTable					
	WHERE chrom='$chrom' 
	AND $pos < txStart
	ORDER BY txStart
	LIMIT 1;
	";
	$logger->debug( "Query: " . $query );
	$out = $dbh->prepare($query) || exit $logger->error($DBI::errstr); #get overlapping genes
	$out->execute || exit $logger->error($DBI::errstr);
	my ($dsGene,$dsDist) = $out->fetchrow_array;
	
	return ("NULL",0) unless defined($usGene) && defined($dsGene);
	return ($dsGene,0) if !defined($usGene) ||  $dsDist<$usDist;
	return ($usGene,0);
	
	
}







=head1 NAME

insertBellerophon.pl

=head1 SYNOPSIS

insertBellerophon.pl -p merged.rmdup.bam -se hg19_plus

=head1 DESCRIPTION

This script inserts translocations called by bellerophon into the database


=head1 OPTIONS

 -p	</path/to/outfile_prefix> usually the path of the bam file bellerophon was run on; required
 -se	settings; required
 -s	samplename; required
 -d	delete entries for given sample before inserting
 -lf	<log file>; default: pipeline.log
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page

=head1 AUTHOR

Thomas Wieland

=cut