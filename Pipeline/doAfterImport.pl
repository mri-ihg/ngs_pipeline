#!/usr/bin/perl

use strict;
use Getopt::Long;
use DBI;
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);
use List::MoreUtils qw/ uniq /;
use Pod::Usage;

my $prog_path = dirname( abs_path($0) );
require $prog_path . "/Utilities.pm";





my $noHuman		= 0;
my $snvQual		= 0;
my $gq			= 0;

my $settings    = "default";
my $logfile     = "SCREEN";
my $loglevel    = "INFO";
my $help	    = 0;
my $man         = 0;

GetOptions(
"s=s"  => \$snvQual,
"g=s"  => \$gq,
"nh"   => \$noHuman,
"se=s" => \$settings,
"lf=s" => \$logfile,
"ll=s" => \$loglevel, 
"man"  => \$man,
"h"    => \$help 
);

pod2usage( {-exitval => 0  ,-verbose => 1} ) if $help;
pod2usage( {-exitval => 0  ,-verbose => 2} ) if $man;
#pod2usage( {-exitval => 1  ,-verbose => 1} ) if $infile eq "" || !$possibleCallers{$caller};

my $sth         = "";
my $rv          = "";
my @row         = ();


Utilities::initLogger( $logfile, $loglevel );
my $logger = Utilities::getLogger();
my $params = Utilities::getParams();

my $dbh = Utilities::connectExomeDB($settings);

my $host                      = $params->{settings}->{$settings}->{exomedb}->{host};
my $user                      = $params->{settings}->{$settings}->{exomedb}->{user};
my $passwd                    = $params->{settings}->{$settings}->{exomedb}->{password};
my $variantdb                 = $params->{settings}->{$settings}->{exomedb}->{database};
my $snvTable                  = $params->{settings}->{$settings}->{exomedb}->{snvtable};
my $snvsampleTable            = $params->{settings}->{$settings}->{exomedb}->{snvsampletable};
my $geneTable                 = $params->{settings}->{$settings}->{exomedb}->{genetable};
my $snvgeneTable              = $params->{settings}->{$settings}->{exomedb}->{snvgenetable};
my $snv2diseaseTable          = $params->{settings}->{$settings}->{exomedb}->{snv2diseasetable};
my $snpeffTable               = $params->{settings}->{$settings}->{exomedb}->{snpefftable};
my $additionalannotationtable = $params->{settings}->{$settings}->{exomedb}->{additionalannotationtable};
my $variantStatTable          = $params->{settings}->{$settings}->{exomedb}->{variantstattable};

my $coredb                    = $params->{coredb}->{database};
my $sampleTable               = $params->{coredb}->{sampletable};
my $diseaseTable              = $params->{coredb}->{diseasetable};
my $disease2sampleTable       = $params->{coredb}->{disease2sampletable};

my $refdb                     = $params->{settings}->{$settings}->{variationdb}->{database};
my $dbsnp                     = $params->{settings}->{$settings}->{variationdb}->{snvtable};
my $knownGene                 = $params->{settings}->{$settings}->{variationdb}->{genetable};
my $knownGenePep              = $params->{settings}->{$settings}->{variationdb}->{maxpeplength};
my $tgenomes                  = $params->{settings}->{$settings}->{variationdb}->{tgenomestable};

my $quality = "";
if ($gq == 0) {
	$quality     = "x.snvqual >= $snvQual";
}
elsif ($snvQual == 0) {
	$quality     = "x.gtqual >= $gq";
}
else {
	$quality     = "x.snvqual >= $snvQual AND x.gtqual >= $gq";
}

&snv2diseasegroup();
&maxpeplength();
&frequency();
&average_heterocygosity();
&origin_hapmap();
unless ($noHuman) {
	&th_genomes_dp();
	&th_genomes_af();
}	
&nonsyn_var();
&lof_var();

$logger->info("finished successfully!");


############### snv2diseasegroup #####################

sub snv2diseasegroup {

$logger->info("snv2diseasegroup");

my $tmpExists = 0;
if(-d "tmp"){
	$tmpExists = 1;
}else{
	system("mkdir tmp 2> /dev/null");
}



&snv2diseasegroupRun();

&executeQuery("DROP TABLE IF EXISTS `$snv2diseaseTable`");

&executeQuery("CREATE TABLE `$snv2diseaseTable` (
  `idsnv2diseasegroup` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `fidsnv` int(11) unsigned NOT NULL,
  `fiddiseasegroup` int(11) unsigned NOT NULL,
  `fsample` int(11) NOT NULL DEFAULT '0',
  `fsampleall` int(11) NOT NULL DEFAULT '0',
  `fpedigree` int(11) NOT NULL DEFAULT '0',
  `fpedigreeall` int(11) NOT NULL DEFAULT '0',
  `samplecontrols` int(11) NOT NULL DEFAULT '0',
  `pedigreecontrols` int(11) NOT NULL DEFAULT '0',
  PRIMARY KEY (`idsnv2diseasegroup`),
  UNIQUE KEY `diseasegroupsnv` (`fiddiseasegroup`,`fidsnv`),
  KEY `samplecontrolsgroup` (`fiddiseasegroup`,`samplecontrols`,`fidsnv`),
  KEY `pedigreecontrolsgroup` (`fiddiseasegroup`,`pedigreecontrols`,`fidsnv`),
  KEY `s2didsnvgroup` (`fidsnv`),
  KEY `s2dcontrols` (`samplecontrols`),
  CONSTRAINT `s2diddiseasegroup` FOREIGN KEY (`fiddiseasegroup`) REFERENCES `$coredb`.`diseasegroup` (`iddiseasegroup`) ON DELETE NO ACTION ON UPDATE NO ACTION,
  CONSTRAINT `s2didsnvgroup` FOREIGN KEY (`fidsnv`) REFERENCES `snv` (`idsnv`) ON DELETE NO ACTION ON UPDATE NO ACTION
) ENGINE=InnoDB DEFAULT CHARSET=latin1;");

if($logfile eq "SCREEN"){
	system("mysqlimport -vv -h $host -u $user -p$passwd -L $variantdb tmp/snv2diseasegroup.txt ");
}else{
	system("mysqlimport -vv -h $host -u $user -p$passwd -L $variantdb tmp/snv2diseasegroup.txt >> $logfile");
}


#delete tmp files
#system("rm tmp/snv2diseasegroup*");
#system("rm -r tmp") unless $tmpExists;

}

############### snv2diseasegroup #####################

sub snv2diseasegroupRun {

$logger->info("snv2diseasegroupRun");


&select1;
&select2;
&combine;



}

#####################################################################################
sub combine {
my $file = "";
my %s2d = ();
my ($idsnv,$fsampleall,$fpedigreeall);
my ($idsnv2diseasegroup,$iddiseasegroup,$fsample,$fpedigree,$samplecontrols,$pedigreecontrols);

$file = "tmp/snv2diseasegroup_part2.txt";
open (IN, "$file");
while (<IN>) {
	chomp;
	($idsnv,$fsampleall,$fpedigreeall)=split(/\s+/);
	$s2d{$idsnv}{fsampleall}=$fsampleall;
	$s2d{$idsnv}{fpedigreeall}=$fpedigreeall;
}
close IN;
$file = "tmp/snv2diseasegroup_part1.txt";
open (IN, "$file");
$file = "tmp/snv2diseasegroup.txt";
open (OUT, ">$file");
while (<IN>) {
	chomp;
	($idsnv2diseasegroup,$idsnv,$iddiseasegroup,$fsample,$fpedigree)=split(/\s+/);
	$fsampleall=$s2d{$idsnv}{fsampleall};
	if ($fsampleall eq "") { $logger->error("$idsnv fsampleall = ''"); exit;}
	$fpedigreeall=$s2d{$idsnv}{fpedigreeall};
	if ($fpedigreeall eq "") { $logger->error("$idsnv fpedigreeall = ''"); exit;}
	$samplecontrols=$fsampleall-$fsample;
	$pedigreecontrols=$fpedigreeall-$fpedigree;
	print OUT "$idsnv2diseasegroup\t$idsnv\t$iddiseasegroup\t$fsample\t$fsampleall\t$fpedigree\t$fpedigreeall\t$samplecontrols\t$pedigreecontrols\n";
}
close IN;
close OUT;

}

#####################################################################################
sub select1 {
my $file = "";

my $sth = &executeQuerySth("
SELECT 0,v.idsnv,d.iddiseasegroup,count(v.idsnv),count(DISTINCT s.pedigree)
FROM $snvTable v
INNER JOIN $snvsampleTable x               ON (v.idsnv = x.idsnv) 
INNER JOIN $coredb.$sampleTable s          ON (s.idsample = x.idsample) 
INNER JOIN $coredb.$disease2sampleTable ds ON (s.idsample = ds.idsample) 
INNER JOIN $coredb.$diseaseTable d         ON (ds.iddisease = d.iddisease)
GROUP BY v.idsnv,d.iddiseasegroup
");

$file = "tmp/snv2diseasegroup_part1.txt";
open (OUT, ">$file");
while (@row = $sth->fetchrow_array) {
	print OUT "@row\n";
}
close OUT;

}


#####################################################################################
sub select2 {
my $file = "";

my $sth = &executeQuerySth("
SELECT v.idsnv,count(v.idsnv),count(DISTINCT s.pedigree)
FROM $snvTable v
INNER JOIN $snvsampleTable x               ON (v.idsnv = x.idsnv) 
INNER JOIN $coredb.$sampleTable s          ON (s.idsample = x.idsample) 
INNER JOIN $coredb.$disease2sampleTable ds ON (s.idsample = ds.idsample) 
INNER JOIN $coredb.$diseaseTable d         ON (ds.iddisease = d.iddisease)
GROUP BY v.idsnv
");

$file = "tmp/snv2diseasegroup_part2.txt";
open (OUT, ">$file");
while (@row = $sth->fetchrow_array) {
	print OUT "@row\n";
}
close OUT;

}

############ update maxpeplength #####################

sub maxpeplength {

$logger->info("maxpeplength");

my $rv = &executeQuery("
UPDATE $geneTable g
SET maxpeplength =
(SELECT 
max(length(seq))
FROM $refdb.$knownGene x 
INNER JOIN $refdb.$knownGenePep p on x.name=p.name
WHERE g.genesymbol=x.geneSymbol)
");

$logger->info("updated $rv");

}

############ update frequency information for every sample #####################

sub frequency {

$logger->info("frequency");

my $rv = &executeQuery("
UPDATE $snvTable v
SET freq = (SELECT count(x.idsnv)
FROM $snvsampleTable x 
WHERE v.idsnv = x.idsnv)
");

$logger->info("updated $rv");
}

############ update frequency information for every pedigree #####################

sub frequency_old {

$logger->info("frequency");

my $rv = &executeQuery("
UPDATE $snvTable v
SET freq = (SELECT count(DISTINCT s.pedigree)
FROM $snvsampleTable x 
INNER JOIN $coredb.$sampleTable s ON x.idsample = s.idsample
WHERE v.idsnv = x.idsnv)
");

$logger->info("updated $rv");
}
############ update average heterocygosity information #####################

sub average_heterocygosity {

$logger->info("average_heterocygosity");

my $rv = &executeQuery("
UPDATE
$snvTable e
SET e.avhet = 
(SELECT
max(h.avHet)
FROM
$refdb.$dbsnp h
WHERE
e.rs = h.name
)
");

$logger->info("updated $rv");

}

############ update origin of SNP information (Hapmap) #####################

sub origin_hapmap {

$logger->info("origin_hapmap");
	
my $rv = &executeQuery("
UPDATE
$snvTable e
SET e.valid = 
(SELECT DISTINCT
h.valid
FROM
$refdb.$dbsnp h
WHERE
e.rs = h.name
)
");

$logger->info("updated $rv");

}

############ update 1000 genomes #####################

sub th_genomes_af {

$logger->info("th_genomes_af");

my $rv = &executeQuery("
UPDATE
$snvTable e
SET e.af = 
(SELECT
max(af)
FROM
$refdb.$tgenomes h
WHERE
e.chrom = h.chrom
AND
e.start = h.pos
AND
e.refallele = h.ref
AND
e.allele = h.alt
)
");

$logger->info("updated $rv");

}

############ update 1000 genomes #####################

sub th_genomes_dp {

$logger->info("th_genomes_dp");

my $rv = &executeQuery("
UPDATE
$snvTable e
SET e.dp = 
(SELECT
max(dp)
FROM
$refdb.$tgenomes h
WHERE
e.chrom = h.chrom
AND
e.start = h.pos
AND
e.refallele = h.ref
AND
e.allele = h.alt
)
");

$logger->info("updated $rv");

}

############ update gene table with sum of nonsynonymous variants #####################

sub nonsyn_var {

$logger->info("nonsyn_var");

my $n        = 0;
my $idgene   = "";
my @idgenes  = ();

my $sth =  &executeQuerySth("SELECT idgene FROM gene");

while (@row = $sth->fetchrow_array) {
	push(@idgenes,@row);
}

my $query ="
UPDATE $geneTable g
SET g.nonsynpergene =
(SELECT
count(DISTINCT v.idsnv)
FROM $snvgeneTable y
INNER JOIN $snvTable       v ON y.idsnv = v.idsnv
INNER JOIN $snvsampleTable x ON v.idsnv = x.idsnv
WHERE
((FIND_IN_SET('missense',v.func)> 0)
OR (FIND_IN_SET('nonsense',v.func)> 0)
OR (FIND_IN_SET('stoploss',v.func)> 0)
OR (FIND_IN_SET('splice',v.func)> 0)
OR (FIND_IN_SET('frameshift',v.func)> 0)
OR (FIND_IN_SET('indel',v.func)> 0))
AND 
(
(v.class='snp')
OR (v.class='indel')
)
AND g.idgene = y.idgene
AND x.filter = 'PASS'
AND $quality
AND x.mapqual >= 50
GROUP BY
y.idgene
)
WHERE g.idgene = ?
";

my $sth = $dbh->prepare($query) || $logger->error("Can't prepare statement: $DBI::errstr");

foreach $idgene (@idgenes) {
	#print  "$idgene\n";
	$sth->execute($idgene);
	$n++;
}

$logger->info("updated $n");

}
################ update gene table with sum of LoF variants ###########################

sub lof_var {

$logger->info("lof_var");

my $rv = &executeQuery("
UPDATE $geneTable g
SET g.delpergene =
(SELECT
count(distinct v.idsnv)
FROM $snvgeneTable y
INNER JOIN $snvTable       v ON y.idsnv = v.idsnv
INNER JOIN $snvsampleTable x ON v.idsnv = x.idsnv
WHERE
(
(FIND_IN_SET('nonsense',v.func)> 0)
OR (FIND_IN_SET('splice',v.func)> 0)
OR (FIND_IN_SET('frameshift',v.func)> 0)
)
AND
(
(v.class='snp')
OR (v.class='indel')
)
AND g.idgene = y.idgene
AND $quality
AND x.mapqual >= 50
AND x.filter = 'PASS')
");


$logger->info("updated = $rv");

}


######################################### execute query ##############################################
sub executeQuery {
	my $sql   = shift;
	my $value = shift;
	$logger->debug($sql);
    my $sth = $dbh->prepare($sql) || $logger->error("Can't prepare statement: $DBI::errstr");
    if($value){
    	return $sth->execute($value) || $logger->error($DBI::errstr);
    }
    return $sth->execute() || $logger->error($DBI::errstr);
}


sub executeQuerySth {
	my $sql = shift;
	$logger->debug($sql);
    my $sth = $dbh->prepare($sql) || $logger->error("Can't prepare statement: $DBI::errstr");
    $sth->execute() || $logger->error($DBI::errstr);
    return $sth;
}


=head1 NAME

doAfterImport.pl

=head1 SYNOPSIS

doAfterImport.pl

=head1 DESCRIPTION

This script must be run after importing a batch of new variants into the database. It calculates
sums and additional annotations that are required for some web-interface queries to run properly.

NOTE: This script creates several temp files in the current working directory, thus write permissions are required.

=head1 OPTIONS

 -s	[snvQuality] minimum SNVquality (SAMtools) or QD (GATK) for variants to be included in sums; default: 0
 -g	[genotypeQuality] minimum genotype quality (GQ) for variants to be included in sums; default: 0
 -nh	[not_human] choose if the database defined in -se does NOT hold human variants --> no human specific
 		databases (e.g. 1000 genomes) are queried
 -se	name of the settings in the current.config.xml file that holds path to reference genome, 
 	to the annotation file and to possible additional annotation files; use default settings if nothing is given
 -lf	log file; default: print to screen
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page

=head1 AUTHOR

Tim M. Strom, Thomas Wieland

=cut
