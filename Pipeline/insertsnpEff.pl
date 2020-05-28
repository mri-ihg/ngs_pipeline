#!/usr/bin/perl

use Getopt::Long qw(:config no_ignore_case);
use Pod::Usage;
use File::Basename;
use Cwd qw(abs_path);

my $prog_path = dirname( abs_path($0) );
require $prog_path."/Utilities.pm";

my $help        = 0;
my $man         = 0;
my $infile      = "";
my $settings    = "";
my $logfile  	= "pipeline.log";
my $loglevel 	= "INFO";

GetOptions(
"i=s" => \$infile,
"se=s"=> \$settings,
"lf=s" => \$logfile,
"ll=s" => \$loglevel,
"h"   => \$help,
"man" => \$man
);

pod2usage( {-exitval => 1  ,-verbose => 1} ) if $infile eq "" || $settings eq "";
pod2usage( {-exitval => 0  ,-verbose => 1} ) if $help;
pod2usage( {-exitval => 0  ,-verbose => 2} ) if $man;

Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();

my $params  	      = Utilities::getParams();
require $params->{programs}->{vcftools}->{pm};

my $exomedb 	      = $params->{settings}->{$settings}->{exomedb}->{database};
my $snvTable          = $params->{settings}->{$settings}->{exomedb}->{snvtable};
my $geneTable         = $params->{settings}->{$settings}->{exomedb}->{genetable};

my $variationdb       = $params->{settings}->{$settings}->{variationdb}->{database};
my $knownGeneSymbol   = $params->{settings}->{$settings}->{variationdb}->{genetable};

exit $logger->error("No snpEff table defined in settings!") unless $params->{settings}->{$settings}->{exomedb}->{snpefftable};
my $snpefftable       = $params->{settings}->{$settings}->{exomedb}->{snpefftable};
my $snpeffeffecttable = $params->{settings}->{$settings}->{exomedb}->{snpeffeffecttable};

my $dbh = Utilities::connectExomeDB($settings);

#open VCF file
my $vcf = Vcf->new( file => $infile );
$vcf->parse_header();


my $snvsNotFound  = 0;
my $genesNotFound = 0;

$logger->info("Inserting snpEff annotations...");
#parse VCF file
while ( my $vcfline = $vcf->next_data_hash() ) {
	
	my $alleleCounter = 0;
	
	#for all alternative alleles
	foreach my $alt_nuc(@{$vcfline->{ALT}}){
		$alleleCounter++;
        next if ( $alt_nuc eq "*");
		
		my $sql = qq{select idsnv from $exomedb.$snvTable where chrom = '$vcfline->{CHROM}' and start = $vcfline->{POS} and allele = '$alt_nuc' and refallele='$vcfline->{REF}'};
		$logger->debug($sql);
		my $sth = $dbh->prepare($sql) || $logger->error("Can't prepare statement: $DBI::errstr");
		$sth->execute() || $logger->error($DBI::errstr);
			
		if(my ($idsnv) = $sth->fetchrow_array()){	#SNV found in database
			my @snpEffEntries = split(",",$vcfline->{INFO}->{EFF});
			foreach(@snpEffEntries){ #parsing single entries
				
				if($_ =~ /(.+)\((.+)\)/ ){
					my $effect = $1;
					my ( $effect_impact, $functional_class, $codon_change, $amino_acid_change, $amino_acid_length, $gene_name, $transcript_biotype, $gene_coding, $transcript_ID, $exon_rank, $genotype_number,@errors) = split(/\|/,$2);
					#skip if not the correct allele
					next if $genotype_number != $alleleCounter;
					
					#get gene id if a gene is defined
					my $sql2 = qq{select g.idgene from $exomedb.$geneTable g inner join $variationdb.$knownGeneSymbol kg on kg.geneSymbol=g.genesymbol where kg.name='$transcript_ID';};
					$logger->debug($sql2);
					my $sth2 = $dbh->prepare($sql2) || $logger->error("Can't prepare statement: $DBI::errstr");
					$sth2->execute() || $logger->error($DBI::errstr);
					my $idgene;	
	
					unless(($idgene) = $sth2->fetchrow_array()){
						$genesNotFound++ if $gene_name ne "";
						$idgene = "NULL";
					}
					
					#insert snpeff entry
					$amino_acid_length = 0 if $amino_acid_length eq "";
					$exon_rank         = 0 if $exon_rank eq "";
					$sql2 = "insert ignore $exomedb.$snpefftable (idsnv,idgene,idsnpeffeffect,transcript,effect_impact,functional_class,codon_change,amino_acid_change,amino_acid_length,exon_rank) values ($idsnv,$idgene,(select idsnpeffeffect from $exomedb.$snpeffeffecttable where name='$effect'),
						                                                                        '$transcript_ID','$effect_impact','$functional_class','$codon_change','$amino_acid_change',$amino_acid_length,$exon_rank);";
					$sth2 = $dbh->prepare($sql2) || $logger->error("Can't prepare statement: $DBI::errstr");
					$logger->debug($sql2);
					$sth2->execute() || $logger->error($DBI::errstr);
				}
			}			
		}else{
			$snvsNotFound++;
		}
	}
}
$logger->info("Finished inserting snpEff entries. $snvsNotFound SNVs not found in $exomedb.$snvTable; $genesNotFound genes not found in $exomedb.$geneTable");



=head1 NAME

 insertnpEff.pl

=head1 SYNOPSIS

 insertsnpEff.pl -i ontarget.snpEff.vcf -se hg19_test

=head1 DESCRIPTION

This script inserts snpEff annotations into the database. As for now it requires the variants and genes to be already
present in their respective tables.

=head1 OPTIONS

 -i	input file; REQUIRED
 -se	settings; REQUIRED
 -lf	<log file>; default: pipeline.log
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	show help text
 -man	show man page

=head1 AUTHOR

Thomas Wieland

=cut
