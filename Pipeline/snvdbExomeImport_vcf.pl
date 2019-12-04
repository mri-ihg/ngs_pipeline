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
require $prog_path . "/BEDRecord.pm";
require $prog_path . "/SequenceDictionary.pm";

my $run = 1;

# database
my $dbh           = "";
my $sql           = "";
my $sth           = "";
my $logfile       = "SCREEN";
my $loglevel      = "INFO";
my $duplicateSNPs = 0;
my $delete        = 0;
my $deleteChrom   = "";

my $infile    = "";
my $line      = "";
my $settings  = "default";
my $help      = 0;
my $man       = 0;
my $rollbacks = 0;

my $preClass = "";
my $caller   = "samtools";

my $minOverlap = 88;    #min overlap in %

GetOptions(
	"i=s"  => \$infile,
	"d"    => \$delete,
	"dc=s" => \$deleteChrom,
	"se=s" => \$settings,
	"l=s"  => \$preClass,
	"c=s"  => \$caller,
	"lf=s" => \$logfile,
	"ll=s" => \$loglevel,
	"h"    => \$help,
	"man"  => \$man
);

pod2usage( { -exitval => 0, -verbose => 1 } ) if $help;
pod2usage( { -exitval => 0, -verbose => 2 } ) if $man;
pod2usage( { -exitval => 1, -verbose => 1 } ) if $infile eq "";

my $params = Utilities::getParams();

require $params->{programs}->{vcftools}->{pm};    #VCFTools library to parse VCF file

my $exomedb                   = $params->{settings}->{$settings}->{exomedb}->{database};
my $snvTable                  = $params->{settings}->{$settings}->{exomedb}->{snvtable};
my $snvsampleTable            = $params->{settings}->{$settings}->{exomedb}->{snvsampletable};
my $geneTable                 = $params->{settings}->{$settings}->{exomedb}->{genetable};
my $snvgeneTable              = $params->{settings}->{$settings}->{exomedb}->{snvgenetable};
my $snv2diseaseTable          = $params->{settings}->{$settings}->{exomedb}->{snv2diseasetable};
my $snpeffTable               = $params->{settings}->{$settings}->{exomedb}->{snpeffeffecttable};
my $additionalannotationtable = $params->{settings}->{$settings}->{exomedb}->{additionalannotationtable};
my $variantStatTable          = $params->{settings}->{$settings}->{exomedb}->{variantstattable};

my $coredb                    = $params->{coredb}->{database};
my $sampleTable               = $params->{coredb}->{sampletable};


Utilities::initLogger( $logfile, $loglevel );
my $logger = Utilities::getLogger();

$dbh = Utilities::connectExomeDB($settings);

my $sampleDir = dirname( $infile );
my $fileName  = basename( $infile );

#open VCF file
my $vcf = Vcf->new( file => $infile );
$vcf->parse_header();

my @samples = $vcf->get_samples(); #get sample ids for all samples  in VCF files
my %sampleIds;

## prepare mysql-statments ##
#checkSampleIds 
my $checkSampleId_sql = qq{select idsample from $coredb.$sampleTable where name = ? };
my $checkSampleId_sth = $dbh->prepare($checkSampleId_sql) || $logger->error("Can't prepare statement: $DBI::errstr");
#delete snvsample-entries for sampleid
my $deleteSnvSample_sql = qq{delete from $exomedb.$snvsampleTable where idsample=? and caller=?};
my $deleteSnvSample_sth = $dbh->prepare($deleteSnvSample_sql) || $logger->error("Can't prepare statement: $DBI::errstr");
#delete variantstat-entriy for sampleid
my $deleteVariantStat_sql = qq{delete from $exomedb.$variantStatTable where idsample=?};
my $deleteVariantStat_sth = $dbh->prepare($deleteVariantStat_sql) || $logger->error("Can't prepare statement: $DBI::errstr");
#delete snvsample-chromosome-entries for sampleid
my $deleteSnvSampleChrom_sql = qq{delete ss.* from $exomedb.$snvsampleTable ss inner join $exomedb.$snvTable v on ss.idsnv=v.idsnv where ss.idsample=? and ss.caller=? and v.chrom=?};
my $deleteSnvSampleChrom_sth = $dbh->prepare($deleteSnvSampleChrom_sql) || $logger->error("Can't prepare statement: $DBI::errstr");
###
###
###
#check if SNV already exists
my $checkSNVexistance_sql = qq{select idsnv from $exomedb.$snvTable where chrom = ? and start = ? and allele = ? and refallele=?};
my $checkSNVexistance_sth = $dbh->prepare($checkSNVexistance_sql) || $logger->error("Can't prepare statement: $DBI::errstr");
#check SNV
my $checkSNV_sql = qq{select chrom,start,end,idsnv from $exomedb.$snvTable where chrom=? AND (end-start)>1 AND find_in_set(?,class) AND ((start>=? and start<=?) OR  (end>=? and end<=?) OR (start<? and end>?)) };
$checkSNV_sql .= " AND allele=?"; #TW 09.10.2015: 
my $checkSNV_sth = $dbh->prepare($checkSNV_sql) || $logger->error("Can't prepare statement: $DBI::errstr");
#get Isoform
my $getIsoform_sql = qq{select transcript,func from $exomedb.$snvTable where idsnv=?};
my $getIsoform_sth = $dbh->prepare($getIsoform_sql) || $logger->error("Can't prepare statement: $DBI::errstr");
#update idsnv					
my $updateIdsnv_sql = qq{update IGNORE $exomedb.$snvsampleTable set idsnv=? where idsnv=?};
my $updateIdsnv_sth = $dbh->prepare($updateIdsnv_sql) || $logger->error("Can't prepare statement: $DBI::errstr");
#updateAdditionalAnnotation
my $updateAdditionalAnnotation_sql = qq{update IGNORE $exomedb.$additionalannotationtable set idsnv=? where idsnv=?};
my $updateAdditionalAnnotation_sth = $dbh->prepare($updateAdditionalAnnotation_sql) || $logger->error("Can't prepare statement: $DBI::errstr");
#delete violentSamples
my $deleteViolentSamples_sql = qq{delete from $exomedb.$snvsampleTable where idsnv=?};
my $deleteViolentSamples_sth = $dbh->prepare($deleteViolentSamples_sql) || $logger->error("Can't prepare statement: $DBI::errstr");
#delete old snv gene
my $deleteOldSnvGene_sql = qq{delete from $exomedb.$snvgeneTable where idsnv=?};
my $deleteOldSnvGene_sth = $dbh->prepare($deleteOldSnvGene_sql) || $logger->error("Can't prepare statement: $DBI::errstr");
#delete old snv2disease
my $deleteOldSnv2disease_sql = qq{delete from $exomedb.$snv2diseaseTable where fidsnv=?};
my $deleteOldSnv2disease_sth = $dbh->prepare($deleteOldSnv2disease_sql) || $logger->error("Can't prepare statement: $DBI::errstr");
#delete old snpeff table
my $deleteOldSnpEff_sql = qq{delete from $exomedb.$snpeffTable where idsnv=?};
my $deleteOldSnpEff_sth = $dbh->prepare($deleteOldSnpEff_sql) || $logger->error("Can't prepare statement: $DBI::errstr");
#delete old snv
my $deleteOldSnv_sql = qq{delete from $exomedb.$snvTable where idsnv=?};
my $deleteOldSnv_sth = $dbh->prepare($deleteOldSnv_sql) || $logger->error("Can't prepare statement: $DBI::errstr");
#change Transcript String
my $changeTranscriptString_sql = qq{update $exomedb.$snvTable set transcript=?,end=?,func=?,length=? where idsnv=?};
my $changeTranscriptString_sth = $dbh->prepare($changeTranscriptString_sql)  || $logger->error("Can't prepare statement: $DBI::errstr");
#insert SNV
my $insertSNV_sql = qq{insert into $exomedb.$snvTable (idsnv,chrom,start,      end ,rs    ,allele ,class    ,func    ,transcript  ,freq,clinical,avhet,valid,dp,af,refallele,length) values (idsnv,?,?,?,?,?,?,?,?,'0','0','0','','0','0',?,?) };
my $insertSNV_sth = $dbh->prepare($insertSNV_sql) || $logger->error("Can't prepare statement: $DBI::errstr");
#check for ref allele
my $checkForRefAllele_sql = qq{select idsnv from $exomedb.$snvTable where chrom = ? and start = ? and end = ? and class = ? and allele = ? AND refallele=? };
my $checkForRefAllele_sth = $dbh->prepare($checkForRefAllele_sql) || $logger->error("Can't prepare statement: $DBI::errstr");
#insert into snv gene
my $insertSnvGene_sql = qq{INSERT IGNORE INTO $exomedb.$snvgeneTable (idsnv,idgene) VALUES (?,(SELECT min(idgene) FROM $geneTable WHERE geneSymbol=? ) )};
my $insertSnvGene_sth = $dbh->prepare($insertSnvGene_sql) || $logger->error("Can't prepare statement: $DBI::errstr");
#add additional annotation
my $insertAdditionalAnnotation_sql = qq{insert into $exomedb.$additionalannotationtable (idsnv,annotationname,idannotation) values (?,?,?); };
my $insertAdditionalAnnotation_sth = $dbh->prepare($insertAdditionalAnnotation_sql) || $logger->error("Can't prepare statement: $DBI::errstr");
##############################

foreach my $sample (@samples) {

	#check if sample id is valid
	$checkSampleId_sth->execute($sample) || $logger->error($DBI::errstr);
	$sampleIds{$sample} = $checkSampleId_sth->fetchrow_array();
	
	if ( !defined( $sampleIds{$sample} ) ) {
		$logger->error("patient name $sample not found in $coredb.$sampleTable - exiting");
		exit(1);
	}

	#delete from snvsample table
	if ($delete) {
		$deleteSnvSample_sth->execute($sampleIds{$sample}, $caller) || $logger->error($DBI::errstr);
		
		#&updateVariantStat($sampleIds{$sample}); #TS 09.02.2017: for now just delete the entire entry and recalculate it afterwards
		$deleteVariantStat_sth->execute($sampleIds{$sample}) || $logger->error($DBI::errstr);
		
	}
	elsif ( $deleteChrom ne "" ) {
		$deleteSnvSampleChrom_sth->execute($sampleIds{$sample},$caller,$deleteChrom) || $logger->error($DBI::errstr);
		#&updateVariantStat($sampleIds{$sample}); #TS 09.02.2017: for now just delete the entire entry and recalculate it afterwards
	}
}
my $counter = 0;
#print "\n";




open MYSQLIMPORTFILE, ">$sampleDir/snvsample.".$samples[0].".$fileName.tsv"|| $logger->error("cannot open outputfile $infile.mysqlimport.csv");
while ( my $vcfline = $vcf->next_data_hash() ) {
	&insertSnv( $vcfline, \%sampleIds );
}
close MYSQLIMPORTFILE;

#foreach my $currId(values %sampleIds){
#	&updateVariantStat($currId);
#}

$logger->info("$duplicateSNPs positions with more than one corresponding rsSNPs found");
$logger->info("$rollbacks SNVs couldn't be inserted because they were inserted by another job in the mean time.");
my $elapsed = ( time - $^T ) / 60;
$logger->info("finished in $elapsed minutes");

#############
#subroutines#
#############
###########################
#insertSnv: fill databases#
###########################
sub insertSnv {

	my $vcfline   = shift;
	my $pointer   = shift;
	my %sampleIds = %$pointer;
	my $pos       = $vcfline->{POS};

	my $class     = "snp";
	my $alt_nuc   = $vcfline->{ALT}[0];
	my $nucChange = "";
	my $length    = 1;
	my $reflength = 0;
	my $altlength = 0;
	
	return if defined $vcfline->{INFO}->{SVTYPE} && $vcfline->{INFO}->{SVTYPE} eq "INV"; #TW 28.10.2015: problems with inversions --> do not insert them for now into the normal snv table

	if ( $caller eq "exomedepth" || ($vcfline->{INFO}->{SVTYPE} && !($vcfline->{REF} =~/[acgt]+/i && $alt_nuc =~/[acgt]+/i)) ) {		#TW 09.10.2015: pindel is a special case, because it has SVTYPE but usually also gives 
		
		if ( $vcfline->{REF} eq "LD" ) {
			$length = $alt_nuc;
			$reflength = $length;
			$altlength = 0;
			
			$vcfline->{REF} = "N";	#TW: 09.10.2015: change VCF format of SVs to actual representation from the VCF specification
			$alt_nuc = "<DEL>";
			
		} elsif ( $vcfline->{REF} eq "LI" ) {
			$length = $alt_nuc;
			$reflength = 0;
			$altlength = $length;
			
			$vcfline->{REF} = "N";	#TW: 09.10.2015: change VCF format of SVs to actual representation from the VCF specification
			$alt_nuc = "<DUP>";
		} else {
			$length         = $vcfline->{INFO}->{END} - $vcfline->{POS};
			$vcfline->{REF} = "N";
			$alt_nuc        = "<".$vcfline->{INFO}->{SVTYPE}.">";
			if($vcfline->{INFO}->{SVTYPE} eq "DEL"){
				$reflength = $length;
				$altlength = 0;
			}elsif($vcfline->{INFO}->{SVTYPE} eq "DUP"){ 		
				$reflength = 0;
				$altlength = $length;
			}elsif($vcfline->{INFO}->{SVTYPE} eq "INS"){
				$reflength = 0;
				$altlength = $vcfline->{INFO}->{SVLEN};
				$length    = 1;
			}
			elsif($vcfline->{INFO}->{SVTYPE} eq "INV"){
				$reflength = $length;
				$altlength = $length;
			}else{
				return;						#currently only the SVTYPES above are supported
			}
		}
	} else {
		( $vcfline->{REF}, $alt_nuc ) =  Utilities::getMinimalIndel( $vcfline->{REF}, $alt_nuc ); #TW: 01.04.2015: minimize Indels before importing them into the database! VCF files are NOT minimized anymore
		$reflength = length $vcfline->{REF};
		$altlength = length $alt_nuc;

		#parse indel
		if ( $reflength != 1 || $altlength != 1 ) {
			$class = "indel";

			#get length for deletions
			if ( $reflength > $altlength ) {
				$length = $reflength - $altlength;
			}
			if ( $reflength > 255 || $altlength > 255 ) {
				if ( $reflength > $altlength ) {
					#$vcfline->{REF} = "LD";    #long deletion that doesn't fit into database
					#$alt_nuc = $reflength - $altlength;
					$vcfline->{REF} = "N";	#TW: 09.10.2015: change VCF format of SVs to actual representation from the VCF specification
					$alt_nuc = "<DEL>";
				} else {
					#$vcfline->{REF} =
					#  "LI";    #long insertion that doesn't fit into database
					#$alt_nuc = $altlength - $reflength;
					$vcfline->{REF} = "N";	#TW: 09.10.2015: change VCF format of SVs to actual representation from the VCF specification
					$alt_nuc = "<INS>";
				}
			}
		}
		$nucChange = ",".$vcfline->{REF}."->".$alt_nuc;
	}
	$class = $preClass if $preClass ne "";

	my $end = $pos + $length;

	#parse rs number from rs string
	my $rs = $vcfline->{ID};
	if ( $rs eq "." ) {
		$rs = "";
	} else {
		$rs =~ /^(rs\d+),*/;
		$rs = $1;
	}

	#parse some of the INFO fields
	my $refScore = $vcfline->{INFO}->{FQ};
	$refScore = $vcfline->{INFO}->{QD} if $caller eq "gatk";
	$refScore = 999 if ( $caller eq "pindel" || $caller eq "exomedepth" );
	$refScore = 0 unless $refScore;
	if ( $refScore < 0 ) {
		$refScore *= -1;
	}
	my $snvScore = $refScore;
	my ( $rf, $rb, $sf, $sb ) = ( 0, 0, 0, 0 );
	( $rf, $rb, $sf, $sb ) = split( ",", $vcfline->{INFO}->{DP4} ) if $vcfline->{INFO}->{DP4};
	my $varPercent = 0;
	$varPercent = ( $sf + $sb ) / ( $sf + $sb + $rf + $rb )        if ( $sf + $sb + $rf + $rb ) != 0;    #note that these values only include high quality reads!
	$varPercent *= 100;

	$varPercent = $vcfline->{INFO}->{EXONS} if $caller eq "exomedepth";
	my $strandPercent = 0;
	if ( $sf + $sb != 0 ) {
		if ( $sf > $sb ) {
			$strandPercent = $sb / ( $sf + $sb );
		} else {
			$strandPercent = $sf / ( $sf + $sb );
		}
	}
	$strandPercent = $vcfline->{INFO}->{SF}	  if $vcfline->{INFO}->{SF} && $caller eq "gatk"; #TODO SF value is calculated by filterSNPqual.pl --> maybe use this one for samtools as well???
	$strandPercent *= 100;
	$strandPercent = $vcfline->{INFO}->{RATIO} if $caller eq "exomedepth";

	my $filter = join( ',', @{ $vcfline->{FILTER} } );
	if ( $filter eq "." ) {
		$filter = "PASS";
	}

	my $medianVarQual = 0;
	if ( $vcfline->{INFO}->{MED} ) {
		$medianVarQual = $vcfline->{INFO}->{MED};
	}

	my $mq = 999;
	$mq = $vcfline->{INFO}->{MQ} if $vcfline->{INFO}->{MQ};

	if ($run) {
		#check if entry in snv already exists
		$checkSNVexistance_sth->execute($vcfline->{CHROM},$pos,$alt_nuc,$vcfline->{REF}) || $logger->error($DBI::errstr);
		my ($idsnv) = $checkSNVexistance_sth->fetchrow_array();

		#if no entry -> insert
		if (!defined($idsnv) ) {
			#build isoform string
			my $isoform = "";
			my $dbindel = "";

			my @allgenes;

			if ( $vcfline->{INFO}->{DBINDEL} ) {

				$dbindel = "," . $vcfline->{INFO}->{DBPOS} . ",";
				if ( length( $vcfline->{INFO}->{DBINDEL} ) < 22 ) {
					$dbindel .= $vcfline->{INFO}->{DBINDEL};
				}
				else {
					$dbindel .= "indellength=" . ( length( $vcfline->{INFO}->{DBINDEL} ) - 1 );
				}
			}
			if ( $reflength != 1 || $altlength != 1 ) {
				if ( $reflength > $altlength ) {
					$dbindel .= ",del";
				} elsif ( $reflength < $altlength ) {
					$dbindel .= ",ins";
				}else {
					$dbindel .= ",inv";
				}
			}

			#CDS
			if ( $vcfline->{INFO}->{CDSTRANS} ) {
				my @cdstrans = split( ",", $vcfline->{INFO}->{CDSTRANS} );
				my @cdsgene  = split( ",", $vcfline->{INFO}->{CDSGENE} );
				push( @allgenes, @cdsgene );
				my @cdsstrand = split( ",", $vcfline->{INFO}->{CDSSTRAND} );
				my @cdsnuc    = split( ",", $vcfline->{INFO}->{CDSNUC} );
				my @cdsframe  = split( ",", $vcfline->{INFO}->{CDSFRAME} );
				my @cdsaa;
				unless ( $vcfline->{INFO}->{INDEL} ) {
					@cdsaa = split( ",", $vcfline->{INFO}->{CDSAA} );
				}
				my @cdsfunc = split( ",", $vcfline->{INFO}->{CDSFUNC} );
				for ( my $i = 0 ; $i < @cdstrans ; $i++ ) {

					$isoform .= "$cdstrans[$i],$cdsgene[$i],$cdsstrand[$i],$cdsnuc[$i],$cdsframe[$i],$cdsfunc[$i]$dbindel";
					if ( @cdsaa && $cdsaa[$i] ) {
						$isoform .= ",$cdsaa[$i]";
					}
					$isoform .= ": ";
				}
			}

			#noncoding
			if ( $vcfline->{INFO}->{NCTRANS} ) {
				my @nctrans = split( ",", $vcfline->{INFO}->{NCTRANS} );
				my @ncgene  = split( ",", $vcfline->{INFO}->{NCGENE} );
				push( @allgenes, @ncgene );
				my @ncstrand = split( ",", $vcfline->{INFO}->{NCSTRAND} );
				for ( my $i = 0 ; $i < @nctrans ; $i++ ) {
					$isoform .= "$nctrans[$i],$ncgene[$i],$ncstrand[$i],noncoding$dbindel$nucChange: ";
				}
			}

			#5utr
			if ( $vcfline->{INFO}->{UTR5TRANS} ) {
				my @utr5trans = split( ",", $vcfline->{INFO}->{UTR5TRANS} );
				my @utr5gene  = split( ",", $vcfline->{INFO}->{UTR5GENE} );
				push( @allgenes, @utr5gene );
				my @utr5strand = split( ",", $vcfline->{INFO}->{UTR5STRAND} );
				for ( my $i = 0 ; $i < @utr5trans ; $i++ ) {
					$isoform .= "$utr5trans[$i],$utr5gene[$i],$utr5strand[$i],5utr$dbindel$nucChange: ";
				}
			}

			#3utr
			if ( $vcfline->{INFO}->{UTR3TRANS} ) {
				my @utr3trans = split( ",", $vcfline->{INFO}->{UTR3TRANS} );
				my @utr3gene  = split( ",", $vcfline->{INFO}->{UTR3GENE} );
				push( @allgenes, @utr3gene );
				my @utr3strand = split( ",", $vcfline->{INFO}->{UTR3STRAND} );
				for ( my $i = 0 ; $i < @utr3trans ; $i++ ) {
					$isoform .= "$utr3trans[$i],$utr3gene[$i],$utr3strand[$i],3utr$dbindel$nucChange: ";
				}
			}

			#direct splice site hit
			if ( $vcfline->{INFO}->{DSTRANS} ) {
				my @dstrans = split( ",", $vcfline->{INFO}->{DSTRANS} );
				my @dsgene  = split( ",", $vcfline->{INFO}->{DSGENE} );
				push( @allgenes, @dsgene );
				my @dsstrand = split( ",", $vcfline->{INFO}->{DSSTRAND} );
				my @dsmut;
				if ( $vcfline->{INFO}->{DSMUT} ) {
					@dsmut = split( ",", $vcfline->{INFO}->{DSMUT} );
				}
				for ( my $i = 0 ; $i < @dstrans ; $i++ ) {
					$isoform .= "$dstrans[$i],$dsgene[$i],$dsstrand[$i]$dbindel";
					if ( @dsmut && $dsmut[$i] ) {
						$isoform .= ",$dsmut[$i]";
					}
					$isoform .= ",direct splicesite: ";
				}
			}

			#near splice site hit
			if ( $vcfline->{INFO}->{NSTRANS} ) {
				my @nstrans = split( ",", $vcfline->{INFO}->{NSTRANS} );
				my @nsgene  = split( ",", $vcfline->{INFO}->{NSGENE} );
				push( @allgenes, @nsgene );
				my @nsstrand = split( ",", $vcfline->{INFO}->{NSSTRAND} );
				for ( my $i = 0 ; $i < @nstrans ; $i++ ) {
					$isoform .= "$nstrans[$i],$nsgene[$i],$nsstrand[$i]$dbindel$nucChange,near splicesite: ";
				}
			}

			#intronic
			if ( $vcfline->{INFO}->{INTTRANS} ) {
				my @inttrans = split( ",", $vcfline->{INFO}->{INTTRANS} );
				my @intgene  = split( ",", $vcfline->{INFO}->{INTGENE} );
				push( @allgenes, @intgene );
				my @intstrand = split( ",", $vcfline->{INFO}->{INTSTRAND} );
				for ( my $i = 0 ; $i < @inttrans ; $i++ ) {
					$isoform .= "$inttrans[$i],$intgene[$i],$intstrand[$i],intronic$nucChange: ";
				}
			}

			#CNV
			if ( $vcfline->{INFO}->{CNVTRANS} ) {
				my @cnvtrans = split( ",", $vcfline->{INFO}->{CNVTRANS} );
				my @cnvgene  = split( ",", $vcfline->{INFO}->{CNVGENE} );
				push( @allgenes, @cnvgene );
				my @cnvstrand = split( ",", $vcfline->{INFO}->{CNVSTRAND} );
				my $cnvclass;
				$cnvclass = $vcfline->{INFO}->{CLASS} if $vcfline->{INFO}->{CLASS};
				$cnvclass = $vcfline->{INFO}->{SVTYPE} if $vcfline->{INFO}->{SVTYPE};
				for ( my $i = 0 ; $i < @cnvtrans ; $i++ ) {
					$isoform .= "$cnvtrans[$i],$cnvgene[$i],$cnvstrand[$i],cnv,$cnvclass: ";
				}
			}

			#upstream & downstream variants
			if ( $vcfline->{INFO}->{UPSTRANS} ) {
				$isoform .= "$vcfline->{INFO}->{UPSTRANS},$vcfline->{INFO}->{UPSGENE},$vcfline->{INFO}->{UPSSTRAND},distance-$vcfline->{INFO}->{UPSDIST}, upstream gene$nucChange: ";
				push( @allgenes, $vcfline->{INFO}->{UPSGENE} );
			}
			if ( $vcfline->{INFO}->{DOSTRANS} ) {
				$isoform .= "$vcfline->{INFO}->{DOSTRANS},$vcfline->{INFO}->{DOSGENE},$vcfline->{INFO}->{DOSSTRAND},distance-$vcfline->{INFO}->{DOSDIST}, downstream gene$nucChange: ";
				push( @allgenes, $vcfline->{INFO}->{DOSGENE} );
			}

			my $insertVariant = 1;

			if ( ( $caller eq "pindel" || $caller eq "exomedepth" || $vcfline->{INFO}->{SVTYPE}) && $length > 1 ) { #if the current variant is a deletion called by pindel (07.08.2013: OR exomeDepth): look in the database for other long deletions in this area
				 #that don't overlap completely but are quite the same and merge them
				
				$checkSNV_sth->execute($vcfline->{CHROM},$class,$pos,$end,$pos,$end,$pos,$end,$alt_nuc) || $logger->error($DBI::errstr);

				my $currDeletion = BEDRecord->new( $vcfline->{CHROM}, $pos, $end );
				my @overlappingDeletions;
				my $minStart = $pos;
				my $maxEnd   = $end;
				while ( my ( $tmp1, $tmp2, $tmp3, $tmp4 ) =	$checkSNV_sth->fetchrow_array() ) {
					my $dbDeletion = BEDRecord->new( $tmp1, $tmp2, $tmp3, $tmp4 );
					
					#print "testoverlap: ".$currDeletion->calcOverlap($dbDeletion)." : ".$currDeletion->toString()." : ".$dbDeletion->toString()."\n";
					if ( $currDeletion->calcOverlap($dbDeletion) >= $minOverlap && $dbDeletion->calcOverlap($currDeletion) >= $minOverlap )	{
						$minStart = $dbDeletion->startpos()  if $minStart > $dbDeletion->startpos();
						$maxEnd = $dbDeletion->endpos()      if $maxEnd < $dbDeletion->endpos();
						push( @overlappingDeletions, $dbDeletion );
						$insertVariant = 0;
					}
				}

				unless ($insertVariant) {
					#if one or more overlapping variants have been found: take the first variant in the array and let it be the one all other variants reference to
					#if there are more than one overlapping variants --> delete all these variants but add the isoform information to the main variant
					my $currfunc = $vcfline->{INFO}->{FUNC};

					#print "functest: ".$vcfline->{INFO}->{FUNC}."\n";
					foreach my $currBed (@overlappingDeletions) {
						my $currId = $currBed->name();

						#get isoform string
						$getIsoform_sth->execute($currId) || $logger->error($DBI::errstr);
						my ( $tmp1, $tmp2 ) = $getIsoform_sth->fetchrow_array();
						$isoform  .= $tmp1;
						$currfunc .= "," . $tmp2;

						if ( !defined($idsnv) )	{    #--> first entry --> main entry
							$idsnv = $currId;
						} else {
							#change all referencing snvsamples
							$updateIdsnv_sth->execute($idsnv,$currId) || $logger->error($DBI::errstr);

							#change all referencing additional annotations
							$updateAdditionalAnnotation_sth->execute($idsnv, $currId) || $logger->error($DBI::errstr);

							#some entries may not be updated because they violent a unique key: delete them
							$deleteViolentSamples_sth->execute($currId) || $logger->error($DBI::errstr);

							#delete old snv
							$deleteOldSnvGene_sth->execute($currId) || $logger->error($DBI::errstr);

							$deleteOldSnv2disease_sth->execute($currId) || $logger->error($DBI::errstr);

							if ($snpeffTable) {
								$deleteOldSnpEff_sth->execute($currId) || $logger->error($DBI::errstr);
							}
							$deleteOldSnv_sth->execute($currId) || $logger->error($DBI::errstr);
						}
					}

					#change transcript string
					$isoform = join( ":", uniq( split( ":", $isoform ) ) );
					my $dbLength = $maxEnd - $minStart;

					$changeTranscriptString_sth->execute($isoform,$maxEnd,$currfunc,$dbLength,$idsnv) || $logger->error($DBI::errstr);
				}
			}
			
			if ($insertVariant) {
				#my $dbLength = $length;
				$alt_nuc     =~ s/^<//;
				$alt_nuc     =~ s/>$//;
				my $dbLength = 1;						#is different than "$length" for Insertions --> length of insertion is 1
				if ( $reflength > $altlength ) {
					$dbLength = $reflength - $altlength;
				} elsif ( $reflength < $altlength ) {
					$dbLength = $altlength - $reflength;
				}

				#snv insert
				$insertSNV_sth->execute($vcfline->{CHROM},$pos,$end,$rs,$alt_nuc,$class,$vcfline->{INFO}->{FUNC},$isoform,$vcfline->{REF},$dbLength) || $rollbacks++;

				#TW: 13.03.2015: added check for REF allelele --> CNV entries for insertion/deletion of same exon were inserted twice

				$sth = 
				$checkForRefAllele_sth->execute($vcfline->{CHROM},$pos,$end,$class,$alt_nuc,$vcfline->{REF}) || $logger->error($DBI::errstr);
				($idsnv) = $checkForRefAllele_sth->fetchrow_array();

				#add link into snvgene
				foreach my $currgene ( uniq(@allgenes) ) {
					$currgene =~ s/\s+/_/g;
					$insertSnvGene_sth->execute($idsnv,$currgene) || $logger->error($DBI::errstr."; affected Gene Name: $currgene");
				}


				#add additional annotation into database table
				foreach my $anno (values %{$params->{settings}->{$settings}->{additionalannotations}}) {
					if ( $anno->{info} && $anno->{info}->{what} eq "bedname" ) {
						if ( $vcfline->{INFO}->{ $anno->{info}->{name} } ) {
							my @addAnnoIds = split( ",", $vcfline->{INFO}->{ $anno->{info}->{name} } );

							foreach my $currId (@addAnnoIds) {
								$insertAdditionalAnnotation_sth->execute($idsnv,$anno->{info}->{name},$currId) || $logger->error($DBI::errstr);
							}
						}
					}
				}

			}
		}

		foreach my $sample ( keys %sampleIds ) {

			#snvsample insert for each sample in VCF file
			my $alleles = 0;

			next if $vcfline->{gtypes}->{$sample}->{GT} eq "./.";   #skip non-calls

			next if $vcfline->{gtypes}->{$sample}->{GT} eq "0/0";    # skip hom ref calls

			if ( $vcfline->{gtypes}->{$sample}->{GT} =~ /(\d)\/(\d)/ ) {
				$alleles = $1 + $2;
			}

			#if ( $alleles != 0 ) {

			my $dp = 0;
			$dp = $vcfline->{gtypes}->{$sample}->{DP}  if $vcfline->{gtypes}->{$sample}->{DP};

			$dp = $strandPercent * 100 if $caller eq "exomedepth";
			$alleles = 2 if $caller eq "exomedepth" && $strandPercent <= 0.1;

			if ( $caller eq "gatk" || $caller eq "pindel" ) {
				$varPercent = 0;

				if ( $vcfline->{gtypes}->{$sample}->{AD} ) { #TW 11.09.2013: may be not there for some variants. don't know why
					my @columns = split( ",", $vcfline->{gtypes}->{$sample}->{AD} );    #get %variants for each sample in a GATK vcf file
					$columns[1] = 0 unless $columns[1];
					$varPercent = $columns[1] / ( $columns[0] + $columns[1] )  if ( $columns[0] + $columns[1] ) != 0;
					$varPercent *= 100;
					$dp = ( $columns[0] + $columns[1] );   # if $caller eq "pindel";
				}
			}

			my $gq = 99;
			$gq = $vcfline->{gtypes}->{$sample}->{GQ} if $vcfline->{gtypes}->{$sample}->{GQ};

			#$sql = qq{insert into $exomedb.$snvsampleTable (idsnvsample,idsnv,idsample,alleles,percentvar,percentfor,medianqual,refqual,snvqual,mapqual,coverage,filter,gtqual,caller) 
			#	                              values       (idsnvsample, '$idsnv', '$sampleIds{$sample}', '$alleles', '$varPercent', '$strandPercent', '$medianVarQual', '$refScore', '$snvScore', '$mq', '$dp','$filter',$gq,'$caller') };
			#$logger->debug($sql);
			
			#$sth = $dbh->prepare($sql) || $logger->error("Can't prepare statement: $DBI::errstr");
			#$sth->execute() || $duplicateSNPs++;
			
			print MYSQLIMPORTFILE "0\t$idsnv\t$sampleIds{$sample}\t$alleles\t$varPercent\t$strandPercent\t$medianVarQual\t$refScore\t$snvScore\t$mq\t$dp\t$filter\t$gq\t$caller\n"; #0 for auto_increment field			

		}
	}
}




=head1 NAME

insertVCF.pl

=head1 SYNOPSIS

insertVCF.pl -i variants.vcf 

=head1 DESCRIPTION

This script inserts annotated variants from a VCF file into the database defined in the config file. The VCF
file must have been annotated with the script annotateVCF.pl.

=head1 OPTIONS

 -i	<input.vcf>; REUIRED
 -se	name of the settings in the current.config.xml file that holds path to reference genome, 
 	to the annotation file and to possible additional annotation files; use default settings if nothing is given
 -d	delete entries for this sample(s) before inserting
 -dc	delete all entries on the specified chromosome for this sample(s) before inserting
 -c	<caller> program that produced the calls; currently available: pindel, samtools, gatk, exomedepth; default: samtools
 -l	<class> class of the variant that should be inserted, determined by length of variant by default;
	useful for e.g. large deletions detected by pindel
 -lf	log file; default: print to screen
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page

=head1 AUTHOR

Thomas Wieland

=cut

