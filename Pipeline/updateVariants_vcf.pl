#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use DBI;
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);
use Data::Dumper;

my $prog_path = dirname( abs_path($0) );
require $prog_path . "/Utilities.pm";

my $run = 1;

# database
my $dbh           = "";
my $sql           = "";
my $sth           = "";
my $logfile       = "pipeline.log";
my $loglevel      = "INFO";
my $duplicateSNPs = 0;
my $column        = "refseq";

my $infile           = "";
my $line             = "";
my $settings         = "";
my $help             = 0;
my $updateFunc       = 0;
my $replaceFunc      = 0;
my $updateAdditional = 0;

GetOptions(
	"i=s"  => \$infile,
	"c=s"  => \$column,
	"f"    => \$updateFunc,
	"a"    => \$updateAdditional,
	"rf"   => \$replaceFunc,
	"se=s" => \$settings,
	"lf=s" => \$logfile,
	"ll=s" => \$loglevel,
	"h"    => \$help,
);

my $params   = Utilities::getParams();
require $params->{programs}->{vcftools}->{pm};

my $exomedb  = $params->{settings}->{$settings}->{exomedb}->{database};
my $snvTable = $params->{settings}->{$settings}->{exomedb}->{snvtable};

Utilities::initLogger( $logfile, $loglevel );
my $logger = Utilities::getLogger();

if ( $help == 1 ) {
	print "
-i	<infile> 
-c	<column> which column in the snv table should be updated; default: $column
-f	concat new function to  \"function\" column
-a	also update additionalannotation entries
-rf replace \"function\" column with new function; only works if -f is NOT chosen
-se	settings
-lf	log file; default: pipeline.log
-ll	log level: ERROR,INFO,DEBUG; default: INFO
-h	this help\n";
	exit(0);
}

if ( $infile eq "" ) {
	$logger->error("no infile given");
	exit(1);
}
if ( $settings eq "" ) {
	$logger->error("no settings given");
	exit(1);
}

$dbh = Utilities::connectExomeDB($settings);
my $additionalannotationtable =
  $params->{settings}->{$settings}->{exomedb}->{additionalannotationtable};

#open VCF file
my $vcf = Vcf->new( file => $infile );
$vcf->parse_header();

while ( my $vcfline = $vcf->next_data_hash() ) {

	&updateSnv($vcfline);

}

#&checkInsert();
#$logger->info(
#	"$duplicateSNPs positions with more than one corresponding rsSNPs found");

my $elapsed = ( time - $^T ) / 60;
$logger->info("finished in $elapsed minutes");

#############
#subroutines#
#############
###########################
#insertSnv: fill databases#
###########################
sub updateSnv {

	my $vcfline = shift;
	
	my $class     = "snp";
	my $alt_nuc   = $vcfline->{ALT}[0];
	my $length    = 1;
	my $reflength = length $vcfline->{REF};
	my $altlength = length $alt_nuc;


    # Support (tolerate) VCF 4.3 - Ignore * allele right now
    my $whichAlt = 1;	# Which alternative allele pick in multi fields 
    my $hasStar = 0;

    if ( $alt_nuc eq "*" )
    {
        $alt_nuc = $vcfline->{ALT}[1];
        $whichAlt=2;
        $hasStar=1;
    }

	if ( $vcfline->{REF} eq "LD" || $vcfline->{REF} eq "LI" ) {			#old
		$length = $alt_nuc;
		if ( $vcfline->{REF} eq "LD" ) {
			$reflength = $length;
			$altlength = 0;
		}
		else {
			$reflength = 0;
			$altlength = $length;
		}
	} elsif ($vcfline->{INFO}->{SVTYPE} && !($vcfline->{REF} =~/[acgt]+/i && $alt_nuc =~/[acgt\*]+/i)) {		#TW 09.10.2015: pindel is a special case, because it has SVTYPE but usually also gives 
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
			$altlength = $length;
			$length    = 1;
		}
		elsif($vcfline->{INFO}->{SVTYPE} eq "INV"){
			$reflength = $length;
			$altlength = $length;
		}
	}
	else {

		#parse indel
		if ( $reflength != 1 || $altlength != 1 ) {
			$class = "indel";

			#get length for deletions
			if ( $reflength > $altlength ) {
				$length = $reflength - $altlength;
			}
			if ( $reflength > 255 || $altlength > 255 ) {
				if ( $reflength > $altlength ) {
					$vcfline->{REF} =
					  "LD";    #long deletion that doesn't fit into database
					$alt_nuc = $reflength - $altlength;
				}
				else {
					$vcfline->{REF} =
					  "LI";    #long insertion that doesn't fit into database
					$alt_nuc = $altlength - $reflength;
				}
			}
		}
	}
	
	

	#build isoform string
	my $isoform = "";

	#CDS
	if ( $vcfline->{INFO}->{CDSTRANS} ) {
		my @cdstrans  = split( ",", $vcfline->{INFO}->{CDSTRANS} );
		my @cdsgene   = split( ",", $vcfline->{INFO}->{CDSGENE} );
		my @cdsstrand = split( ",", $vcfline->{INFO}->{CDSSTRAND} );
		my @cdsnuc    = split( ",", $vcfline->{INFO}->{CDSNUC} );
		my @cdsframe  = split( ",", $vcfline->{INFO}->{CDSFRAME} );
		my @cdsaa;
		unless ( $vcfline->{INFO}->{INDEL} ) {
			@cdsaa = split( ",", $vcfline->{INFO}->{CDSAA} );
		}
		my @cdsfunc = split( ",", $vcfline->{INFO}->{CDSFUNC} );
		for ( my $i = 0 ; $i < @cdstrans ; $i++ ) {

#print "$cdstrans[$i],$cdsgene[$i],$cdsstrand[$i],$cdsnuc[$i],$cdsframe[$i],$cdsaa[$i],$cdsfunc[$i]: \n";
			$isoform .=
"$cdstrans[$i],$cdsgene[$i],$cdsstrand[$i],$cdsnuc[$i],$cdsframe[$i],$cdsfunc[$i]";
			if ( @cdsaa && $cdsaa[$i] ) {
				$isoform .= ",$cdsaa[$i]";
			}
			$isoform .= ": ";
		}
	}

	#noncoding
	if ( $vcfline->{INFO}->{NCTRANS} ) {
		my @nctrans  = split( ",", $vcfline->{INFO}->{NCTRANS} );
		my @ncgene   = split( ",", $vcfline->{INFO}->{NCGENE} );
		my @ncstrand = split( ",", $vcfline->{INFO}->{NCSTRAND} );
		for ( my $i = 0 ; $i < @nctrans ; $i++ ) {
			$nctrans[$i]  = "" unless $nctrans[$i];
			$ncgene[$i]   = "" unless $ncgene[$i];
			$ncstrand[$i] = "" unless $ncstrand[$i];
			$isoform .= "$nctrans[$i],$ncgene[$i],$ncstrand[$i],noncoding: ";
		}
	}

	#5utr
	if ( $vcfline->{INFO}->{UTR5TRANS} ) {
		my @utr5trans  = split( ",", $vcfline->{INFO}->{UTR5TRANS} );
		my @utr5gene   = split( ",", $vcfline->{INFO}->{UTR5GENE} );
		my @utr5strand = split( ",", $vcfline->{INFO}->{UTR5STRAND} );
		for ( my $i = 0 ; $i < @utr5trans ; $i++ ) {
			$utr5trans[$i]  = "" unless $utr5trans[$i];
			$utr5gene[$i]   = "" unless $utr5gene[$i];
			$utr5strand[$i] = "" unless $utr5strand[$i];
			$isoform .= "$utr5trans[$i],$utr5gene[$i],$utr5strand[$i],5utr: ";
		}
	}

	#3utr
	if ( $vcfline->{INFO}->{UTR3TRANS} ) {
		my @utr3trans  = split( ",", $vcfline->{INFO}->{UTR3TRANS} );
		my @utr3gene   = split( ",", $vcfline->{INFO}->{UTR3GENE} );
		my @utr3strand = split( ",", $vcfline->{INFO}->{UTR3STRAND} );
		for ( my $i = 0 ; $i < @utr3trans ; $i++ ) {
			$utr3trans[$i]  = "" unless $utr3trans[$i];
			$utr3gene[$i]   = "" unless $utr3gene[$i];
			$utr3strand[$i] = "" unless $utr3strand[$i];
			$isoform .= "$utr3trans[$i],$utr3gene[$i],$utr3strand[$i],3utr: ";
		}
	}

	#direct splice site hit
	if ( $vcfline->{INFO}->{DSTRANS} ) {
		my @dstrans  = split( ",", $vcfline->{INFO}->{DSTRANS} );
		my @dsgene   = split( ",", $vcfline->{INFO}->{DSGENE} );
		my @dsstrand = split( ",", $vcfline->{INFO}->{DSSTRAND} );
		my @dsmut;
		if ( $vcfline->{INFO}->{DSMUT} ) {
			@dsmut = split( ",", $vcfline->{INFO}->{DSMUT} );
		}
		for ( my $i = 0 ; $i < @dstrans ; $i++ ) {
			$dstrans[$i]  = "" unless $dstrans[$i];
			$dsgene[$i]   = "" unless $dsgene[$i];
			$dsstrand[$i] = "" unless $dsstrand[$i];
			$isoform .= "$dstrans[$i],$dsgene[$i],$dsstrand[$i]";
			if ( @dsmut && $dsmut[$i] ) {
				$isoform .= ",$dsmut[$i]";
			}
			$isoform .= ",direct splicesite: ";
		}
	}

	#near splice site hit
	if ( $vcfline->{INFO}->{NSTRANS} ) {
		my @nstrans  = split( ",", $vcfline->{INFO}->{NSTRANS} );
		my @nsgene   = split( ",", $vcfline->{INFO}->{NSGENE} );
		my @nsstrand = split( ",", $vcfline->{INFO}->{NSSTRAND} );
		for ( my $i = 0 ; $i < @nstrans ; $i++ ) {
			$nstrans[$i]  = "" unless $nstrans[$i];
			$nsgene[$i]   = "" unless $nsgene[$i];
			$nsstrand[$i] = "" unless $nsstrand[$i];
			$isoform .=
			  "$nstrans[$i],$nsgene[$i],$nsstrand[$i],near splicesite: ";
		}
	}

	#intronic
	if ( $vcfline->{INFO}->{INTTRANS} ) {
		my @inttrans  = split( ",", $vcfline->{INFO}->{INTTRANS} );
		my @intgene   = split( ",", $vcfline->{INFO}->{INTGENE} );
		my @intstrand = split( ",", $vcfline->{INFO}->{INTSTRAND} );
		for ( my $i = 0 ; $i < @inttrans ; $i++ ) {
			$inttrans[$i]  = "" unless $inttrans[$i];
			$intgene[$i]   = "" unless $intgene[$i];
			$intstrand[$i] = "" unless $intstrand[$i];
			$isoform .= "$inttrans[$i],$intgene[$i],$intstrand[$i],intronic: ";
		}
	}

	#upstream & downstream variants
	if ( $vcfline->{INFO}->{UPSTRANS} ) {
		$isoform .=
"$vcfline->{INFO}->{UPSTRANS},$vcfline->{INFO}->{UPSGENE},$vcfline->{INFO}->{UPSSTRAND},distance-$vcfline->{INFO}->{UPSDIST}, upstream gene: ";
	}
	if ( $vcfline->{INFO}->{DOSTRANS} ) {
		$isoform .=
"$vcfline->{INFO}->{DOSTRANS},$vcfline->{INFO}->{DOSGENE},$vcfline->{INFO}->{DOSSTRAND},distance-$vcfline->{INFO}->{DOSDIST}, downstream gene: ";
	}

	my $funcstring = "";
	if ( $vcfline->{INFO}->{FUNC} && $updateFunc ) {
		$funcstring = "func=concat_ws(',',func,'$vcfline->{INFO}->{FUNC}') ";
	}
	elsif ( $vcfline->{INFO}->{FUNC} && $replaceFunc ) {
		$funcstring = "func='$vcfline->{INFO}->{FUNC}' ";
	}

	if ($run) {
		if ($column eq "mirna") {
			if ($funcstring ne "") {
				$sql = qq{update $exomedb.$snvTable set $funcstring where chrom = '$vcfline->{CHROM}' and start = '$vcfline->{POS}' and allele = '$alt_nuc' and refallele='$vcfline->{REF}'};
				$sth = $dbh->prepare($sql) || $logger->error("Can't prepare statement: $DBI::errstr");
				$sth->execute() || $logger->error($DBI::errstr);
			}
		} else {
			my $comma="";
			if ( $funcstring ne "" ) { $comma=",";}
			$sql = qq{update $exomedb.$snvTable set $column='$isoform' $comma  $funcstring where chrom = '$vcfline->{CHROM}' and start = '$vcfline->{POS}' and allele = '$alt_nuc' and refallele='$vcfline->{REF}'};
			$sth = $dbh->prepare($sql) || $logger->error("Can't prepare statement: $DBI::errstr");
			$sth->execute() || $logger->error($DBI::errstr);
			#check if entry in snv already exists
			#print $sql."\n";
			#qq{select idsnv from $exomedb.$snvTable where chrom = '$vcfline->{CHROM}' and start = '$pos' and allele = '$alt_nuc' and refallele='$vcfline->{REF}'};
		}

		if ($updateAdditional) {

			#delete old entries
			$sql = qq{delete from $exomedb.$additionalannotationtable where idsnv=(select idsnv from $exomedb.$snvTable where chrom = '$vcfline->{CHROM}' and start = '$vcfline->{POS}' and allele = '$alt_nuc' and refallele='$vcfline->{REF}'); };
			$logger->debug($sql);
			$sth = $dbh->prepare($sql) || $logger->error("Can't prepare statement: $DBI::errstr");
			$sth->execute() || $logger->error($DBI::errstr);

			foreach my $anno (values %{ $params->{settings}->{$settings}->{additionalannotations} } ) {
				if ( $anno->{info} && $anno->{info}->{what} eq "bedname" ) {
					if ( $vcfline->{INFO}->{ $anno->{info}->{name} } ) {
						my @addAnnoIds = split( ",", $vcfline->{INFO}->{ $anno->{info}->{name} } );

						foreach my $currId (@addAnnoIds) {
							$sql = qq{insert into $exomedb.$additionalannotationtable (idsnv,annotationname,idannotation) values ((select idsnv from $exomedb.$snvTable where chrom = '$vcfline->{CHROM}' and start = '$vcfline->{POS}' and allele = '$alt_nuc' and refallele='$vcfline->{REF}'),'$anno->{info}->{name}',$currId); };
							$logger->debug($sql);
							$sth = $dbh->prepare($sql) || $logger->error("Can't prepare statement: $DBI::errstr");
							$sth->execute() || $logger->error($DBI::errstr);
						}
					}
				}
			}
		}
	}
}
