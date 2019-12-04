#!/usr/bin/perl 


use strict;
use Getopt::Long;
use DBI;
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);
use Data::Dumper;
use IO::File;
use Pod::Usage;


my $prog_path = dirname( abs_path($0) );
require $prog_path . "/Utilities.pm";
require $prog_path . "/BEDRecord.pm";
require $prog_path . "/SequenceDictionary.pm";



my $dbh = "";  

my $settings = "default";

my $infile = "";

my $outfile      = "";
my $refSeqAli    = 0;
my $lincRNA      = 0;
my $miRNA		 = 0;
my $spliceWindow = 5;
my $help         = 0;
my $man			 = 0;
my $logfile      = "SCREEN";
my $loglevel     = "INFO";
my $sequenceDict;
my $table    = "";
my $cdstable = "";
my $canonicaltable = "";
my $canonicaljoin = "";
my $addFunc  = "";
my $annoMode     = "auto";
my $canonical = 0;


GetOptions(
	"i=s"  => \$infile,
	"o=s"  => \$outfile,
	"r"    => \$refSeqAli,
	"l"    => \$lincRNA,
	"mir"   => \$miRNA,
	"f=s"  => \$addFunc,
	"w=s"  => \$spliceWindow,
	"t=s"  => \$table,
	"c=s"  => \$cdstable,
	"se=s" => \$settings,
	"m=s"  => \$annoMode,
	"canonical" => \$canonical,
	"lf=s" => \$logfile,
	"ll=s" => \$loglevel,
	"man"  => \$man,
	"h"    => \$help
);

pod2usage( { -exitval => 0, -verbose => 1 } ) if $help;
pod2usage( { -exitval => 0, -verbose => 2 } ) if $man;
pod2usage( { -exitval => 1, -verbose => 1 } ) if $infile eq "";


Utilities::initLogger( $logfile, $loglevel );
my $logger = Utilities::getLogger();




if ( $outfile eq "" ) {
	$outfile = $infile;
	$outfile =~ s/.vcf$/.plus.vcf/;
}


$dbh = Utilities::connectVariationDB($settings);
my $params = Utilities::getParams();


require $params->{programs}->{vcftools}->{pm};    #VCFTools library to parse VCF file

my $database = $params->{settings}->{$settings}->{variationdb}->{database};

# Canonical is supported ONLY with knownGene and knownCanonical UCSC tables
# Check it and exit with error if canonical option activated improperly
if ( $table eq "" ) {
	if ($refSeqAli) {
		$table = $params->{settings}->{$settings}->{variationdb}->{refseqgene};
		if ($canonical == 1)
		{
			$logger->error("Canonical annotation not supported in this annotation mode");
			exit -1;
		}
	}
	elsif ($lincRNA) {
		$table = $params->{settings}->{$settings}->{variationdb}->{lincrna};
		if ($canonical == 1)
		{
			$logger->error("Canonical annotation not supported in this annotation mode");
			exit -1;
		}
	}
	elsif ($miRNA) {
		$table = $params->{settings}->{$settings}->{variationdb}->{mirna};
		if ($canonical == 1)
		{
			$logger->error("Canonical annotation not supported in this annotation mode");
			exit -1;
		}
	}
	else {
		$table = $params->{settings}->{$settings}->{variationdb}->{genetable};		
		if ($canonical == 1)
		{
			$canonicaltable = $params->{settings}->{$settings}->{variationdb}->{canonicaltable};
			$canonicaljoin= " INNER JOIN $database.$canonicaltable CAN ON CAN.transcript=T.name ";
			$logger->debug("Annotating just with canonical transcripts");
		}
	}
}
else
{
	if ($canonical == 1)
	{
		$logger->error("Canonical annotation not supported in this annotation mode");
		exit -1;
	}
}

if ( $cdstable eq "" ) {
	if ($refSeqAli) {
		$cdstable =
		  $params->{settings}->{$settings}->{variationdb}->{refseqcds};
	}
	elsif ($lincRNA) {
		$cdstable =
		  $params->{settings}->{$settings}->{variationdb}->{lincrnacds};
	}
	elsif ($miRNA) {
		$cdstable = 
		  $params->{settings}->{$settings}->{variationdb}->{mirnacds};
	}
	else {
		$cdstable = $params->{settings}->{$settings}->{variationdb}->{cds};
	}
}

#global hitcounters and variables
my $codingHits      = 0;
my $intronHits      = 0;
my $synHits         = 0;
my $nonSynHits      = 0;
my $spliceHits      = 0;
my $knownSNPs       = 0;
my $UTR5hits        = 0;
my $UTR3hits        = 0;
my $UTRhits         = 0;
my %snpHash         = ();
my %indelHash       = ();
my %tgenomeHash     = ();
my %hapMapExomeHash = ();
my %snvdbExomeHash  = ();
my %hitCountHash    = ();    #record number of cds hits for each gene, filled in get triplet

#my %geneSymbol  = ();
my $currSNPName = ".";

#&initGeneSymbol();
&annotateVariant;


#$logger->info(
#"found $codingHits variants in coding regions, $synHits synonymous substitutions and $nonSynHits non-synonymous substitutions"
#);
#$logger->info("found $spliceHits hits in splice sites");
#$logger->info("$UTRhits variants hit (one or more) UTR regions");
#$logger->info("$intronHits variants hit one (or more) introns");
#$logger->info(
#	"total of $UTR5hits 5'UTRs and $UTR3hits 3'UTRs are affected by variants");

my $elapsed = ( time - $^T ) / 60;
$logger->info("finished in $elapsed minutes");


#############
#subroutines#
#############

sub annotateVariant {

	#vcsf variables
	my ( $chr, $pos, $alt_nuc );
	my $vcfline;

	my (
		$name,         $chrom,      $strand,     $txStart,
		$txEnd,        $cdsStart,   $cdsEnd,     $exonCount,
		$exonStarts,   $exonEnds,   $blockSizes, $spliceCount,
		$spliceStarts, $spliceEnds, $geneSymbol
	);
	my $chr_old        = "";
	my $query          = "";
	my $out            = "";
	my @exonStarts     = ();
	my @exonEnds       = ();
	my @blockSizes     = ();
	my @spliceStarts   = ();
	my @spliceEnds     = ();
	my $tmpExonStart   = 0;
	my $i              = 0;
	my $spos           = 0;
	my $aa             = "";
	my $codingprint    = "";
	my $spliceprint    = "";
	my $spliceMutation = "";
	my $knownSNPprint  = "";
	my $UTRprint       = "";
	my $searchkey      = "";
	my $outprint       = "";
	my $occurence      = 0;
	my @line           = ();
	my @knownAnnot     = ();

	open( OUT, ">$outfile" ) || exit $logger->error("Cannot open $outfile");

	my $vcf = Vcf->new( file => $infile );
	$vcf->parse_header();							#add header entries
	
	$vcf->add_header_line(
		{
			key         => 'INFO',
			ID          => 'FUNC',
			Number      => '.',
			Type        => 'String',
			Description => 'Putative function of this variant, according to all isoforms.'
		}
	);	
	
	$vcf->add_header_line(
		{
			key         => 'INFO',
			ID          => 'UPSTRANS',
			Number      => '1',
			Type        => 'String',
			Description => 'Next upstream transcript if variant is intergenic'
		}
	);
	$vcf->add_header_line(
		{
			key         => 'INFO',
			ID          => 'UPSGENE',
			Number      => '1',
			Type        => 'String',
			Description => 'Next upstream gene if variant is intergenic'
		}
	);
	$vcf->add_header_line(
		{
			key         => 'INFO',
			ID          => 'UPSSTRAND',
			Number      => '1',
			Type        => 'String',
			Description => 'Strand of next upstream gene if variant is intergenic'
		}
	);
	$vcf->add_header_line(
		{
			key         => 'INFO',
			ID          => 'UPSDIST',
			Number      => '1',
			Type        => 'Integer',
			Description => 'Distance to next upstream gene if variant is intergenic'
		}
	);
	$vcf->add_header_line(
		{
			key         => 'INFO',
			ID          => 'DOSTRANS',
			Number      => '1',
			Type        => 'String',
			Description => 'Next downstream transcript if variant is intergenic'
		}
	);
	$vcf->add_header_line(
		{
			key         => 'INFO',
			ID          => 'DOSGENE',
			Number      => '1',
			Type        => 'String',
			Description => 'Next downstream gene if variant is intergenic'
		}
	);
	$vcf->add_header_line(
		{
			key         => 'INFO',
			ID          => 'DOSSTRAND',
			Number      => '1',
			Type        => 'String',
			Description => 'Strand of next downstream gene if variant is intergenic'
		}
	);
	$vcf->add_header_line(
		{
			key         => 'INFO',
			ID          => 'DOSDIST',
			Number      => '1',
			Type        => 'Integer',
			Description => 'Distance to next downstream gene if variant is intergenic'
		}
	);
	
	if($annoMode eq "cnv" || $annoMode eq "auto"){
		
		$vcf->add_header_line(
			{
				key         => 'INFO',
				ID          => 'CNVTRANS',
				Number      => '.',
				Type        => 'String',
				Description => 'Transcript(s) this CNV overlaps with'
			}
		);
		$vcf->add_header_line(
			{
				key         => 'INFO',
				ID          => 'CNVGENE',
				Number      => '.',
				Type        => 'String',
				Description => 'Gene(s) this CNV overlaps with'
			}
		);
		$vcf->add_header_line(
			{
				key         => 'INFO',
				ID          => 'CNVSTRAND',
				Number      => '.',
				Type        => 'String',
				Description => 'Strand(s) of the Transcript(s)/Gene(s) this CNV overlaps with'
			}
		);
	}
	if($annoMode eq "short" || $annoMode eq "auto"){

		#add lines to header
		$vcf->add_header_line(
			{
				key         => 'INFO',
				ID          => 'TG',
				Number      => 0,
				Type        => 'Flag',
				Description => 'The variant was also found in the 1000G project'
			}
		);
		$vcf->add_header_line(
			{
				key         => 'INFO',
				ID          => 'UTR5TRANS',
				Number      => '.',
				Type        => 'String',
				Description => 'Transcript(s) in whose 5\'UTR this variant lies'
			}
		);
		$vcf->add_header_line(
			{
				key         => 'INFO',
				ID          => 'UTR5GENE',
				Number      => '.',
				Type        => 'String',
				Description => 'Gene(s) in whose 5\'UTR this variant lies'
			}
		);
		$vcf->add_header_line(
			{
				key         => 'INFO',
				ID          => 'UTR5STRAND',
				Number      => '.',
				Type        => 'String',
				Description => 'Strand(s) of the Transcript(s)/Gene(s) in whose 5\'UTR this variant lies'
			}
		);
		$vcf->add_header_line(
			{
				key         => 'INFO',
				ID          => 'UTR3TRANS',
				Number      => '.',
				Type        => 'String',
				Description => 'Transcript(s) in whose 3\'UTR this variant lies'
			}
		);
		$vcf->add_header_line(
			{
				key         => 'INFO',
				ID          => 'UTR3GENE',
				Number      => '.',
				Type        => 'String',
				Description => 'Gene(s) in whose 3\'UTR this variant lies'
			}
		);
		$vcf->add_header_line(
			{
				key         => 'INFO',
				ID          => 'UTR3STRAND',
				Number      => '.',
				Type        => 'String',
				Description => 'Strand(s) of the Transcript(s)/Gene(s) in whose 3\'UTR this variant lies'
			}
		);
		$vcf->add_header_line(
			{
				key         => 'INFO',
				ID          => 'NCTRANS',
				Number      => '.',
				Type        => 'String',
				Description => 'Transcript(s) in whose noncoding area this variant lies'
			}
		);
		$vcf->add_header_line(
			{
				key         => 'INFO',
				ID          => 'NCGENE',
				Number      => '.',
				Type        => 'String',
				Description => 'Gene(s) in whose noncoding area this variant lies'
			}
		);
		$vcf->add_header_line(
			{
				key         => 'INFO',
				ID          => 'NCSTRAND',
				Number      => '.',
				Type        => 'String',
				Description => 'Strand(s) of the Transcript(s)/Gene(s) in whose noncoding area this variant lies'
			}
		);
		$vcf->add_header_line(
			{
				key         => 'INFO',
				ID          => 'CDSTRANS',
				Number      => '.',
				Type        => 'String',
				Description => 'Transcript(s) in whose coding sequence this variant lies'
			}
		);
		$vcf->add_header_line(
			{
				key         => 'INFO',
				ID          => 'CDSGENE',
				Number      => '.',
				Type        => 'String',
				Description => 'Gene(s) in whose coding sequence this variant lies'
			}
		);
		$vcf->add_header_line(
			{
				key         => 'INFO',
				ID          => 'CDSSTRAND',
				Number      => '.',
				Type        => 'String',
				Description => 'Strand(s) of the Transcript(s)/Gene(s) in whose coding sequence this variant lies'
			}
		);
		$vcf->add_header_line(
			{
				key         => 'INFO',
				ID          => 'CDSNUC',
				Number      => '.',
				Type        => 'String',
				Description => 'Nucleotide change and its position in the current transcript'
			}
		);
		$vcf->add_header_line(
			{
				key         => 'INFO',
				ID          => 'CDSFRAME',
				Number      => '.',
				Type        => 'Integer',
				Description => 'Position of the variant within its frame'
			}
		);
		$vcf->add_header_line(
			{
				key         => 'INFO',
				ID          => 'CDSAA',
				Number      => '.',
				Type        => 'String',
				Description => 'Amino acid change triggered by the variant'
			}
		);
		$vcf->add_header_line(
			{
				key         => 'INFO',
				ID          => 'CDSFUNC',
				Number      => '.',
				Type        => 'String',
				Description => 'Putative function of the variant'
			}
		);
		$vcf->add_header_line(
			{
				key         => 'INFO',
				ID          => 'DSTRANS',
				Number      => '.',
				Type        => 'String',
				Description => 'Transcript(s) in direct splice site this variant lies'
			}
		);
		$vcf->add_header_line(
			{
				key         => 'INFO',
				ID          => 'DSGENE',
				Number      => '.',
				Type        => 'String',
				Description => 'Gene(s) in whose direct splice site this variant lies'
			}
		);
		$vcf->add_header_line(
			{
				key         => 'INFO',
				ID          => 'DSSTRAND',
				Number      => '.',
				Type        => 'String',
				Description => 'Strand(s) of the Transcript(s)/Gene(s) in whose direct splice site this variant lies'
			}
		);
		$vcf->add_header_line(
			{
				key         => 'INFO',
				ID          => 'DSMUT',
				Number      => '.',
				Type        => 'String',
				Description => 'Mutation at the current direct splice site.'
			}
		);
		$vcf->add_header_line(
			{
				key         => 'INFO',
				ID          => 'NSTRANS',
				Number      => '.',
				Type        => 'String',
				Description => 'Transcript(s) near whose splice site this variant lies'
			}
		);
		$vcf->add_header_line(
			{
				key         => 'INFO',
				ID          => 'NSGENE',
				Number      => '.',
				Type        => 'String',
				Description => 'Gene(s) near whose splice site this variant lies'
			}
		);
		$vcf->add_header_line(
			{
				key         => 'INFO',
				ID          => 'NSSTRAND',
				Number      => '.',
				Type        => 'String',
				Description => 'Strand(s) of the Transcript(s)/Gene(s) near whose splice site this variant lies'
			}
		);
		$vcf->add_header_line(
			{
				key         => 'INFO',
				ID          => 'INTTRANS',
				Number      => '.',
				Type        => 'String',
				Description => 'Transcript(s) in whose intron this variant lies'
			}
		);
		$vcf->add_header_line(
			{
				key         => 'INFO',
				ID          => 'INTGENE',
				Number      => '.',
				Type        => 'String',
				Description => 'Gene(s) in whose intron this variant lies'
			}
		);
		$vcf->add_header_line(
			{
				key         => 'INFO',
				ID          => 'INTSTRAND',
				Number      => '.',
				Type        => 'String',
				Description => 'Strand(s) of the Transcript(s)/Gene(s) in whose intron this variant lies'
			}
		);
		$vcf->add_header_line(
			{
				key         => 'INFO',
				ID          => 'DBPOS',
				Number      => '.',
				Type        => 'Integer',
				Description => 'Position of an indel in database (last possible position), for comparison with e.g. HGMD.'
			}
		);
		$vcf->add_header_line(
			{
				key         => 'INFO',
				ID          => 'DBINDEL',
				Number      => '.',
				Type        => 'String',
				Description => 'Nominclature of an indel for database, for comparison with e.g. HGMD.'
			}
		);


		#add header for additional annotations
		#and open file pointers
		my $readSD = 0;
		foreach my $anno (values %{ $params->{settings}->{$settings}->{additionalannotations} }) {
			$readSD = 1;
			$anno->{filepointer} = new IO::File $anno->{file}  || exit $logger->error("Can't open $anno->{file}!");
			$anno->{readbuffer} = undef;
			if ( $anno->{info} ) {
				my $type   = 'Flag';
				my $number = 0;
				if ( $anno->{info}->{what} eq "bedname" ) {
					$type   = 'String';
					$number = ".";
				}
				$vcf->add_header_line(
					{
						key         => 'INFO',
						ID          => $anno->{info}->{name},
						Number      => $number,
						Type        => $type,
						Description => $anno->{info}->{description}
					}
				);
			}
			elsif ( $anno->{filter} ) {
				$vcf->add_header_line(
					{
						key         => 'FILTER',
						ID          => $anno->{info}->{name},
						Description => $anno->{info}->{description}
					}
				);
			}
			else {
				next;
			}
		}
	
		#Read sequence dictionary if necessary and bamfile given --> otherwise generic sort order of chromosomes will be used
		if ($readSD) {
			$sequenceDict = SequenceDictionary->new($params->{settings}->{$settings}->{reference}.".fai");
		}
	}

	print OUT $vcf->format_header();

	while ( $vcfline = $vcf->next_data_hash() ) {
		
		if($annoMode eq "cnv" || ($annoMode eq "auto" && $vcfline->{INFO}->{END} && !($vcfline->{REF} =~/[acgt]+/i && $vcfline->{ALT}[0] =~/[acgt]+/i) ) ){			#CNV annotation mode
			$vcfline->{INFO}->{FUNC} = "unknown";
		
			
	
			my $hitsGene = 0;
	 
			my $query = "Select distinct T.name,T.strand,T.geneSymbol FROM $database.$table T $canonicaljoin
					 WHERE T.chrom='$chr' AND (($vcfline->{INFO}->{END}>=T.txStart AND $vcfline->{INFO}->{END}<=T.txEnd) OR (T.txEnd>=$vcfline->{POS} AND T.txEnd<=$vcfline->{INFO}->{END}))";
			$logger->debug( "query: " . $query );
			my $out = $dbh->prepare($query) || exit $logger->error($DBI::errstr);
			$out->execute || exit $logger->error($DBI::errstr);
	
			
			#
			# go through every intersecting sequence
			#
			while (my ($name,$strand,$geneSymbol) = $out->fetchrow_array)
			{

				$hitsGene = 1;
		

				if ( $vcfline->{INFO}->{CNVTRANS} ) {
					$vcfline->{INFO}->{CNVTRANS}  .= "," . $name;
					$vcfline->{INFO}->{CNVGENE}   .= "," . $geneSymbol;
					$vcfline->{INFO}->{CNVSTRAND} .= "," . $strand;
				}
				else {
					$vcfline->{INFO}->{CNVTRANS}  = $name;
					$vcfline->{INFO}->{CNVGENE}   = $geneSymbol;
					$vcfline->{INFO}->{CNVSTRAND} = $strand;
				}


			}
			
			unless($hitsGene){									#---> if no transcript is hit --> get upstream & downstream genes
				$vcfline->{INFO}->{FUNC} = "intergenic";
				&getUpstreamDownstream($vcfline->{CHROM},$vcfline->{POS},$vcfline->{INFO}->{END},$vcfline);
			}
		
			
		}elsif ($annoMode eq "short" || ($annoMode eq "auto" && $vcfline->{REF} =~/[acgt]+/i && $vcfline->{ALT}[0] =~/[acgt]+/i ) ) {
			
			$codingprint   = "";
			$spliceprint   = "";
			$UTRprint      = "";
			$knownSNPprint = "";
			$outprint      = "";
			my $func = $vcfline->{INFO}->{FUNC} . ",";
			$occurence = 0;
			my $is1000G = 0;
			chomp;
	
			#get info from VCF line
			$chr     = $vcfline->{CHROM};
			$pos     = $vcfline->{POS};
			$alt_nuc = $vcfline->{ALT}[0];
	
			#splice read handling
			if ( $chr =~ /^uc/ ) {
				$chr =~ /(\d+)$/;
				$pos = $1 + $pos;
				if ( $chr =~ /(chr\d+)/ ) {
					$chr = $1;
				}
				elsif ( $chr =~ /chrX/ ) {
					$chr = "chrX";
				}
				elsif ( $chr =~ /chrY/ ) {
					$chr = "chrY";
				}
			}
	
			$spos = $pos;
			$pos--;    # database starts with 0
	
			my $topIndel = "";
			my $length   = 1;
	
			#IF INDEL
			if ( length $vcfline->{REF} != 1 || length $alt_nuc != 1 ) {
	
	#parse indels from VCF to "old" format
	#this is quite difficult, because the position of an indel is nonambiguous
	#for example:
	#   Pos: 123456789
	#		 TACCCCCGT
	#          |   |
	#
	#If a single "C" is deleted in this example it is not clear which "C" (Pos 3-7) is really deleted.
	#Actually it doesn't matter which "C" is deleted, but there must be a convention which position to save,
	#especially to compare them to known indels. At the moment we look for our indels in dbSNP and HGMD, which
	#use different annotations:
	# -) dbSNP uses the the first possible position (i.e. Pos 3 in the example)
	# -) HGMD uses the last possible position, but dependent on the strand the corresponding gene is on
	#    (i.e. if the gene at this position is on the + strand Pos 7 is used, if it is on the - strand Pos 3 is used)
	#Currently we save the position according to HGMD in our database (VCF field DBPOS), because the comparison with
	#HGMD is done on database level. But since we can only save a single position per variant, this might lead to problems
	#at positions with genes annotated on BOTH strands. For now this means that the position according to the last gene
	#that comes from the database query. Within this script both positions (and deletions) are calculated and used as
	#$posFirst, $indelFirst and $posLast, $indelLast, respectively (for the example above: $posFirst=3,$indelFirst=-G; $posLast=7,$indelLast=-C).
	#
	# UPDATE:
	# It might be smarter to just take variant databases and convert their format into VCF format. --> see Utilities::indel2VCF
	#
	
				if ( length $vcfline->{REF} < length $alt_nuc ) {    #insertion
					$topIndel = "+";
					my $offset =
					  length $vcfline
					  ->{REF}; #the number of bases not to be taken from the alternative
	
					$topIndel .=
					  uc substr( $alt_nuc, 1, ( length $alt_nuc ) - $offset );
	
				}
				else {         #deletion
					$topIndel = "-";
					my $offset =
					  length $alt_nuc
					  ;    #the number of bases not to be taken from the reference
					$length = ( length $topIndel ) - 1;
	
					$topIndel .= uc substr( $vcfline->{REF}, 1,
						( length $vcfline->{REF} ) - $offset );
				}
				
				$pos++;		#TW, 14.04.2012: instead of $pos += $offset; --> don't shift it to the rightmost position, but rather to the first possible position --> should fit to the position that is shown in the webinterface as genomic position
				
				$vcfline->{INFO}->{DBPOS}   = $pos;
				$vcfline->{INFO}->{DBINDEL} = $topIndel;

			}
	
			$chr_old = $chr;
	
			# get all intersecting coding sequences
			$query = "Select T.name,T.chrom,T.strand,T.txStart,T.txEnd,T.cdsStart,T.cdsEnd,
		   T.exonCount,T.exonStarts,T.exonEnds,T.geneSymbol FROM $database.$table T $canonicaljoin
		WHERE T.chrom='$chr' 
		AND $pos >= T.txStart  
		AND $pos <= T.txEnd 
		";
			$logger->debug( "query: " . $query );
			$out = $dbh->prepare($query) || exit $logger->error($DBI::errstr);
			$out->execute || exit $logger->error($DBI::errstr);
	
			my (
				$isoforms,   $stopCount, $gainCount, $nonSyn, $syn,
				$frameshift, $indel,     $directSS,  $nearSS, $utr5,
				$utr3,       $noncoding, $introns,   $missense
			  )
			  = 0;    #counters for overall function annotation
	
			my $intronic;
	
			while (
				(
					$name,       $chrom,    $strand, $txStart,
					$txEnd,      $cdsStart, $cdsEnd, $exonCount,
					$exonStarts, $exonEnds, $geneSymbol
				)
				= $out->fetchrow_array
			  )
			{
				$isoforms++;
				$intronic = 1;
	
				(@exonStarts) = split( /\,/, $exonStarts );
				(@exonEnds)   = split( /\,/, $exonEnds );
	
				$spliceCount = $exonCount;
				(@spliceStarts) = split( /\,/, $exonStarts );
				(@spliceEnds)   = split( /\,/, $exonEnds );
	
				#(2) check hits in 5' and 3' UTRs or noncoding
				if ( $cdsStart == $cdsEnd ) {
					my ( $ncTrans, $ncGene ) = checkNc(
						$chrom,     $pos,         $name,      $strand,
						$txStart,   $txEnd,       $cdsStart,  $cdsEnd,
						$exonCount, \@exonStarts, \@exonEnds, $geneSymbol
					);
	
					if ( $ncTrans && $ncTrans ne "" ) {
						$intronic = 0;
						if ( $vcfline->{INFO}->{NCTRANS} ) {
							$vcfline->{INFO}->{NCTRANS}  .= "," . $ncTrans;
							$vcfline->{INFO}->{NCGENE}   .= "," . $ncGene;
							$vcfline->{INFO}->{NCSTRAND} .= "," . $strand;
						}
						else {
							$vcfline->{INFO}->{NCTRANS}  = $ncTrans;
							$vcfline->{INFO}->{NCGENE}   = $ncGene;
							$vcfline->{INFO}->{NCSTRAND} = $strand;
						}
						$noncoding++;
					}
				}
				else {
					my (
						$utr5trans, $utr5gene, $utr3trans, $utr3gene
	
					  ) = checkUTR(
						$chrom,     $pos,         $name,      $strand,
						$txStart,   $txEnd,       $cdsStart,  $cdsEnd,
						$exonCount, \@exonStarts, \@exonEnds, $geneSymbol
					  );
					if ($utr5trans) {
						$intronic = 0;
						if ( $vcfline->{INFO}->{UTR5TRANS} ) {
							$vcfline->{INFO}->{UTR5TRANS}  .= "," . $utr5trans;
							$vcfline->{INFO}->{UTR5GENE}   .= "," . $utr5gene;
							$vcfline->{INFO}->{UTR5STRAND} .= "," . $strand;
						}
						else {
							$vcfline->{INFO}->{UTR5TRANS}  = $utr5trans;
							$vcfline->{INFO}->{UTR5GENE}   = $utr5gene;
							$vcfline->{INFO}->{UTR5STRAND} = $strand;
						}
						$utr5++;
					}
					if ($utr3trans) {
						$intronic = 0;
						if ( $vcfline->{INFO}->{UTR3TRANS} ) {
							$vcfline->{INFO}->{UTR3TRANS}  .= "," . $utr3trans;
							$vcfline->{INFO}->{UTR3GENE}   .= "," . $utr3gene;
							$vcfline->{INFO}->{UTR3STRAND} .= "," . $strand;
						}
						else {
							$vcfline->{INFO}->{UTR3TRANS}  = $utr3trans;
							$vcfline->{INFO}->{UTR3GENE}   = $utr3gene;
							$vcfline->{INFO}->{UTR3STRAND} = $strand;
						}
						$utr3++;
					}
				}
	
				#(3) check for hits in coding sequence
				if ( $topIndel eq "" ) {
	
					my ( $csTrans, $csGene, $csNuc, $csFrame, $csFunc, $csAA ) =
					  checkCds(
						$chrom,       $pos,       $vcfline->{REF},
						$alt_nuc,     $name,      $strand,
						$cdsStart,    $cdsEnd,    $exonCount,
						\@exonStarts, \@exonEnds, $geneSymbol
					  );
	
					if ($csTrans) {
						$codingHits++;
						$intronic = 0;
						if ( $vcfline->{INFO}->{CDSTRANS} ) {
							$vcfline->{INFO}->{CDSTRANS}  .= "," . $csTrans;
							$vcfline->{INFO}->{CDSGENE}   .= "," . $csGene;
							$vcfline->{INFO}->{CDSSTRAND} .= "," . $strand;
							$vcfline->{INFO}->{CDSNUC}    .= "," . $csNuc;
							$vcfline->{INFO}->{CDSFRAME}  .= "," . $csFrame;
							$vcfline->{INFO}->{CDSAA}     .= "," . $csAA;
							$vcfline->{INFO}->{CDSFUNC}   .= "," . $csFunc;
						}
						else {
							$vcfline->{INFO}->{CDSTRANS}  = $csTrans;
							$vcfline->{INFO}->{CDSGENE}   = $csGene;
							$vcfline->{INFO}->{CDSSTRAND} = $strand;
							$vcfline->{INFO}->{CDSNUC}    = $csNuc;
							$vcfline->{INFO}->{CDSFRAME}  = $csFrame;
							$vcfline->{INFO}->{CDSAA}     = $csAA;
							$vcfline->{INFO}->{CDSFUNC}   = $csFunc;
						}
						if ( $csAA =~ /^Stop/ ) {
							$gainCount++;
							$nonSyn++;
						}
						elsif ( $csAA =~ /Stop$/ ) {
							$stopCount++;
							$nonSyn++;
						}
						elsif ( $csFunc eq "non-syn" ) {
							$missense++;
							$nonSyn++;
						}
						else {
							$syn++;
						}
					}
				}
				else {
					my (
						$indelType,    $tripletPos, $aminoAcids,
						$nomenclature, $gsymbol
					  )
					  = checkCodingIndel(
						$chrom,       $pos,       $topIndel, $name,
						$strand,      $cdsStart,  $cdsEnd,   $exonCount,
						\@exonStarts, \@exonEnds, $geneSymbol
					  );
					if ($indelType) {
						$codingHits++;
						$intronic = 0;
						if ( $vcfline->{INFO}->{CDSTRANS} ) {
							$vcfline->{INFO}->{CDSTRANS}  .= "," . $name;
							$vcfline->{INFO}->{CDSGENE}   .= "," . $gsymbol;
							$vcfline->{INFO}->{CDSSTRAND} .= "," . $strand;
							$vcfline->{INFO}->{CDSNUC}    .= "," . $nomenclature;
							$vcfline->{INFO}->{CDSFRAME}  .= "," . $tripletPos;
							$vcfline->{INFO}->{CDSAA}     .= "," . $aminoAcids;
							$vcfline->{INFO}->{CDSFUNC}   .= "," . $indelType;
						}
						else {
							$vcfline->{INFO}->{CDSTRANS}  = $name;
							$vcfline->{INFO}->{CDSGENE}   = $gsymbol;
							$vcfline->{INFO}->{CDSSTRAND} = $strand;
							$vcfline->{INFO}->{CDSNUC}    = $nomenclature;
							$vcfline->{INFO}->{CDSFRAME}  = $tripletPos;
							$vcfline->{INFO}->{CDSAA}     = $aminoAcids;
							$vcfline->{INFO}->{CDSFUNC}   = $indelType;
						}
						if ( $indelType eq "frame-shift" ) {
							$frameshift++;
						}
						else {
							$indel++;
						}
					}
				}
	
				#(4) check splice site hits
	
				#if single exon gene -> no splice sites, skip
				if ( $spliceCount == 1 ) { next; }
	
				if ( $topIndel eq "" ) {
					my (
						$dsName, $dsGene, $dsMut, $dsStrand,
						$nsName, $nsGene, $nsStrand
					  )
					  = checkSpliceSiteSNP(
						$chrom,         $spos,        $vcfline->{REF},
						$alt_nuc,       $name,        $strand,
						$cdsStart,      $cdsEnd,      $spliceCount,
						\@spliceStarts, \@spliceEnds, $geneSymbol
					  );
					if ($dsName) {
						$intronic = 0;
						if ( $vcfline->{INFO}->{DSTRANS} ) {
							$vcfline->{INFO}->{DSTRANS}  .= "," . $dsName;
							$vcfline->{INFO}->{DSGENE}   .= "," . $dsGene;
							$vcfline->{INFO}->{DSSTRAND} .= "," . $dsStrand;
							$vcfline->{INFO}->{DSMUT}    .= "," . $dsMut;
						}
						else {
							$vcfline->{INFO}->{DSTRANS}  = $dsName;
							$vcfline->{INFO}->{DSGENE}   = $dsGene;
							$vcfline->{INFO}->{DSSTRAND} = $dsStrand;
							$vcfline->{INFO}->{DSMUT}    = $dsMut;
						}
						$directSS++;
					}
	
					if ($nsName) {
						$intronic = 0;
						if ( $vcfline->{INFO}->{NSTRANS} ) {
							$vcfline->{INFO}->{NSTRANS}  .= "," . $nsName;
							$vcfline->{INFO}->{NSGENE}   .= "," . $nsGene;
							$vcfline->{INFO}->{NSSTRAND} .= "," . $nsStrand;
						}
						else {
							$vcfline->{INFO}->{NSTRANS}  = $nsName;
							$vcfline->{INFO}->{NSGENE}   = $nsGene;
							$vcfline->{INFO}->{NSSTRAND} = $nsStrand;
						}
						$nearSS++;
					}
				}
				else {
					my $spGene = checkSpliceSiteIndel(
						$chrom,       $spos,      $topIndel, $name,
						$strand,      $cdsStart,  $cdsEnd,   $exonCount,
						\@exonStarts, \@exonEnds, $geneSymbol
					);
	
					if ($spGene) {
						$intronic = 0;
						if ( $vcfline->{INFO}->{DSTRANS} ) {
							$vcfline->{INFO}->{DSTRANS}  .= "," . $name;
							$vcfline->{INFO}->{DSGENE}   .= "," . $spGene;
							$vcfline->{INFO}->{DSSTRAND} .= "," . $strand;
						}
						else {
							$vcfline->{INFO}->{DSTRANS}  = $name;
							$vcfline->{INFO}->{DSGENE}   = $spGene;
							$vcfline->{INFO}->{DSSTRAND} = $strand;
						}
						$directSS++;
	
					}
				}
	
	#if it is within a gene but not in a region specified above, the variant is intronic
				if ( $intronic == 1 ) {
					$intronHits++;
					$introns++;
					if ( $vcfline->{INFO}->{INTTRANS} ) {
						$vcfline->{INFO}->{INTTRANS}  .= "," . $name;
						$vcfline->{INFO}->{INTGENE}   .= "," . $geneSymbol;
						$vcfline->{INFO}->{INTSTRAND} .= "," . $strand;
					}
					else {
						$vcfline->{INFO}->{INTTRANS}  = $name;
						$vcfline->{INFO}->{INTGENE}   = $geneSymbol;
						$vcfline->{INFO}->{INTSTRAND} = $strand;
					}
				}
	
			}
	
			if ( $stopCount > 0 ) { 
				$func .= "nonsense,";
			}
			if ( $gainCount > 0 ) { 
				$func .= "stoploss,";
			}
	
			if ( $missense > 0 ) {  
				$func .= "missense,";
			}
			if ( $syn > 0 ) {
				$func .= "syn,";
			}
			if ( $frameshift > 0 ) {
				$func .= "frameshift,";
			}
			if ( $indel > 0 ) {
				$func .= "indel,";
			}
			if ( $directSS > 0 ) {
				$func .= "splice,";
			}
			if ( $nearSS > 0 ) {
				$func .= "nearsplice,";
			}
			if ( $utr5 > 0 ) {
				$func .= "5utr,";
			}
			if ( $utr3 > 0 ) {
				$func .= "3utr,";
			}
			if ( $noncoding > 0 ) {
				$func .= "noncoding,";
			}
			if ( $introns > 0 ) {
				$func .= "intronic,";
			}
	
			if ( $func eq "," ) {			#if no function is annotated until now the $func field contains ","
				if ( $isoforms == 0 ) {
					$func = "intergenic,";
					&getUpstreamDownstream($chr,$pos,$pos,$vcfline);
					
					
				}
				else {
					$func = "unknown,";
				}
			}
			if ( $addFunc ne "" ) {
				if ( $isoforms > 0 ) {
					$func = $addFunc . ",";
				}
				else {
					$func = "";
				}
			}
	
			#check additional annotations
			my $currentBed = BEDRecord->new(
				$chr,
				$vcfline->{POS} - 1,
				( $vcfline->{POS} - 1 + length $vcfline->{REF} )
			  )
			  ; #add the length of the reference --> for indels this assures that the whole
			    #area where the indel may occure is annotated.
			 #and take -1 since bedfiles are index origin 0 (as the databases), but VCF files are index origin 1
			foreach my $anno (
				values
				%{ $params->{settings}->{$settings}->{additionalannotations} } )
			{
				my $ovPointer;
	
				( $anno->{readbuffer}, $ovPointer ) =
				  $currentBed->readOverlappingSequences( $anno->{filepointer},
					$sequenceDict, $anno->{readbuffer} );
	
				if ( defined($ovPointer) ) {
					my @overlaps = @$ovPointer;
					if ( $anno->{info} ) {
						if ( $anno->{info}->{what} eq "bedname" ) {
							my $names = "";
							foreach my $overlap (@overlaps) {
								$names .= $overlap->name() . ",";
							}
							$names =~ s/,$//;
							$vcfline->{INFO}->{ $anno->{info}->{name} } = $names;
						}
						else {
							$vcfline->{INFO}->{ $anno->{info}->{name} } = undef;
						}
					}
					if ( $anno->{filter} ) {
						if ( $vcfline->{FILTER}[0] && $vcfline->{FILTER}[0] eq "." )
						{
							delete $vcfline->{FILTER}[0];
						}
						push( @{ $vcfline->{FILTER} }, $anno->{info}->{name} );
					}
					if ( $anno->{function} ) {
						$func .= $anno->{function} . ",";
					}
				}
			}
	
			$func =~ s/,$//;
			$func =~ s/^,//;
			$vcfline->{INFO}->{FUNC} = $func;
		}
		print OUT $vcf->format_line($vcfline);
	}
	close OUT;

}


#######################################################
#(2) getUpstreamDownstream: get closest transcripts   #
#######################################################
sub getUpstreamDownstream {
	my $chr     = shift;
	my $start   = shift;
	my $end		= shift;
	my $vcfline = shift;
	
	
	#TW 14.07.2014: if a variant is intergenic --> give the upstream and the downstream gene
	my $query = "Select T.name,T.chrom,T.strand,T.txStart,T.txEnd,T.cdsStart,T.cdsEnd,
		   T.exonCount,T.exonStarts,T.exonEnds,T.geneSymbol FROM $database.$table T $canonicaljoin
		WHERE T.chrom='$chr' 
		AND $start > T.txEnd
		ORDER BY T.txEnd DESC
		LIMIT 1;
	";
	$logger->debug( "query: " . $query );
	my $out = $dbh->prepare($query) || exit $logger->error($DBI::errstr);
	$out->execute || exit $logger->error($DBI::errstr);

			
	if (
		my (
			$name,       $chrom,    $strand, $txStart,
			$txEnd,      $cdsStart, $cdsEnd, $exonCount,
			$exonStarts, $exonEnds, $geneSymbol
		)
		= $out->fetchrow_array
	  )
	{
		$vcfline->{INFO}->{UPSTRANS}  = $name;
		$vcfline->{INFO}->{UPSGENE}   = $geneSymbol;
		$vcfline->{INFO}->{UPSSTRAND} = $strand;
		$vcfline->{INFO}->{UPSDIST}   = $start-$txEnd;
	}
					
					
	$query = "Select T.name,T.chrom,T.strand,T.txStart,T.txEnd,T.cdsStart,T.cdsEnd,
		   T.exonCount,T.exonStarts,T.exonEnds,T.geneSymbol FROM $database.$table T $canonicaljoin
		WHERE T.chrom='$chr' 
		AND $end < T.txStart
		ORDER BY T.txStart
		LIMIT 1;
	";
	$logger->debug( "query: " . $query );
	$out = $dbh->prepare($query) || exit $logger->error($DBI::errstr);
	$out->execute || exit $logger->error($DBI::errstr);
			

			
	if (
		my (
			$name,       $chrom,    $strand, $txStart,
			$txEnd,      $cdsStart, $cdsEnd, $exonCount,
			$exonStarts, $exonEnds, $geneSymbol
		)
		= $out->fetchrow_array
	  )
	{
		$vcfline->{INFO}->{DOSTRANS}  = $name;
		$vcfline->{INFO}->{DOSGENE}   = $geneSymbol;
		$vcfline->{INFO}->{DOSSTRAND} = $strand;
		$vcfline->{INFO}->{DOSDIST}   = $txStart - $end;
	}

}


#######################################################
#(2) checkNc : check for hits in noncoding transcripts#
#######################################################

sub checkNc {

	my (
		$chr,       $pos,        $name,     $strand,
		$txStart,   $txEnd,      $cdsStart, $cdsEnd,
		$exonCount, $exonStarts, $exonEnds, $geneSymbol
	  )
	  = @_;
	my $UTRprint = "";

#-> $pos between txStart and cdsStart or between cdsEnd and txEnd and within exon
#5'UTR

	for ( my $k = 0 ; $k < $exonCount ; $k++ ) {
		if ( $strand eq "+" ) {
			if ( $pos >= $$exonStarts[$k] and $pos < $$exonEnds[$k] ) {
				$logger->debug(
					"noncoding at $name,$geneSymbol,$strand,$chr,$pos");
				$UTRprint .= "noncoding at $name,$strand,$chr,$pos: ";

				return ( $name, $geneSymbol );
			}
		}
		if ( $strand eq "-" ) {
			if ( $pos >= $$exonStarts[$k] and $pos < $$exonEnds[$k] ) {
				$logger->debug(
					"noncoding at $name,$geneSymbol,$strand,$chr,$pos");
				$UTRprint .= "noncoding at $name,$strand,$chr,$pos: ";

				return ( $name, $geneSymbol );
			}
		}
	}

	return ( "", "" );

	#end sub
}

#################################################
#(2) checkUTR : check for hits in 5' and 3' UTRs#
#################################################

sub checkUTR {

	my (
		$chr,       $pos,        $name,     $strand,
		$txStart,   $txEnd,      $cdsStart, $cdsEnd,
		$exonCount, $exonStarts, $exonEnds, $geneSymbol
	  )
	  = @_;
	my $UTRprint = "";

	my ( $FutrName, $FutrGene, $TutrName, $TutrGene, );

#-> $pos between txStart and cdsStart or between cdsEnd and txEnd and within exon
#5'UTR
	for ( my $k = 0 ; $k < $exonCount ; $k++ ) {
		if ( $strand eq "+" ) {
			if (    ( $pos >= $txStart )
				and ( $pos < $cdsStart )
				and ( $pos >= $$exonStarts[$k] )
				and ( $pos < $$exonEnds[$k] ) )
			{
				$logger->debug("5UTR at $name,$geneSymbol,$strand,$chr,$pos");

				if ($FutrName) {
					$FutrName .= "," . $name;
					$FutrGene .= "," . $geneSymbol;
				}
				else {
					$FutrName = $name;
					$FutrGene = $geneSymbol;
				}
				$UTR5hits++;
				last;
			}

			#break loop after 5'UTR
			if ( $pos >= $cdsStart ) {
				last;
			}
		}
		if ( $strand eq "-" ) {
			if (    ( $pos >= $txStart )
				and ( $pos < $cdsStart )
				and ( $pos >= $$exonStarts[$k] )
				and ( $pos < $$exonEnds[$k] ) )
			{
				$logger->debug("3UTR at $name,$geneSymbol,$strand,$chr,$pos");

				#$UTRprint .= "3UTR at $name,$strand,$chr,$pos: ";
				if ($TutrName) {
					$TutrName .= "," . $name;
					$TutrGene .= "," . $geneSymbol;
				}
				else {
					$TutrName = $name;
					$TutrGene = $geneSymbol;
				}
				$UTR3hits++;
				last;
			}

			#break loop after 5'UTR
			if ( $pos >= $cdsStart ) {
				last;
			}
		}
	}

	#3'UTR -> start loop with last exon and go backwards
	for ( my $o = $exonCount - 1 ; $o >= 0 ; $o-- ) {
		if ( $strand eq "+" ) {
			if (    ( $pos >= $cdsEnd )
				and ( $pos < $txEnd )
				and ( $pos >= $$exonStarts[$o] )
				and ( $pos < $$exonEnds[$o] ) )
			{
				$logger->debug("3UTR at $name,$geneSymbol,$strand,$chr,$pos");

				#$UTRprint .= "3UTR at $name,$strand,$chr,$pos: ";
				if ($TutrName) {
					$TutrName .= "," . $name;
					$TutrGene .= "," . $geneSymbol;
				}
				else {
					$TutrName = $name;
					$TutrGene = $geneSymbol;
				}
				$UTR3hits++;
				last;
			}

			#break loop after 3'UTR
			if ( $pos < $cdsEnd ) {
				last;
			}
		}
		if ( $strand eq "-" ) {
			if (    ( $pos >= $cdsEnd )
				and ( $pos < $txEnd )
				and ( $pos >= $$exonStarts[$o] )
				and ( $pos < $$exonEnds[$o] ) )
			{
				$logger->debug("5UTR at $name,$geneSymbol,$strand,$chr,$pos");

				#$UTRprint .= "5UTR at $name,$strand,$chr,$pos: ";
				if ($FutrName) {
					$FutrName .= "," . $name;
					$FutrGene .= "," . $geneSymbol;
				}
				else {
					$FutrName = $name;
					$FutrGene = $geneSymbol;
				}
				$UTR5hits++;
				last;
			}

			#break loop after 3'UTR
			if ( $pos < $cdsEnd ) {
				last;
			}
		}
	}
	return ( $FutrName, $FutrGene, $TutrName, $TutrGene );

	#end sub
}

################################################
#(3) checkCoding: identify hits in cds sequence#
################################################
sub checkCodingIndel {

	my (
		$chr,        $pos,      $indel,  $name,
		$strand,     $cdsStart, $cdsEnd, $exonCount,
		$exonStarts, $exonEnds, $geneSymbol
	  )
	  = @_;
	my $codingprint  = "";
	my $indelType    = "";
	my $aminoAcids   = "";
	my $nomenclature = "";
	my $gsymbol      = "";
	my $tripletPos;
	
	my $indelLength = length($indel) - 1;
	my $end = $pos;
	$end += length($indel) - 2 if $indel =~ /^-/;		#add endposition for deletion

	if ( (( $pos >= $cdsStart ) and ( $pos < $cdsEnd )) || (( $end >= $cdsStart ) and ( $end < $cdsEnd )) ) {
		for ( my $i = 0 ; $i < $exonCount ; $i++ ) {
			if (   (( $pos >= $$exonStarts[$i] )
				and ( $pos < $$exonEnds[$i] )) || (( $end >= $$exonStarts[$i] ) and ( $end < $$exonEnds[$i] )) )
			{
				
				#TW 26.11.2012: if not the whole deletion is within an exon (or coding sequence), extract the part of the indel that lies within the exon/cds
				if($pos < $$exonStarts[$i]){
					$indel   = "-".substr($indel, ($$exonStarts[$i] - $pos)+1, length($indel)-(($$exonStarts[$i] - $pos)+1));
					$pos = $$exonStarts[$i];
				}
				if($pos < $cdsStart){
					$indel   = "-".substr($indel, ($cdsStart - $pos)+1, length($indel)-(($cdsStart - $pos)+1));
					$pos = $cdsStart;
				}
				if( ($pos + length($indel) - 1) > $$exonEnds[$i]){
					$indel = substr($indel,0,($$exonEnds[$i]-$pos + 1 ) )
				}
				if( ($pos + length($indel) - 1) > $cdsEnd){
					$indel = substr($indel,0,($cdsEnd-$pos + 1) )
				}
				
				( $indelType, $tripletPos, $aminoAcids, $nomenclature ) =
				  &classifyCodingHit(
					$indel,      $chr,      $pos,    $name,
					$strand,     $cdsStart, $cdsEnd, $exonCount,
					$exonStarts, $exonEnds
				  );

				if ($geneSymbol) {
					$gsymbol = $geneSymbol;
				}
				$logger->debug(
"$name,$strand,$nomenclature,$tripletPos,$indel,$indelType,$aminoAcids"
				);
				$nomenclature = "size_".$indelLength if $indelLength > 20;	#add size of indel instead of sequence if the indel is long
				return (
					$indelType,    $tripletPos, $aminoAcids,
					$nomenclature, $gsymbol
				);
			}
		}
	}
	return ( $indelType, $tripletPos, $aminoAcids, $nomenclature, $gsymbol );

	#end sub
}

########################################################################################
#(3) classifyCodingHit: determine exact consequences of coding indel (frame shift etc.)#
########################################################################################
sub classifyCodingHit {

	my (
		$indel,    $chrom,  $pos,       $name,       $strand,
		$cdsStart, $cdsEnd, $exonCount, $exonStarts, $exonEnds
	  )
	  = @_;
	my $tmpExonStart  = 0;
	my $i             = 0;
	my @start         = ();
	my @end           = ();
	my $num           = 0;
	my $cdspos        = 0;
	my $frame         = 0;
	my $aminoAcids    = "";
	my $cds           = "";
	my %code          = ();
	my $lengthIndel   = 0;
	my $indelSeq      = "";
	my $indelClass    = "";
	my $del           = 0;
	my $cds_location  = "";
	my $cds_start_pos = 0;
	my $cds_end_pos   = 0;

	if ( $chrom eq "chrM" )
	{
		&Utilities::init_code_mtdna( \%code );
	}
	else
	{
		&Utilities::init_code( \%code );
	}


	#parse indel length and sequence
	$indel =~ s/\+//;

	#if deletion remove leading "-"
	if ( $indel =~ /^-/ ) {
		$del = 1;
		$indel =~ s/\-//;
	}
	$lengthIndel = length($indel);
	$indelSeq    = $indel;

	#get coding sequence
	my $query = "Select T.cds FROM $database.$cdstable T $canonicaljoin WHERE T.name='$name'";
	my $out = $dbh->prepare($query) || exit $logger->error($DBI::errstr);
	$out->execute || exit $logger->error($DBI::errstr);
	$cds = uc $out->fetchrow_array();

	# an welcher Codon-Position ist pos
	if ( $strand eq "+" ) {
		for ( $i = 0 ; $i < $exonCount ; $i++ ) {
			if (    ( $cdsStart >= $$exonStarts[$i] )
				and ( $pos < $$exonEnds[$i] ) )
			{    #pos in start exon
				push( @start, $cdsStart );
				push( @end,   $pos );
			}
			elsif ( ( $cdsStart >= $$exonStarts[$i] )
				and ( $cdsStart < $$exonEnds[$i] ) )
			{    #start exon
				push( @start, $cdsStart );
				push( @end,   $$exonEnds[$i] );
			}
			elsif ( ( $$exonStarts[$i] > $cdsStart )
				and ( $pos > $$exonEnds[$i] - 1 ) )
			{    #other exons
				push( @start, $$exonStarts[$i] );
				push( @end,   $$exonEnds[$i] );
			}
			elsif ( ( $pos >= $$exonStarts[$i] )
				and ( $pos < $$exonEnds[$i] ) )
			{    #pos in other exon
				push( @start, $$exonStarts[$i] );
				push( @end,   $pos );
			}
		}
		$num    = @start;
		$cdspos = 0;
		for ( $i = 0 ; $i < $num ; $i++ ) {
			$cdspos = $cdspos + $end[$i] - $start[$i];    #exonEnd is +1 in UCSC
		}

		$cdspos++;                                        #pos ist exact end
		$frame = ( $cdspos - 1 ) % 3;
	}
	if ( $strand eq "-" ) {
		for ( $i = 0 ; $i < $exonCount ; $i++ ) {
			if (    ( $pos >= $$exonStarts[$i] )
				and ( $cdsEnd <= $$exonEnds[$i] ) )
			{    #first exon with pos and cdsEnd
				push( @start, $pos );
				push( @end,   $cdsEnd );
			}
			elsif ( ( $pos >= $$exonStarts[$i] )
				and ( $pos < $$exonEnds[$i] ) )
			{    #first exon with pos
				push( @start, $pos );
				push( @end,   $$exonEnds[$i] );
			}
			elsif ( ( $$exonStarts[$i] > $pos )
				and ( $$exonEnds[$i] < $cdsEnd ) )
			{    #other exons
				push( @start, $$exonStarts[$i] );
				push( @end,   $$exonEnds[$i] );
			}
			elsif ( ( $cdsEnd >= $$exonStarts[$i] )
				and ( $cdsEnd <= $$exonEnds[$i] ) )
			{    #end exon
				push( @start, $$exonStarts[$i] );
				push( @end,   $cdsEnd );
			}
		}
		$num    = @start;
		$cdspos = 0;
		for ( $i = 0 ; $i < $num ; $i++ ) {
			$cdspos = $cdspos + $end[$i] - $start[$i];    #exonEnd is +1 in UCSC
		}

		$frame = ( $cdspos - 1 ) % 3;
	}

#if indel starts at first pos of a triplett and has a length of factor 3 -> no frameshift, insertion/deletion of specific amino acids
	if ( ( $frame == 0 ) && ( $lengthIndel % 3 == 0 ) ) {
		my $aas = ( $lengthIndel / 3 );

		if ($del) {
			$indelClass = "del_" . $aas . "aa";
		}
		else {
			$indelClass = "ins_" . $aas . "aa";
		}

		#translate inserted / deleted tripletts in amino acids
		for ( my $j = 0 ; $j < $aas ; $j++ ) {
			my $triplett = substr( $indelSeq, $j * 3, 3 );

			$aminoAcids .= $code{ uc $triplett };
		}
	}
	elsif ( ( $lengthIndel % 3 == 0 ) ) {
		$indelClass = "indel";
	}
	else {
		$indelClass = "frame-shift";
	}

	#calculate cds_location for cariation nomenclature standard:
	if ($del) {

	  #if deletion: cds_location = pos of the first and last deleted nucleotides
		
		if($strand eq "-"){
			$cdspos -= ($lengthIndel -1);
		}
		$cds_start_pos = $cdspos;
		$cds_end_pos   = $cdspos + $lengthIndel - 1;
		$cds_location = $cds_start_pos . "_" . $cds_end_pos . "del" . $indelSeq;
	}
	else {

		#if insertion: cds_location = pos of the flanking nucleotides
		if($strand eq "-"){
			$cds_start_pos = $cdspos ;
			$cds_end_pos   = $cdspos + 1;
		}
		else{
			$cds_start_pos = $cdspos - 1;
			$cds_end_pos   = $cdspos;
		}
		$cds_location = $cds_start_pos . "_" . $cds_end_pos . "ins" . $indelSeq;
	}

	return ( $indelClass, $frame, $aminoAcids, $cds_location );
}

#############################################
#(3) checkCds: identify hits in cds sequence#
#############################################
sub checkCds {

	my (
		$chrom,     $pos,        $nuc,      $alt_nuc,
		$name,      $strand,     $cdsStart, $cdsEnd,
		$exonCount, $exonStarts, $exonEnds, $geneSymbol
	  )
	  = @_;
	my ( $frame, $out_nuc, $out_aa );
	my $codingprint = "";
	my $indelType   = "";
	my $aminoAcids  = "";
	my $tripletPos;
	my ( $retTrans, $retGene, $retNuc, $retFrame, $retAA );

	if ( ( $pos >= $cdsStart ) and ( $pos < $cdsEnd ) ) {
		for ( my $i = 0 ; $i < $exonCount ; $i++ ) {
			if (    ( $pos >= $$exonStarts[$i] )
				and ( $pos < $$exonEnds[$i] ) )
			{

				( $frame, $out_nuc, $out_aa ) = &checkCoding(
					$chrom,     $pos,        $nuc,      $alt_nuc,
					$name,      $strand,     $cdsStart, $cdsEnd,
					$exonCount, $exonStarts, $exonEnds
				);

				my ( $func, $aa ) = split( ",", $out_aa );
				return ( $name, $geneSymbol, $out_nuc, $frame, $func, $aa );
			}
		}
	}

	#return $codingprint;

	#end sub
}

##################################################################################################
#(3) checkCoding: check for SNPs in cds and determine synonymous or non-synonymous substituitions#
##################################################################################################
sub checkCoding {
	my $chrom        = shift;
	my $pos          = shift;
	my $nuc          = shift;
	my $alt_nuc      = shift;
	my $name         = shift;
	my $strand       = shift;
	my $cdsStart     = shift;
	my $cdsEnd       = shift;
	my $exonCount    = shift;
	my $exonStarts   = shift;
	my $exonEnds     = shift;
	my $tmpExonStart = 0;
	my $i            = 0;
	my @start        = ();
	my @end          = ();
	my $num          = 0;
	my $cdspos       = 0;       # base position after coding start
	my $frame        = 0;       # 0=first base in triplet
	my $out_nuc      = "";
	my $out_aa       = "";
	my $cds          = "";

	#get coding sequence
	#TODO : improper usage of prepare  
	my $query = "Select T.cds FROM $database.$cdstable T $canonicaljoin WHERE T.name='$name'";
	my $out = $dbh->prepare($query) || exit $logger->error($DBI::errstr);
	$out->execute || exit $logger->error($DBI::errstr);

	#print "$name $cds\n";
	$cds = uc $out->fetchrow_array();

	if ( $strand eq "+" ) {    # at which Codon-Position is pos
		for ( $i = 0 ; $i < $exonCount ; $i++ ) {
			if (    ( $cdsStart >= $$exonStarts[$i] )
				and ( $pos < $$exonEnds[$i] ) )
			{                  #pos in start exon
				push( @start, $cdsStart );
				push( @end,   $pos );
			}
			elsif ( ( $cdsStart >= $$exonStarts[$i] )
				and ( $cdsStart < $$exonEnds[$i] ) )
			{                  #start exon
				push( @start, $cdsStart );
				push( @end,   $$exonEnds[$i] );
			}
			elsif ( ( $$exonStarts[$i] > $cdsStart )
				and ( $pos > $$exonEnds[$i] - 1 ) )
			{                  #other exons
				push( @start, $$exonStarts[$i] );
				push( @end,   $$exonEnds[$i] );
			}
			elsif ( ( $pos >= $$exonStarts[$i] )
				and ( $pos < $$exonEnds[$i] ) )
			{                  #pos in other exon
				push( @start, $$exonStarts[$i] );
				push( @end,   $pos );
			}
		}
		$num    = @start;
		$cdspos = 0;
		for ( $i = 0 ; $i < $num ; $i++ ) {
			$cdspos =
			  $cdspos + $end[$i] - $start[$i]
			  ;                #exonEnd ist immer 1 groesser in UCSC
			                   #print " $cdspos $end[$i] $start[$i] \n";

		}
		$cdspos++;             #pos ist exact end
		$frame = ( $cdspos - 1 ) % 3;
		( $out_nuc, $out_aa ) =
		  &get_triplet( $cds, $chrom, $cdspos, $strand, $frame, $nuc, $alt_nuc,
			$name );
	}
	if ( $strand eq "-" ) {
		for ( $i = 0 ; $i < $exonCount ; $i++ ) {
			if (    ( $pos >= $$exonStarts[$i] )
				and ( $cdsEnd <= $$exonEnds[$i] ) )
			{    #first exon with pos and cdsEnd
				push( @start, $pos );
				push( @end,   $cdsEnd );
			}
			elsif ( ( $pos >= $$exonStarts[$i] )
				and ( $pos < $$exonEnds[$i] ) )
			{    #first exon with pos
				push( @start, $pos );
				push( @end,   $$exonEnds[$i] );
			}
			elsif ( ( $$exonStarts[$i] > $pos )
				and ( $$exonEnds[$i] < $cdsEnd ) )
			{    #other exons
				push( @start, $$exonStarts[$i] );
				push( @end,   $$exonEnds[$i] );
			}
			elsif ( ( $cdsEnd >= $$exonStarts[$i] )
				and ( $cdsEnd <= $$exonEnds[$i] ) )
			{    #end exon
				push( @start, $$exonStarts[$i] );
				push( @end,   $cdsEnd );
			}
		}
		$num    = @start;
		$cdspos = 0;
		for ( $i = 0 ; $i < $num ; $i++ ) {
			$cdspos =
			  $cdspos + $end[$i] - $start[$i]
			  ;    #exonEnd ist immer 1 groesser in UCSC
		}

		$frame = ( $cdspos - 1 ) % 3;
		( $out_nuc, $out_aa ) =
		  &get_triplet( $cds, $chrom, $cdspos, $strand, $frame, $nuc, $alt_nuc,
			$name );

	}

	return ( $frame, $out_nuc, $out_aa );
}

############################################################################################
#(3) get_triplet: translate codons to amino acids and detemine syn or non-syn substitutions#
############################################################################################
sub get_triplet {
	my $cds         = shift;
	my $chrom       = shift;
	my $pos         = shift;
	my $strand      = shift;
	my $frame       = shift;
	my $nuc         = shift;
	my $alt_nuc     = shift;
	my $name        = shift;
	my %code        = ();
	my $aa          = "";
	my $alt_aa      = "";
	my $out_nuc     = "";
	my $out_aa      = "";
	my $pos_aa      = "";      # aa position after start codon
	my $start       = "";
	my $end         = "";
	my $i           = 0;
	my $triplet     = "";
	my $alt_triplet = "";
	my $out_syn     = "syn";
	
	if ( $chrom eq "chrM" )
	{
		&Utilities::init_code_mtdna( \%code );
	}
	else
	{
		&Utilities::init_code( \%code );
	}

	# aa position after start codon
	$pos_aa = int( $pos / 3 );
	if ( ( $frame == 0 ) or ( $frame == 1 ) ) {
		$pos_aa++;
	}

	# define postion of triplet
	if ( $frame == 0 ) {
		$start = $pos;
		$end   = $pos + 3;
	}
	if ( $frame == 1 ) {
		$start = $pos - 1;
		$end   = $pos + 2;
	}
	if ( $frame == 2 ) {
		$start = $pos - 2;
		$end   = $pos + 1;
	}

	$triplet = substr( $cds, $start - 1, 3 );

	# get aminio acids for triplet and alt_triplet
	if ( $triplet ne "" ) {
		$alt_triplet = $triplet;
		$alt_nuc     = &Utilities::iub($alt_nuc);
		$alt_nuc =~ s/$nuc//;
		$aa      = $code{$triplet};
		$out_nuc = $pos . $nuc . "->" . $alt_nuc;

		for ( $i = 0 ; $i < length($alt_nuc) ; $i++ ) {
			if ( $strand eq "+" ) {
				substr( $alt_triplet, $frame, 1 ) = substr( $alt_nuc, $i, 1 );
			}
			elsif ( $strand eq "-" ) {
				$alt_nuc = &Utilities::comp_rev($alt_nuc);
				substr( $alt_triplet, $frame, 1 ) = substr( $alt_nuc, $i, 1 );
			}

			$alt_aa = $code{$alt_triplet};
			if ( $aa ne $alt_aa ) {
				$out_syn = "non-syn";
			}
			if ( $i == 0 ) {
				$out_aa = $aa . "$pos_aa" . $alt_aa;
			}
			else {
				$out_aa = $out_aa . "/" . $alt_aa;
			}
		}
		$out_aa = $out_syn . "," . $out_aa;
	}

	$logger->debug(
		"$chrom $frame $strand $triplet $alt_triplet $out_nuc $out_aa");

	if ( $out_syn eq "syn" ) {
		$synHits++;
	}
	if ( $out_syn eq "non-syn" ) {
		$nonSynHits++;

		#for each gene record number of non syn cds hits
		if ( !exists( $hitCountHash{$name} ) ) { $hitCountHash{$name} = 1; }
		else { $hitCountHash{$name}++; }
	}

	return ( $out_nuc, $out_aa );
}

############################################################
#(4) checkSplice: identify indels in canonical splice sites#
############################################################
sub checkSpliceSiteIndel {

	my (
		$chr,        $spos,     $indel,  $name,
		$strand,     $cdsStart, $cdsEnd, $exonCount,
		$exonStarts, $exonEnds, $geneSymbol
	  )
	  = @_;
	my $spliceprint = "";
	
	
	
	$spos++;		#TW 26.11.2012: required, because actual indel starts one position after the position in VCF
	my $end = $spos;
	$end += length($indel) - 2 if $indel =~ /^-/;		#add endposition for deletion

	#splice sites between exonEnds[i] and exonStarts[i+1]
	for ( my $j = 0 ; $j < $exonCount - 1 ; $j++ ) {
		if ( $spos <= $$exonEnds[$j] ) { last; }
		if (   ( $spos == $$exonEnds[$j] + 1 )
			or ( $spos == $$exonEnds[$j] + 2 )
			or ( $end  == $$exonEnds[$j] + 1 )
			or ( $end  == $$exonEnds[$j] + 2 )
			or ( $spos < $$exonEnds[$j] + 1 && $end > $$exonEnds[$j] + 2)
			or ( $spos == $$exonStarts[ $j + 1 ] )
			or ( $spos == $$exonStarts[ $j + 1 ] - 1 ) 
			or ( $end  == $$exonStarts[ $j + 1 ] - 1 )
			or ( $end  == $$exonStarts[ $j + 1 ] )
			or ( $spos < $$exonStarts[ $j + 1 ] - 1 && $end > $$exonStarts[ $j + 1 ] )
			
			)
		{

#call checkSplice
			$logger->debug("splice site: $name,$chr,$spos,$strand,$indel");
			my $gsymbol = "";
			if ($geneSymbol) {
				$gsymbol = $geneSymbol;
			}
			return $gsymbol;
		}
	}

	return;

	#end sub
}

##############################################################
#(4) checkSpliceSite: identify SNPs at canonical splice sites#
##############################################################
sub checkSpliceSiteSNP {

	my (
		$chr,       $spos,       $nuc,      $alt_nuc,
		$name,      $strand,     $cdsStart, $cdsEnd,
		$exonCount, $exonStarts, $exonEnds, $geneSymbol
	  )
	  = @_;
	my $spliceprint    = "";
	my $spliceMutation = "";
	my ( $dsName, $dsGene, $dsMut, $dsStrand );
	my ( $nsName, $nsGene, $nsMut, $nsStrand );

	#splice sites between exonEnds[i] and exonStarts[i+1]
	for ( my $j = 0 ; $j < $exonCount - 1 ; $j++ ) {
		if ( $spos <= ($$exonEnds[$j]-3) ) { last; }		#changed
		if (   ( $spos == $$exonEnds[$j] + 1 )
			or ( $spos == $$exonEnds[$j] + 2 )
			or ( $spos == $$exonStarts[ $j + 1 ] )
			or ( $spos == $$exonStarts[ $j + 1 ] - 1 ) )
		{

			#call checkSplice
			$spliceMutation =
			  &checkSplice( $chr, $spos, $nuc, $alt_nuc, $name, $strand,
				$exonCount, $$exonStarts[ $j + 1 ],
				$$exonEnds[$j] );

			$spliceprint .=
"direct splice site at $name,$geneSymbol,$chr,$spos,$strand,$spliceMutation: ";
			if ($dsGene) {
				$dsName   .= "," . $name;
				$dsGene   .= "," . $geneSymbol;
				$dsMut    .= "," . $spliceMutation;
				$dsStrand .= "," . $strand;
			}
			else {
				$dsName   = $name;
				$dsGene   = $geneSymbol;
				$dsMut    = $spliceMutation;
				$dsStrand = $strand;
			}
		}

		#SNPs in close poximity to canonical splice sites
		elsif (($spos > $$exonEnds[$j] + 2 && $spos <= $$exonEnds[$j] + $spliceWindow + 3)
			|| ($spos < $$exonStarts[ $j + 1 ] - 1	&& $spos >= $$exonStarts[ $j + 1 ] - $spliceWindow - 2 )
			|| ($spos > $$exonEnds[$j] - 3 && $spos <= $$exonEnds[$j])
			|| ($spos < $$exonStarts[ $j +1 ] + 4 && $spos >= $$exonStarts[ $j + 1 ] +1))
		{
			$spliceprint .= "SNP near splice Site at $name,$geneSymbol,$chr,$spos,$strand: ";
			if ($nsGene) {
				$nsName   .= "," . $name;
				$nsGene   .= "," . $geneSymbol;
				$nsStrand .= "," . $strand;
			} else {
				$nsName   = $name;
				$nsGene   = $geneSymbol;
				$nsStrand = $strand;
			}
		}
	}
	return ( $dsName, $dsGene, $dsMut, $dsStrand, $nsName, $nsGene, $nsStrand );

	#end sub
}
########################################################################
#(4) checkSplice : determine nucleotide substitution at splice site SNP#
########################################################################
sub checkSplice {
	my $chrom     = shift;
	my $pos       = shift;
	my $nuc       = shift;
	my $alt_nuc   = shift;
	my $name      = shift;
	my $strand    = shift;
	my $exonCount = shift;
	my $exonStart = shift;
	my $exonEnd   = shift;

	my $start_splice = 0;
	my $end_splice   = 0;
	my $num          = 0;
	my $coding       = 0;
	my $coding2      = 0;

	my $aa             = "";
	my $spliceSite     = "";
	my $spliceMutation = "";

	#heterozygote snp call translation
	$alt_nuc = &Utilities::iub($alt_nuc);
	$alt_nuc =~ s/$nuc//;

	#forward strand
	if ( $strand eq "+" ) {

#splice sites are 2 nucleotides after exon end i and 2 nucleotides before exon start i+1
		if ( ( $pos == $exonEnd + 1 ) or ( $pos == $exonEnd + 2 ) ) {
			$logger->debug("splice site,$name,$chrom,$strand,$pos, ");

			$start_splice = $exonEnd + 1;
			$end_splice   = $exonEnd + 2;    # +2;

			$spliceSite =
			  Utilities::twoBit( $chrom, $start_splice, $end_splice, $strand,
				$settings );
			chomp($spliceSite);
			$spliceSite = uc($spliceSite);
		}

		if ( ( $pos == $exonStart ) or ( $pos == $exonStart - 1 ) ) {
			$logger->debug("splice site,$name,$chrom,$strand, $pos, ");

			$start_splice = $exonStart - 1;    # -2;
			$end_splice   = $exonStart;

			$spliceSite =
			  Utilities::twoBit( $chrom, $start_splice, $end_splice, $strand,
				$settings );
			chomp($spliceSite);
			$spliceSite = uc($spliceSite);
			
		}
		#generate output
		if ( $pos == $exonEnd + 1 ) {
			$spliceMutation =
			  $spliceSite . "->" . $alt_nuc . substr( $spliceSite, 1, 1 );
			$logger->debug("$spliceMutation");
			return $spliceMutation;
		}
		if ( $pos == $exonEnd + 2 ) {
			$spliceMutation =
			  $spliceSite . "->" . substr( $spliceSite, 0, 1 ) . $alt_nuc;
			$logger->debug("$spliceMutation");
			return $spliceMutation;
		}
		if ( $pos == $exonStart ) {
			$spliceMutation =
			  $spliceSite . "->" . substr( $spliceSite, 0, 1 ) . $alt_nuc;
			$logger->debug("$spliceMutation");
			return $spliceMutation;
		}
		if ( $pos == $exonStart - 1 ) {
			$spliceMutation =
			  $spliceSite . "->" . $alt_nuc . substr( $spliceSite, 1, 1 );
			$logger->debug("$spliceMutation");
			return $spliceMutation;
		}
	}

	#reverse strand
	if ( $strand eq "-" ) {

#splice sites are 2 nucleotides after exon end i and 2 nucleotides before exon start i+1
		if ( ( $pos == $exonEnd + 1 ) or ( $pos == $exonEnd + 2 ) ) {
			$logger->debug("splice site,$name,$chrom,$strand,$pos,");

			$start_splice = $exonEnd + 1;
			$end_splice   = $exonEnd + 2;    # +2;

			$spliceSite =
			  Utilities::twoBit( $chrom, $start_splice, $end_splice, $strand,
				$settings );
			chomp($spliceSite);
			$spliceSite = uc($spliceSite);
		}

		if ( ( $pos == $exonStart ) or ( $pos == $exonStart - 1 ) ) {
			$logger->debug("splice site,$name,$chrom,$strand,$pos,");

			$start_splice = $exonStart - 1;    # -2;
			$end_splice   = $exonStart;

			$spliceSite =
			  Utilities::twoBit( $chrom, $start_splice, $end_splice, $strand,
				$settings );
			chomp($spliceSite);
			$spliceSite = uc($spliceSite);
		}

		#generate output
		if ( $pos == $exonEnd + 1 ) {
			$spliceMutation =
			    $spliceSite . "->"
			  . substr( $spliceSite, 0, 1 )
			  . Utilities::comp_rev($alt_nuc);
			$logger->debug("$spliceMutation");
			return $spliceMutation;
		}
		if ( $pos == $exonEnd + 2 ) {
			$spliceMutation =
			    $spliceSite . "->"
			  . Utilities::comp_rev($alt_nuc)
			  . substr( $spliceSite, 1, 1 );
			$logger->debug("$spliceMutation");
			return $spliceMutation;
		}
		if ( $pos == $exonStart ) {
			$spliceMutation =
			    $spliceSite . "->"
			  . Utilities::comp_rev($alt_nuc)
			  . substr( $spliceSite, 1, 1 );
			$logger->debug("$spliceMutation");
			return $spliceMutation;
		}
		if ( $pos == $exonStart - 1 ) {
			$spliceMutation =
			    $spliceSite . "->"
			  . substr( $spliceSite, 0, 1 )
			  . Utilities::comp_rev($alt_nuc);
			$logger->debug("$spliceMutation");
			return $spliceMutation;
		}
	}

	#end subroutine
}

###########################################################
#(5) calculateCdsLength: calculates coding sequence length#
###########################################################
sub calculateCdsLength {
	my ( $cdsStart, $cdsEnd, $exonStarts, $exonEnds, $exonCount, $strand ) = @_;
	my $cdsLength = 0;

	if ( $strand eq "+" ) {    #get coding sequence length
		for ( my $i = 0 ; $i < $exonCount ; $i++ ) {

			#single exon gene
			if (    ( $cdsStart >= $$exonStarts[$i] )
				and ( $cdsEnd < $$exonEnds[$i] ) )
			{
				$cdsLength += $cdsEnd - $cdsStart;
			}

			#exon with cdsStart
			elsif ( ( $cdsStart >= $$exonStarts[$i] )
				and ( $cdsStart < $$exonEnds[$i] ) )
			{
				$cdsLength += $$exonEnds[$i] - $cdsStart;
			}

			#other coding exons
			elsif ( ( $cdsStart < $$exonStarts[$i] )
				and ( $cdsEnd > $$exonEnds[$i] ) )
			{
				$cdsLength += $$exonEnds[$i] - $$exonStarts[$i];
			}

			#exon with cdsEnd
			elsif ( ( $cdsEnd > $$exonStarts[$i] )
				and ( $cdsEnd <= $$exonEnds[$i] ) )
			{
				$cdsLength += $cdsEnd - $$exonStarts[$i];
			}
		}
	}
	if ( $strand eq "-" ) {    #get coding sequence length
		for ( my $i = $exonCount - 1 ; $i >= 0 ; $i-- ) {

			#single exon gene
			if (    ( $cdsStart >= $$exonStarts[$i] )
				and ( $cdsEnd < $$exonEnds[$i] ) )
			{
				$cdsLength += $cdsEnd - $cdsStart;
			}

			#exon with cdsStart
			elsif ( ( $cdsStart >= $$exonStarts[$i] )
				and ( $cdsStart < $$exonEnds[$i] ) )
			{
				$cdsLength += $$exonEnds[$i] - $cdsStart;
			}

			#internal coding exons
			elsif ( ( $cdsStart < $$exonStarts[$i] )
				and ( $cdsEnd > $$exonEnds[$i] ) )
			{
				$cdsLength += $$exonEnds[$i] - $$exonStarts[$i];
			}

			#exon with cdsEnd
			elsif ( ( $cdsEnd > $$exonStarts[$i] )
				and ( $cdsEnd <= $$exonEnds[$i] ) )
			{
				$cdsLength += $cdsEnd - $$exonStarts[$i];
			}
		}
	}
	return $cdsLength;

	#end sub
}


=head1 NAME

annotateVCF.pl

=head1 SYNOPSIS

annotateVCF.pl -i variants.dbSNP.vcf

=head1 DESCRIPTION

This script annotates a VCF file with information such as in which gene a variant lies, which amino acids (if any)
are changed, which presumable function a variant has, ...

=head1 OPTIONS

 -i	<infile.vcf> file containing the variants to annotate; REQUIRED
 -o	<outfile.vcf> output file; default: infile.plus.vcf
 -m	<annotation mode> currently there are three possible values:
 		"short": all variants are treated as short (i.e. short indels and SNVs) and the exact annotation (e.g. amino acid change, ...) are produced
 		"cnv": only the overlap with genes is reported (useful for large SVs/CNVs); for this mode the "END"" field must be defined in the INFO part of the VCF
 				line (see definition of SV representation in VCF from the 1000G project)
 		"auto" (DEFAULT): "short" annotations are produced for all normal VCF variants, i.e. variants where the REF and the ALT field contain only A,C,G or T and 
 				"cnv" annotations are produced for the other variants. This mode is useful for instance for variants called by Pindel where both short and long Indels
 				can be present in a single file and longer variants are encoded as in the 1000G definition (e.g. a long deletion has <DEL> in the alt field).
 -se	name of the settings in the current.config.xml file that holds path to reference genome, 
 	to the annotation file and to possible additional annotation files; use default settings if nothing is given
 -t	<knownGeneTable> is taken from settings file by default
 -c	<codingSequenceTable> is taken from settings file by default
 -r	enables use of RefSeqAli table as annotation basis (from settings file)
 -l	enables use of lincRNA table as annotation basis (from settings file)
 -mir enables use of miRNA table as annotation basis (from settings file)
 -f	function that will be written in function field if a variant lies within a gene (coding,noncoding,intronic), useful e.g. for lincRNAs, miRNAs;
	note that if this option is chosen the annotated function will be only the here specified or no function!
 -w	<window size> size of window around splice sites [20]
 -canonical		annotate with canonical transcripts (only with UCSC tables)
 -lf	log file; default: print to screen
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page
 

=head1 AUTHOR

Thomas Wieland
Riccardo Berutti

=cut