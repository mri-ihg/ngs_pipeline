#!/usr/bin/perl

use strict;

#use warnings;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);
use Bio::DB::Sam;
use List::Util qw(sum);

my $prog_path = dirname( abs_path($0) );
require $prog_path . "/Utilities.pm";

#my $vcf = "test.all.bcf.out.coding";
#my $qual = "test.all.bcf.out.coding.qual";
#$qual = "test.bcf.out.coding.part";
#my $cutoff = 2;

#my $out = "$vcf.filtered$cutoff";
#my $lowQual = "$vcf.filtered$cutoff.lowQual";

my $logfile  = "pipeline.log";
my $loglevel = "INFO";
my $vcffile  = "";
my $help     = 0;
my $settings = "";
my $out      = "";
my @columns  = ();
my $window     = 10;
my $bamfile  = "";

my $helptext = "-i	<vcf file to check> 
-b	<bamfile> to extract read depth; default: merged.rmdup.bam in folder of the infile
-o	outfile
-w	how many bases in direction of read end to check; default: $window 
-se	settings
-lf	log file; default: pipeline.log
-ll	log level: ERROR,INFO,DEBUG; default: INFO
-h	this help\n";

GetOptions(
	"i=s"  => \$vcffile,
	"b=s"  => \$bamfile,
	"o=s"  => \$out,
	"w=s"  => \$window,
	"lf=s" => \$logfile,
	"ll=s" => \$loglevel,
	"se=s" => \$settings,
	"h"    => \$help
);

if ( $help == 1 ) {
	print $helptext;
	exit(1);
}

Utilities::initLogger( $logfile, $loglevel );
my $logger = Utilities::getLogger();
my $params = Utilities::getParams();

require $params->{programs}->{vcftools}->{pm};


if ( $vcffile eq "" || !-e $vcffile ) {
	$logger->error("VCF file not given or it doesn't exist");
	exit(1);
}

my (
	$vcfline,    $qualline, $pos,        $depth, $altdepth,
	$medianqual, $medread,  $strandfreq, $varpercent
);

my $outdir    = dirname($vcffile);
my $refgenome = $params->{settings}->{$settings}->{reference};

$bamfile = "$outdir/merged.rmdup.bam" if $bamfile eq "";

my $sam = Bio::DB::Sam->new(
	-fasta => $refgenome,
	-bam   => $bamfile
);

#open(VCF,"$vcffile")|| die "Cannot open $vcffile\n";
#open(QUA,"$prog_path/checkQual_api.pl -se $settings -i $vcffile -lf $logfile -ll $loglevel|")|| die "Cannot open stream on checkQual_api.pl\n";
open( OUT, ">$out" ) || die "Cannot open $out\n";

my $vcf = Vcf->new( file => $vcffile );
$vcf->parse_header();

$vcf->add_header_line(
	{
		key         => 'FILTER',
		ID          => 'Q20',
		Description =>
		  'Filtered out by filterSNPqual.pl because median quality < 20'
	}
);    #add lines to header
$vcf->add_header_line(
	{
		key         => 'FILTER',
		ID          => 'Q15',
		Description =>
		  'Filtered out by filterSNPqual.pl because median quality < 15'
	}
);
$vcf->add_header_line(
	{
		key         => 'FILTER',
		ID          => 'Q10',
		Description =>
		  'Filtered out by filterSNPqual.pl because median quality < 10'
	}
);
$vcf->add_header_line(
	{
		key         => 'FILTER',
		ID          => 'Q5',
		Description =>
		  'Filtered out by filterSNPqual.pl because median quality < 5'
	}
);
$vcf->add_header_line(
	{
		key         => 'FILTER',
		ID          => 'Q3',
		Description =>
		  'Filtered out by filterSNPqual.pl because median quality < 3'
	}
);

$vcf->add_header_line(
	{
		key         => 'INFO',
		ID          => 'MED',
		Number      => 1,
		Type        => 'Integer',
		Description => 'Median read quality at this position'
	}
);
$vcf->add_header_line(
	{
		key         => 'INFO',
		ID          => 'SF',
		Number      => 1,
		Type        => 'Float',
		Description => 'Lower strand frequency (<=0.5) of variant read'
	}
);
$vcf->add_header_line(
	{
		key         => 'INFO',
		ID          => 'VP',
		Number      => 1,
		Type        => 'Float',
		Description => 'Percentage of reads containing the (best) variant'
	}
);

print OUT $vcf->format_header();

while ( $vcfline = $vcf->next_data_hash() ) {

#print "test:".$vcfline->{CHROM}.",".$vcfline->{POS}.",".$vcfline->{ALT}[0]."\n";
	( $medianqual, $medread, $strandfreq, $varpercent ) = &calcMedian(
		$vcfline->{CHROM},  $vcfline->{POS},
		$vcfline->{ALT}[0], $vcfline->{REF}
	);

	my $filter = "";

	if ( $medianqual < 20 ) {
		$filter = "Q20";
	}
	if ( $medianqual < 15 ) {
		$filter = "Q15";
	}
	if ( $medianqual < 10 ) {
		$filter = "Q10";
	}
	if ( $medianqual < 5 ) {
		$filter = "Q5";
	}
	if ( $medianqual < 3 ) {
		$filter = "Q3";
	}
	if ( $filter ne "" ) {
		if($vcfline->{FILTER}[0] && ($vcfline->{FILTER}[0] eq "." || $vcfline->{FILTER}[0] eq "PASS") ){
			delete $vcfline->{FILTER}[0];
		}
		push( @{ $vcfline->{FILTER} }, $filter );
	}

	#	if($filter ne ""){
	#		if($columns[6] eq "."){
	#			$columns[6] = $filter;
	#		}else{
	#			$columns[6] .= ",".$filter;
	#		}
	#	}

	$strandfreq = sprintf( "%.2f", $strandfreq );
	$varpercent = sprintf( "%.2f", $varpercent );

	$vcfline->{INFO}->{MED} = $medread;
	$vcfline->{INFO}->{SF}  = $strandfreq;
	$vcfline->{INFO}->{VP}  = $varpercent;

	#print OUT join("\t",@columns),"\n";
	print OUT $vcf->format_line($vcfline);

}

close OUT;

#close LOW;
my $elapsed = ( time - $^T ) / 60;
$logger->info("finished in $elapsed minutes");

#################################################################################################################
sub calcMedian {
	my @medianqual = ();
	my @basemedian = ();
	my $readpos;
	my $read;
	my $strand;
	my $cap;
	my @quality;
	my $pos = 0;
	my $base;
	my $snvchr    = shift;
	my $snvpos    = shift;
	my $snv       = shift;
	my $ref       = shift;
	my $plStrand  = 0;
	my $minStrand = 0;

	#($snv)=split(/\,/,$snv);

	my $segment =
	  $sam->segment( -seq_id => $snvchr, -start => $snvpos, -end => $snvpos );
	my @alignments = $segment->features;

#open(IN, "$samtools view $outdir/merged.rmdup.bam $snvchr:$snvpos-$snvpos|")|| exit $logger->error("Cannot open samtools view stream on $outdir/merged.rmdup.bam");
	for my $a (@alignments) {

		$readpos = $a->start;
		$read    = $a->query->dna;
		$strand  = $a->strand;

		@quality = $a->qscore;

		my $cigar = $a->cigar_str;
		my $isD = 0;
		($isD,$pos) = &getPosition( $snvpos, $readpos, $cigar, ((length $ref) - (length $snv)) );
		$base = substr( $read, $pos, length $snv );
		my $refbase = substr( $read, $pos, length $ref );
#		if($snvpos == 144923838 ){
#			print $ref.":".$refbase." - ".$snv.":".$base." - ".$cigar." ".$isD." ".(length $ref)."-".(length $snv)." ";
#		}


		if (
			(uc $base eq uc $snv
			&&  length $ref <= length $snv)
				|| $isD#uc $refbase ne uc $ref )
		  )
		{ #check if base(s) in read equal variant base(s); for deletions it should also be different to the reference

#			if($snvpos == 144923838 ){
#				print " - in if";
#			}
			
			if ( $strand > 0 ) {
				$plStrand++;
			}
			else {
				$minStrand++;
			}

			push( @basemedian, $quality[$pos] );
			push( @medianqual, &medianqual( $strand, \@quality, $pos, $read ) );

		}

		my $algnanz = @alignments;

	}    #end IN
	
	if(($plStrand+$minStrand) == 0){ #for some complex positions (e.g. where deletions, insertions and SNV accumulate) there may be no reads that match the reported variant
		
		$logger->debug("Couldn't find a matching read for reported variant at $snvchr:$snvpos!");
		return (0,0,0,0);
	}
	
	my $medianqual = &median( \@medianqual );
	my $basemed    = &median( \@basemedian );
	my $strandfreq = 0;
	if ( $plStrand > $minStrand ) {
		$strandfreq = $minStrand / ( $plStrand + $minStrand );
	}
	else {
		$strandfreq = $plStrand / ( $plStrand + $minStrand );
	}
	my $varpercent = ( $plStrand + $minStrand ) / @alignments;
	return ( $medianqual, $basemed, $strandfreq, $varpercent );

}

sub getPosition {
	my $snvpos  = shift;
	my $readpos = shift;
	my $cigar   = shift;
	my $snvLength = shift;

	if ( $cigar =~ /S/ || $cigar =~ /D/ || $cigar =~ /I/ ) {
		my $bases  = 0;
		my $offset = 0;
		my $diff   = $snvpos - $readpos;
		my $type   = "";
		my $curr   = 0;
		my $cigarPart = "";

		while ( $cigar =~ m/(\d+\D)/g ) {
			$cigarPart .= $1;
			$type = substr( $1, ( length $1 ) - 1, 1 );
			$curr = substr( $1, 0, ( length $1 ) - 1 );
			if ( $type eq "S" || $type eq "I" ) {
				$bases  += $curr;
				$offset += $curr;
			}
			elsif ( $type eq "D" ) {
				$offset -= $curr;
			}
			elsif ( $type eq "M" ) {
				$bases += $curr;
			}

			if ( $bases >= ( $diff + $offset ) ) {
				my $d = 0;
				my $regex = $cigarPart.$snvLength."D";
#				if($snvpos == 144923838 ){
#					print "\n - $cigarPart - $regex ";
#				}
				
				if( $cigar =~ m/$regex/){
					$d = 1;
				}
				return ($d, ($diff + $offset ));

			}
		}
		print "no return: $cigar,$snvpos, $readpos, $diff, $offset, $bases\n";
	}
	else {
		return (0,( $snvpos - $readpos ));
	}
}

##############################################################################
sub medianqual {

	#print "\nmedianqual\n";
	my $strand     = shift;
	my $ref        = shift;
	my @quality    = @{$ref};
	my $pos        = shift;
	my $read       = shift;
	my $i          = 0;
	my $qual       = 0;
	my @qual       = ();
	my $medianqual = 0;

	if ($strand > 0) {
		for ($i=$pos;$i<&min(length($read),$pos+$window);$i++) {
			push(@qual,$quality[$i]);
		}
	}
	else {
		for ($i=&max(0,$pos-$window);$i<=$pos;$i++) {
			push(@qual,$quality[$i]);
		}
	}
	

	$medianqual = &median( \@qual );
	return ($medianqual);
}
##############################################################################
sub median {

	#print "\nmedian\n";
	my $ref = shift;
	my $median;
	my @array = ();
	my $n     = 0;
	my $t     = 0;

	$n = ( @{$ref} );
	@{$ref} = sort { $a <=> $b } @{$ref};

	#print "median: @{$ref}\n";
	if ( $n % 2 != 0 ) {
		$median = ${$ref}[ int( $n / 2 ) ];
	}
	else {
		$median =
		  ( ${$ref}[ int( $n / 2 - 1 ) ] + ${$ref}[ int( $n / 2 ) ] ) / 2;
	}

	return ($median);
}

##############################################################################
sub min ($$) { $_[$_[0] > $_[1]] }
sub max ($$) { $_[$_[0] < $_[1]] }
