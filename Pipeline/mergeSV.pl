#!/usr/bin/perl 

use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);
use Pod::Usage;
use DBI;
use Data::Dumper;

my $prog_path = dirname( abs_path($0) );
require $prog_path . "/Utilities.pm";
require $prog_path . "/BEDRecord.pm";
require $prog_path . "/SequenceDictionary.pm";


my $settings = "default";

my $outfile      = "";

my $cnvnator     = "";
my $lumpy        = "";
my $pindel       = "";
my $breakdancer  = "";
my $manta		 = "";

my $overlap      = 0.8;
my $pindelLength = 50;

my $help         = 0;
my $man          = 0;
my $logfile      = "SCREEN";
my $loglevel     = "INFO";

GetOptions(
	"o=s"  => \$outfile,
	"c=s"  => \$cnvnator,
	"l=s"  => \$lumpy,
	"p=s"  => \$pindel,
	"pl=s" => \$pindelLength,
	"b=s"  => \$breakdancer,
	"m=s"  => \$manta,
	"v=s"  => \$overlap,
	"se=s" => \$settings,
	"lf=s" => \$logfile,
	"ll=s" => \$loglevel,
	"h"    => \$help,
	"man"  => \$man
);

pod2usage( { -exitval => 0, -verbose => 1 } ) if $help;
pod2usage( { -exitval => 0, -verbose => 2 } ) if $man;
pod2usage( { -exitval => 1, -verbose => 1 } ) if $outfile eq "" || ( $cnvnator eq "" && $lumpy eq "" && $pindel eq  "" && $breakdancer eq "" && $manta eq "" );

my $params = Utilities::getParams();
Utilities::initLogger( $logfile, $loglevel );		#initialize logging system
my $logger = Utilities::getLogger();

require $params->{programs}->{vcftools}->{pm};    	#VCFTools library to parse VCF file


#this is the top level hash that will contain hashes for every variant type which in turn will contain arrays for every chromosome
my %variants;


#read cnvnator
if ($cnvnator ne "") {
	
	$logger->info("Adding CNVnator variants...");
	open IN, "$cnvnator" or exit $logger->error("Can't open $cnvnator!");
	while(<IN>) {
		chomp;
		my @columns = split("\t");
		my $svtype  = "DEL";
		$svtype     = "DUP" if $columns[0] eq "duplication";
		my ($chrom,$rest) = split(":",$columns[1]);
		my ($start,$end)  = split("-",$rest);
		
		# CNVnator returns always positive length, make it negative for deletions		
		my $svlen = abs($columns[2]) * ( $svtype eq "DEL" ? -1 : 1 );
		
		my @values;
		push(@values,"CNDOSAGE=$columns[3]");
		push(@values,"CNSCORE1=$columns[4]");
		push(@values,"CNSCORE2=$columns[5]");
		push(@values,"CNSCORE3=$columns[6]");
		push(@values,"CNSCORE4=$columns[7]");
		push(@values,"CNUNIQUE=$columns[8]");

		&mergeVariant($chrom,$start,$end,$svtype,$svlen,"cnvnator",\@values);
	}
	
	close IN;
	
}


#read lumpy
if ($lumpy ne "") {
	
	$logger->info("Adding LUMPY-SV variants...");
	my $vcf = Vcf->new( file => $lumpy );
	$vcf->parse_header();
	while ( my $vcfline = $vcf->next_data_hash() ) {
		
		next unless ($vcfline->{INFO}->{SVTYPE} eq "DEL" or $vcfline->{INFO}->{SVTYPE} eq "INS" or $vcfline->{INFO}->{SVTYPE} eq "DUP" or $vcfline->{INFO}->{SVTYPE} eq "INV");
		
		my @values;
		push(@values,"LPPE=$vcfline->{INFO}->{PE}");
		push(@values,"LPSR=$vcfline->{INFO}->{SR}");
		
		&mergeVariant($vcfline->{CHROM},$vcfline->{POS},$vcfline->{INFO}->{END},$vcfline->{INFO}->{SVTYPE},$vcfline->{INFO}->{SVLEN},"lumpy-sv",\@values );
	}
	
	$vcf->close();
	
}


#read pindel
if ($pindel ne "") {
	
	$logger->info("Adding Pindel variants...");
	my $vcf = Vcf->new( file => $pindel );
	$vcf->parse_header();
	my @samples = $vcf->get_samples();
	my $sample  = $samples[0];
	
	while ( my $vcfline = $vcf->next_data_hash() ) {
		
		next unless ($vcfline->{INFO}->{SVTYPE} eq "DEL" or $vcfline->{INFO}->{SVTYPE} eq "INS" or $vcfline->{INFO}->{SVTYPE} eq "DUP" or $vcfline->{INFO}->{SVTYPE} eq "INV");
		next if abs($vcfline->{INFO}->{SVLEN}) < $pindelLength;
		
		unless($vcfline->{gtypes}->{$sample}->{GT}){
			print STDERR Dumper($vcfline);
		}
		
		my $alleles;
		if ( $vcfline->{gtypes}->{$sample}->{GT} =~ /(\d)\/(\d)/ ) {
			$alleles = $1 + $2;
		}
		
		my ($ref,$alt) = split(",",$vcfline->{gtypes}->{$sample}->{AD});
		
		
		my @values;
		push(@values,"PIDP=".($ref+$alt));
		push(@values,"PIPERCENTVAR=".($alt/($ref+$alt)));
		push(@values,"PIALLELES=$alleles");
			
		&mergeVariant($vcfline->{CHROM},$vcfline->{POS},$vcfline->{INFO}->{END},$vcfline->{INFO}->{SVTYPE},$vcfline->{INFO}->{SVLEN},"pindel",\@values );
	}
	
	$vcf->close();
	
}


#read breakdancer
if ($breakdancer ne "") {
	
	$logger->info("Adding breakdancer variants...");
	my $vcf = Vcf->new( file => $breakdancer );
	$vcf->parse_header();
	while ( my $vcfline = $vcf->next_data_hash() ) {
		
		next unless ($vcfline->{INFO}->{SVTYPE} eq "DEL" or $vcfline->{INFO}->{SVTYPE} eq "INS" or $vcfline->{INFO}->{SVTYPE} eq "DUP" or $vcfline->{INFO}->{SVTYPE} eq "INV");
				
		# SVLEN Patch for old breakdancer VCFs run with buggy breakdancer2VCF [DEL] must have SVLEN<0 , the others >0 )
		$vcfline->{INFO}->{SVLEN} = abs( $vcfline->{INFO}->{SVLEN} ) * ( $vcfline->{INFO}->{SVTYPE} eq "DEL" ? -1 : 1 ); 
				
		my @values;
		push(@values,"BDDP=".$vcfline->{INFO}->{DP});
		push(@values,"BDOR1=".$vcfline->{INFO}->{OR1});
		push(@values,"BDOR2=".$vcfline->{INFO}->{OR2});
		
		&mergeVariant($vcfline->{CHROM},$vcfline->{POS},$vcfline->{INFO}->{END},$vcfline->{INFO}->{SVTYPE},$vcfline->{INFO}->{SVLEN},"breakdancer",\@values );
	}
	
	$vcf->close();
	
}


#read manta
if ($manta ne "") {
	
	$logger->info("Adding manta variants...");
	if($manta =~ /\.gz$/){
		my $bgzip = $params->{programs}->{bgzip}->{path};
		my $mantagz = $manta;
		$manta =~ s/\.gz$//;
		system("$bgzip -cd $mantagz > $manta");
	}
	my $vcf = Vcf->new( file => $manta );
	$vcf->parse_header();
	my @samples = $vcf->get_samples();
	while ( my $vcfline = $vcf->next_data_hash() ) {
		
		next unless ($vcfline->{INFO}->{SVTYPE} eq "DEL" or $vcfline->{INFO}->{SVTYPE} eq "INS" or $vcfline->{INFO}->{SVTYPE} eq "DUP" or $vcfline->{INFO}->{SVTYPE} eq "INV");
		next unless $vcfline->{FILTER}[0] eq "PASS";
		
		my $alleles = 0;
		if ( $vcfline->{gtypes}->{$samples[0]}->{GT} =~ /(\d)\/(\d)/ ) {
			$alleles = $1 + $2;
		}
		
		my @values;
		push(@values,"MTGT=".$alleles);
		push(@values,"MTGQ=".$vcfline->{gtypes}->{$samples[0]}->{GQ});

		if (!($vcfline->{INFO}->{SVLEN})) {		
			if($vcfline->{INFO}->{SVTYPE} eq "INS"){
				if($vcfline->{INFO}->{SVINSLEN}){
					$vcfline->{INFO}->{SVLEN} = $vcfline->{INFO}->{SVINSLEN};
				}else{
					my $leftlength = 0;
					$leftlength = length($vcfline->{INFO}->{LEFT_SVINSSEQ}) if $vcfline->{INFO}->{LEFT_SVINSSEQ};
					
					my $rightlength = 0;
					$rightlength = length($vcfline->{INFO}->{RIGHT_SVINSSEQ}) if $vcfline->{INFO}->{RIGHT_SVINSSEQ};
					
					$leftlength++ if ($leftlength + $rightlength) == 0;
					
					
					$vcfline->{INFO}->{SVLEN} = $leftlength + $rightlength;
				}
			}
		}
		&mergeVariant($vcfline->{CHROM},$vcfline->{POS},$vcfline->{INFO}->{END},$vcfline->{INFO}->{SVTYPE},$vcfline->{INFO}->{SVLEN},"manta",\@values );
	}
	
	$vcf->close();
	
}


#print merged SVS
open OUT, "| perl $prog_path/vcfsorter.pl -o $outfile -se $settings -lf $logfile -ll $loglevel" or exit $logger->error("Can't open $outfile!");
$logger->info("Sorting and printing variants...");

#print header
print OUT "##fileformat=VCFv4.2
##source=mergeSV.pl
##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">
##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">
##INFO=<ID=CNSCORE1,Number=1,Type=Float,Description=\"CNVnator score 1\">
##INFO=<ID=CNSCORE2,Number=1,Type=Float,Description=\"CNVnator score 2\">
##INFO=<ID=CNSCORE3,Number=1,Type=Float,Description=\"CNVnator score 3\">
##INFO=<ID=CNSCORE4,Number=1,Type=Float,Description=\"CNVnator score 4\">
##INFO=<ID=CNUNIQUE,Number=1,Type=Float,Description=\"CNVnator unique\">
##INFO=<ID=LPPE,Number=1,Type=Integer,Description=\"LUMPY-SV: Number of paired-end reads supporting the variant across all samples\">
##INFO=<ID=LPSR,Number=1,Type=Integer,Description=\"LUMPY-SV: Number of split reads supporting the variant across all samples\">
##INFO=<ID=PIDP,Number=1,Type=Integer,Description=\"Pindel: Read depth\">
##INFO=<ID=PIPERCENTVAR,Number=1,Type=Float,Description=\"Pindel: Percent of reads showing the variant\">
##INFO=<ID=PIALLELES,Number=1,Type=Integer,Description=\"Pindel: Number of variant alleles.\">
##INFO=<ID=BDDP,Number=1,Type=Integer,Description=\"Breakdancer: Read depth\">
##INFO=<ID=BDOR1,Number=1,Type=String,Description=\"Breakdancer: Orientation of reads at breakpoint 1\">
##INFO=<ID=BDOR2,Number=1,Type=String,Description=\"Breakdancer: Orientation of reads at breakpoint 2\">
##INFO=<ID=MTGT,Number=1,Type=Integer,Description=\"Manta: alleles of variant: 1 - heterozygous, 2 - homozygous\">
##INFO=<ID=MTGQ,Number=1,Type=Integer,Description=\"Manta: genotype quality of variant\">
##ALT=<ID=DEL,Description=\"Deletion\">
##ALT=<ID=DUP,Description=\"Duplication\">
##ALT=<ID=INV,Description=\"Inversion\">
##ALT=<ID=INS,Description=\"Insertion of novel sequence\">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
";

foreach my $svtype(keys %variants){
	foreach my $chrom (keys %{$variants{$svtype}}){
		foreach my $sv (@{$variants{$svtype}->{$chrom}}){
			
			#print Dumper($sv);
			print OUT "$chrom\t".$sv->startpos()."\t.\tN\t<$svtype>\t.\t.\tSVLEN=".$sv->overlap().";SVTYPE=$svtype;END=".$sv->endpos().";CALLER=".$sv->name().";";
			print OUT join(";",@{$sv->rest()})."\n";
		}
	}
}
close OUT;

##################### merge variants ########################
sub mergeVariant {
	
	my $chrom   = shift;
	my $start   = shift;
	my $end     = shift;
	my $svtype  = shift;
	my $svlen   = shift;
	my $caller  = shift;
	my $pointer = shift;
	
	#create BEDRecord for current
	my $currSV = BEDRecord->new($chrom,$start,$end,$caller);
	$currSV->rest($pointer);
	$currSV->overlap($svlen);			#store the svlen in the overlap field of the BED entry
	
	if(defined $variants{$svtype}){			#check if variants for the current type are already available
		
		my @overlapInd = $currSV->getOverlappingIndices($variants{$svtype},$overlap);			#get list of reciprocal overlapping indices
		foreach(@overlapInd){
			my $ovSV = $variants{$svtype}->{$chrom}->[$_];
			$currSV->startpos($currSV->startpos()+$ovSV->startpos());				#merge overlapping SVs
			$currSV->endpos($currSV->endpos()+$ovSV->endpos());
			$currSV->overlap($currSV->overlap()+$ovSV->overlap());
			next if $ovSV->name() =~ /$caller/;						#skip entries of same variant caller
			$currSV->name($currSV->name().",".$ovSV->name());
			my @values = @{$currSV->rest()};
			push(@values,@{$ovSV->rest()});
			$currSV->rest(\@values);
		}
		$currSV->startpos(int ( $currSV->startpos() / (@overlapInd+1) ) );			#normalize startpos and endpos
		$currSV->endpos(int ( $currSV->endpos() / (@overlapInd+1) ) );
		$currSV->overlap(int ( $currSV->overlap() / (@overlapInd+1) ) );
		
		BEDRecord::removeBEDRecords($variants{$svtype}->{$chrom},\@overlapInd);		#remove old entries
		
	}else{
		my %empty;
		$variants{$svtype} = \%empty;
	}
	
	#insert current entry into hash
	$currSV->insertBEDRecord($variants{$svtype});
	
}

=head1 NAME

mergeSV.pl

=head1 SYNOPSIS

 mergeSV.pl -c cnvnator.out -l lumpy.vcf -p pindel.vcf -b breakdancer.vcf -o merged.vcf 

=head1 DESCRIPTION

This script merges structural variants (SVs) from (currently) for different callers. If two SVs of the same type
overlap to a large extend (defined by -v) the two variants are merged and the mean of their start and end coordinates is used
for the merged variant.

NOTE: There will be no genotype 

=head1 OPTIONS

 -o	<outfile.vcf> output VCF file; REQUIRED
 -c	<cnvnator.out> cnvnator output file
 -l	<lumpy.vcf> lumpy-sv output file
 -p	<pindel.vcf> pindel output file
 -pl minimum length of a pindel variant to merge; default: 50
 -b	<breakdancer.vcf> breakdancer output file
 -m	<manta.vcf> manta (diploid) ouput file
 -se	name of the settings in the current.config.xml file that holds path to reference genome, 
 	to the annotation file and to possible additional annotation files; use default settings if nothing is given
 -v	required reciprocal overlap to join two variants of the same sample; default: 0.8
 -lf	log file; default: print to screen
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page

=head1 AUTHOR

Thomas Wieland

=cut
