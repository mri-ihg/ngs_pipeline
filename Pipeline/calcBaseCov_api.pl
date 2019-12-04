#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use DBI;
use POSIX;
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);
use Bio::DB::Sam;
use Pod::Usage;


my $prog_path = dirname( abs_path($0) );
require $prog_path . "/Utilities.pm";


my $infile  = "";
my $outfile = "";

#my $refFa = "";
my $settings         = "";
my $bam              = "";
my $COV              = 5;
my $logfile          = "SCREEN";
my $loglevel         = "INFO";
my $bed              = 0;
my $sample           = "coverage";
my $zeros            = 0;
my $regFile          = "";
my $maxCov           = 0;
my $calcGC           = 0;
my $wigFile          = "";
my $covPerTargetFile = "";
my $covPerTargetSeq  = "";
my $calcMQ           = 0;
my $targetType       = "file";
my $maxDepth		 = 99999;
my $coveredTarget    = "";

#program flags
my $help = 0;
my $man  = 0;

GetOptions(
	"t=s"  => \$infile,
	"tt=s" => \$targetType,
	"se=s" => \$settings,
	"b=s"  => \$bam,
	"c=s"  => \$COV,
	"gc"   => \$calcGC,
	"m=s"  => \$maxCov,
	"lf=s" => \$logfile,
	"ll=s" => \$loglevel,
	"r=s"  => \$regFile,
	"bed"  => \$bed,
	"q"    => \$calcMQ,
	"o=s"  => \$outfile,
	"w=s"  => \$wigFile,
	"p=s"  => \$covPerTargetFile,
	"ps=s" => \$covPerTargetSeq,
	"ct=s" => \$coveredTarget,
	"z"    => \$zeros,
	"man"  => \$man,
	"h"    => \$help,
	"d=s"  => \$maxDepth
);


pod2usage( {-exitval => 0  ,-verbose => 1} ) if $help;
pod2usage( {-exitval => 0  ,-verbose => 2} ) if $man;
pod2usage( {-exitval => 1  ,-verbose => 1} ) if $bam eq "" || $outfile eq "";

Utilities::initLogger( $logfile, $loglevel );
my $logger = Utilities::getLogger();

my $params = Utilities::getParams();
my $ref      = $params->{settings}->{$settings}->{reference};
my $samtools = $params->{programs}->{samtools}->{path};



#get sample name from BAM file
open BAM, "$samtools view -H $bam |" || exit $logger->error("Could not open $samtools view -H $bam");
while(<BAM>){
	if($_ =~ /^\@/){
		if($_ =~ /^\@RG/){
			my @columns = split();
			foreach my $field(@columns){
				if($field =~ /SM:(.+)/){
					$sample = $1;
					last;
				}
			}
		}
	}
}
close BAM;


#open Bio-Samtools object
my $sam = Bio::DB::Sam->new(
	-fasta => $ref,
	-bam   => $bam,
	-expand_flags  => 1
);
#Original command in Bio::DB::Sam, gives erratic results
#$sam->max_pileup_cnt([$maxDepth]);		#set max depth
$sam->max_pileup_cnt($maxDepth);		#set max depth


####### open output files and print headers ########
my $outdir = dirname($outfile);
if ($zeros) {
	open( ZEROS, ">$outdir/pileup.zeros" )
	  || exit $logger->error("Could not open pileup.zeros")
	  ; # for deletion_hom --> report all positions with zero coverage (if traget file is specified)
}

if($coveredTarget ne ""){
	my $bedtools = $params->{programs}->{bedtools}->{path};
	open COVT, "| $bedtools/mergeBed -i stdin > $coveredTarget" or  exit $logger->error("Could not open pipe: | $bedtools/mergeBed -i stdin > $coveredTarget");
}

open( OUT, ">$outfile" ) || exit $logger->error("Could not open $outfile");

$wigFile = "$outdir/coverage.wig" if $wigFile eq "";
open( WIG, ">$wigFile" )
  || exit $logger->error("Could not open $wigFile");

$covPerTargetFile = "$outdir/coverage_per_target.txt"
  if $covPerTargetFile eq "";
open( AVG, ">$covPerTargetFile" )
  || exit $logger->error("Could not open $covPerTargetFile");

$covPerTargetSeq = "$outdir/coverage_per_target.seq"
  if $covPerTargetSeq eq "";
open( AVGSEG, ">$covPerTargetSeq" )
  || exit $logger->error("Could not open $covPerTargetSeq");

if ( $bed && $regFile ne "" ) {
	open( REG, ">$regFile" )
	  || exit $logger->error("Could not open $regFile!");
	if ($calcMQ) {
		print REG
"#Region\tRegion chr\tRegion starts\tRegion ends\tCoverage >=20 (\%)\tavg. coverage\tavg. mapping quality\tsingle coverage >=20 (\%)\tsingle avg. coverages\tsingle avg. mapping qualities\n";
	}
	else {
		print REG
"#Region\tRegion chr\tRegion starts\tRegion ends\tCoverage >=20 (\%)\tavg. coverage\tsingle coverage >=20 (\%)\tsingle avg. coverages\n";
	}
}
print WIG
"track type=wiggle_0 name=\"$sample\" description=\"Coverage of $sample\" visibility=full color=0,0,0 altColor=255,0,0\n";
print AVG "Target	total_coverage	average_coverage";
if ($calcGC) {
	print AVG "\tGC Content (%)";
}
print AVG "\n";

print AVGSEG "#type=GENE_EXPRESSION
#track  graphtype=points color=0,0,0 altColor=0,0,0 viewLimits=-1.2:75 maxHeightPixels=120:120:120
feature	chrom	start	end	value
";



my $uncovered   = 0;
my $coveredOne  = 0;
my $coveredCOV  = 0;
my $covered4    = 0;
my $covered8    = 0;
my $coveredSum  = 0;
my $targetCount = 0;
my @median;
my ( $chr,   $pos );
my ( $start, $end );
my $region;
my $old_region;
my $old_region_chr;

my @percCoverages;

my $autoCov   = 0;
my $autoRange = 0;

my $regionCov    = 0;
my $regionRange  = 0;
my $region20x    = 0;
my $regionMQual  = 0;
my $regionReads  = 0;
my $singleCov    = "";
my $single20x    = "";
my $singleMQual  = "";
my $singleStarts = "";
my $singleEnds   = "";
my $currMQ       = 0;



if ( $infile ne "" && $targetType ne "chr") {

	if ( $targetType eq "db" ) {
		my $dbh = Utilities::connectCoreDB();
		my $sth =
		  $dbh->prepare(
"select chrom,exonStarts,exonEnds,name,idtranscript from $infile"
		  )
		  || $logger->error("Can't prepare statement: $DBI::errstr");
		$sth->execute()
		  || $logger->error("Can't execute statement: $DBI::errstr");

		my ( $exonStarts, $exonEnds, $idtranscript );
		while ( ( $chr, $exonStarts, $exonEnds, $region, $idtranscript ) =
			$sth->fetchrow_array() )
		{    #get regions from database
			$exonStarts =~ s/,$//;
			$exonEnds   =~ s/,$//;
			$region .= ",$idtranscript";

			if (   $regFile ne ""
				&& $region
				&& $old_region
				&& $region ne $old_region )
			{ #calculate the average coverage for more than one regions, e.g. for all exons of a gene
				&calcRegionCov();
			}

			my @exonStarts = split( ",", $exonStarts );
			my @exonEnds   = split( ",", $exonEnds );

			for ( my $i = 0 ; $i < @exonStarts ; $i++ ) {
				$currMQ = &calcMapQForRegion( $chr, $exonStarts[$i], $exonEnds[$i] )
				  if $calcMQ;
				&calcCovForRegion( $chr, $exonStarts[$i], $exonEnds[$i], $currMQ );

			}
			$old_region     = $region;
			$old_region_chr = $chr;
		}

	}
	else {

		$infile = $params->{settings}->{$settings}->{$infile}
		  if $targetType eq "settings";

		open( IN, "$infile" ) || exit $logger->error("Could not open $infile");
		while (<IN>) {
			chomp;
			if ($bed) {
				( $chr, $start, $end, $region ) = split();
				if ($region) {
					$region =~ s/\s/_/g;
				}
			}
			else {
				( $chr, $pos ) = split(":");
				( $start, $end ) = split( "-", $pos );
			}

			if (   $regFile ne ""
				&& $region
				&& $old_region
				&& $region ne $old_region )
			{ #calculate the average coverage for more than one regions, e.g. for all exons of a gene

				&calcRegionCov();
			}
			$currMQ = &calcMapQForRegion( $chr, $start, $end ) if $calcMQ;
			&calcCovForRegion( $chr, $start, $end, $currMQ );

			$old_region     = $region;
			$old_region_chr = $chr;
		}
	}

	if ( $regFile ne "" && $old_region ) {
		&calcRegionCov();
	}

}
else {
	&calcWholeGenomeCov($infile);
}

my $median = 0;
if (@median) {
	my $sum = 0;
	for(my $i = 0; $i < @median; $i++){			#calc median from historgram
		next unless $median[$i];
		$sum += $median[$i];
		if($sum >= ($targetCount/2) ){
			$median = $i;
			last;
		}
	}
}

my $uncoveredPer  = ( $uncovered * 100 ) / $targetCount;
my $coveredOnePer = ( $coveredOne * 100 ) / $targetCount;
my $coveredCOVPer = ( $coveredCOV * 100 ) / $targetCount;
my $covered4Per   = ( $covered4 * 100 ) / $targetCount;
my $covered8Per   = ( $covered8 * 100 ) / $targetCount;
my $coveredAvg    = $coveredSum / $targetCount;
my $autoAvg       = 0;
unless ( $autoRange == 0 ) {
	$autoAvg = $autoCov / $autoRange;
}

$uncoveredPer  = sprintf( "%.2f", $uncoveredPer );
$coveredOnePer = sprintf( "%.2f", $coveredOnePer );
$coveredCOVPer = sprintf( "%.2f", $coveredCOVPer );
$covered4Per   = sprintf( "%.2f", $covered4Per );
$covered8Per   = sprintf( "%.2f", $covered8Per );
$coveredAvg    = sprintf( "%.2f", $coveredAvg );
$autoAvg       = sprintf( "%.2f", $autoAvg );

my $variance       = 0;
my $medianVariance = 0;
my $n              = 0;


for(my $i = 0; $i < @median; $i++){ #calc variance from histogram
	next unless $median[$i];
	$variance       += (( $i - $coveredAvg )**2)*$median[$i];
	$medianVariance += (( $i - $median )**2)*$median[$i];
	$n 				+= $median[$i];
}

#foreach (@median) {
#	$variance       += ( $_ - $coveredAvg )**2;
#	$medianVariance += ( $_ - $median )**2;
#	$n++;
#}

my $stdDev = sqrt( $variance /       ( $n - 1 ) );
my $medDev = sqrt( $medianVariance / ( $n - 1 ) );

$stdDev = sprintf( "%.2f", $stdDev );
$medDev = sprintf( "%.2f", $medDev );

#print coverage distribution file, if required
if ( $maxCov != 0 ) {
	my $outfilename = basename($outfile);
	open DIST,
	  ">$outdir/distribution.$outfilename"
	  || exit $logger->error(
		"Could not open $outdir/distribution.$outfilename!");
	print DIST "#Depth of Coverage\t\% covered\n";
	for ( my $j = 1 ; $j <= $maxCov ; $j++ ) {
		print DIST $j;
		if ( $percCoverages[$j] ) {
			my $perc =
			  sprintf( "%.2f",
				( ( $percCoverages[$j] / $targetCount ) * 100 ) );
			print DIST "\t$perc\n";
		}
		else {
			print DIST "\t0\n";
		}
	}
	close DIST;
}

print OUT "\n$targetCount bases targeted\n";
print OUT "$uncovered ($uncoveredPer\%) bases not covered\n";
print OUT "$coveredOne ($coveredOnePer\%) bases covered at least once\n";
print OUT "$covered4 ($covered4Per\%) bases covered at 4 X or more\n";
print OUT "$covered8 ($covered8Per\%) bases covered at 8 X or more\n";
print OUT "$coveredCOV ($coveredCOVPer\%) bases covered at $COV X or more\n";
print OUT "$coveredAvg average coverage ($stdDev standard deviation)\n";
print OUT "$median median coverage ($medDev)\n";
print OUT "$autoAvg average coverage of autosomes\n\n";

my $elapsed = ( time - $^T ) / 60;
print OUT "finished in $elapsed minutes\n";

close OUT;
close IN;
close ZEROS;
close COVT if $coveredTarget ne "";
close WIG;
close AVG;
close AVGSEG;
system("gzip -f $wigFile");

##############################################################################
sub calcRegionCov {

	my $avgCov = 0;
	$avgCov = sprintf( "%.2f", ( $regionCov / $regionRange ) )
	  if $regionRange > 0;
	my $perc20x = 0;
	$perc20x = sprintf( "%.2f", ( ( $region20x / $regionRange ) * 100 ) )
	  if $regionRange > 0;

	$singleCov    =~ s/,$//;
	$single20x    =~ s/,$//;
	$singleMQual  =~ s/,$//;
	$singleStarts =~ s/,$//;
	$singleEnds   =~ s/,$//;
	if ($calcMQ) {
		my $mq = 0;
		$mq = sprintf( "%.2f", ( $regionMQual / $regionReads ) )
		  if $regionReads > 0;
		print REG
"$old_region\t$old_region_chr\t$singleStarts\t$singleEnds\t$perc20x\t$avgCov\t$mq\t$single20x\t$singleCov\t$singleMQual\n";
	}
	else {
		print REG
"$old_region\t$old_region_chr\t$singleStarts\t$singleEnds\t$perc20x\t$avgCov\t$single20x\t$singleCov\n";
	}

	$regionCov    = 0;
	$regionRange  = 0;
	$region20x    = 0;
	$singleCov    = "";
	$single20x    = "";
	$singleMQual  = "";
	$singleStarts = "";
	$singleEnds   = "";
	$regionMQual  = 0;
	$regionReads  = 0;
}

##############################################################################
sub calcMapQForRegion {
	my $chr   = shift;
	my $start = shift;
	my $end   = shift;

	my @alignments = $sam->get_features_by_location(
		-seq_id => $chr,
		-start  => $start,
		-end    => $end,
		-flags=>{DUPLICATE=>0}
	);

	unless(@alignments){
		if($regFile ne ""){
			$singleMQual .= "0.00,";
		}
		return 0;
	}

	my $quals = 0;

	foreach my $align (@alignments) {
		$quals += $align->qual;
	}

	if ( $regFile ne "" ) {
		if(@alignments > 0){
			$singleMQual .= sprintf( "%.2f", ( $quals / @alignments ) ) . ",";
		}else{
			$singleMQual .= "0,";
		}
		
		$regionMQual += $quals;
		$regionReads += @alignments;
	}

	my $retval = 0;
	$retval    = ( $quals / @alignments ) if @alignments>0;

	undef @alignments;
	return $retval;

}

##############################################################################
sub calcCovForRegion {
	my $chr   = shift;
	my $start = shift;
	my $end   = shift;
	my $currMQ= shift;

	my $coverage;


	($coverage) = $sam->features(
		-type   => 'coverage',
		-seq_id => $chr,
		-start  => $start,
		-end    => $end,
		-flags=>{DUPLICATE=>0}
	);
	
	if ( $regFile ne "" ) {
		$singleStarts .= $start . ",";
		$singleEnds   .= $end . ",";
	}

	unless($coverage){
		if ( $regFile ne "" ) {
			$singleCov .= "0.00,";
			$single20x .= "0.00,";
		}
		return; #skip region if $coverage is not defined --> may be region that is not in ref genome
	}
	my @cov = $coverage->coverage;

	my $range    = @cov;
	my $wigstart = $start + 1;

	print WIG "\nfixedStep chrom=$chr start=$wigstart step=1\n";

	my $summedCov      = 0;
	my $single20xCount = 0;

	if ( $range == 0 ) {    #if no reads overlap this region @cov would be empty
		$range = $end - $start + 1;
		@cov[ 0 .. ( $range - 1 ) ] = 0;
	}

	for ( my $i = 0 ; $i < $range ; $i++ ) {
		print WIG "$cov[$i]\n";

		if ( $chr =~ m/^(c|C)hr\d+$/ ) {    # only for chromosomes 1-22
			$autoCov += $cov[$i];
		}
		if ( $cov[$i] == 0 ) {
			if ($zeros) {
				$pos = $start + $i;
				print ZEROS "$chr\t$pos\tn\t0\n";
			}

			$uncovered++;
		}
		else {
			if($coveredTarget ne ""){
				$pos = $start + $i;
				print COVT "$chr\t".($pos-1)."\t".$pos."\n";
			}
			#push( @median, $cov[$i] );
			if($median[$cov[$i]]){
				$median[$cov[$i]]++;
			}else{
				$median[$cov[$i]] = 1;
			}

			$coveredOne++;
			$coveredSum += $cov[$i];
			$summedCov  += $cov[$i];
			$regionCov  += $cov[$i];

			if ( $cov[$i] >= $COV ) { $coveredCOV++; }
			if ( $cov[$i] >= 4 )    { $covered4++; }
			if ( $cov[$i] >= 8 )    { $covered8++; }
			if ( $cov[$i] >= 20 )   {
				$region20x++;
				$single20xCount++;
			}
			if ( $maxCov != 0 ) {
				for ( my $j = 1 ; ( $j <= $cov[$i] && $j <= $maxCov ) ; $j++ ) {
					if ( $percCoverages[$j] ) {
						$percCoverages[$j]++;
					}
					else {
						$percCoverages[$j] = 1;
					}
				}
			}
		}
	}
	my $avgCov = $summedCov / $range;
	$targetCount += $range;
	$regionRange += $range;
	if ( $chr =~ m/^(c|C)hr\d+$/ ) {    # only for chromosomes 1-22
		$autoRange += $range;
	}
	$avgCov = sprintf( "%.2f", $avgCov );
	print AVG "$chr:$start-$end\t$summedCov\t$avgCov";

	if ( $regFile ne "" ) {
		$singleCov .= $avgCov . ",";
		$single20x .=
		  sprintf( "%.2f", ( ( $single20xCount / $range ) * 100 ) ) . ",";
	}

	if ($calcGC) {
		my $seq =
		  uc $sam->seq( $chr, $start, $end )
		  ;    #calculate GC content of this target
		my $all = 0;
		my $gc  = 0;
		while ( $seq =~ /(.)/g ) {    # for each base
			$all++;
			$gc++ if ( $1 eq "G" || $1 eq "C" );
		}
		my $gcCont = sprintf( "%.2f", ( ( $gc / $all ) * 100 ) );
		print AVG "\t$gcCont";
	}
	print AVG "\t" . sprintf( "%.2f", $currMQ ) if $calcMQ;
	print AVG "\n";
	if($calcMQ == 0 || $currMQ >= 45){
		print AVGSEG "$sample\t$chr\t$start\t$end\t".($avgCov+2)."\n";
	}
	
	undef $coverage;
	undef @cov;
	
}

##############################################################################
#Calculating coverage for the whole genome takes too long with the api
sub calcWholeGenomeCov() {
	my $targetChrom = shift;
	my %chromLength;
	my @chroms = ();

	
	
	#get chromosome information from sequence dictionary
	my $dict = $ref;
	$dict    =~ s/fa$/dict/;
	$dict    =~ s/fasta$/dict/;
	
	open DICT, "grep \"\@SQ\" $dict |" || exit $logger->error("Can't open $dict!");
	while (<DICT>) {
		chomp;
		if ( $_ =~ /\@SQ/ ) {
			$_ =~ /SN:(.*)\sLN/;
			my $chr = $1;
			next if $targetChrom ne "" && $chr ne $targetChrom;
			$_ =~ /LN:(\d*)/;
			$chromLength{$chr} = $1;
			push( @chroms, $chr );
		}
	}
	close DICT;
	
	
	

	foreach my $currChr (@chroms) {
		open( PILEUP, "$samtools mpileup -r $currChr $bam|" )
		  ;    #check all chromosomes found in header
		my $posCounter = 1;
		my $summedCov  = 0;

		print "Processing $currChr...\n";

		print WIG "\nfixedStep chrom=$currChr start=1 step=1\n";

		while (<PILEUP>) {
			my ( $dummy, $pos, $base, $cov ) = split();
			while ( $posCounter < $pos ) {
				print WIG "0\n";
				if ($zeros) {
					print ZEROS "$currChr\t$posCounter\tn\t0\n";
				}

				$uncovered++;
				$posCounter++;
			}

			if ( $currChr =~ m/^(c|C)hr\d+$/ ) {    # only for chromosomes 1-22
				$autoCov += $cov;
			}

			print WIG "$cov\n";
			$coveredOne++;
			$coveredSum += $cov;
			$summedCov  += $cov;

			if ( $cov >= $COV ) { $coveredCOV++; }
			if ( $cov >= 4 )    { $covered4++; }
			if ( $cov >= 8 )    { $covered8++; }
			
			if ( $maxCov != 0 ) {
				for ( my $j = 1 ; ( $j <= $cov && $j <= $maxCov ) ; $j++ ) {
					if ( $percCoverages[$j] ) {
						$percCoverages[$j]++;
					}
					else {
						$percCoverages[$j] = 1;
					}
				}
			}

			$posCounter++;
		}

		while ( $posCounter < $chromLength{$currChr} )
		{    #fill uncovered end of chromosome with zeros
			print WIG "0\n";
			if ($zeros) {
				print ZEROS "$currChr\t$posCounter\tn\t0\n";
			}

			$uncovered++;
			$posCounter++;
		}

		my $avgCov = $summedCov / $chromLength{$currChr};
		$targetCount += $chromLength{$currChr};
		if ( $currChr =~ m/^(c|C)hr\d+$/ ) {    # only for chromosomes 1-22
			$autoRange += $chromLength{$currChr};
		}
		$avgCov = sprintf( "%.2f", $avgCov );
		print AVG "$currChr:1-$chromLength{$currChr}\t$summedCov\t$avgCov\n";
		print AVGSEG "$sample\t$currChr\t1\t$chromLength{$currChr}\t".($avgCov+2)."\n";
		close PILEUP;
	}

}

##############################################################################
sub median {
	my $ref = shift;
	my $median;
	my @array = ();
	my $n     = 0;
	my $t     = 0;

	$n = ( @{$ref} );
	@{$ref} = sort { $a <=> $b } @{$ref};
	if ( $n % 2 != 0 ) {
		$median = ${$ref}[ int( $n / 2 ) ];
	}
	else {
		$median =
		  ( ${$ref}[ int( $n / 2 - 1 ) ] + ${$ref}[ int( $n / 2 ) ] ) / 2;
	}

	return ($median);
}

=head1 NAME

calcBaseCov.pl

=head1 SYNOPSIS

calcBaseCov.pl -b mrkdup.bam -o testCov.out -t transcript -tt db -m 200 

=head1 DESCRIPTION

This script calculates the coverage for given target regions (or the whole genome) from
a BAM file. It uses the SAMtools Perl API and outputs a variety of different files (see options).

=head1 OPTIONS

 -b	<bamfile> BAM file to calculate coverage from; REQUIRED
 -o	<outfile> main output file; REQUIRED
	The directory of the outfile will be used as default outdir for other output files.
 -t	<targetfile> targeted regions; if empty: calc coverage for whole genome
 -tt	[target type] specifies where to take the target regions from. Currently three possibilities:
		- "file" (default): directly use filename given in -t
		- "settings"      : use a file specified in the config XML file under {settings(-se)}->{targetfile(-t)}
		- "db"            : get targets from the database. Then -t should give a table (or view) that contains the
			                fields 'chrom','exonStarts','exonEnds','name','idtranscript'.
			                The fields should be formated as the 'refGene' table from UCSC.
			                Note that the database that this script connects to is the 'coredb' specified
			                in the config XML file.
		- "chr"           : The entry in -t is the name of a chromosome
 -bed	target file is in BED format
 -w	<wiggle_output.wig>; default: outdir/coverage.wig
 -p	<coverage_per_target_file.txt>; default: outdir/coverage_per_target.txt
 -ps	<coverage_per_target_file.seg>; default: outdir/coverage_per_target.seg
 -c	<coverage> threshhold for which % coverage is given in the outputfile
 -d	<max readdepth> that is reported; normally the samtools api cuts the max coverage at 8000 --> might be to low for e.g. RNA-seq; default: 99999
 -m	<max coverage> if this value is given, a file distribution.<outfile> is created that gives the % of bases covered n-times
	for n=1-<max-coverage>
 -q	calculate mean mapping quality for each region
 -z	save regions with no coverage into a file for homologous deletion detection
 -r	<region_output_file> if this name is given, the script calculates the coverage per region, i.e. for all regions in a BED file
	that have the same name. Useful for calculating e.g. the coverage of a gene from a list of its exons; requires -bed
 -gc	calculate the GC content of each target (i.e. number of G and C bases divided by total number of bases) and write them in coverage
	per target file
- ct	<covered.target.bed> create a BED file that includes the regions from the given target region that are covered by at least one red
	mergeBed from the bedtools is used to create this file
 -se	name of the settings in the current.config.xml file that holds path to reference genome, 
 	to the annotation file and to possible additional annotation files; use default settings if nothing is given
 -lf	log file; default: print to screen
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page

=head1 AUTHOR

Thomas Wieland

=cut

