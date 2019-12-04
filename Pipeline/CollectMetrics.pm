#!/usr/bin/perl
package CollectMetrics;

use strict;
use File::Basename;
use Cwd qw(abs_path);


my $prog_path = dirname( abs_path($0) );
require $prog_path."/Utilities.pm";

##################################################################################################
# Calls picard tools "CollectAlignmentSummaryMetrics"; doesn't return anything (for now)
##################################################################################################
sub calcAlignmentMetrics{
	my $bamfile = shift;
	my $outdir  = shift;
	my $settings= shift;
	
	my $params  = Utilities::getParams();
	my $java    = $params->{programs}->{java}->{path};
	my $picard  = $params->{programs}->{picard}->{path};
	my $ref     = $params->{settings}->{$settings}->{reference};
	
	my $output  = "$outdir/alignment.metric";
	
	# Returns error status of command
	return system("$java -Xmx6g -XX:ParallelGCThreads=1 -jar $picard CollectAlignmentSummaryMetrics METRIC_ACCUMULATION_LEVEL=LIBRARY INPUT=$bamfile OUTPUT=$output REFERENCE_SEQUENCE=$ref VALIDATION_STRINGENCY=SILENT");
}

##################################################################################################
# Calls picard tools "CollectBaseDistributionByCycle"; doesn't return anything (for now)
##################################################################################################
sub calcBaseDistributionByCycle{
	my $bamfile = shift;
	my $outdir  = shift;
	my $settings= shift;
	
	my $params  = Utilities::getParams();
	my $java    = $params->{programs}->{java}->{path};
	my $picard  = $params->{programs}->{picard}->{path};
	my $convert = $params->{programs}->{convert}->{path};
	
	my $output  = "$outdir/basedistribution.metric";
	my $chart_output  = "$outdir/basedistribution.pdf";
		my $chart_output_png = 	"$outdir/basedistribution.png";
	my $aligned_reads_only 	= "false";
	my $pf_reads_only		= "false";

	# Returns error status of command
	my $ret = system("$java -Xmx6g -XX:ParallelGCThreads=1 -jar $picard CollectBaseDistributionByCycle INPUT=$bamfile OUTPUT=$output CHART_OUTPUT=$chart_output ALIGNED_READS_ONLY=$aligned_reads_only PF_READS_ONLY=$pf_reads_only VALIDATION_STRINGENCY=SILENT");
	system("$convert $chart_output $chart_output_png");
	return $ret;
}


##################################################################################################
# Calls picard tools "CollectGcBiasMetrics"; doesn't return anything (for now)
##################################################################################################
sub calcGCMetrics{
	my $bamfile = shift;
	my $outdir  = shift;
	my $settings= shift;
	
	my $params  = Utilities::getParams();
	my $java    = $params->{programs}->{java}->{path};
	my $picard  = $params->{programs}->{picard}->{path};
	my $ref     = $params->{settings}->{$settings}->{reference};
	my $convert = $params->{programs}->{convert}->{path};
	
	my $output	= "$outdir/gc.metric";
	my $chart_output = "$outdir/gc.metric.hist.pdf";
		my $chart_output_png = "$outdir/gc.metric.hist.png";
	my $summary_output = "$outdir/gc.metric.summary";
	
	# Returns error status of command
	my $ret = system("$java -Xmx6g -XX:ParallelGCThreads=1 -jar $picard CollectGcBiasMetrics INPUT=$bamfile OUTPUT=$output CHART_OUTPUT=$chart_output SUMMARY_OUTPUT=$summary_output REFERENCE_SEQUENCE=$ref VALIDATION_STRINGENCY=SILENT");
	system("$convert $chart_output $chart_output_png");
	return $ret;
}

##################################################################################################
# Calls picard tools "CollectInsertSizeMetrics" 
##################################################################################################
sub calcInsertSizeMetrics{
	my $bamfile = shift;
	my $outdir  = shift;
	
	my $params  = Utilities::getParams();
	my $java    = $params->{programs}->{java}->{path};
	my $picard  = $params->{programs}->{picard}->{path};
	my $convert = $params->{programs}->{convert}->{path};
	
	my $output	= "$outdir/insert.metric";
	my $output_hist = "$outdir/insert.hist.pdf";
		my $output_hist_png = "$outdir/insert.hist.png";
		
	# Returns error status of command
	my $ret = system("$java -Xmx6g -XX:ParallelGCThreads=1 -jar $picard CollectInsertSizeMetrics VERBOSITY=ERROR INPUT=$bamfile OUTPUT=$output HISTOGRAM_FILE=$output_hist METRIC_ACCUMULATION_LEVEL=LIBRARY VALIDATION_STRINGENCY=SILENT");
	system("$convert $output_hist $output_hist_png");
	return $ret;
}

##################################################################################################
# Calls picard tools "MeanQualityByCycle" 
##################################################################################################
sub calcMeanQualityByCycle{
	my $bamfile = shift;
	my $outdir  = shift;
	
	my $params  = Utilities::getParams();
	my $java    = $params->{programs}->{java}->{path};
	my $picard  = $params->{programs}->{picard}->{path};
	my $convert = $params->{programs}->{convert}->{path};
	
	my $output	= "$outdir/meanqualitybycycle.metric";
	my $output_chart = "$outdir/meanqualitybycycle.chart.pdf";
		my $output_chart_png = "$outdir/meanqualitybycycle.chart.png";
		
	# Returns error status of command
	my $ret = system("$java -Xmx6g -XX:ParallelGCThreads=1 -jar $picard MeanQualityByCycle VERBOSITY=ERROR INPUT=$bamfile OUTPUT=$output CHART=$output_chart VALIDATION_STRINGENCY=SILENT");
	system("$convert $output_chart $output_chart_png");
	return $ret;
}

##################################################################################################
# Calls picard tools "CollectMultipleMetrics" 
##################################################################################################
sub calcMultipleMetrics{
	my $bamfile = shift;
	my $outdir  = shift;
	my $settings= shift;
	
	my $params  = Utilities::getParams();
	my $java    = $params->{programs}->{java}->{path};
	my $picard  = $params->{programs}->{picard}->{path};
	my $ref     = $params->{settings}->{$settings}->{reference};
	
	my $output	= "$outdir/multiple.metric";
	
	#New versions should have more modules
	#my @modules=("CollectAlignmentSummaryMetrics", "CollectInsertSizeMetrics", "QualityScoreDistribution", "MeanQualityByCycle", "CollectBaseDistributionByCycle", "CollectGcBiasMetrics", "RnaSeqMetrics", "CollectSequencingArtifactMetrics", "CollectQualityYieldMetrics");
	my @modules=("CollectAlignmentSummaryMetrics", "CollectInsertSizeMetrics", "QualityScoreDistribution", "MeanQualityByCycle", "CollectBaseDistributionByCycle", "CollectGcBiasMetrics");
	my $numthreads=scalar(@modules);
	my $PROG_COMMAND="PROGRAM=".join(" PROGRAM=", @modules);	
		
	# Returns error status of command
	return system("$java -Xmx6g -XX:ParallelGCThreads=$numthreads -jar $picard CollectMultipleMetrics VERBOSITY=ERROR INPUT=$bamfile OUTPUT=$output REFERENCE_SEQUENCE=$ref $PROG_COMMAND VALIDATION_STRINGENCY=SILENT");
}

##################################################################################################
# Calls picard tools "QualityScoreDistribution" 
##################################################################################################
sub calcQualityScoreDistribution{
	my $bamfile = shift;
	my $outdir  = shift;
	
	my $params  = Utilities::getParams();
	my $java    = $params->{programs}->{java}->{path};
	my $picard  = $params->{programs}->{picard}->{path};
	my $convert = $params->{programs}->{convert}->{path};
	
	my $output	= "$outdir/qualityscoredistribution.metric";
	my $output_chart = "$outdir/qualityscoredistribution.hist.pdf";
		my $output_chart_png = "$outdir/qualityscoredistribution.hist.png";
		
	# Returns error status of command
	my $ret = system("$java -Xmx6g -XX:ParallelGCThreads=1 -jar $picard QualityScoreDistribution VERBOSITY=ERROR INPUT=$bamfile OUTPUT=$output CHART=$output_chart VALIDATION_STRINGENCY=SILENT");
	system("$convert $output_chart $output_chart_png");
	return $ret;
}

##################################################################################################
# Calls picard tools "CollectWgsMetrics" 
##################################################################################################
sub calcWgsMetrics{
	my $bamfile = shift;
	my $outdir  = shift;
	my $settings= shift;
	
	
	my $params  = Utilities::getParams();
	my $java    = $params->{programs}->{java}->{path};
	my $picard  = $params->{programs}->{picard}->{path};
	my $ref     = $params->{settings}->{$settings}->{reference};
	
	my $output	= "$outdir/wgs.metric";
	my $minimum_mapping_quality = 20;
	my $minimum_base_quality	= 20;
	my $coverage_cap = 250; # defaults for WG #TODO: CHECK
	my $include_bq_histogram="true";
	my $count_unpaired="false";			# Include unpaired reads in the stats - don't really want them
		
	# Returns error status of command
	return system("$java -Xmx6g -XX:ParallelGCThreads=1 -jar $picard CollectWgsMetrics VERBOSITY=ERROR INPUT=$bamfile OUTPUT=$output REFERENCE_SEQUENCE=$ref MINIMUM_MAPPING_QUALITY=$minimum_mapping_quality MINIMUM_BASE_QUALITY=$minimum_base_quality COVERAGE_CAP=$coverage_cap INCLUDE_BQ_HISTOGRAM=$include_bq_histogram COUNT_UNPAIRED=$count_unpaired VALIDATION_STRINGENCY=SILENT");
}

##################################################################################################
# Calls picard tools "EstimateLibraryComplexity"; returns library size
##################################################################################################
sub calcEstimateLibraryComplexity{
	my $bamfile = shift;
	my $outdir  = shift;
	
	my $params  = Utilities::getParams();
	my $java    = $params->{programs}->{java}->{path};
	my $picard  = $params->{programs}->{picard}->{path};
	
	my $output	= "$outdir/libcomplexity.metric";

	# Returns error status of command
	return system("$java -Xmx32g -XX:ParallelGCThreads=1 -jar $picard EstimateLibraryComplexity VERBOSITY=ERROR INPUT=$bamfile OUTPUT=$output VALIDATION_STRINGENCY=SILENT TMP_DIR=$outdir");
}



##################################################################################################
# Retrieval commands:
##################################################################################################


##################################################################################################
# Return mean insert sizes and sds per library. If not calculated calls calcInsertSizemetrics
##################################################################################################
sub getInsertSize{
	my $bamfile = shift;
	my $outdir  = shift;
	
	my $params  = Utilities::getParams();
	my $java    = $params->{programs}->{java}->{path};
	my $picard  = $params->{programs}->{picard}->{path};
	
	my $infile  = "$outdir/insert.metric";
	
	# Calc metric if not yet performed
	if(!(-e $infile)){
		calcInsertSizeMetrics( $bamfile, $outdir );
	}
	
	# Open metric	
	open(IN, $infile)|| die print "Can't open $infile\n";
	
	#Skip first lines
	while(<IN>){
		chomp;
		last if $_ =~ /MEDIAN_INSERT_SIZE/;
	}
	
	my $line = <IN>;
		my @columns = split(' ',$line);
			$columns[4] =~ s/,/\./g;
			$columns[5] =~ s/,/\./g;
			my ($mean,$dummy) = split("\\.",$columns[4]);
			my ($std,$dummy2) = split("\\.",$columns[5]);
	
	my $libs;
	while(<IN>){
		@columns = split(' ');
		last if !($columns[0] =~ /\d+/);
		my $tmp;
		($tmp->[0],$dummy) = split("\\.",$columns[4]);
		($tmp->[1],$dummy2) = split("\\.",$columns[5]);

		$libs->{$columns[-1]} = $tmp;
	}
	
	close IN;
	
	return ($mean,$std,$libs);
}

##################################################################################################
# Returns library size estimated by Picard EstimateLibraryComplexity 
##################################################################################################
sub estimateLibraryComplexity{
	my $bamfile = shift;
	my $outdir  = shift;
	
	my $params  = Utilities::getParams();
	my $java    = $params->{programs}->{java}->{path};
	my $picard  = $params->{programs}->{picard}->{path};
	
	my $infile	= "$outdir/libcomplexity.metric";

	# Calc metric if not yet performed	
	if(!(-e "$infile")){
		return 0;
		calcEstimateLibraryComplexity($bamfile, $outdir);
	}
	
	# Open metric
	open(IN, "$infile")|| die print "Can't open $infile\n";

	while(<IN>){
		chomp;
		last if $_ =~ /ESTIMATED_LIBRARY_SIZE/;	#skip first lines
	}
	
	my $line = <IN>;
	chomp $line;
	my @columns = split("\t",$line);
	return $columns[-1];
}

##################################################################################################
# Returns library size estimated by Picard EstimateLibraryComplexity normalised to genome size 
##################################################################################################
sub estimateLibraryComplexityNormalised{
	my $bamfile = shift;
	my $outdir  = shift;
	my $settings= shift;
	
	return estimateLibraryComplexity($bamfile, $outdir) / Utilities::getReferenceLen($settings);
}


##################################################################################################
# Returns fraction of bases with quals (>Q10,>Q20,>Q30) as estimated by QualityScoreDistribution
##################################################################################################
sub getQFraction{
		my $bamfile = shift;
		my $outdir  = shift;
		
		my $params  = Utilities::getParams();
		my $infile  = "$outdir/qualityscoredistribution.metric";
		
		# Q10/20/30 counters
			my $qlt10 = 0;	#10
			my $qgt10 = 0;
			my $qlt20 = 0;	#20
			my $qgt20 = 0;
			my $qlt30 = 0;	#30
			my $qgt30 = 0;
		
		# Calc metric if not yet performed
		if (!(-e $infile )){
			calcQualityScoreDistribution($bamfile,$outdir);
		}

		# Open metric		
		open(IN, $infile) || die print "Can't open $infile\n";
		
		#Skipe first lines
		while(<IN>){
			chomp;
			last if $_ =~ /^QUALITY/;
		}
		
		# Reads histogram
		# format:
		# 		Q [TAB] COUNTBASES
		while(<IN>){
			chomp;
			my @entry=split("\t", $_);

			if ( defined $entry[1] ){			
				#10
				if ($entry[0]<10){  
					$qlt10 += $entry[1];
				}else{
					$qgt10 += $entry[1];
				}
				#20
				if ($entry[0]<20){  
					$qlt20 += $entry[1];
				}else{
					$qgt20 += $entry[1];
				}
				#30
				if ($entry[0]<30){  
					$qlt30 += $entry[1];
				}else{
					$qgt30 += $entry[1];
				}
			}
		}
		
		# Returns % (1-based) of q>(10,20,30) over total
		my $q10frac=( $qgt10 / ( $qlt10 + $qgt10 ));
		my $q20frac=( $qgt20 / ( $qlt20 + $qgt20 ));
		my $q30frac=( $qgt30 / ( $qlt30 + $qgt30 ));
		
		return ($q10frac, $q20frac, $q30frac);
}

##################################################################################################
# Returns fraction of bases with quals >Q30 as estimated by QualityScoreDistribution
##################################################################################################
sub getQ30Fraction{
		my $bamfile = shift;
		my $outdir  = shift;	
		my @Q=getQFraction($bamfile, $outdir);
		
		return $Q[2];
}

##################################################################################################
# Parse ALIGNMENT METRIC and return hash with statistics OR alignment metric value if specified name as 4th argument
##################################################################################################
sub getAlignmentMetric{
		my $bamfile = shift;
		my $outdir  = shift;
		my $settings= shift;
		my $which_metric = shift;

		my $params  = Utilities::getParams();
		my $infile  = "$outdir/alignment.metric";
		
		# Calc metric if not yet performed
		if (!(-e $infile )){
			calcAlignmentMetrics($bamfile,$outdir, $settings);
		}

		# Open metric
		open(IN, $infile) || die print "Can't open $infile\n";
		
		my $head="";
		my $data="";
		
		#Get HEAD
		while(<IN>){
			chomp;
			if ($_ =~ /^CATEGORY/){
				$head=$_;
			}
			
			$data=$_;
			
			last if $_ =~ /^PAIR/;
		}
		
		#my %alignment_metric{split("\t", $head) } = split("\t", $data);
		my %alignment_metric;
		    @alignment_metric{split("\t", $head) } = split("\t", $data);
		    
		return ( $which_metric eq "" ? 
											%alignment_metric
											:
											$alignment_metric{$which_metric});
			
}

##################################################################################################
# Parse BASE DISTRIBUTION and return BASE PERCENTAGE AVERAGE PLUS STDEV WITHIN CYCLES
##################################################################################################
sub getBaseDistribution{
		my $bamfile = shift;
		my $outdir  = shift;
		my $settings= shift;

		my $params  = Utilities::getParams();
		my $infile  = "$outdir/basedistribution.metric";
		
		# Calc metric if not yet performed
		if (!(-e $infile )){
			calcBaseDistributionByCycle($bamfile,$outdir, $settings);
		}

		# Open metric
		open(IN, $infile) || die print "Can't open $infile\n";
	
		# Vars
		my $data="";
		my $total_cycles=0;
		# ACGT Counts and SD
			my $A=0; my $ASD=0;
			my $C=0; my $CSD=0;
			my $G=0; my $GSD=0;
			my $T=0; my $TSD=0;
			my $N=0; my $NSD=0;  
		
		#Skip HEAD
		while(<IN>){
			chomp;
			last if $_ =~ /^READ_END/
		}
		
		#Read ACGTN content entries and calculate STDEV
		while(<IN>){
			chomp;
			$data=$_;
			my @datarr=split("\t", $data);
			if ( defined $datarr[6])
			{
				$A=$A+$datarr[2];	$ASD=$ASD+$datarr[2]*$datarr[2];
				$C=$C+$datarr[3];	$CSD=$CSD+$datarr[3]*$datarr[3];
				$G=$G+$datarr[4];	$GSD=$GSD+$datarr[4]*$datarr[4];
				$T=$T+$datarr[5];	$TSD=$TSD+$datarr[5]*$datarr[5];
				$N=$N+$datarr[6];	$NSD=$NSD+$datarr[6]*$datarr[6];
				$total_cycles++;
			}
		}
		
		#A
		$ASD = sqrt (  (($total_cycles*$ASD) - ($A*$A))/($total_cycles*($total_cycles-1)) );
			$A = $A/$total_cycles;
		#C		
		$CSD = sqrt (  (($total_cycles*$CSD) - ($C*$C))/($total_cycles*($total_cycles-1)) );
			$C = $C/$total_cycles;
		#G			
		$GSD = sqrt (  (($total_cycles*$GSD) - ($G*$G))/($total_cycles*($total_cycles-1)) );
			$G = $G/$total_cycles;
		#T			
		$TSD = sqrt (  (($total_cycles*$TSD) - ($T*$T))/($total_cycles*($total_cycles-1)) );
			$T = $T/$total_cycles;	
		#N		
		$NSD = sqrt (  (($total_cycles*$NSD) - ($N*$N))/($total_cycles*($total_cycles-1)) );
			$N = $N/$total_cycles;

		#Build result hash
		my %baseDistribution;
			$baseDistribution{"A"}=$A;		$baseDistribution{"ASD"}=$ASD;
			$baseDistribution{"C"}=$C;		$baseDistribution{"CSD"}=$CSD;
			$baseDistribution{"G"}=$G;		$baseDistribution{"GSD"}=$GSD;
			$baseDistribution{"T"}=$T;		$baseDistribution{"TSD"}=$TSD;
			$baseDistribution{"N"}=$N;		$baseDistribution{"NSD"}=$NSD;
						
		return ( %baseDistribution );
}

##################################################################################################
# Parse BASE DISTRIBUTION and return count of cycles exceeding n (or default 4 ) sigmas from normal base composition
##################################################################################################
sub getBaseDistributionBadCycles{
		my $bamfile = shift;
		my $outdir  = shift;
		my $settings= shift;
		my $req_sigmas = shift;
		
		if ( $req_sigmas eq "" )
		{
			$req_sigmas=4;
		}

		my $params  = Utilities::getParams();
		my $infile  = "$outdir/basedistribution.metric";
		
		# Calc metric if not yet performed
		if (!(-e $infile )){
			calcBaseDistributionByCycle($bamfile,$outdir, $settings);
		}

		# Vars
		my $data="";
		my $badCycles=0;
		
		# Get Computed base composition and stdev
		my %baseDistribution = getBaseDistribution($bamfile, $outdir, $settings);
			
		# Open metric
		open(IN, $infile) || die print "Can't open $infile\n";
	
		#Skip HEAD
		while(<IN>){
			chomp;
			last if $_ =~ /^READ_END/
		}
		
		#Read ACGTN content entries and evaulate number deviation from average
		while(<IN>){
			chomp;
			$data=$_;
			my @datarr=split("\t", $data);
			if ( defined $datarr[6])
			{
				if ( 
						$datarr[2] > ( $baseDistribution{"A"}+$req_sigmas*$baseDistribution{"ASD"} )	||
						$datarr[2] < ( $baseDistribution{"A"}-$req_sigmas*$baseDistribution{"ASD"} )	||
						$datarr[3] > ( $baseDistribution{"C"}+$req_sigmas*$baseDistribution{"CSD"} )	||
						$datarr[3] < ( $baseDistribution{"C"}-$req_sigmas*$baseDistribution{"CSD"} )	||
						$datarr[4] > ( $baseDistribution{"G"}+$req_sigmas*$baseDistribution{"GSD"} )	||
						$datarr[4] < ( $baseDistribution{"G"}-$req_sigmas*$baseDistribution{"GSD"} )	||
						$datarr[5] > ( $baseDistribution{"T"}+$req_sigmas*$baseDistribution{"TSD"} )	||
						$datarr[5] < ( $baseDistribution{"T"}-$req_sigmas*$baseDistribution{"TSD"} )	
					){
						$badCycles++;
					}
			}
		}
						
		return ( $badCycles );
}
	
##################################################################################################
# Parse GC Bias and returns top val
##################################################################################################
sub getGCMetric{
		my $bamfile = shift;
		my $outdir  = shift;
		my $settings= shift;
		my $which_metric=shift;

		my $params  = Utilities::getParams();
		my $infile  = "$outdir/gc.metric.summary";
		
		# Calc metric if not yet performed
		if (!(-e $infile )){
			calcGCMetrics($bamfile,$outdir, $settings);
		}

		# Open metric
		open(IN, $infile) || die print "Can't open $infile\n";
	
			my $head="";
		my $data="";
		
		#Get HEAD
		while(<IN>){
			chomp;
			if ($_ =~ /^ACCUMULATION_LEVEL/){
				$head=$_;
			}
			
			$data=$_;
			
			last if $_ =~ /^All/;
		}
		
		#my %alignment_metric{split("\t", $head) } = split("\t", $data);
		my %metric;
		    @metric{split("\t", $head) } = split("\t", $data);
		    
		return ( $which_metric eq "" ? 
											%metric
											:
											$metric{$which_metric});

}

##################################################################################################
# Parse meanqualitybycycle and GET average quality
##################################################################################################
sub getMeanQuality{
		my $bamfile = shift;
		my $outdir  = shift;
		my $settings= shift;

		my $params  = Utilities::getParams();
		my $infile  = "$outdir/meanqualitybycycle.metric";
		
		# Calc metric if not yet performed
		if (!(-e $infile )){
			calcMeanQualityByCycle($bamfile,$outdir, $settings);
		}

		# Open metric
		open(IN, $infile) || die print "Can't open $infile\n";
	
		# Vars
		my $data="";
		my $total_cycles=0;
		my $Q=0;
		my $QSD=0;
		
		#Skip HEAD
		while(<IN>){
			chomp;
			last if $_ =~ /^CYCLE/;
		}
		
		#Read Qs
		while(<IN>){
			chomp;
			$data=$_;
			my @datarr=split("\t", $data);
			if ( defined $datarr[1])
			{
				$Q+=$datarr[1];
				$QSD+=$datarr[1]*$datarr[1];
				$total_cycles++;
			}
		}
		
		$QSD = sqrt (  (($total_cycles*$QSD) - ($Q*$Q))/($total_cycles*($total_cycles-1)) );
		$Q=($Q/$total_cycles);
		
		my %metric;
			$metric{"Q"}=$Q;
			$metric{"QSD"}=$QSD;
			
		return %metric;
}

##################################################################################################
# Parse meanqualitybycycle and GET average quality for the last N cycles 4th argument 
##################################################################################################
sub getMeanQualityLastN{
		my $bamfile = shift;
		my $outdir  = shift;
		my $settings= shift;
		
		my $lastN = shift;
		
			if ( $lastN eq "" )
			{
				$lastN=5;
			}

		my $params  = Utilities::getParams();
		my $infile  = "$outdir/meanqualitybycycle.metric";
		
		# Calc metric if not yet performed
		if (!(-e $infile )){
			calcMeanQualityByCycle($bamfile,$outdir, $settings);
		}

		# Open metric
		open(IN, $infile) || die print "Can't open $infile\n";
	
		# Vars
		my $data="";
		my $total_cycles=0;
		
		my @Qs;
		my $QTot=0;
		
		my $Q=0;
		my $QSD=0;
		
		#Skip HEAD
		while(<IN>){
			chomp;
			last if $_ =~ /^CYCLE/;
		}
		
		# Accum		
		while(<IN>){
			chomp;
			$data=$_;
			my @datarr=split("\t", $data);
			if ( defined $datarr[1])
			{
				$Qs[$QTot]=$datarr[1];
			 	$QTot++;
			}
		}
		
		#Prevent index < 0 
		if (( $QTot-1-$lastN )<0)
		{
			$lastN=$QTot-1;
		}
		
		# Get it
		for ( my $Loop=$QTot-1; $Loop>$QTot-1-$lastN; $Loop--)
		{
			$Q+=$Qs[$Loop];
			$QSD+=$Qs[$Loop]*$Qs[$Loop];
		}

		
		$QSD = $lastN>1 ? sqrt (  (($lastN*$QSD) - ($Q*$Q))/($lastN*($lastN-1)) ) : 0 ;
		$Q=($Q/$lastN);
		
		my %metric;
			$metric{"Q"}=$Q;
			$metric{"QSD"}=$QSD;
			
		return %metric;
}

##################################################################################################
# Parse meanqualitybycycle and return count of cycles below n (or default 4) sigma from normal qscore
##################################################################################################
sub getMeanQualityBadCycles{
		my $bamfile = shift;
		my $outdir  = shift;
		my $settings= shift;
		my $req_sigmas = shift;
		
		if ( $req_sigmas eq "" )
		{
			$req_sigmas=4;
		}

		my $params  = Utilities::getParams();
		my $infile  = "$outdir/meanqualitybycycle.metric";
		
		# Calc metric if not yet performed
		if (!(-e $infile )){
			calcMeanQualityByCycle($bamfile,$outdir, $settings);
		}

		# Vars
		my $data="";
		my $badCycles=0;
		
		# Get Computed base composition and stdev
		my %baseQual = getMeanQuality($bamfile, $outdir, $settings);
			
		# Open metric
		open(IN, $infile) || die print "Can't open $infile\n";
	
		#Skip HEAD
		while(<IN>){
			chomp;
			last if $_ =~ /^CYCLE/
		}
		
		#Read Qs
		while(<IN>){
			chomp;
			$data=$_;
			my @datarr=split("\t", $data);
			if ( defined $datarr[1])
			{
				if ( 
						$datarr[1] < ( $baseQual{"Q"}-$req_sigmas*$baseQual{"QSD"} )
					){
						$badCycles++;
					}
			}
		}
						
		return ( $badCycles );
}

##################################################################################################
# Parse rmdup.metric
##################################################################################################
sub getDuplicates{
		my $bamfile = shift;
		my $outdir  = shift;
		my $settings= shift;

		my $params  = Utilities::getParams();
		my $infile  = "$outdir/rmdup.metric";
		
		# Return NULL if metric not calculated
		if (!(-e $infile )){
			return 'NULL';
		}

		# Open metric
		open(IN, $infile) || die print "Can't open $infile\n";
	
			my $head="";
		my $data="";
		
		#Get HEAD
		while(<IN>){
			chomp;
			last if ($_ =~ /^LIBRARY/)
		}			
		
		my $total_reads_examined=0;
		my $total_duplicates=0;
		my $total_optical_duplicates=0;

		# Sum over libraries
		
		while(<IN>){
			chomp;
			$data=$_;

			if ($data ne "")
			{
				#0     		1						2					3				4							5						6								7					8
				#LIBRARY	UNPAIRED_READS_EXAMINED	READ_PAIRS_EXAMINED	UNMAPPED_READS	UNPAIRED_READ_DUPLICATES	READ_PAIR_DUPLICATES	READ_PAIR_OPTICAL_DUPLICATES	PERCENT_DUPLICATION	ESTIMATED_LIBRARY_SIZE
				$data =~ /\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t([0-9]*\.*[0-9]*)\t(\d+)/;
				$total_reads_examined     += $1 + (2 * $2);
				$total_duplicates         += $4 + (2 * $5);
				$total_optical_duplicates +=      (2 * $6);				
			}
		}
		
		
		my %metric;
			$metric{"DUPLICATION_RATE"} 		=	$total_duplicates 			/ 	$total_reads_examined;
		    $metric{"OPTICAL_DUPLICATION_RATE"} =	$total_optical_duplicates 	/ 	$total_reads_examined;
		    
		    
		return ( %metric );
		

}
