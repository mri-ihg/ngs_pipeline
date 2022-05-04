#!/usr/bin/perl

##############################################################
# Tim M Strom May 2010
# generates a config file for pipeline.pl
# Sebastian Eck May 2010
# added subroutines:
# -exome stats (locateReads_capture, calcBaseCov & parseStats)
# -annotateIndelpe
# -filterSNP
# -checkPileup
# -depth (add uncovered bases to pileup)
# -del_hom (homozygote deletions)
##############################################################

use strict;
use Getopt::Long;
use File::Basename;
use File::Path qw(make_path);
use Cwd qw(abs_path);
umask(002);

my $prog_path = dirname( abs_path($0) );
require $prog_path . "/Utilities.pm";

my $inifile = "pipeline.ini";

GetOptions( "i=s" => \$inifile );
my $outfile = $inifile;
$outfile =~ s/\.ini/\.cfg/;
my $indir    = "";
my $lane     = "";
my $prefix   = "";
my @indirs   = ();
my @lanes    = ();
my @prefixes = ();

my $outdir     = "";
my $folder     = "";
my $subfolder  = "";
my $fastqdir   = "";
my $project    = "";
my $sample     = "";
my $target     = "";
my $type       = "";
my $item1      = "";
my $slot       = "";
my $line       = "";
my $settings   = "";
my $organism   = "";
my $aligner    = "";
my %runs       = ();
my $params     = Utilities::getParams();
my $assay      = "";
my $targetBED  = "";
my $maxInssize = -1;
my $vcf        = 0;
my @removes;
my @libraries       = ();
my @infiletype      = ();    #possible types: ILLUMINA, FASTQ, BAM
my $targetWin       = 0;
my $doTrim			= "";
my $targetMarginBED = "";
my $targetMarginPileup;
my %validFiles;
my @failedLanes;
my $mergeFast	= "";
my $noRmdup     = 0;
my $pindelreg   = "";
my %slots;				#SGE slots/threads that should be given to the script
my $starInfile1 = "";
my $starInfile2 = "";
my $starReadgroup = "";
my $tophatReadgroup = "";
my $tophatInfile1 = "";
my $tophatInfile2 = "";
my $isStrandedRNA = 0;
my $tmpArgument ="";

my $removeAdapters=0;
my $removePhix=0;

#read ini file
open( IN,  "$inifile" )  || die "Cannot open $inifile\n";
open( OUT, ">$outfile" ) || die "Cannot open $outfile\n";
while (<IN>) {
	if (/^\s+/) { next; }
	if (/^#/)   { next; }
	chomp;
	s/\s+//g;
	$line = $_;
	( $type, $item1,$slot ) = split( /\:/, $line );
	if ( $type eq "indir" ) {    #gerald dir
		                         #$indir = $item1;
		push( @indirs, $item1 );
	}
	elsif ( $type eq "lane" ) {
		push( @lanes, $item1 );
	}
	elsif ( $type eq "infiletype" ) {
		push( @infiletype, $item1 );
	}
	elsif ( $type eq "prefix" ) {

		#$prefix = $item1;
		push( @prefixes, $item1 );
	}
	elsif ( $type eq "outdir" ) {
		$outdir = $item1;
	}
	elsif ( $type eq "fastqdir" ) {
		$fastqdir = $item1;
	}
	elsif ( $type eq "folder" ) {
		$folder = $item1;
	}
	elsif ( $type eq "failedlanes" ) {
		@failedLanes = split( ",", $item1 );
	}
	elsif ( $type eq "subfolder" ) {
		$subfolder = $item1;
	}
	elsif ( $type eq "project" ) {
		$project = $item1;
	}
	elsif ( $type eq "sample" ) {
		$sample = $item1;
	}
	elsif ( $type eq "dotrim" ) {
		if ($item1 eq "TRUE") {
			$doTrim = "true";
		}
	}
	elsif ( $type eq "fastmerge") {
		if ($item1 eq "TRUE") {
			$doTrim = "true";
		}
	}
	elsif ( $type eq "target" ) {
		$target = $item1;
	}
	elsif ( $type eq "maxinssize" ) {
		$maxInssize = $item1;
	}
	elsif ( $type eq "targetWin" ) {
		$targetWin = $item1;
	}
	elsif ( $type eq "settings" ) {
		$settings = $item1;
	}
	elsif ( $type eq "organism" ) {
		$organism = $item1;
	}
	elsif ( $type eq "run" ) {
		
		$runs{$item1} = 1;
		if($slot && $slot ne ""){
			$slots{$item1} = $slot;
		}
		
	}
	
	elsif ( $type eq "assay" ) {
		$assay = $item1;
	}
	elsif ( $type eq "rm" ) {
		push( @removes, $item1 );
	}
	elsif ( $type eq "version" ) {
		if ( lc $item1 eq "vcf" ) {
			$vcf = 1;
		}
		elsif ( lc $item1 eq "gatk" ) {
			$vcf = 2;
		}
	}
	elsif ( $type eq "library" ) {
		push( @libraries, $item1 );
	}
	elsif ( $type eq "aligner" ) {
		$aligner = $item1;
	}
	elsif ( $type eq "normdup") {
		$noRmdup = 1 if $item1 eq "TRUE";
	}
	elsif ( $type eq "pindelreg" ) {
		$pindelreg = $item1;
	} elsif ($type eq "strandedRNA") {
		$isStrandedRNA = $item1;
	} elsif ($type eq "tmpArg") {
		$tmpArgument = $item1;
	} elsif ($type eq "removeAdapters") {
		$removeAdapters=1 if $item1 eq "TRUE";
	} elsif ($type eq "removePhix") {
		$removePhix=1 if $item1 eq "TRUE";
	}
}


my $debug = 0;

unless ($debug) {
	#make directories
	make_path("$outdir/$project/$sample/$folder/$subfolder", {
		mode => 0775
	});
	
	if ( $params->{settings}->{$settings} eq "" ) {
		die "settings not found in config file\n";
	}

	my $refgenome = $params->{settings}->{$settings}->{reference};

}
my $projectdatabase =
  $params->{settings}->{$settings}->{variationdb}->{database};

if ( $target ne "" ) {

	#build name to target bed file
	my @columns = split( /\./, $target );
	$columns[-1] = "bed";
	$targetBED = join( ".", @columns );
	unless ( -e $targetBED ) {
		system("cat $target | $prog_path/exomeTargets2BED.pl > $targetBED")
		  ;    #build target file in bed format
	}

	#build name of target+window bed
	$targetMarginBED = $targetBED;
	if ( $targetWin != 0 ) {
		my $marginname = "_" . $targetWin . "b_margin.bed";
		$targetMarginBED =~ s/\.bed/$marginname/;
		unless ( -e $targetMarginBED ) {
			system(
"cat $targetBED | $prog_path/addMarginToBED.pl -m $targetWin > $targetMarginBED"
			);
		}
	}

	#build name to pileup file
	if ( $vcf == 0 ) {
		$targetMarginPileup = $targetMarginBED;
		$targetMarginPileup =~ s/\.bed/\.pileup\.list/;
		unless ( -e $targetMarginPileup ) {
			system(
"cat $targetMarginBED | $prog_path/bed2pileup.pl > $targetMarginPileup"
			);
		}
	}
}

#get library type from outfolder
my $libtype = "";
my $libpair = "";
if ( $folder =~ /(.*)out$/ ) {
	$libtype = $1;
}
if ( $subfolder =~ /(.*)out$/ ) {
	$libpair = $1;
}

if ($isStrandedRNA == 0 && $assay && $assay ne "") {
	if ($assay =~ m/stranded/gi) {
		$isStrandedRNA = 1;
	}
}

print OUT "
##############################################################
# general information for the pipeline execution script
##############################################################
";
print OUT "outdir	: $outdir/$project/$sample/$folder/$subfolder\n";
print OUT "project  : $project\n";
print OUT "sample   : $sample\n";
print OUT "settings : $settings\n";
print OUT "organism : $organism\n";
print OUT "libtype  : $libtype\n";
print OUT "libpair  : $libpair\n\n";


&alignment;
&bam_merge;
&varfilter;

&deepVariant;

&checkContamination;
&transcriptStats;


&runSissrs;
&runMacs;

&cnvnator;
&annotateSNPs_sam;
&checkPileup;

&snvdbInsertExome;

&callChrM;
&annotateChrM;
&insertChrM;

&exomeDepth;
&snpEff;
&homozygosity;
&updateRefSeq;
&updatelncRNA;
&updatemiRNA;
&runPindel;
&runBreakdancer;
&runLumpy;
&runManta;
&runWhamg;
&insertSV;

&bellerophon;
&IGVtools;


#stats
&stats;
&transversion;



#rnaseq specific stuff
&htseqCount;
&calcFPKM;
&insertRNAresults2DB;
&calcRNAstats;
&insertRNAStats2DB;


sub alignment {
	for ( my $i = 0 ; $i < @indirs ; $i++ ) {
		my @lane = split( //, $lanes[$i] );

		my $library;
		if ( $libraries[$i] ) {
			$library = $libraries[$i];
		}
		else {
			$library = "$sample\_LIB1";
		}
		
		if((-e $indirs[$i]) && !(-d $indirs[$i])){		#file name to align is given directly
			&printAln(
				$indirs[$i],
				$indirs[$i],
				$infiletype[$i],
				$prefixes[$i],
				$library,
				"0"
			);
		}elsif($indirs[$i] =~ /,/){								#two files for alignment are given --> read1 & read2 files
			my ($file1,$file2,@dummy) = split(",",$indirs[$i]);
			&printAln(
				$file1,
				$file2,
				$infiletype[$i],
				$prefixes[$i],
				$library,
				"0"
			);
		}else {
			foreach (@lane) {
	
				#find input files
				
				if ( -e "$indirs[$i]/s\_$_\_1\_sequence.txt" )
				{    #old directory structure
					&printAln(
						"$indirs[$i]/s\_$_\_1\_sequence.txt",
						"$indirs[$i]/s\_$_\_2\_sequence.txt",
						$infiletype[$i],
						$prefixes[$i],
						$library,
						$_,
						""
					);
				}
				elsif ( -e "$indirs[$i]/s\_$_\_1\_sequence.fastq.gz" )
				{    #old directory structure, but zipped fastq files
					&printAln(
						"$indirs[$i]/s\_$_\_1\_sequence.fastq.gz",
						"$indirs[$i]/s\_$_\_2\_sequence.fastq.gz",
						"STANDARD",
						$prefixes[$i],
						$library,
						$_,
						""
					);
	
				}
				elsif ( $infiletype[$i] eq "BAM"
					&& -e "$indirs[$i]/s\_$_\_sequence.bam" )
				{    #bam file to align
					&printAln(
						"$indirs[$i]/s\_$_\_sequence.bam",
						"$indirs[$i]/s\_$_\_sequence.bam",
						$infiletype[$i],
						$prefixes[$i],
						$library,
						$_,
						""
					);
	
				}
				else {    #new structure bam file
					my $formatedLane = sprintf( "%03d", $_ );
	
					my @files;
					if ( $infiletype[$i] eq "BAM" ) {
						@files = glob( $indirs[$i] . "/$sample*L$formatedLane*.bam" );
					}
					else {
						@files =
						  glob( $indirs[$i] . "/$sample*L$formatedLane\_R1*.fastq.gz" );
						if ( scalar(@files) eq 0  )
							{
								my $sampleUnderscoreCured=$sample; $sampleUnderscoreCured =~ s/\_/\-/g;
								@files = glob( $indirs[$i] . "/$sampleUnderscoreCured*L$formatedLane\_R1*.fastq.gz" );
							}
					}
					foreach my $file (@files) {
						chomp $file;
						my $series;
						if ( $file =~ /_(\d+).fastq/ ) {
							$series = "_" . $1;
						}
						my $file1 = $file;
						$file =~ s/_R1_/_R2_/;
	
						&printAln( $file1, $file, $infiletype[$i], $prefixes[$i],
							$library, $_, $series );
						if ( $doTrim ) {
							if ($mergeFast) {
								#get output prefix
								my $outprefix = &Utilities::prepareOutputPrefix($prefix, "$outdir/$project/$sample/$folder/$subfolder/", $file1);
	
	
								#TODO: add functionality for fast merge
								#2nd step: do not realign the entire trimmed fastq file but only those reads not in the first bam file
								#3rd step: just trim reads that weren't aligned in 1st and 2nd step and align those
								
								#1st trim: default 0,20
								&prepareUnmappedFastq($outprefix . ".gem.bam", $file1, $file, $infiletype[$i], $prefixes[$i], $library, $_, $series, "0,20");
								#result file1 & file2 of prepareUnmappedFastq must be input for printAln() below
								&printAln(&Utilities::prepareTrimmingFileName($outprefix.".trim", "0,20") . ".fastq", &Utilities::prepareTrimmingFileName($outprefix.".trim", "0,20") . ".fastq", $infiletype[$i], $prefixes[$i], $library, $_, $series, "0,0");  #TODO: pass the correct infiles (modified fastq)
								
								#2nd trim: default 5,20
								&prepareUnmappedFastq($outprefix . ".gem.bam", $file1, $file, $infiletype[$i], $prefixes[$i], $library, $_, $series, "5,20");
								#result file1 & file2 of prepareUnmappedFastq must be input for printAln() below
								&printAln(&Utilities::prepareTrimmingFileName($outprefix.".trim", "5,20") . "._1.fastq", &Utilities::prepareTrimmingFileName($outprefix.".trim", "5,20") . "_2.fastq", $infiletype[$i], $prefixes[$i], $library, $_, $series, "0,0");  #TODO: pass the correct infiles (modified fastq)
							} else {
								&printAln( $file1, $file, $infiletype[$i], $prefixes[$i], $library, $_, $series, "0,20" );
								&printAln( $file1, $file, $infiletype[$i], $prefixes[$i], $library, $_, $series, "5,20" );
								&performGEMMerge($file1, $file, $infiletype[$i], $prefixes[$i], $library, $_, $series, "0,20", "5,20");
							}
						} else {
							#sort the single bam file
							#&performGEMSingleSorting($file1, $file, $infiletype[$i], $prefixes[$i], $library, $_, $series) if $aligner eq "gem";
						}
					}
				}
			}
		}
	}

	#remove files that are not valid anymore (e.g. the library has been set to "failed")
	if (@failedLanes) {
		
		foreach (@failedLanes) {
			
			my @availableFiles =
		  		glob("$outdir/$project/$sample/$folder/$subfolder/*$_*");
	
			foreach my $currFile (@availableFiles) {
				print "remove $currFile\n";
				unlink($currFile);
			}
		}
	}
	
	if ($aligner eq "star") {
		&printSTARAln($starInfile1, $starInfile2, $starReadgroup);
	} elsif ($aligner eq "tophat") {
		&printTophatAln($tophatInfile1, $tophatInfile2, $tophatReadgroup);
	}
}

sub prepareUnmappedFastq() {
	my $bamfile  = shift;
	my $file1    = shift;
	my $file2    = shift;
	my $filetype = shift;
	my $prefix   = shift;
	my $library  = shift;
	my $lane     = shift;
	my $series   = shift;
	my $trim     = shift;
	
	#get output prefix
	my $outprefix = &Utilities::prepareOutputPrefix($prefix, "$outdir/$project/$sample/$folder/$subfolder/", $file1);
	
	print OUT "
##############################################################
# prepare fastq file for gem mapping $lane (trimming $trim)
##############################################################
		";
	print OUT "pgr   : gem_prepare_fastq.pl\n";
	print OUT "# settings from config file\n";
	print OUT "param : se : $settings\n";
	print OUT "# infile1 \n";
	print OUT "param : f : $file1\n";
	if ( -e $file2 ) {
		print OUT "# infile2 \n";
		print OUT "param : F : $file2\n";
	}
	print OUT "# bamfile\n";
	print OUT "param : b : $bamfile\n";
	print OUT "# output prefix\n";
	print OUT "param : o : " . &Utilities::prepareTrimmingFileName($outprefix.".trim", $trim) . "\n";
	print OUT "#Read Group Sample\n";
	print OUT "param : s : $sample\n";
	print OUT "run\n";
}


# Print Alignment SLOT
sub printAln {

	my $file1    = shift;
	my $file2    = shift;
	my $filetype = shift;
	my $prefix   = shift;
	my $library  = shift;
	my $lane     = shift;
	my $series   = shift;
	my $trim     = shift;
	$trim = "0,0" unless $trim;

	my $external_bam = 0;
	
	# External BAM is marked as BAMX
	if ( $filetype eq "BAMX" )
	{
		$external_bam = 1;
		$filetype = "BAM";
	}

	#get output prefix
	my $outprefix =  &Utilities::prepareOutputPrefix($prefix, "$outdir/$project/$sample/$folder/$subfolder/", $file1);
	my %trimP = &Utilities::parseTrimmingParameter($trim);
	if (($trimP{"x"} > 0) || ($trimP{"y"} > 0)) {
		$outprefix = &Utilities::prepareTrimmingFileName($outprefix.".trim", $trim);
	}

	if ( $aligner eq "gem" ) {
		print OUT "
##############################################################
# perform gem mapping lane $lane (trimming $trim)
##############################################################
";
		print OUT "pgr   : gem_aln.pl\n";
		print OUT "# settings from config file\n";
		print OUT "param : se : $settings\n";
		print OUT "# infile1 \n";
		print OUT "param : f : $file1\n";
		if ( -e $file2 ) {
			print OUT "# infile2 \n";
			print OUT "param : F : $file2\n";
		}
		print OUT "# output prefix\n";
		print OUT "param : o : $outprefix.gem\n";

		print OUT "#Read Group Sample\n";
		print OUT "param : s : $sample\n";
		print OUT "#Read Group Library Name\n";
		print OUT "param : l : $library\n";

		if ( $filetype eq "BAM" ) {
			print OUT "#input is in BAM format\n";
			print OUT "param : b\n";
		}
		
		print OUT "# trim sequence file\n";
		print OUT "param : T : $trim\n";
		&printSlot("alignment");
		
		if ($runs{alignment} == 1 && !(-e "$outprefix.gem.sort.bam")) {
			print OUT "run\n";
			&performGEMSingleSorting($file1, $file2, $filetype, $prefix, $library, $lane, $series) if $aligner eq "gem";
		} else {
			print OUT "#run\n";
		}
	} elsif ($aligner eq "star") {
		#build the infile strings for the star alignment
		$starInfile1 .= "$file1,";
		$starInfile2 .= "$file2," if (-e $file2);
		$starReadgroup .= basename($outprefix) . ",";
	} elsif ($aligner eq "tophat") {
		#build the infile strings for the tophat alignment
		$tophatInfile1 .= "$file1,";
		$tophatInfile2 .= "$file2," if (-e $file2);
		$tophatReadgroup = basename($outprefix);
	} else {
		print OUT "
##############################################################
# perform bwa mapping lane $lane
##############################################################
";
		print OUT "pgr   : bwa_aln_all_Formats.pl\n";
		print OUT "# settings from config file\n";
		print OUT "param : se : $settings\n";
		print OUT "# infile1 \n";

		print OUT "param : f : $file1\n";
		if ( -e $file2 ) {
			print OUT "# infile2 \n";
			print OUT "param : F : $file2\n";
		}
		print OUT "# output prefix\n";
		print OUT "param : o : $outprefix\n";
		if ( $filetype eq "ILLUMINA" ) {
			print OUT "#input is in Illumina format\n";
			print OUT "param : I\n";
		}
		elsif ( $filetype eq "BAM" ) {
			print OUT "#input is in BAM format\n";
			print OUT "param : b\n";
			# If external supplied bam then tell the aligner
			print OUT "param : externalBAM\n" if $external_bam;
		}
   
   	
   		if ( $project eq "S0201" || $removePhix )
   		{
   			print OUT "#DZHKomics samples have Phix spiked in and adapters are not trimmed\n" if $project eq "S0201";
			print OUT "param : removePhix\n";
   		}
   		
   		if ( $project eq "S0201" || $removeAdapters )
   		{
   			print OUT "#DZHKomics samples have Phix spiked in and adapters are not trimmed\n" if $project eq "S0201";
   			print OUT "param : removeAdapters\n";
   		}
		
		print OUT "#Read Group sample\n";
		print OUT "param : s : $sample\n";
		print OUT "#Read Group Library Name\n";
		print OUT "param : l : $library\n";
		if($libtype eq "genomic"){
			print OUT "#For genomic library: create BAM index\n";
			print OUT "param : n\n";
		}
		if ( $maxInssize != -1 ) {
			print OUT "#Maximum insert size for bwa sampe\n";
			print OUT "param : a : $maxInssize\n";
		}

		if ( $aligner eq "bwamem" ) {
			print OUT "#use \"bwa mem\" algorithm for alignment instead of old algorithm\n";
			print OUT "param : m\n";
		}
		&printSlot("alignment");
		
		if ($runs{alignment} == 1 && !(-e "$outprefix.sort.bam")) {
			print OUT "run\n";
		}
		else {
			print OUT "#run\n";
		}
	}
}

sub printSTARAln() {
	my $file1 = shift;
	my $file2 = shift;
	my $srg = shift;
	
	my $filetype = "FASTQ";
	$filetype = "BAM" if ($file1 =~ m/\.bam/g);
	
	my $outprefix =  &Utilities::prepareOutputPrefix($sample, "$outdir/$project/$sample/$folder/$subfolder/", "");
	
	$file1 =~ s/,$//;
	$file2 =~ s/,$//;
	$srg =~ s/,$//;
	
	print OUT "
##############################################################
# perform STAR mapping lane
##############################################################
";
	print OUT "pgr   : star_aln.pl\n";
	print OUT "# settings from config file\n";
	print OUT "param : se : $settings\n";
	print OUT "# infile1 \n";
	print OUT "param : f : $file1\n";
	if ($file2 ne "") {
		print OUT "# infile2 \n";
		print OUT "param : F : $file2\n";
	}
	print OUT "# output prefix\n";
	print OUT "param : o : $outprefix\n";
	
	if ( $filetype eq "BAM" ) {
		print OUT "#input is in BAM format\n";
		print OUT "param : b\n";
	}
		
	print OUT "#Read Group Sample\n";
	print OUT "param : s : $sample\n";
	
	print OUT "#Read Group\n";
	print OUT "param : rg : $srg\n";
	
	&printSlot("alignment");

	if ($runs{alignment} == 1 && !(-e "" . $outprefix . "Aligned.out.bam")) {
		print OUT "run\n";
	} else {
		print OUT "#run\n";
	}
}

sub printTophatAln() {
	my $file1 = shift;
	my $file2 = shift;
	my $trg = shift;
	
	my $filetype = "FASTQ";
	$filetype = "BAM" if ($file1 =~ m/\.bam/g);
	
	my $outprefix =  &Utilities::prepareOutputPrefix($sample, "$outdir/$project/$sample/$folder/$subfolder/TophatAlignment/", "");
	
	$file1 =~ s/,$//;
	$file2 =~ s/,$//;
	
	print OUT "
##############################################################
# perform Tophat mapping lane
##############################################################
";
	print OUT "pgr   : tophat_aln.pl\n";
	print OUT "# settings from config file\n";
	print OUT "param : se : $settings\n";
	print OUT "# infile1 \n";
	print OUT "param : f : $file1\n";
	if ($file2 ne "") {
		print OUT "# infile2 \n";
		print OUT "param : F : $file2\n";
	}
	print OUT "# output prefix\n";
	print OUT "param : o : $outprefix\n";
	
	if ( $filetype eq "BAM" ) {
		print OUT "#input is in BAM format\n";
		print OUT "param : b\n";
	}
		
	print OUT "#Read Group Sample\n";
	print OUT "param : s : $sample\n";
	
	if ($isStrandedRNA) {
		print OUT "#Stranded RNA-seq data\n";
		print OUT "param : sr"
	}
	
	print OUT "#Read Group\n";
	print OUT "param : rg : $trg\n";

	&printSlot("alignment");

	if ($runs{alignment} == 1 && !(-e "" . dirname($outprefix) . "/accepted_hits.bam")) {
		print OUT "run\n";
	} else {
		print OUT "#run\n";
	}
}

sub performGEMMerge() {
	my $file1    = shift;
	my $file2    = shift;
	my $filetype = shift;
	my $prefix   = shift;
	my $library  = shift;
	my $lane     = shift;
	my $series   = shift;
	my $trim1    = shift;
	my $trim2	 = shift;
	
	#get output prefix
	my $outprefix = &Utilities::prepareOutputPrefix($prefix, "$outdir/$project/$sample/$folder/$subfolder/", $file1);
	
	print OUT "
##############################################################
# perform gem merge (of full length, 1st trimmed and 2nd trimmed reads)
##############################################################
";
	print OUT "pgr	:	gem_merge.pl\n";
	print OUT "# output directory\n";
	print OUT "param	: o  : $outprefix\n";
	print OUT "# input files\n";
	print OUT "param	: i1 : $outprefix.gem.bam\n";
	print OUT "param	: i2 : " . &Utilities::prepareTrimmingFileName($outprefix.".trim", $trim1) . ".gem.bam\n";
	print OUT "param	: i3 : " . &Utilities::prepareTrimmingFileName($outprefix.".trim", $trim2) . ".gem.bam\n";
	print OUT "# reference file (usually the initial fastq file\n";
	print OUT "param	: r	: $file1\n";
	print OUT "#sort the resulting bam file by position\n";
	print OUT "param	: s\n";
	print OUT "#Read Group sample\n";
	print OUT "param : n : $sample\n";
	print OUT "#Read Group Library Name\n";
	print OUT "param : l : $library\n";
	print OUT "dependson : gem_aln.pl\n";
	&printSlot("alignment");
	
	if ($runs{alignment} == 1 && !(-e "$outprefix.gem.sort.bam") && $aligner eq "gem") {
		print OUT "run\n";
	} else {
		print OUT "#run\n";
	}
}
sub performGEMSingleSorting() {
	my $file1    = shift;
	my $file2    = shift;
	my $filetype = shift;
	my $prefix   = shift;
	my $library  = shift;
	my $lane     = shift;
	my $series   = shift;
	
	#get output prefix
	my $outprefix = &Utilities::prepareOutputPrefix($prefix, "$outdir/$project/$sample/$folder/$subfolder/", $file1);
	print OUT "
##############################################################
# perform gem single sorting (just full length reads)
##############################################################
";
	print OUT "pgr	:	gem_single_sort.pl\n";
	print OUT "# output directory\n";
	print OUT "param	: o  : $outprefix\n";
	print OUT "# input files\n";
	print OUT "param	: i : $outprefix.gem.bam\n";
	print OUT "#Read Group sample\n";
	print OUT "param : n : $sample\n";
	print OUT "#Read Group Library Name\n";
	print OUT "param : l : $library\n";
	print OUT "dependson : gem_aln.pl\n";
	&printSlot("alignment");
	
	if ($runs{alignment} == 1 && !(-e "$outprefix.gem.sort.bam") && $aligner eq "gem") {
		print OUT "run\n";
	} else {
		print OUT "#run\n";
	}
}

sub bam_merge {
		print OUT "
##############################################################
# perform bam merge
##############################################################
";
	print OUT "pgr   : bwa_merge_gatk.pl\n";
	print OUT "# settings from config file\n";
	print OUT "param : se : $settings\n";
	print OUT "# output directory\n";
	print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder\n";
	print OUT "# aligner\n";
	print OUT "param : a : $aligner\n";
	if($noRmdup || $libtype eq "MIP"){
		print OUT "# don't remove duplicates\n";
		print OUT "param : d\n";
	}
	if($libtype eq "MIP"){
		print OUT "# clip MIP primers instead of removing PCR duplicates\n";
		print OUT "param : c\n";
	}
	if($libtype eq "genomic" || $libtype eq "exomic" || $libtype eq "mtDNA" || $libtype eq "RNA")
	{
		print OUT "# use Picard Markdup\n";
		print OUT "param : picard\n";
	}


	#($libtype eq "genomic"  && $organism eq "human")						
	if( $libtype ne "mtDNA" && $organism eq "human"  && $libtype ne "RNA" ) #Only NovaSeq Samples && human (we need known sites to recalibrate)    #If variant calling for RNAseq should be done, perform recal but perform GATK's SplitNCigarReads
	{
		print OUT "# recal\n";
		print OUT "param : recal\n";
	}	
	
	print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
	print OUT "dependson : bwa_aln_all_Formats.pl,gem_aln.pl,gem_merge.pl,gem_single_sort.pl,star_aln.pl,tophat_aln.pl\n";
	&printSlot("bam_merge");
	
	if ( $runs{bam_merge} == 1 ) {
		print OUT "run\n";
	}
	else {
		print OUT "#run\n";
	}
}





sub varfilter {
	print OUT "
##############################################################
# perform SAMtools/GATK varfilter
##############################################################
";
	if ($libtype eq "RNA"){
		print OUT "pgr   : varfilter_rna.pl\n";
		
		print OUT "#infile\n";
		print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam\n";
		
		print OUT "# output file\n";
		print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/rna.gatk\n";
		
		print OUT "# settings\n";
		print OUT "param : se : $settings\n";

		print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
		print OUT "dependson : bwa_merge_gatk.pl,gem_aln.pl,gem_merge.pl,gem_single_sort.pl\n";
			
		&printSlot("varfilter");
		if ( $runs{varfilter} == 1 ) {
			print OUT "run\n";
		}
		else {
			print OUT "#run\n";
		}
	} else {
		if ( $vcf == 1 ) {
			print OUT "pgr   : varfilter_vcf.pl\n";
			print OUT "#infile\n";
			if($libtype eq "genomic" && !(-e "$outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam.bai") ){
				print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/\n";
			}else{
				print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam\n";
			}
			print OUT "# output directory\n";
			print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder\n";
			
			if ( $target ne "" ) {
				print OUT "# target region\n";
				print OUT "param : l : $targetMarginBED\n";
			}
		}
		elsif ( $vcf == 2 ) {
			print OUT "pgr   : varfilter_gatk.pl\n";
			print OUT "#infile\n";
			print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam\n";
			
			if($libtype eq "MIP")
			{
				print OUT "param : hcminpruning : 10\n";
			}
			print OUT "# output file\n";
			if($libtype eq "genomic" ){
				make_path("$outdir/$project/$sample/$folder/$subfolder/HaplotypeCaller", {				#create directory for single files
					mode => 0775
				});
				
				
				print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/HaplotypeCaller/haplotypecaller.vcf\n";
				print OUT "# start this job as an array job, i.e. one job will be started by split BED file\n";
				print OUT "param : ajb\n";
				print OUT "#signal for parallelpipeline.pl to start this job as an array job on Grid Engine\n";
				print OUT "slots : perBEDArray\n";
				print OUT "param : pcrfree\n" if ($organism eq "human");
				print OUT "param : gatk4\n" if ($organism eq "human");
			}else{
				print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/gatk.ontarget.haplotypecaller.vcf\n";
				print OUT "# target region\n";
				print OUT "param : l : $targetMarginBED\n";
				&printSlot("varfilter");
			}
		
		}
	
		print OUT "# settings\n";
		print OUT "param : se : $settings\n";
		print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
		print OUT "dependson : bwa_merge_gatk.pl,gem_aln.pl,gem_merge.pl,gem_single_sort.pl\n";
	
		if ( $runs{varfilter} == 1 ) {
			print OUT "run\n";
		}
		else {
			print OUT "#run\n";
		}
		
		if($vcf == 1){
			print OUT "
##############################################################
# resort samtools VCF file
##############################################################
";
			print OUT "pgr   : vcfsorter.pl\n";
			print OUT "#infile\n";
			print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/ontarget.varfilter.vcf\n";
			print OUT "# settings from config file\n";
			print OUT "param : se : $settings\n";
			print OUT "# replace infile with sorted file\n";
			print OUT "param : r\n";
			print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
			print OUT "dependson : varfilter_vcf.pl\n";
			&printSlot("varfilter");
			
			if ( $runs{varfilter} == 1 ) {
				print OUT "run\n";
			}
			else {
				print OUT "#run\n";
			}
		}
		
		
		if($vcf == 2 && $libtype eq "genomic"){
			
			
			
			print OUT "
##############################################################
# concat single chromosome VCF & gVCF files
##############################################################
";
			print OUT "pgr   : concatVCF.pl\n";
			print OUT "#indir\n";
			print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/HaplotypeCaller\n";
			print OUT "#outfile\n";
			print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/gatk.ontarget.haplotypecaller.vcf\n";
			print OUT "# files are numbered\n";
			print OUT "param : n\n";
			print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
			print OUT "dependson : varfilter_gatk.pl\n";
			&printSlot("varfilter");
			
			if ( $runs{varfilter} == 1 ) {
				print OUT "run\n";
			}
			else {
				print OUT "#run\n";
			}
			
			print OUT "\npgr   : concatVCF.pl\n";
			print OUT "#indir\n";
			print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/HaplotypeCaller\n";
			print OUT "#outfile\n";
			print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/gatk.ontarget.haplotypecaller.gvcf.gz\n";
			print OUT "#file ending\n";
			print OUT "param : e : gvcf.gz\n";
			print OUT "# files are numbered\n";
			print OUT "param : n\n";
			print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
			print OUT "dependson : varfilter_gatk.pl\n";
			&printSlot("varfilter");
			
			if ( $runs{varfilter} == 1 ) {
				print OUT "run\n";
			}
			else {
				print OUT "#run\n";
			}
		}
		
		
		
		
		if($vcf == 2){
			print OUT "
##############################################################
# filter GATK variants
##############################################################
";
			print OUT "pgr   : filterGATK.pl\n";
			print OUT "#infile\n";
			print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/gatk.ontarget.haplotypecaller.vcf\n";
			print OUT "# output file\n";
			print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/gatk.ontarget.haplotypecaller.filtered.vcf\n";
			print OUT "# GATK HaplotypeCaller was run instead of UnifiedGenotyper\n";
			print OUT "param : gh\n";
			if($libtype eq "exomic" && $organism eq "human"){
				print OUT "# single exomic sample --> too less indels for variant qualit score recalibration --> static filters\n";
				print OUT "param : if\n";
			}elsif($libtype eq "MIP"){
				print OUT "# MIP sample --> use static filters instead of variant recalibration\n";
				print OUT "param : if\n";
				print OUT "param : sf\n";
			}elsif($libtype eq "genomic"){
				print OUT "# use DP for recalibration\n";
				print OUT "param : w\n";
				#print OUT "param : if\n";
				#print OUT "param : sf\n";
			}
			
			if($organism ne "human"){
				print OUT "# no human sample --> use static filters instead of variant recalibration\n";
				print OUT "param : if\n";
				print OUT "param : sf\n";
			}
			print OUT "# settings\n";
			print OUT "param : se : $settings\n";
			print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
			print OUT "dependson : filterVCFforRegion.pl,varfilter_gatk.pl,concatVCF.pl\n";
			&printSlot("varfilter");
			
			if ( $runs{varfilter} == 1 ) {
				print OUT "run\n";
			}
			else {
				print OUT "#run\n";
			}
		}
	}
}


sub deepVariant {
	
	if ( $vcf == 2 && $libtype eq "genomic" && $organism eq "human" )
	{
		print OUT "
##############################################################
# calls variant with DeepVariant
##############################################################
";
		print OUT "pgr   : varfilter_deepvariant.pl\n";
		print OUT "#infile\n";
		print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam\n";
		print OUT "# output directory\n";
		print OUT "param : outdir : $outdir/$project/$sample/$folder/$subfolder/deepvariant\n";
		print OUT "# Sample name\n";
		print OUT "param : s : $sample\n";
		print OUT "# Sample type\n";
		print OUT "param : type : $libtype\n";
	
		print OUT "# settings\n";
		print OUT "param : se : $settings\n";
		print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
		print OUT "dependson : bwa_merge_gatk.pl\n";
		print OUT "slots : 8\n";
		&printSlot("deepvariant");
		
		if ( $runs{deepvariant} == 1 ) {
			print OUT "run\n";
		}
		else {
			print OUT "#run\n";
		}
	} 
}



sub checkContamination {
	print OUT "
##############################################################
# checks contamination of a BAM file
##############################################################
";
	print OUT "pgr   : checkContamination.pl\n";
	print OUT "# bam file to check\n";
	print OUT "param : b : $outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam\n";
	print OUT "# settings\n";
	print OUT "param : se : $settings\n";
	print OUT "# output directory\n";
	print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/verifyBAM\n";
	print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
	print OUT "dependson : bwa_merge_gatk.pl,gem_aln.pl,gem_merge.pl,gem_single_sort.pl\n";
	&printSlot("checkContamination");

	if ( $runs{checkContamination} == 1 && $organism eq "human" )
	{ #check contamination currently only works for humans because we don't have a hapmap file for other species
		print OUT "run\n";
	}
	else {
		print OUT "#run\n";
	}
}

sub transcriptStats {

	print OUT "
##############################################################
# calculate transcript stats
##############################################################
";

	print OUT "pgr   : calcBaseCov_api.pl\n";
	print OUT "# target file name in config.xml\n";
	print OUT "param : t : transcript\n";
	print OUT "#get filename from database\n";
	print OUT "param : tt : db\n";
	print OUT "# outfile\n";
	print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/refSeqCoverage.txt\n";
	print OUT "# bamfile\n";
	print OUT "param : b : $outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam\n";
	print OUT "# settings\n";
	print OUT "param : se : $settings\n";
	print OUT "# wigfile\n";
	print OUT "param : w : $outdir/$project/$sample/$folder/$subfolder/refSeqCoverage.wig\n";
	print OUT "# coverage threshold\n";
	print OUT "param : c : 20\n";
	print OUT "#Sample name\n";
	print OUT "param : s : $sample\n";
	print OUT "#coverage per target file\n";
	print OUT "param : p : $outdir/$project/$sample/$folder/$subfolder/refSeqCoveragePerTarget.txt\n";
	print OUT "#coverage per target seg file\n";
	print OUT "param : ps : $outdir/$project/$sample/$folder/$subfolder/refSeqCoveragePerTarget.seg\n";
	print OUT "#coverage per transcript file\n";
	print OUT "param : r : $outdir/$project/$sample/$folder/$subfolder/refSeqCoveragePerTranscript.txt\n";
	print OUT "#target file is in bed format\n";
	print OUT "param : bed\n";
	print OUT "#calculate mean mapping quality for regions\n";
	print OUT "param : q\n";
	print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
	print OUT "dependson : bwa_merge_gatk.pl,gem_aln.pl,gem_merge.pl,gem_single_sort.pl\n";
	&printSlot("transcriptStats");

	if ( $runs{transcriptStats} == 1 && $organism eq "human" )
	{ 
		print OUT "run\n\n";
	}
	else {
		print OUT "#run\n\n";
	}

	print OUT "pgr   : insertTranscriptStats.pl\n";
	print OUT "#inputfile - coverage per transcript file\n";
	print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/refSeqCoveragePerTranscript.txt\n";
	print OUT "#Sample name\n";
	print OUT "param : s : $sample\n";
	print OUT "#delete old entries before inserting\n";
	print OUT "param : d\n";
	print OUT "#libtype and libpair\n";
	print OUT "param : lt : $libtype\n";
	print OUT "param : lp : $libpair\n";

	if ( $assay ne "" ) {
		print OUT "# assay\n";    #assay is empty for whole genome experiments
		print OUT "param : a : $assay\n";
	}
	print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
	print OUT "dependson : calcBaseCov_api.pl,gem_aln.pl,gem_merge.pl,gem_single_sort.pl\n";
	&printSlot("transcriptStats");
	
	if ( $runs{transcriptStats} == 1 && $organism eq "human" )
	{ 
		print OUT "run\n\n";
	}
	else {
		print OUT "#run\n\n";
	}

}





sub runSissrs {

		print OUT "
##############################################################
# run sissrs
##############################################################
";
	print OUT "pgr   : runSissrs.pl\n";
	print OUT "# infile\n";
	print OUT "param : c : $outdir/$project/$sample/$folder/$subfolder/merged.bam\n";
	print OUT "# settings\n";
	print OUT "param : se : $settings\n";
	print OUT "# output directory\n";
	print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/\n";
	print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
	print OUT "dependson : bwa_merge_gatk.pl\n";
	&printSlot("runSissrs");
	
	if ( $runs{runSissrs} == 1) {
		print OUT "run\n";
	}
	else {
		print OUT "#run\n";
	}

}


sub runMacs {

	print OUT "
##############################################################
# run MACS2
##############################################################
";
	print OUT "pgr   : runMACS2.pl\n";
	print OUT "# infile\n";
	print OUT "param : c : $outdir/$project/$sample/$folder/$subfolder/merged.bam\n";
	print OUT "# settings\n";
	print OUT "param : se : $settings\n";
	print OUT "# output directory\n";
	print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/\n";
	print OUT "# organism\n";
	print OUT "param : r : $organism\n";
	#print OUT "# use --broad peak option of MACS2\n";
	#print OUT "param : b\n";
	print OUT "# samplename\n";
	print OUT "param : sn : $sample\n";
	print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
	print OUT "dependson : bwa_merge_gatk.pl\n";
	&printSlot("runMacs");
	
	if ( $runs{runMacs} == 1) {
		print OUT "run\n";
	}
	else {
		print OUT "#run\n";
	}
	
	
	print OUT "
##############################################################
# calculate ChIP-Seq stats
##############################################################
";
	print OUT "pgr   : calcChIPstats.pl\n";
	print OUT "# infile\n";
	print OUT "param : c : $outdir/$project/$sample/$folder/$subfolder/merged.bam\n";
	print OUT "# settings\n";
	print OUT "param : se : $settings\n";
	print OUT "# output file\n";
	print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/ChIP.macs.stats\n";
	print OUT "# peak file\n";
	print OUT "param : p : $outdir/$project/$sample/$folder/$subfolder/ChIP.macs_peaks.narrowPeak\n";
	print OUT "# spp phantom output file\n";
	print OUT "param : s : $outdir/$project/$sample/$folder/$subfolder/spp_phantom_qc.txt\n";
	print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
	print OUT "dependson : runMACS2.pl\n";
	&printSlot("runMacs");
	
	if ( $runs{runMacs} == 1) {
		print OUT "run\n";
	}
	else {
		print OUT "#run\n";
	}
	
	print OUT "
##############################################################
# insert ChIP-Seq stats into database
##############################################################
";
	print OUT "pgr   : insertChIPstats.pl\n";
	print OUT "# stats file\n";
	print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/ChIP.macs.stats\n";
	print OUT "# use --broad peak option of MACS2\n";
	print OUT "param : b\n";
	print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
	print OUT "dependson : calcChIPstats.pl\n";
	&printSlot("runMacs");
	
	if ( $runs{runMacs} == 1) {
		print OUT "run\n";
	}
	else {
		print OUT "#run\n";
	}
	

}


sub cnvnator {

print OUT "
##############################################################
# run CNVnator
##############################################################
";
print OUT "pgr   : runCNVnator.pl\n";
print OUT "# infile\n";
print OUT "param : b : $outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam\n";
print OUT "# settings\n";
print OUT "param : se : $settings\n";
print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
print OUT "dependson : bwa_merge_gatk.pl\n";
&printSlot("cnvnator");

if ( $runs{cnvnator} == 1) {
	print OUT "run\n";
}
else {
	print OUT "#run\n";
}

		
}


sub annotateSNPs_sam {
	print OUT "
##############################################################
# annotating SNPs
##############################################################
";
	if ( $vcf > 0 ) {
		print OUT "#First, annotate known SNPs from dbSNP\n";
		print OUT "pgr   : annotatedbSNP.pl\n"
		  ;    #TW,Dec 2011: annotation of known dbSNPs now in extra script
		if($vcf == 1){
			print OUT "# infile\n";
			print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/ontarget.varfilter.vcf\n";
		}elsif($vcf == 2){
			print OUT "# infile\n";
			print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/gatk.ontarget.haplotypecaller.filtered.vcf\n";
		}
		print OUT "# settings\n";
		print OUT "param : se : $settings\n\n";
		print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
		print OUT "dependson : vcfsorter.pl,filterVCFforRegion.pl,varfilter_vcf.pl,varfilter_gatk.pl,filterGATK.pl,varfilter_rna.pl\n";
		&printSlot("annotateSNPs");

		if ( $runs{annotateSNPs} == 1 ) {
			print OUT "run\n";
		}
		else {
			print OUT "#run\n";
		}

		print OUT "pgr   : annotateVCF.pl\n"
		  ;    #TW,Oct 2011: annotation of snps and indels is now in one script
		if($vcf == 1){
			print OUT "# infile\n";
			print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/ontarget.varfilter.dbSNP.vcf\n";
			print OUT "# outfile \n";
			print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/ontarget.varfilter.dbSNP.plus.vcf\n";
		}elsif($vcf == 2){
			print OUT "# infile\n";
			print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/gatk.ontarget.haplotypecaller.filtered.dbSNP.vcf\n";
			print OUT "# outfile \n";
			print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/gatk.ontarget.haplotypecaller.filtered.dbSNP.plus.vcf\n";
		}
	}
	print OUT "# settings\n";
	print OUT "param : se : $settings\n";
	print OUT "# size of window around splice sites [5]\n";
	print OUT "param : w : 5\n";
	print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
	print OUT "dependson : vcfsorter.pl,filterVCFforRegion.pl,varfilter_vcf.pl,varfilter_gatk.pl,annotatedbSNP.pl,filterGATK.pl,varfilter_rna.pl\n";
	&printSlot("annotateSNPs");

	if ( $runs{annotateSNPs} == 1 ) {
		print OUT "run\n";
	}
	else {
		print OUT "#run\n";
	}
}


sub checkPileup {

	print OUT "
##############################################################
# check SNP quality by looking at the surrounding median quality
##############################################################
";
	print OUT "pgr   : filterSNPqual.pl\n";
	if($vcf == 1){
		print OUT "# input file\n";
		print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/ontarget.varfilter.dbSNP.plus.vcf\n";
		print OUT "# output file\n";
		print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/ontarget.varfilter.dbSNP.plus.checked.vcf\n";
	}elsif($vcf == 2){
		print OUT "# input file\n";
		print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/gatk.ontarget.haplotypecaller.filtered.dbSNP.plus.vcf\n";
		print OUT "# output file\n";
		print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/gatk.ontarget.haplotypecaller.filtered.dbSNP.plus.checked.vcf\n";
	}
	print OUT "# settings\n";
	print OUT "param : se : $settings\n";
	print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
	print OUT "dependson : annotateVCF.pl\n";
	&printSlot("checkPileup");
	
	if ( $runs{checkPileup} == 1 ) {
		print OUT "run\n";
	}
	else {
		print OUT "#run\n";
	}

}




sub snvdbInsertExome {
	
	if( $vcf == 2 && $libtype eq "genomic" && $target ne ""){
		
		print OUT "
##############################################################
# filter Whole Genome variants for target region
##############################################################
";
		print OUT "\n\npgr   : filterVCFforRegion.pl\n";

		print OUT "#infile\n";
		print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/gatk.ontarget.haplotypecaller.filtered.dbSNP.plus.checked.vcf\n";
		print OUT "# output file\n";
		print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/exome.gatk.ontarget.haplotypecaller.filtered.dbSNP.plus.checked.vcf\n";
		print OUT "# target file\n";
		print OUT "param : t : $targetMarginBED\n";
		print OUT "# settings\n";
		print OUT "param : se : $settings\n";
		print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
		print OUT "dependson : filterSNPqual.pl\n";
		&printSlot("varfilter");
		
		if ( $runs{snvdbInsertExome} == 1 ) {
			print OUT "run\n";
		}
		else {
			print OUT "#run\n";
		}
		
		
	}
	
	
	print OUT "
##############################################################
# insert variants and stats in database
##############################################################
";
	print OUT "#delete old entries before inserting new entries\n";
	print OUT "pgr   : snvdbExomeDelete.pl\n";
	print OUT "# sample name\n";
	print OUT "param : p : $sample\n";
	print OUT "# settings\n";
	print OUT "param : se : $settings\n";

	if ( $vcf == 1 ) {
		print OUT "# delete only entries generated by samtools\n";
		print OUT "param : c : samtools\n";
	}elsif ($vcf == 2 ) {
		print OUT "# delete only entries generated by gatk\n";
		print OUT "param : c : gatk\n";
	}
	print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n"
	  ; #start deleting when indel annotation is finished --> afterwards start inserting
	print OUT "dependson : annotateVCF.pl\n";
	&printSlot("snvdbInsertExome");
	
	if ( $runs{snvdbInsertExome} == 1 ) {
		print OUT "run\n\n";
	}
	else {
		print OUT "#run\n\n";
	}


	if ($libtype eq "genomic" && 
		$vcf == 2 && 
		$target eq "" && 
		(($organism eq "human" || $organism eq "mouse")) ) {				#execute snvdbExomeImport_vcf.pl instead of snvdbExomeInsert_vcf.pl only if sample is a genome and is inserted into genomegatk db
		print OUT "pgr   : snvdbExomeImport_vcf.pl\n";
		print OUT "# input file\n";
		print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/gatk.ontarget.haplotypecaller.filtered.dbSNP.plus.checked.vcf\n";
			
		print OUT "# caller\n";
		print OUT "param : c : gatk\n";
		
		print OUT "# settings\n";
		print OUT "param : se : $settings\n";
		print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
		print OUT "dependson : snvdbExomeDelete.pl,filterSNPqual.pl,filterVCFforRegion.pl\n";
		&printSlot("snvdbInsertExome");
		if ( $runs{snvdbInsertExome} == 1 ) {
			print OUT "run\n\n";
		}
		else {
			print OUT "#run\n\n";
		}
	} else {
		print OUT "pgr   : snvdbExomeInsert_vcf.pl\n";
		if ($vcf == 1){
			print OUT "# input file\n";
			print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/ontarget.varfilter.dbSNP.plus.checked.vcf\n";
		}elsif ($vcf == 2){
			print OUT "# input file\n";
			if($libtype eq "genomic" && $target ne ""){
				print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/exome.gatk.ontarget.haplotypecaller.filtered.dbSNP.plus.checked.vcf\n";
			}else{
				print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/gatk.ontarget.haplotypecaller.filtered.dbSNP.plus.checked.vcf\n";
			}
			
			print OUT "# caller\n";
			print OUT "param : c : gatk\n";
		}
		
		print OUT "# settings\n";
		print OUT "param : se : $settings\n";
		print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
		print OUT "dependson : snvdbExomeDelete.pl,filterSNPqual.pl,filterVCFforRegion.pl\n";
		&printSlot("snvdbInsertExome");
		if ( $runs{snvdbInsertExome} == 1 ) {
			print OUT "run\n\n";
		}
		else {
			print OUT "#run\n\n";
		}
	}

	print OUT "pgr   : insertSizeIntoDB.pl\n";
	print OUT "# bamfile to calculate insert size\n";
	print OUT "param : b : $outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam\n";
	print OUT "#sample name\n";
	print OUT "param : s : $sample\n";
	print OUT "# libtype\n";
	print OUT "param : l : $libtype\n";
	print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
	print OUT "dependson : bwa_merge_gatk.pl,gem_aln.pl,gem_merge.pl,gem_single_sort.pl\n";
	&printSlot("insertSizeIntoDB");

	if ( $runs{insertSizeIntoDB} ) {
		print OUT "run\n\n";
	}
	else {
		print OUT "#run\n\n";
	}

}







sub callChrM {
	if ( $vcf == 2 ) {
		print OUT "
##############################################################
# run GATK HaplotypeCaller on chrM
##############################################################
";

		print OUT "pgr   : varfilter_gatk.pl\n";
		print OUT "#infile\n";
		if($libtype eq "genomic" && !(-e "$outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam.bai") ){
			print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/\n";
		}else{
			print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam\n";
		}
		print OUT "# output file\n";
    	print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/gatk.chrM.haplotypecaller.vcf\n";
		print OUT "# target region\n";
		print OUT "param : r : chrM\n";
		#print OUT "param : fa\n";
		

		print OUT "# settings\n";
		print OUT "param : se : $settings\n";
		&printSlot("callChrM");
		if($libtype eq "genomic" && !(-e "$outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam.bai") ){
			print OUT "dependson : bwa_aln_all_Formats.pl,gem_aln.pl,gem_merge.pl,gem_single_sort.pl\n";
		}else{
			print OUT "dependson : bwa_merge_gatk.pl,gem_aln.pl,gem_merge.pl,gem_single_sort.pl\n";
		}

		if ( $runs{callChrM} == 1 ) {
			print OUT "run\n";
		}
		else {
			print OUT "#run\n";
		}
	}
	
	
	
	
}
sub annotateChrM {
	if ( $vcf == 2 ) {
		print OUT "
##############################################################
# annotate chrM VCF file
##############################################################
";

	print OUT "pgr   : annotateVCF.pl\n";
	print OUT "# infile\n";
	print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/gatk.chrM.haplotypecaller.vcf\n";
	print OUT "# outfile \n";
	print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/gatk.chrM.haplotypecaller.plus.vcf\n";
				

	print OUT "# settings\n";
	print OUT "param : se : $settings\n";
	print OUT "# use chrM gencode table\n";
	print OUT "param : t : gencodeV20chrM\n";
	print OUT "# use chrM gencode cds table\n";
	print OUT "param : c : gencodeV20chrM_cds\n";
	print OUT "# size of window around splice sites [5]\n";
	print OUT "param : w : 5\n";
	print OUT "# compare with additional SNP tables, default hg19\n";
	print OUT "param : p : $projectdatabase\n";
	print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
	print OUT "dependson : varfilter_vcf.pl,varfilter_gatk.pl,annotatedbSNP.pl,vcfsorter.pl,filterGATK.pl,varfilter_rna.pl\n";
	&printSlot("annotateChrM");

	if ( $runs{annotateChrM} == 1 ) {
		print OUT "run\n";
	}
	else {
		print OUT "#run\n";
	}
	
	}

	
	
	
	
}

sub insertChrM {
	print OUT "
##############################################################
# insert chrM VCF file into DB
##############################################################
";
	if ($libtype eq "genomic" && 
		$vcf == 2 && 
		$target eq "" && 
		(($organism eq "human" || $organism eq "mouse")) ) {				#execute snvdbExomeImport_vcf.pl instead of snvdbExomeInsert_vcf.pl only if sample is a genome and is inserted into genomegatk db
		print OUT "pgr   : snvdbExomeImport_vcf.pl\n";
		print OUT "# input file\n";
		print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/gatk.chrM.haplotypecaller.plus.vcf\n";
		print OUT "# caller\n";
		print OUT "param : c : gatk\n";
		print OUT "# settings\n";
		print OUT "param : se : $settings\n";
		print OUT "# delete old entries before inserting\n";
		print OUT "param : dc : chrM\n";
		print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
		print OUT "dependson : snvdbExomeDelete.pl,annotateVCF.pl\n";
		&printSlot("insertChrM");
		if ( $runs{insertChrM} == 1 ) {
			print OUT "run\n\n";
		}
		else {
			print OUT "#run\n\n";
		}
	} elsif ( $vcf == 2 ) {
		print OUT "pgr   : snvdbExomeInsert_vcf.pl\n";
		print OUT "# input file\n";
		print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/gatk.chrM.haplotypecaller.plus.vcf\n";
		print OUT "# caller\n";
		print OUT "param : c : gatk\n";
		print OUT "# settings\n";
		print OUT "param : se : $settings\n";
		print OUT "# delete old entries before inserting\n";
		print OUT "param : dc : chrM\n";
		print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
		print OUT "dependson : snvdbExomeDelete.pl,annotateVCF.pl\n";
		&printSlot("insertChrM");
	
		if ( $runs{insertChrM} == 1 ) {
			print OUT "run\n";
		}
		else {
			print OUT "#run\n";
		}
	}
}



sub exomeDepth {
	if ( $vcf >= 1  && ($assay ne "" || $libtype eq "genomic" )) {			#TW 18.07.2013 --> should be run only for VCF pipeline, because the CNV table is only available in the exomevcf database # RB 05.10.2016 BAD HARDCODING, genomic / genomicmouse Exomedepth configuration should go into settings
		print OUT "
##############################################################
# run exomeDepth
##############################################################
";
		print OUT "pgr   : exomeDepth.pl\n";
		if($libtype eq "genomic" && !(-e "$outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam.bai") ){
			print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/\n";
		}else{
			print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam\n";
		}
		print OUT "# settings\n";
		print OUT "param : se : $settings\n";
		print OUT "# output directory\n";
		print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/\n";
		print OUT "# sample name\n";
		print OUT "param : s : $sample\n";
		print OUT "# assay\n";
		if($assay ne ""){
			print OUT "param : a : $assay\n";
		}else{
			my $genomicassay=( $organism eq "mouse" ? "genomicmouse" : $libtype );
			print OUT "param : a : $genomicassay\n";
		}
		
		if($libtype eq "genomic" && !(-e "$outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam.bai") ){
			print OUT "dependson : bwa_aln_all_Formats.pl\n";
		}else{
			print OUT "dependson : bwa_merge_gatk.pl\n";
		}
		&printSlot("exomeDepth");
		
		if ( $runs{exomeDepth} == 1) {
			print OUT "run\n";
		}
		else {
			print OUT "#run\n";
		}
	
		
		print OUT "
##############################################################
# insert exomeDepth in snv table of database
##############################################################
";
		print OUT "#delete old entries before inserting new entries\n";
		print OUT "pgr   : snvdbExomeDelete.pl\n";
		print OUT "# sample name\n";
		print OUT "param : p : $sample\n";
		print OUT "# settings\n";
		print OUT "param : se : $settings\n";
		print OUT "# delete only entries generated by exomedepth\n";
		print OUT "param : c : exomedepth\n";
		print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
		print OUT "dependson : exomeDepth.pl\n";
		&printSlot("exomeDepthInsert");
		if ( $runs{exomeDepthInsert} == 1) {					
			print OUT "run\n\n";
		}
		else {
			print OUT "#run\n\n";
		}
		
		if ($libtype eq "genomic" && 
			$vcf == 2 && 
			$target eq "" && 
			(($organism eq "human" || $organism eq "mouse")) ) {				#execute snvdbExomeImport_vcf.pl instead of snvdbExomeInsert_vcf.pl only if sample is a genome and is inserted into genomegatk db
			print OUT "pgr   : snvdbExomeImport_vcf.pl\n";
			print OUT "# input file\n";
			print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/exomeDepth.all.vcf\n";
			print OUT "# caller\n";
			print OUT "param : c : exomedepth\n";
			print OUT "# class of inserted variants\n";
			print OUT "param : l : cnv\n";
			print OUT "# settings\n";
			print OUT "param : se : $settings\n";
			print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
			print OUT "dependson : exomeDepth.pl,snvdbExomeDelete.pl\n";
			&printSlot("exomeDepthInsert");
			if ( $runs{exomeDepthInsert} == 1 ) {
				print OUT "run\n\n";
			}
			else {
				print OUT "#run\n\n";
			}
		} else {	
			print OUT "\npgr   : snvdbExomeInsert_vcf.pl\n";
			print OUT "# input file\n";
			print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/exomeDepth.all.vcf\n";
			print OUT "# settings\n";
			print OUT "param : se : $settings\n";
			print OUT "# class of inserted variants\n";
			print OUT "param : l : cnv\n";
			print OUT "# delete only entries generated by exomedepth\n";
			print OUT "param : c : exomedepth\n";
			print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
			print OUT "dependson : exomeDepth.pl,snvdbExomeDelete.pl\n";
			&printSlot("exomeDepthInsert");
			if ( $runs{exomeDepthInsert} == 1) {
				print OUT "run\n\n";
			}
			else {
				print OUT "#run\n\n";
			}
		}
		
		if ($libtype eq "genomic" && 
			$vcf == 2 && 
			$target eq "" && 
			(($organism eq "human" || $organism eq "mouse")) ) {				#execute snvdbExomeImport_vcf.pl instead of snvdbExomeInsert_vcf.pl only if sample is a genome and is inserted into genomegatk db
			print OUT "pgr   : snvdbExomeImport_vcf.pl\n";
			print OUT "# input file\n";
			print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/exomeDepth.chrX.all.vcf\n";
			print OUT "# caller\n";
			print OUT "param : c : exomedepth\n";
			print OUT "# class of inserted variants\n";
			print OUT "param : l : cnv\n";
			print OUT "# settings\n";
			print OUT "param : se : $settings\n";
			print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
			print OUT "dependson : exomeDepth.pl,snvdbExomeDelete.pl\n";
			&printSlot("exomeDepthInsert");
			if ( $runs{exomeDepthInsert} == 1 ) {
				print OUT "run\n\n";
			}
			else {
				print OUT "#run\n\n";
			}
		} else {	
			print OUT "\npgr   : snvdbExomeInsert_vcf.pl\n";
			print OUT "# input file\n";
			print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/exomeDepth.chrX.all.vcf\n";
			print OUT "# settings\n";
			print OUT "param : se : $settings\n";
			print OUT "# class of inserted variants\n";
			print OUT "param : l : cnv\n";
			print OUT "# delete only entries generated by exomedepth\n";
			print OUT "param : c : exomedepth\n";
			print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
			print OUT "dependson : exomeDepth.pl,snvdbExomeDelete.pl\n";
			&printSlot("exomeDepthInsert");
			if ( $runs{exomeDepthInsert} == 1) {
				print OUT "run\n\n";
			}
			else {
				print OUT "#run\n\n";
			}
		}
	}
}




sub snpEff {
if ( $vcf > 0 ) {
	print OUT "
########################################################################
# annotate SNVs/ small Indels with snpEff and insert annotations into DB
########################################################################
";
	
	print OUT "pgr   : runsnpEff.pl\n";
	if ($vcf == 1){
		print OUT "# input file\n";
		print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/ontarget.varfilter.dbSNP.plus.checked.vcf\n";
	}elsif ($vcf == 2){
		print OUT "# input file\n";
		if($libtype eq "genomic" && $target ne ""){
			print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/exome.gatk.ontarget.haplotypecaller.filtered.dbSNP.plus.checked.vcf\n";
		}else{
			print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/gatk.ontarget.haplotypecaller.filtered.dbSNP.plus.checked.vcf\n";
		}
	}
	
	print OUT "# settings\n";
	print OUT "param : se : $settings\n";
	print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
	print OUT "dependson : snvdbExomeInsert_vcf.pl,snvdbExomeImport_vcf.pl,filterVCFforRegion.pl\n";
	&printSlot("snpEff");

	if ( $runs{snpEff} == 1 ) # 13.02.2014: now also for mice; && $organism eq "human" ) #currently only run for humans --> other species don't have the refSeq tables (yet)
	{ 
		print OUT "run\n";
	}
	else {
		print OUT "#run\n";
	}
	
	print OUT "\n\npgr   : insertsnpEff.pl\n";
	if ($vcf == 1){
		print OUT "# input file\n";
		print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/ontarget.varfilter.dbSNP.plus.checked.snpEff.vcf\n";
	}elsif ($vcf == 2){
		print OUT "# input file\n";
		if($libtype eq "genomic" && $target ne ""){
			print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/exome.gatk.ontarget.haplotypecaller.filtered.dbSNP.plus.checked.snpEff.vcf\n";
		}else{
			print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/gatk.ontarget.haplotypecaller.filtered.dbSNP.plus.checked.snpEff.vcf\n";
		}
	}
	print OUT "# settings\n";
	print OUT "param : se : $settings\n";
	print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
	print OUT "dependson : runsnpEff.pl\n";
	&printSlot("snpEffInsert");

	if ( $runs{snpEffInsert} == 1 ) # 13.02.2014: now also for mice; && $organism eq "human" ) #currently only run for humans --> other species don't have the refSeq tables (yet)
	{ 
		print OUT "run\n";
	}
	else {
		print OUT "#run\n";
	}
}
}

#################################################################################
#
# TW, 06.11.2013; Note that this script is part of the "insertDB" module of runs
#                 since it requires the variants in the database to work properly
#
#################################################################################
sub homozygosity {

	print OUT "
##############################################################
# calculate stretches of homozyogosity and insert them into DB
##############################################################
";
	print OUT "pgr   : homozygosity.pl\n";
	print OUT "# sample name \n";
	print OUT "param : s : $sample\n";
	print OUT "# settings\n";
	print OUT "param : se : $settings\n";
	if ($vcf == 1){
		print OUT "# input file\n";
		print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/ontarget.varfilter.dbSNP.plus.checked.vcf\n";
	}elsif ($vcf == 2){
		print OUT "# input file\n";
		print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/gatk.ontarget.haplotypecaller.filtered.dbSNP.plus.checked.vcf\n";
	}
	print OUT "# outfile \n";
	print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/homozygosity.out\n";
	print OUT "# insert into db \n";
	print OUT "param : insert\n"; 
	print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
	print OUT "dependson : snvdbExomeInsert_vcf.pl,snvdbExomeImport_vcf.pl\n";
	&printSlot("homozygosity");

	if ( $runs{homozygosity} == 1 )
	{
		print OUT "run\n";
	}
	else {
		print OUT "#run\n";
	}
	
	
}


#################################################################################
sub updateRefSeq {
	print OUT "
##############################################################
# Annotate SNPs and Indels again with refSeq tables instead
# of UCSC tables and update SNVs with it.
##############################################################
";

	print OUT "
##############################################################
# annotating SNPs
##############################################################
";
	if ($vcf == 1){
		print OUT "pgr   : annotateVCF.pl\n"
		  ;    #TW,Oct 2011: annotation of snps and indels is now in one script
		print OUT "# infile\n";
		print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/ontarget.varfilter.dbSNP.vcf\n";
		print OUT "# outfile \n";
		print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/ontarget.varfilter.dbSNP.refSeq.vcf\n";
	}elsif ($vcf == 2){
		print OUT "pgr   : annotateVCF.pl\n"		  ;    
		print OUT "# infile\n";
		
		print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/gatk.ontarget.haplotypecaller.filtered.dbSNP.vcf\n";
		print OUT "# outfile \n";
		print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/gatk.ontarget.haplotypecaller.filtered.dbSNP.refSeq.vcf\n";
				
	}
	print OUT "# settings\n";
	print OUT "param : se : $settings\n";
	print OUT "# use refSeq tables\n";
	print OUT "param : r\n";
	print OUT "# size of window around splice sites [5]\n";
	print OUT "param : w : 5\n";
	print OUT "# compare with additional SNP tables, default hg19\n";
	print OUT "param : p : $projectdatabase\n";
	print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
	print OUT "dependson : varfilter_vcf.pl,varfilter_gatk.pl,annotatedbSNP.pl,vcfsorter.pl,filterGATK.pl,varfilter_rna.pl\n";
	&printSlot("RefSeq");

	if ( $runs{RefSeq} == 1 && $organism eq "human" )
	{ #currently only run for humans --> other species don't have the refSeq tables (yet)
		print OUT "run\n";
	}
	else {
		print OUT "#run\n";
	}


	print OUT "
##############################################################
# update snv table
##############################################################
";

	
		
	if( $vcf == 2 && $libtype eq "genomic" && $target ne ""){
	
		print OUT "
		##############################################################
		# filter Whole Genome variants for target region
		##############################################################
		";
			print OUT "\n\npgr   : filterVCFforRegion.pl\n";

			print OUT "#infile\n";
			print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/gatk.ontarget.haplotypecaller.filtered.dbSNP.refSeq.vcf\n";
			print OUT "# output file\n";
			print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/exome.gatk.ontarget.haplotypecaller.filtered.dbSNP.refSeq.vcf\n";
			print OUT "# target file\n";
			print OUT "param : t : $targetMarginBED\n";
			print OUT "# settings\n";
			print OUT "param : se : $settings\n";
			print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
			print OUT "dependson : annotateVCF.pl,snvdbExomeInsert_vcf.pl,snvdbExomeImport_vcf.pl\n";
			&printSlot("varfilter");
		
			if ( $runs{RefSeq} == 1 ) {
				print OUT "run\n";
			}
			else {
				print OUT "#run\n";
			}	
		
	}
	

	print OUT "pgr   : updateVariants_vcf.pl\n";
	if ($vcf == 1){
		print OUT "# infile\n";
		print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/ontarget.varfilter.dbSNP.refSeq.vcf\n";
	}elsif ($vcf == 2){
		print OUT "# infile\n";
		if($libtype eq "genomic" && $target ne ""){
			print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/exome.gatk.ontarget.haplotypecaller.filtered.dbSNP.refSeq.vcf\n";
		}else{
			print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/gatk.ontarget.haplotypecaller.filtered.dbSNP.refSeq.vcf\n";
		}
							
	}
	print OUT "# column to update\n";
	print OUT "param : c : refseq\n";
	print OUT "# settings\n";
	print OUT "param : se : $settings\n";
	print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
	print OUT "dependson : annotateVCF.pl,snvdbExomeInsert_vcf.pl,snvdbExomeImport_vcf.pl\n";
	&printSlot("updateRefSeq");

	if ( $runs{updateRefSeq} == 1 && $organism eq "human" )
	{ #currently only run for humans --> other species don't have the refSeq tables (yet)
		print OUT "run\n";
	}
	else {
		print OUT "#run\n";
	}


}
#########################################################################################
sub updatelncRNA {
	if ( $vcf != 0 ) {
		print OUT "
##############################################################
# Annotate SNPs and Indels again with lncRNA tables instead
# of UCSC tables and update SNVs with it.
##############################################################
";

		print OUT "
##############################################################
# annotating SNPs and Indels
##############################################################
";

		print OUT "pgr   : annotateVCF.pl\n"
		  ;    #TW,Oct 2011: annotation of snps and indels is now in one script
		if ($vcf == 1){
			print OUT "# infile\n";
			print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/ontarget.varfilter.dbSNP.vcf\n";
			print OUT "# outfile \n";
			print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/ontarget.varfilter.dbSNP.lncRNA.vcf\n";
		}elsif ($vcf == 2){
			print OUT "# infile\n";
			print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/gatk.ontarget.haplotypecaller.filtered.dbSNP.vcf\n";
			print OUT "# outfile \n";
			print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/gatk.ontarget.haplotypecaller.filtered.dbSNP.lncRNA.vcf\n";
					
		}
		print OUT "# settings\n";
		print OUT "param : se : $settings\n";
		print OUT "# use lncRNA tables\n";
		print OUT "param : l\n";
		print OUT "# annotate lncrna function\n";
		print OUT "param : f  : lincrna\n";
		print OUT "# size of window around splice sites [5]\n";
		print OUT "param : w : 5\n";
		print OUT "# compare with additional SNP tables, default hg19\n";
		print OUT "param : p : $projectdatabase\n";
		print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
		print OUT "dependson : varfilter_vcf.pl,varfilter_gatk.pl,annotatedbSNP.pl,filterGATK.pl,varfilter_rna.pl\n";
		&printSlot("lncRNA");

		if ( $runs{lncRNA} == 1 && $organism eq "human" )
		{ #currently only run for humans --> other species don't have the lncrna tables (yet)
			print OUT "run\n";
		}
		else {
			print OUT "#run\n";
		}
		
		if( $vcf == 2 && $libtype eq "genomic" && $target ne ""){
		
			print OUT "
##############################################################
# filter Whole Genome variants for target region
##############################################################
	";
			print OUT "\n\npgr   : filterVCFforRegion.pl\n";
	
			print OUT "#infile\n";
			print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/gatk.ontarget.haplotypecaller.filtered.dbSNP.lncRNA.vcf\n";
			print OUT "# output file\n";
			print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/exome.gatk.ontarget.haplotypecaller.filtered.dbSNP.lncRNA.vcf\n";
			print OUT "# target file\n";
			print OUT "param : t : $targetMarginBED\n";
			print OUT "# settings\n";
			print OUT "param : se : $settings\n";
			print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
			print OUT "dependson : annotateVCF.pl,snvdbExomeInsert_vcf.pl,snvdbExomeImport_vcf.pl\n";
			&printSlot("varfilter");
			
			if ( $runs{lncRNA} == 1 ) {
				print OUT "run\n";
			}
			else {
				print OUT "#run\n";
			}
			
			
		}

		print OUT "
##############################################################
# update snv table
##############################################################
";
		print OUT "pgr   : updateVariants_vcf.pl\n";
		if ($vcf == 1){
			print OUT "# infile\n";
			print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/ontarget.varfilter.dbSNP.lncRNA.vcf\n";
		}elsif ($vcf == 2){
			print OUT "# infile\n";
			if($libtype eq "genomic" && $target ne ""){
				print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/exome.gatk.ontarget.haplotypecaller.filtered.dbSNP.lncRNA.vcf\n";
			}else{
				print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/gatk.ontarget.haplotypecaller.filtered.dbSNP.lncRNA.vcf\n";
			}
					
		}
		print OUT "# column to update\n";
		print OUT "param : c : lincrna\n";
		print OUT "# update function column\n";
		print OUT "param : f\n";
		print OUT "# settings\n";
		print OUT "param : se : $settings\n";
		print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
		print OUT "dependson : annotateVCF.pl,snvdbExomeInsert_vcf.pl,snvdbExomeImport_vcf.pl,filterVCFforRegion.pl\n";
		&printSlot("updatelncRNA");

		if ( $runs{updatelncRNA} == 1 && $organism eq "human" )
		{ #currently only run for humans --> other species don't have the refSeq tables (yet)
			print OUT "run\n";
		}
		else {
			print OUT "#run\n";
		}
	}

}

sub updatemiRNA {
		if ( $vcf != 0 ) {
		print OUT "
##############################################################
# Annotate SNPs and Indels again with miRNA tables instead
# of UCSC tables and update SNVs with it.
##############################################################
";

		print OUT "
##############################################################
# annotating SNPs and Indels
##############################################################
";

		print OUT "pgr   : annotateVCF.pl\n";
		if ($vcf == 1){
			print OUT "# infile\n";
			print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/ontarget.varfilter.dbSNP.vcf\n";
			print OUT "# outfile \n";
			print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/ontarget.varfilter.dbSNP.miRNA.vcf\n";
		}elsif ($vcf == 2){
			print OUT "# infile\n";
			print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/gatk.ontarget.haplotypecaller.filtered.dbSNP.vcf\n";
			print OUT "# outfile \n";
			print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/gatk.ontarget.haplotypecaller.filtered.dbSNP.miRNA.vcf\n";
					
		}
		print OUT "# settings\n";
		print OUT "param : se : $settings\n";
		print OUT "# use miRNA tables\n";
		print OUT "param : mir\n";
		print OUT "# annotate miRNA function\n";
		print OUT "param : f  : mirna\n";
		print OUT "# size of window around splice sites [0]\n";
		print OUT "param : w : 0\n";
		print OUT "# compare with additional SNP tables, default hg19\n";
		print OUT "param : p : $projectdatabase\n";
		print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
		print OUT "dependson : varfilter_vcf.pl,varfilter_gatk.pl,annotatedbSNP.pl,filterGATK.pl,varfilter_rna.pl\n";
		&printSlot("miRNA");

		if ( $runs{miRNA} == 1 && $organism eq "human" )
		{ #currently only run for humans --> other species don't have the miRNA tables (yet)
			print OUT "run\n";
		}
		else {
			print OUT "#run\n";
		}
		
		if( $vcf == 2 && $libtype eq "genomic" && $target ne ""){
		
			print OUT "
##############################################################
# filter Whole Genome variants for target region
##############################################################
	";
			print OUT "\n\npgr   : filterVCFforRegion.pl\n";
	
			print OUT "#infile\n";
			print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/gatk.ontarget.haplotypecaller.filtered.dbSNP.miRNA.vcf\n";
			print OUT "# output file\n";
			print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/exome.gatk.ontarget.haplotypecaller.filtered.dbSNP.miRNA.vcf\n";
			print OUT "# target file\n";
			print OUT "param : t : $targetMarginBED\n";
			print OUT "# settings\n";
			print OUT "param : se : $settings\n";
			print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
			print OUT "dependson : annotateVCF.pl,snvdbExomeInsert_vcf.pl,snvdbExomeImport_vcf.pl\n";
			&printSlot("varfilter");
			
			if ( $runs{miRNA} == 1 && $organism eq "human" ) {
				print OUT "run\n";
			}
			else {
				print OUT "#run\n";
			}
			
			
		}

		print OUT "
##############################################################
# update snv table
##############################################################
";
		print OUT "pgr   : updateVariants_vcf.pl\n";
		if ($vcf == 1){
			print OUT "# infile\n";
			print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/ontarget.varfilter.dbSNP.miRNA.vcf\n";
		}elsif ($vcf == 2){
			print OUT "# infile\n";
			if($libtype eq "genomic" && $target ne ""){
				print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/exome.gatk.ontarget.haplotypecaller.filtered.dbSNP.miRNA.vcf\n";
			}else{
				print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/gatk.ontarget.haplotypecaller.filtered.dbSNP.miRNA.vcf\n";
			}
					
		}
		print OUT "# column to update\n";
		print OUT "param : c : mirna\n";
		print OUT "# update function column\n";
		print OUT "param : f\n";
		print OUT "# settings\n";
		print OUT "param : se : $settings\n";
		print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
		print OUT "dependson : annotateVCF.pl,snvdbExomeInsert_vcf.pl,snvdbExomeImport_vcf.pl,filterVCFforRegion.pl\n";
		&printSlot("updatemiRNA");

		if ( $runs{updatemiRNA} == 1 && $organism eq "human" )
		{ #currently only run for humans --> other species don't have the miRNA tables (yet)
			print OUT "run\n";
		}
		else {
			print OUT "#run\n";
		}
	}
}




sub runPindel {
	if ( $vcf > 0 ) {
		
		print OUT "
##############################################################
# create Pindel config file
##############################################################
";

		print OUT "pgr   : createPindelCfg.pl\n";
		print OUT "# bam file to check\n";
		print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam\n";
		print OUT "# output file\n";
		print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/pindel.cfg\n";
		print OUT "# sample name\n";
		print OUT "param : s : $sample\n";
		print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
		print OUT "dependson : bwa_merge_gatk.pl,collectMetrics.pl,gem_aln.pl,gem_merge.pl,gem_single_sort.pl\n";
		
		&printSlot("runPindel");
		
		if ( $runs{runPindel} == 1  ) { #&& $vcf == 1) {
			print OUT "run\n";
		}
		else {
			print OUT "#run\n";
		}
		
		
		
		print OUT "
##############################################################
# run Pindel
##############################################################
";
		make_path("$outdir/$project/$sample/$folder/$subfolder/Pindel", {				#create directory for single files
				mode => 0775
			}) if ( $runs{runPindel} == 1  ); #{ && $vcf == 1);
		#print STDERR "make_path\n";
		print OUT "pgr   : runPindel.pl\n";
		print OUT "# settings\n";
		print OUT "param : se : $settings\n";
		print OUT "# output directory\n";
		print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/Pindel\n";
		print OUT "# config file\n";
		print OUT "param : c :  $outdir/$project/$sample/$folder/$subfolder/pindel.cfg\n";
		
		print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
		if($libtype eq "genomic" && !(-e "$outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam.bai") ){
			print OUT "dependson : createPindelCfg.pl,gem_aln.pl,gem_merge.pl,gem_single_sort.pl\n";
		}else{
			print OUT "dependson : createPindelCfg.pl,collectMetrics.pl,gem_aln.pl,gem_merge.pl,gem_single_sort.pl\n";
		}
		
		if($libtype eq "mtDNA" )
		{
			print OUT "param : r  : chrM\n";
			print OUT "param : he : 0.01\n";
			print OUT "param : ho : 0.99\n";
			print OUT "#signal for parallelpipeline.pl to start this job as an array job on Grid Engine\n";
			print OUT "slots : 1\n";	
		}
		else
		{
			print OUT "# start this job as an array job, i.e. one job will be started by chromosome\n";
			print OUT "param : aj\n";
			print OUT "#signal for parallelpipeline.pl to start this job as an array job on Grid Engine\n";
			print OUT "slots : perChrArray\n";
		}
		
		
		if ( $runs{runPindel} == 1  ) { #&& $vcf == 1) {
			print OUT "run\n";
		}
		else {
			print OUT "#run\n";
		}
		
print OUT "
##############################################################
# concat single chromosome pindel VCFs
##############################################################
";
		print OUT "pgr   : concatVCF.pl\n";
		print OUT "#indir\n";
		print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/Pindel\n";
		print OUT "#outfile\n";
		print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/all.pindel.vcf\n";
		print OUT "#files to look for\n";
		print OUT "param : e : pindel.vcf\n";		
		print OUT "# settings from config file\n";
		print OUT "param : se : $settings\n";
		print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
		print OUT "dependson : runPindel.pl\n";
		&printSlot("runPindel");
		
		if ( $runs{runPindel} == 1 ) { #&& $vcf == 1) {
			print OUT "run\n";
		}
		else {
			print OUT "#run\n";
		}
		
		
		
				print OUT "
##############################################################
# annotate Pindel
##############################################################
";
		print OUT "#First, annotate known SNPs from dbSNP\n";
		print OUT "pgr   : annotatedbSNP.pl\n";
		print OUT "#infile\n";
		print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/all.pindel.vcf\n";
		print OUT "# settings\n";
		print OUT "param : se : $settings\n\n";
		print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
		print OUT "dependson : concatVCF.pl,runPindel.pl\n";
		&printSlot("annotatePindel");

		if ( $runs{annotatePindel} == 1 ) {
			print OUT "run\n";
		}
		else {
			print OUT "#run\n";
		}

		print OUT "\npgr   : annotateVCF.pl\n";
		
		print OUT "#infile\n";
		print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/all.pindel.dbSNP.vcf\n";
		print OUT "# output file\n";
		print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/all.pindel.dbSNP.plus.vcf\n";
		print OUT "# settings\n";
		print OUT "param : se : $settings\n";
		print OUT "# size of window around splice sites [5]\n";
		print OUT "param : w : 5\n";
		print OUT "# compare with additional SNP tables, default hg19\n";
		print OUT "param : p : $projectdatabase\n";
		print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
		print OUT "dependson : concatVCF.pl,runPindel.pl,annotatedbSNP.pl\n";
		&printSlot("annotatePindel");

		if ( $runs{annotatePindel} == 1 ) {
			print OUT "run\n";
		}
		else {
			print OUT "#run\n";
		}
		
		
		print OUT "
##############################################################
# filter Pindel
##############################################################
";		

		
		print OUT "pgr   : filterVCFforRegion.pl\n";
		print OUT "#infile\n";
		print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/all.pindel.dbSNP.plus.vcf\n";


		if ( $libtype eq "genomic" || $libtype eq "mtDNA" ) {
			print OUT "# output file\n";
			print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/genome.all.pindel.dbSNP.plus.vcf\n";
		}
		elsif ( $target ne "" ) {
			print OUT "# target region\n";
			print OUT "param : t : $targetMarginBED\n";
			print OUT "# output file\n";
			print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/exome.all.pindel.dbSNP.plus.vcf\n";
		} else {
			print OUT "# output file\n";
			print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/exome.all.pindel.dbSNP.plus.vcf\n";
		}
		print OUT "#print overlapping variants only once\n";
		print OUT "param : u\n";
		print OUT "# settings\n";
		print OUT "param : se : $settings\n";
		
		print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
		print OUT "dependson : concatVCF.pl,annotateVCF.pl\n";
		&printSlot("filterPindel");
		if ( $runs{filterPindel} == 1 ) {
			print OUT "run\n";
		}
		else {
			print OUT "#run\n";
		}
		


		print OUT "\npgr   : filterVCFforRegion.pl\n";
		if($vcf == 1){
			print OUT "# samtools vcf file to filter\n";
			print OUT "param : t : $outdir/$project/$sample/$folder/$subfolder/ontarget.varfilter.vcf\n";
		}else{
			print OUT "# samtools vcf file to filter\n";
			if ( $libtype eq "mtDNA" ){
				print OUT "param : t : $outdir/$project/$sample/$folder/$subfolder/gatk.chrM.haplotypecaller.vcf\n";
			}else{
				print OUT "param : t : $outdir/$project/$sample/$folder/$subfolder/gatk.ontarget.haplotypecaller.vcf\n";
			}
			print OUT "#minimize Pindel indels before filtering\n";
			print OUT "param : m\n";
		}
		if($libtype eq "genomic" || $libtype eq "mtDNA"){
			print OUT "#infile\n";
			print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/genome.all.pindel.dbSNP.plus.vcf\n";
			print OUT "# output file\n";
			print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/filtered.genome.all.pindel.dbSNP.plus.vcf\n";
		}else{
			print OUT "#infile\n";
			print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/exome.all.pindel.dbSNP.plus.vcf\n";
			print OUT "# output file\n";
			print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/filtered.exome.all.pindel.dbSNP.plus.vcf\n";
			
		}
		print OUT "#fraction that needs to overlap\n";
		print OUT "param : f : 0.8\n";
		print OUT "#variants in the variant file should be filtered out\n";
		print OUT "param : v\n";
		print OUT "# settings\n";
		print OUT "param : se : $settings\n";
		print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
		print OUT "dependson : concatVCF.pl,filterVCFforRegion.pl,varfilter_gatk.pl,varfilter_vcf.pl,runPindel.pl,vcfsorter.pl,varfilter_rna.pl\n";
		&printSlot("filterPindel");
		if ( $runs{filterPindel} == 1 ) {
			print OUT "run\n";
		}
		else {
			print OUT "#run\n";
		}
	
		
		



		print OUT "
##############################################################
# insert variants in database
##############################################################
";
		print OUT "#delete old entries before inserting new entries\n";
		print OUT "pgr   : snvdbExomeDelete.pl\n";
		print OUT "# sample name\n";
		print OUT "param : p : $sample\n";
		print OUT "# settings\n";
		print OUT "param : se : $settings\n";
		print OUT "# delete only entries generated by pindel\n";
		print OUT "param : c : pindel\n";
		print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n"
		  ; #start deleting when indel annotation is finished --> afterwards start inserting
		print OUT "dependson : filterVCFforRegion.pl\n";
		&printSlot("insertPindel");

		if ( $runs{insertPindel} == 1 ) {
			print OUT "run\n\n";
		}
		else {
			print OUT "#run\n\n";
		}


		if ($libtype eq "genomic" && 
			$vcf == 2 && 
			$target eq "" && 
			(($organism eq "human" || $organism eq "mouse")) ) {				#execute snvdbExomeImport_vcf.pl instead of snvdbExomeInsert_vcf.pl only if sample is a genome and is inserted into genomegatk db
			print OUT "pgr   : snvdbExomeImport_vcf.pl\n";
			print OUT "# input file\n";
			print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/filtered.genome.all.pindel.dbSNP.plus.vcf\n";
				
			print OUT "# caller\n";
			print OUT "param : c : pindel\n";
			print OUT "# class of inserted variants\n";
			print OUT "param : l : deletion\n";
			print OUT "# settings\n";
			print OUT "param : se : $settings\n";
			print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
			print OUT "dependson : snvdbExomeDelete.pl\n";
			&printSlot("insertPindel");
			if ( $runs{insertPindel} == 1 ) {
				print OUT "run\n\n";
			}
			else {
				print OUT "#run\n\n";
			}
		} else {
			print OUT "\npgr   : snvdbExomeInsert_vcf.pl\n";
			if($libtype eq "genomic" || $libtype eq "mtDNA"){
				print OUT "#infile\n";
				print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/filtered.genome.all.pindel.dbSNP.plus.vcf\n";
				
				
			}else{
				print OUT "#infile\n";
				print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/filtered.exome.all.pindel.dbSNP.plus.vcf\n";	
			}
			print OUT "# settings\n";
			print OUT "param : se : $settings\n";
			print OUT "# class of inserted variants\n";
			print OUT "param : l : deletion\n";
			print OUT "# delete only entries generated by pindel\n";
			print OUT "param : c : pindel\n";
			print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
			print OUT "dependson : snvdbExomeDelete.pl\n";
			&printSlot("insertPindel");
	
			if ( $runs{insertPindel} == 1 ) {
				print OUT "run\n\n";
			}
			else {
				print OUT "#run\n\n";
			}
		}
	}
}


sub runBreakdancer {
	if ( $vcf > 0 && $libtype eq "genomic") {
		
		print OUT "
##############################################################
# run breakdancer
##############################################################
";
		make_path("$outdir/$project/$sample/$folder/$subfolder/Breakdancer", {				#create directory for single files
				mode => 0775
			}) if ( $runs{runBreakdancer} == 1  ); #{ && $vcf == 1);
		print OUT "pgr   : runBreakdancer.pl\n";
		print OUT "# bamfile\n";
		print OUT "param : b : $outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam\n";
		print OUT "# settings\n";
		print OUT "param : se : $settings\n";
		print OUT "# start this job as an array job, i.e. one job will be started by chromosome\n";
		print OUT "param : aj\n";
		print OUT "# output directory\n";
		print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/Breakdancer\n";
		print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
		print OUT "dependson : bwa_merge_gatk.pl\n";
		print OUT "#signal for parallelpipeline.pl to start this job as an array job on Grid Engine\n";
		print OUT "slots : perChrArray\n";
		
		if ( $runs{runBreakdancer} == 1  ) { #&& $vcf == 1) {
			print OUT "run\n";
		}
		else {
			print OUT "#run\n";
		}
		
print OUT "
##############################################################
# concat single chromosome breakdancer VCFs
##############################################################
";
		print OUT "pgr   : concatVCF.pl\n";
		print OUT "#indir\n";
		print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/Breakdancer\n";
		print OUT "#outfile\n";
		print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/all.breakdancer.vcf\n";
		print OUT "#files to look for\n";
		print OUT "param : e : vcf\n";		
		print OUT "# settings from config file\n";
		print OUT "param : se : $settings\n";
		print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
		print OUT "dependson : runBreakdancer.pl\n";
		&printSlot("runBreakdancer");
		
		if ( $runs{runBreakdancer} == 1 ) { #&& $vcf == 1) {
			print OUT "run\n";
		}
		else {
			print OUT "#run\n";
		}
	}
}


sub runLumpy {
	if ( $vcf > 0 && $libtype eq "genomic") {
		
		print OUT "
##############################################################
# run lumpy-sv
##############################################################
";
	
		print OUT "pgr   : runLumpy.pl\n";
		print OUT "# bamfile\n";
		print OUT "param : b : $outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam\n";
		print OUT "# settings\n";
		print OUT "param : se : $settings\n";
		print OUT "# output directory\n";
		print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/\n";
		print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
		print OUT "dependson : bwa_merge_gatk.pl\n";
		&printSlot("runLumpy");
		
		if ( $runs{runLumpy} == 1  ) { #&& $vcf == 1) {
			print OUT "run\n";
		}
		else {
			print OUT "#run\n";
		}

	}
}

sub runManta {
	if ( $vcf > 0 && $libtype eq "genomic") {
		
		print OUT "
##############################################################
# run manta
##############################################################
";
		make_path("$outdir/$project/$sample/$folder/$subfolder/manta", {				#create directory for single files
				mode => 0775
			}) if ( $runs{runManta} == 1  );
	
		print OUT "pgr   : runManta.pl\n";
		print OUT "# bamfile\n";
		print OUT "param : b : $outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam\n";
		print OUT "# settings\n";
		print OUT "param : se : $settings\n";
		print OUT "# output directory\n";
		print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/manta/\n";
		print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
		print OUT "dependson : bwa_merge_gatk.pl\n";
		&printSlot("runManta");
		
		if ( $runs{runManta} == 1  ) { #&& $vcf == 1) {
			print OUT "run\n";
		}
		else {
			print OUT "#run\n";
		}

	}
}

sub runWhamg {
        if ( $vcf > 0 && $libtype eq "genomic") {

                print OUT "
##############################################################
# run whamg
##############################################################
";
                print OUT "pgr   : runWhamg.pl\n";
                print OUT "# bamfile\n";
                print OUT "param : b : $outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam\n";
                print OUT "# settings\n";
                print OUT "param : se : $settings\n";
                print OUT "# output directory\n";
                print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/\n";
                print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
                print OUT "dependson : bwa_merge_gatk.pl\n";
                &printSlot("runWhamg");

                if ( $runs{runWhamg} == 1  ) { #&& $vcf == 1) {
                        print OUT "run\n";
                }
                else {
                        print OUT "#run\n";
                }

        }
}


sub insertSV {
	if ( $vcf > 0 && $libtype eq "genomic") {
		
		print OUT "
##############################################################
# merge & insert structural variants
##############################################################
";
		
		print OUT "pgr   : mergeSV.pl\n";
		print OUT "# out file\n";
		print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/merged.sv.vcf\n";
		print OUT "# cnvnator file\n";
		print OUT "param : c : $outdir/$project/$sample/$folder/$subfolder/merged.rmdup.cnvnator.out\n";
		print OUT "# lumpy-sv file\n";
		print OUT "param : l : $outdir/$project/$sample/$folder/$subfolder/lumpy.vcf\n";
		print OUT "# pindel file\n";
		print OUT "param : p : $outdir/$project/$sample/$folder/$subfolder/all.pindel.vcf\n";
		print OUT "# breakdancer file\n";
		print OUT "param : b : $outdir/$project/$sample/$folder/$subfolder/all.breakdancer.vcf\n";
		print OUT "# manta file\n";
		print OUT "param : m : $outdir/$project/$sample/$folder/$subfolder/manta/results/variants/diploidSV.vcf.gz\n";
		print OUT "# whamg file\n";
		print OUT "param : w : $outdir/$project/$sample/$folder/$subfolder/whamg.filtered.annotated.vcf\n";
		print OUT "# settings\n";
		print OUT "param : se : $settings\n";

		print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
		print OUT "dependson : runLumpy.pl,concatVCF.pl,runPindel.pl,runWhamg.pl,runCNVnator.pl,runBreakdancer.pl,runManta.pl\n";
		&printSlot("insertSV");
		
		if ( $runs{mergeSV} == 1  ) {
			print OUT "run\n";
		}
		else {
			print OUT "#run\n";
		}
		
		
		print OUT "pgr   : insertSV.pl\n";
		print OUT "# vcf file from mergeSV.pl file\n";
		print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/merged.sv.vcf\n";
		print OUT "# sample name\n";
		print OUT "param : s : $sample\n";
		print OUT "# settings\n";
		print OUT "param : se : $settings\n";
		print OUT "# delete old entries\n";
		print OUT "param : d\n";
		print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
		print OUT "dependson : mergeSV.pl\n";
		&printSlot("insertSV");
		
		if ( $runs{insertSV} == 1  ) { 
			print OUT "run\n";
		}
		else {
			print OUT "#run\n";
		}
	}
}



sub bellerophon {
	print OUT "
##############################################################
# run Bellerophon
##############################################################
";

	print OUT "pgr   : runBellerophon.pl\n";
	print OUT "# settings\n";
	print OUT "param : se : $settings\n";
	print OUT "# bam file\n";
	print OUT "param : b : $outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam\n";
	print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
	print OUT "dependson :  bwa_merge_gatk.pl\n";
	&printSlot("runBellerophon");
	
	if ( $runs{runBellerophon} == 1  ) { 
		print OUT "run\n";
	}
	else {
		print OUT "#run\n";
	}
	
	
	print OUT "
##############################################################
# insert Bellerophon into database
##############################################################
";

	print OUT "pgr   : insertBellerophon.pl\n";
	print OUT "# settings\n";
	print OUT "param : se : $settings\n";
	print OUT "# bam file as outfile prefix\n";
	print OUT "param : p : $outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam\n";
	print OUT "# sample name\n";
	print OUT "param : s : $sample\n";
	print OUT "# delete entries before inserting\n";
	print OUT "param : d\n";
	print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
	print OUT "dependson :  runBellerophon.pl\n";
	&printSlot("insertBellerophon");
	
	if ( $runs{runBellerophon} == 1  ) { 
		print OUT "run\n";
	}
	else {
		print OUT "#run\n";
	}
		
}


sub IGVtools {
	print OUT "
##############################################################
# run IGVtools
##############################################################
";

	print OUT "pgr   : runIGVtools.pl\n";
	print OUT "# settings\n";
	print OUT "param : se : $settings\n";
	print OUT "# bam file\n";
	print OUT "param : b : $outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam\n";
	print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
	print OUT "dependson :  bwa_merge_gatk.pl\n";
	&printSlot("runIGVtools");
	
	if ( $runs{runIGVtools} == 1  ) { 
		print OUT "run\n";
	}
	else {
		print OUT "#run\n";
	}
	

		
}




sub stats {
	print OUT "
##############################################################
# calculate statistics of covered Bases and % reads on target
##############################################################
";
	print OUT "pgr   : calcBaseCov_api.pl\n";
	if ( $target ne "" ) {
		print OUT "# target file\n";
		print OUT "param : t : $target\n";
		print OUT "#Output a file that contains 0 coverage regions\n";
		print OUT "param : z\n"
	}

	print OUT "# outfile\n";
	print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/calcCov.out\n";
	print OUT "# settings\n";
	print OUT "param : se : $settings\n";
	print OUT "# bam file\n";
	print OUT "param : b : $outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam\n";
	print OUT "# coverage threshold\n";
	print OUT "param : c : 20\n";
	print OUT "#Sample name\n";
	print OUT "param : s : $sample\n";
	print OUT "#generate a file that calculates a distribution of % coverage values up to the given number\n";
	print OUT "param : m : 200\n";
	print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
	print OUT "dependson : bwa_merge_gatk.pl,gem_aln.pl,gem_merge.pl,gem_single_sort.pl\n";
	&printSlot("stats");

	if ( $runs{stats} == 1 && $target ne "" )
	{ # && (!(-e "$outdir/$project/$sample/$folder/$subfolder/calcCov.out") || ($assay ne "SureSelect38Mb" && !(-e "$outdir/$project/$sample/$folder/$subfolder/coverage_per_target.txt")) )){
		print OUT "run\n";
	}
	else {
		print OUT "#run\n";
	}
	
	# Parallel metrics calculation on whole genomes:
	if ( $libtype eq "genomic" || $libtype eq "mtDNA" ){
		
		# Parallel steps:	
		#			-a	calculate alignment metrics
		#			-bd calculate base distribution by cycle
		#			-g	calculate gc metrics
		#			-i  calculate insert size distribution
		#			-qcy calculate mean quality by cycle
		#			-mm	multiple metrics ( don't use it ! ) 
		#			-qd  calculate quality score distribution
		#			-w  calculate wgs metrics
		#			-e	estimate library complexity
		
		# Alignment Metrics
			print OUT "pgr   : collectMetrics.pl\n";
			print OUT "# settings\n";
			print OUT "param : se : $settings\n";
			print OUT "# bam file\n";
			print OUT "param : b : $outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam\n";
			print OUT "#calculate alignment metrics\n";
			print OUT "param : a\n";
			print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
			print OUT "dependson : bwa_merge_gatk.pl,gem_aln.pl,gem_merge.pl,gem_single_sort.pl\n";
			&printSlot("stats");
			if ( $runs{stats} == 1 )
			{ 
				print OUT "run\n";
			}
			else {
				print OUT "#run\n";
			}
			
		# Base distribution by cycle
			print OUT "pgr   : collectMetrics.pl\n";
			print OUT "# settings\n";
			print OUT "param : se : $settings\n";
			print OUT "# bam file\n";
			print OUT "param : b : $outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam\n";
			print OUT "#calculate bd base distribution by cycle\n";
			print OUT "param : bd\n";
			print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
			print OUT "dependson : bwa_merge_gatk.pl,gem_aln.pl,gem_merge.pl,gem_single_sort.pl\n";
			&printSlot("stats");
			if ( $runs{stats} == 1 )
			{ 
				print OUT "run\n";
			}
			else {
				print OUT "#run\n";
			}

		# GC Metrics			
			print OUT "pgr   : collectMetrics.pl\n";
			print OUT "# settings\n";
			print OUT "param : se : $settings\n";
			print OUT "# bam file\n";
			print OUT "param : b : $outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam\n";
			print OUT "#calculate gc metrics\n";
			print OUT "param : g\n";
			print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
			print OUT "dependson : bwa_merge_gatk.pl,gem_aln.pl,gem_merge.pl,gem_single_sort.pl\n";
			&printSlot("stats");
			if ( $runs{stats} == 1 )
			{ 
				print OUT "run\n";
			}
			else {
				print OUT "#run\n";
			}
			
		# Insert Size Distribution			
			print OUT "pgr   : collectMetrics.pl\n";
			print OUT "# settings\n";
			print OUT "param : se : $settings\n";
			print OUT "# bam file\n";
			print OUT "param : b : $outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam\n";
			print OUT "#calculate insert size distribution\n";
			print OUT "param : i\n";
			print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
			print OUT "dependson : bwa_merge_gatk.pl,gem_aln.pl,gem_merge.pl,gem_single_sort.pl\n";
			&printSlot("stats");
			if ( $runs{stats} == 1 )
			{ 
				print OUT "run\n";
			}
			else {
				print OUT "#run\n";
			}
			
		# Mean Quality By Cycle		
			print OUT "pgr   : collectMetrics.pl\n";
			print OUT "# settings\n";
			print OUT "param : se : $settings\n";
			print OUT "# bam file\n";
			print OUT "param : b : $outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam\n";
			print OUT "#calculate mean quality by cycle\n";
			print OUT "param : qcy\n";
			print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
			print OUT "dependson : bwa_merge_gatk.pl,gem_aln.pl,gem_merge.pl,gem_single_sort.pl\n";
			&printSlot("stats");
			if ( $runs{stats} == 1 )
			{ 
				print OUT "run\n";
			}
			else {
				print OUT "#run\n";
			}		
			
		# Quality Score Distribution		
			print OUT "pgr   : collectMetrics.pl\n";
			print OUT "# settings\n";
			print OUT "param : se : $settings\n";
			print OUT "# bam file\n";
			print OUT "param : b : $outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam\n";
			print OUT "#calculate quality score distribution\n";
			print OUT "param : qd\n";
			print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
			print OUT "dependson : bwa_merge_gatk.pl,gem_aln.pl,gem_merge.pl,gem_single_sort.pl\n";
			&printSlot("stats");
			if ( $runs{stats} == 1 )
			{ 
				print OUT "run\n";
			}
			else {
				print OUT "#run\n";
			}
			
		# WGS Metrics			
			print OUT "pgr   : collectMetrics.pl\n";
			print OUT "# settings\n";
			print OUT "param : se : $settings\n";
			print OUT "# bam file\n";
			print OUT "param : b : $outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam\n";
			print OUT "#calculate wgs metrics\n";
			print OUT "param : w\n";
			print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
			print OUT "dependson : bwa_merge_gatk.pl,gem_aln.pl,gem_merge.pl,gem_single_sort.pl\n";
			&printSlot("stats");
			if ( $runs{stats} == 1 )
			{ 
				print OUT "run\n";
			}
			else {
				print OUT "#run\n";
			}
			
		# Estimate library complexity			
			print OUT "pgr   : collectMetrics.pl\n";
			print OUT "# settings\n";
			print OUT "param : se : $settings\n";
			print OUT "# bam file\n";
			print OUT "param : b : $outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam\n";
			print OUT "#estimate library complexity\n";
			print OUT "param : e\n";
			print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
			print OUT "dependson : bwa_merge_gatk.pl,gem_aln.pl,gem_merge.pl,gem_single_sort.pl\n";
			#&printSlot("stats");
			print OUT "slots : 5\n";
			if ( $runs{stats} == 1 )
			{ 
				print OUT "run\n";
			}
			else {
				print OUT "#run\n";
			}
			
	}
	else
	{
		# Calculate in one step
		print OUT "pgr   : collectMetrics.pl\n";
		print OUT "# settings\n";
		print OUT "param : se : $settings\n";
		print OUT "# bam file\n";
		print OUT "param : b : $outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam\n";
		print OUT "#calculate alignment metrics\n";
		print OUT "param : a\n";
		print OUT "#calculate base distribution by cycle\n";
		print OUT "param : bd\n";
		print OUT "#calculate gc metrics\n";
		print OUT "param : g\n";
		print OUT "#calculate insert size distribution\n";
		print OUT "param : i\n";
		print OUT "#calculate mean quality by bycle\n";
		print OUT "param : qcy\n";
		print OUT "#calculate quality distribution\n";
		print OUT "param : qd\n";
		print OUT "#estimate library complexity\n";
		print OUT "param : e\n";
		print OUT "#calculate wgs metrics\n";
		print OUT "param : w\n";		
		print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
		print OUT "dependson : bwa_merge_gatk.pl,gem_aln.pl,gem_merge.pl,gem_single_sort.pl\n";
		&printSlot("stats");
		if ( $runs{stats} == 1 && $target ne "" )
		{ 
			print OUT "run\n";
		}
		else {
			print OUT "#run\n";
		}
	}


	# OnOffTargetProfile - Stats for Whole Genome
	print OUT "\npgr   : calcOnOffTargetCoverageProfile.pl\n";
	print OUT "# settings\n";
	print OUT "param : se : $settings\n";
	print OUT "# bam file\n";
	print OUT "param : b : $outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam\n";
	print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
	print OUT "dependson : bwa_merge_gatk.pl,gem_aln.pl,gem_merge.pl,gem_single_sort.pl\n";
	&printSlot("stats");
	if ( $runs{stats} == 1 && $libtype eq "genomic" && ( $organism eq "human" || $organism eq "mouse"))
	{ 
		print OUT "run\n";
	} 
	else 
	{
		print OUT "#run\n";
	}


	print OUT "\npgr   : calcOnTarget.pl\n";
	print OUT "# target file\n";
	print OUT "param : t : $targetBED\n";
	print OUT "# infile\n";
	print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam\n";
	print OUT "# outfile\n";
	print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/locateReads_capture.out\n";
	print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
	print OUT "dependson : bwa_merge_gatk.pl,gem_aln.pl,gem_merge.pl,gem_single_sort.pl\n";
	&printSlot("stats");
	if ( $runs{stats} == 1 && $target ne "" )
	{
		print OUT "run\n";
	}
	else {
		print OUT "#run\n";
	}

	print OUT "\npgr   : parseStats.pl\n";
	print OUT "# bamfile\n";
	print OUT "param : b : $outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam\n";
	print OUT "# number of fields that should be taken from bam file path\n";
	print OUT "param : n : 5\n";
	print OUT "# outfile\n";
	print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/summary.stats.tsv\n";
	print OUT "# flagstat file\n";
	print OUT "param : f : $outdir/$project/$sample/$folder/$subfolder/merged.rmdup.flagstat.out\n";
	print OUT "# settings\n";
	print OUT "param : se : $settings\n";
	print OUT "# picard metrics\n";
	print OUT "param : p\n";
	unless($noRmdup || $libtype eq "genomic" || $libtype eq "MIP"){
		print OUT "# rmdup file\n";
		print OUT "param : r : $outdir/$project/$sample/$folder/$subfolder/merged.rmdup.log\n";
	}
	
	if ( $libtype eq "genomic" && ( $organism eq "human" || $organism eq "mouse")){
		print OUT "#OnOffTargetCoverageProfile\n";
		print OUT "param : ct : $outdir/$project/$sample/$folder/$subfolder/coverageprofile\n";
	}
	
	if ( $organism eq "human" ) {
		print OUT "# check gender\n";
		print OUT "param : g\n";
	}

	if ( $target ne "" ) {
		print OUT "# calcBaseCov file\n";
		print OUT "param : c : $outdir/$project/$sample/$folder/$subfolder/calcCov.out\n";
		print OUT "# locateReads_capture file\n";
		print OUT "param : l : $outdir/$project/$sample/$folder/$subfolder/locateReads_capture.out\n";
		print OUT "# exome depth Rsd\n";
		print OUT "param : e : $outdir/$project/$sample/$folder/$subfolder/exomeDepth.stats\n";
	}
	elsif ( $organism eq "human" ) {
		print OUT "# calcBaseCov file\n";
		print OUT "param : c : $outdir/$project/$sample/$folder/$subfolder/refSeqCoverage.txt\n"
		  ; #for  RNA-seq experiments --> use coverage of refSeq genes instead of target region
	}
	if (   ( -e "$outdir/$project/$sample/$folder/$subfolder/verifyBAM.selfSM" )
		|| ( ( $runs{checkContamination} == 1 ) && $organism eq "human" ) )
	{  #if verifyBamID has been run (or is running) insert contamination into DB
		print OUT "#contamination file\n";
		print OUT "param : v : $outdir/$project/$sample/$folder/$subfolder/verifyBAM.selfSM\n";
	}
	print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
	print OUT "dependson : exomeDepth.pl,checkContamination.pl,calcOnTarget.pl,collectMetrics.pl,calcOnOffTargetCoverageProfile.pl,calcBaseCov_api.pl,bwa_merge_gatk.pl,gem_aln.pl,gem_merge.pl,gem_single_sort.pl\n";
	&printSlot("parsestats");
	
	if ( $runs{parsestats} == 1 || $runs{stats} == 1) {	#OR syntax to preserve backwards compatibility 
		print OUT "run\n";
	}
	else {
		print OUT "#run\n";
	}

	# Library Stats duplicates and optical duplicates per lane 
	print OUT "pgr   : opticalDuplicates.pl\n";
	print OUT "# settings\n";
	print OUT "param : se : $settings\n";
	print OUT "# folder\n";
	print OUT "param : folder : $outdir/$project/$sample/$folder/$subfolder/\n";
	print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
	print OUT "dependson : bwa_merge_gatk.pl,gem_aln.pl,gem_merge.pl,gem_single_sort.pl\n";
	&printSlot("stats");
	if ( $runs{stats} == 1 )
	{ 
		print OUT "run\n";
	}
	else {
		print OUT "#run\n";
	}
	
	# chrM Stats
	print OUT "pgr   : chrMstats.pl\n";
	print OUT "# settings\n";
	print OUT "param : se : $settings\n";
	print OUT "# bam file\n";
	print OUT "param : b : $outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam\n";
	print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
	print OUT "dependson : bwa_merge_gatk.pl,gem_aln.pl,gem_merge.pl,gem_single_sort.pl\n";
	&printSlot("stats");
	if ( $runs{stats} == 1 && $runs{annotateChrM} )
	{ 
		print OUT "run\n";
	}
	else {
		print OUT "#run\n";
	}
	
	# StatsDB
	print OUT "pgr   : statsdb.pl\n";
	print OUT "# input file\n";
	print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/summary.stats.tsv\n";
	print OUT "# settings\n";
	print OUT "param : se : $settings\n";
	if ( $assay ne "" ) {
		print OUT "# assay\n";    #assay is empty for whole genome experiments
		print OUT "param : a : $assay\n";
	}
	print OUT "#libtype and libpair\n";
	print OUT "param : lt : $libtype\n";
	print OUT "param : lp : $libpair\n";
	print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
	print OUT "dependson : parseStats.pl,checkContamination.pl\n";
	&printSlot("statsdb");
	
	if ( $runs{statsdb} == 1 ) {
		print OUT "run\n\n";
	}
	else {
		print OUT "#run\n\n";
	}
	

}




sub transversion {
	print OUT "
#####################################################
# calculate number of transitions to transversion
#####################################################
";
	print OUT "pgr   : tstv.pl\n";
	print OUT "#sample\n";
	print OUT "param : s : $sample\n";
	print OUT "# settings\n";
	print OUT "param : se : $settings\n";
	print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
	print OUT "dependson : statsdb.pl,snvdbExomeInsert_vcf.pl,snvdbExomeImport_vcf.pl\n";
	&printSlot("transversion");
	if ( $runs{transversion} == 1 ) {
		print OUT "run\n";
	}
	else {
		print OUT "#run\n";
	}
}





sub htseqCount() {
	print OUT "
##############################################################
# calculate how many reads are located in a specific feature
# features are specified by a gtf file
##############################################################
";
	print OUT "pgr   : htseqCount.pl\n";
	print OUT "# infile\n";
	print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam\n";
	if ($isStrandedRNA) {
		print OUT "#Stranded RNA-seq data\n";
		print OUT "param : sr\n"
	}
	print OUT "# aligner\n";
	print OUT "param : a : $aligner\n";
	print OUT "# sample name\n";
	print OUT "param : s : $sample\n";
	print OUT "# settings\n";
	print OUT "param : se : $settings\n";
	print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
	print OUT "dependson : bwa_merge_gatk.pl,gem_aln.pl,gem_merge.pl,gem_single_sort.pl\n";
	&printSlot("htseqcount");
	
	if ( $runs{htseqcount} == 1)
	{
		print OUT "run\n";
	}
	else {
		print OUT "#run\n";
	}
	
	
	print OUT "
##############################################################
# calculate how many reads are located in a specific feature
# features are specified by a gtf file
##############################################################
";
	print OUT "pgr   : htseqCount.pl\n";
	print OUT "# infile\n";
	print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/merged.bam\n";
	if ($isStrandedRNA) {
		print OUT "#Stranded RNA-seq data\n";
		print OUT "param : sr\n"
	}
	print OUT "# aligner\n";
	print OUT "param : a : $aligner\n";
	print OUT "# sample name\n";
	print OUT "param : s : $sample\n";
	print OUT "# settings\n";
	print OUT "param : se : $settings\n";
	print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
	print OUT "dependson : bwa_merge_gatk.pl,gem_aln.pl,gem_merge.pl,gem_single_sort.pl\n";
	&printSlot("htseqcount");
	
	if ( $runs{htseqcount} == 1)
	{
		print OUT "run\n";
	}
	else {
		print OUT "#run\n";
	}
}

sub calcFPKM() {
	print OUT "
##############################################################
# calculate the fpkm values
##############################################################
";
	print OUT "pgr   : calcFPKM.pl\n";
	print OUT "# infile\n";
	print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam.htseqcounts\n";
	print OUT "# outprefix\n";
	print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/merged.rmdup.bam\n";
	print OUT "# settings\n";
	print OUT "param : se : $settings\n";
	print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
	print OUT "dependson : htseqCount.pl\n";
	&printSlot("calcFPKM");
	
	if ( $runs{calcFPKM} == 1)
	{ #currently only run for mouse and humans --> other species don't have the lncrna tables (yet)
		print OUT "run\n";
	}
	else {
		print OUT "#run\n";
	}
	
	print OUT "
##############################################################
# calculate the fpkm values
##############################################################
";
	print OUT "pgr   : calcFPKM.pl\n";
	print OUT "# infile\n";
	print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/merged.bam.htseqcounts\n";
	print OUT "# outprefix\n";
	print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/merged.bam\n";
	print OUT "# settings\n";
	print OUT "param : se : $settings\n";
	print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
	print OUT "dependson : htseqCount.pl\n";
	&printSlot("calcFPKM");
	
	if ( $runs{calcFPKM} == 1)
	{ #currently only run for mouse and humans --> other species don't have the lncrna tables (yet)
		print OUT "run\n";
	}
	else {
		print OUT "#run\n";
	}
}

sub insertRNAresults2DB() {
	print OUT "
##############################################################
# insert count and fpkm values into DB
##############################################################
";
	print OUT "pgr   : insertDBrnaSeqGeneWise.pl\n";
	print OUT "# countfile\n";
	print OUT "param : c : $outdir/$project/$sample/$folder/$subfolder/merged.bam.htseqcounts\n";
	print OUT "# fpkm file\n";
	print OUT "param : f : $outdir/$project/$sample/$folder/$subfolder/merged.bam.fpkm\n";
	
	print OUT "# sample name\n";
	print OUT "param : s : $sample\n";
	
	print OUT "# delete entries for this sample in DB before inserting new ones\n";
	print OUT "param : d\n";
	
	print OUT "# settings\n";
	print OUT "param : se : $settings\n";
	print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
	print OUT "dependson : calcFPKM.pl\n";
	&printSlot("RNAdbInsert");
	
	if ( $runs{RNAdbInsert} == 1)
	{
		print OUT "run\n";
	}
	else {
		print OUT "#run\n";
	}
	
}

sub calcRNAstats() {
	print OUT "
##############################################################
# calculate QC metrics for the sample
##############################################################
";
	print OUT "pgr   : calcRNASeQCstats.pl\n";
	print OUT "# infile\n";
	print OUT "param : i : $outdir/$project/$sample/$folder/$subfolder/merged.bam\n";
	print OUT "# sample name\n";
	print OUT "param : s : $sample\n";
	print OUT "# outputprefix\n";
	print OUT "param : o : $outdir/$project/$sample/$folder/$subfolder/merged\n";
	
	print OUT "# settings\n";
	print OUT "param : se : $settings\n";
	
	print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
	print OUT "dependson : bwa_merge_gatk.pl,gem_aln.pl,gem_merge.pl,gem_single_sort.pl\n";
	&printSlot("calcRNAstats");
	
	if ( $runs{calcRNAstats} == 1)
	{
		print OUT "run\n";
	}
	else {
		print OUT "#run\n";
	}
	
	
	#TODO:RSeQC-stuff here
}

sub insertRNAStats2DB() {
	print OUT "
##############################################################
# insert QC metrics for the sample into dB
##############################################################
";
	print OUT "pgr   : insertDBRNASeQCstat.pl\n";
	print OUT "# folder\n";
	print OUT "param : f : $outdir/$project/$sample/$folder/$subfolder/\n";
	print OUT "# sample name\n";
	print OUT "param : s : $sample\n";
	
	print OUT "# aligner\n";
	print OUT "param : m : $aligner\n";
	print OUT "# file from which RNASeQC created results form\n";
	print OUT "param : n : $outdir/$project/$sample/$folder/$subfolder/merged.bam\n";
	
	print OUT "# settings\n";
	print OUT "param : se : $settings\n";
	print OUT "#Dependencies on other scripts --> jobs must have finished before this script can run\n";
	print OUT "dependson : calcRNASeQCstats.pl\n";
	
	&printSlot("RNAStatsDBInsert");
	if ( $runs{RNAStatsDBInsert} == 1)	{
		print OUT "run\n";
	} else {
		print OUT "#run\n";
	}
}

###############################################################################
sub printSlot {
	my $run = shift;
	print OUT "#Number of slots that should be used in Grid Engine\n";
	if($slots{$run}){
		print OUT "slots : $slots{$run}\n";
	}elsif( $run =~ /^[0-9]+$/ ){
		print OUT "slots : $run\n";
	}else{
		print OUT "slots : 1\n";
	}
	
}
