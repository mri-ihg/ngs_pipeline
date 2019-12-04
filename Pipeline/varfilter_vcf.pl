#!/usr/bin/perl

#use strict;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);


my $prog_path = dirname(abs_path($0));
require $prog_path."/Utilities.pm";

my $outdir   = "";
my $help     = 0;
my $params = Utilities::getParams();
my $logfile  = "pipeline.log";
my $loglevel = "INFO";
my $settings = "";
my $infiles  = "";

my $bamutil	= $params->{programs}->{bamutil}->{path};

my $sam		= $params->{programs}->{samtools}->{path};
my $bcftools= $params->{programs}->{samtools}->{bcftools};
my $vcfutils= $params->{programs}->{samtools}->{vcfutils};
my $vcftoolslib = $params->{programs}->{vcftools}->{lib};
my $vcftoolssort= $params->{programs}->{vcftools}->{"sort"};

my $mpC		= $params->{programs}->{samtools}->{mpileupparams}->{C};
my $mpm		= $params->{programs}->{samtools}->{mpileupparams}->{m};
my $mpF		= $params->{programs}->{samtools}->{mpileupparams}->{F};
my $mpd		= $params->{programs}->{samtools}->{mpileupparams}->{d};
my $mpq		= $params->{programs}->{samtools}->{mpileupparams}->{q};
my $mpL		= $params->{programs}->{samtools}->{mpileupparams}->{L};
my $mE      = 0;
if($params->{programs}->{samtools}->{mpileupparams}->{E}){
	$mE = 1;
}

my $vfQ		= $params->{programs}->{samtools}->{varfilterparams}->{Q};
my $vfd		= $params->{programs}->{samtools}->{varfilterparams}->{d};
my $vfD		= $params->{programs}->{samtools}->{varfilterparams}->{D};
my $vfa		= $params->{programs}->{samtools}->{varfilterparams}->{a};
my $vfw		= $params->{programs}->{samtools}->{varfilterparams}->{w};
my $vfW		= $params->{programs}->{samtools}->{varfilterparams}->{W};
my $vfp1	= $params->{programs}->{samtools}->{varfilterparams}->{p1};
my $vfp2	= $params->{programs}->{samtools}->{varfilterparams}->{p2};
my $vfp3	= $params->{programs}->{samtools}->{varfilterparams}->{p3};
my $vfp4	= $params->{programs}->{samtools}->{varfilterparams}->{p4};
my $region  = "";
my $triosamples = "";
my $outputAll   = 0;
my $optionv	    = "v";
my $ignoreRG    = 0;
my $regionFile  = "";
my $clipOverlap = 0;
my $afs			= "";
my $multiplyD   = 1;
my $filePrefix  = "";
my $resetFilter = 0;
my $noVarfilter = 0;


my $helptext      = 
"
-i	<list of infiles> or text file that includes bam files; if empty outdir/merged.rmdup.bam will be taken as infile; if a VCF file
	is given: just run varFilter on the file; if a directory is given: use all *.sort.bam files in this directory
-T	<Trio> Sample Ids for trio analysis; must be in the order \"child father mother\"; .bam file list must be in same order 
-l	<region.bed> BED file with regions for which variants should be called
-r	<region> where variants should be called; outputfile will be <region>ontarget.varfilter.vcf; default: whole genome
-c	clip overlapping reads with \"bam clipOverlap\"
-o	<outdir>
-se	<settings>
-R	ignore readgroup entries --> take filenames as samplenames in VCF file
-n	don't run varFilter
-a	output VCF variant entries for all positions in <region> even if no alternative genotype is called; useful for comparing samples
	implies -n!
-P	initial allele frequency (afs) file from previous runs
-mD	multiply the maximum allowed readdepth (varFilter -D) by this number (should be the number of samples); default: $multiplyD
-fp	outfile prefix string
-rf	if -i specifies a VCF file, reset the filter column to \".\"
-lf	log file; default: pipeline.log
-ll	log level: ERROR,INFO,DEBUG; default: INFO
-h	this help\n";


GetOptions(
"i=s" => \$infiles,
"r=s" => \$region,
"T=s" => \$triosamples,
"l=s" => \$regionFile,
"a"	  => \$outputAll,
"R"		=> \$ignoreRG,
"c"		=> \$clipOverlap,
"P=s"		=> \$afs,
"n"		=> \$noVarfilter,	
"mD=s"	=> \$multiplyD,
"fp=s"	=> \$filePrefix,
"rf"	=> \$resetFilter,
#"ap"  => \$append,
"o=s" => \$outdir, 
"se=s" => \$settings,
"lf=s" => \$logfile,
"ll=s" => \$loglevel,
"h" => \$help);

Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();

if ($help == 1) {
	print $helptext;exit(1);
}

my $BamFileListPrefix = "";
my $clipOverlapString = "";

if($infiles eq ""){
	$infiles = "$outdir/merged.rmdup.bam";
}elsif(-d $infiles){
	system("ls $infiles/*.sort.bam > $infiles/sort.bam.list"); 	#create list with all sort.bam files to call variants from
	$infiles .= "/sort.bam.list";
}




if($region ne "" && $regionFile eq ""){
	my ($chr,$rest)  = split(":",$region);
	my ($start,$end) = split("-",$rest);
	
	$regionFile = "$outdir/$region.region.bed";
	open BED, ">$regionFile" or exit $logger->error("Can't open $regionFile!");
	print BED "$chr\t$start\t$end";
	close BED;
}

if($afs ne ""){
	$afs = "-P ".$afs;
}

if($clipOverlap){
	
	#if infiles contains a file with bamfiles --> read the bam paths
	if(!($infiles =~ /\.bam$"*/)){
		open BAMS, $infiles or exit $logger->error("Can't open $infiles!");
		$infiles = "";
		while(<BAMS>){
			chomp;
			$infiles .= $_." ";
		}
		$infiles =~ s/\s$//;
	}
	
	#more than one bam file
	if($infiles =~ /\s/){
		#construct header
		my @files = split(' ',$infiles);
		system("$sam view -H $files[0] > $outdir/".$filePrefix."$region".".tmpHeader.sam");
		#print "test1: $sam view -H $files[0] > $outdir/$region".".tmpHeader.sam\n";
		for(my $i=1;$i<@files;$i++){
			#print "test2: $sam view -H $files[$i] | grep \"\@RG\" >> $outdir/$region".".tmpHeader.sam\n";
			system("$sam view -H $files[$i] | grep \"\@RG\" >> $outdir/".$filePrefix."$region".".tmpHeader.sam");
		}
		my $regionString = "";
		if($region ne ""){
			$regionString = "-R $region";
		}
		$clipOverlapString = "$sam merge $regionString -h $outdir/".$filePrefix."$region".".tmpHeader.sam -u - $infiles | $bamutil clipOverlap --in -.ubam --out -.ubam |";
		
	}else{
		if($regionFile ne ""){
			$clipOverlapString = "$sam view -h -L $regionFile -u $infiles | $bamutil clipOverlap --in -.ubam --out -.ubam |";
		}else{
			$clipOverlapString = "$bamutil clipOverlap --in $infiles --out -.ubam |";
		}
		
		
	}
	$infiles = "-";
	
}else{

	if(!($infiles =~ /\.bam"*$/)){
		$BamFileListPrefix = "-b";
	}
}


if($outputAll){
	$optionv = "";
}

my $ref = $params->{settings}->{$settings}->{reference};

$vfD *= $multiplyD;





my $filtered    = $outdir."/".$filePrefix.$region."filtered.vcf";
my $notFiltered = $outdir."/".$filePrefix.$region."not_filtered.vcf";
my $outFile     = $outdir."/".$filePrefix.$region."ontarget.varfilter.vcf";

my $command;

if(!($infiles =~ /\.vcf$/) ){
	$command="$clipOverlapString $sam \\
mpileup --ff 1024 -C$mpC -m$mpm -F$mpF -L $mpL \\\n";
	if($mE){
		$command .= "-E \\\n";
	}
	if($ignoreRG){
		$command .= "-R \\\n";
	}
	
	
	if($regionFile ne ""){
		$command .= "-l $regionFile \\\n";
	}
	$command .= "-d$mpd \\
-q$mpq \\
-DSuf $ref \\
$BamFileListPrefix $infiles \\
| $bcftools view -cg$optionv $afs \\
";
	if($regionFile ne ""){
		$command .= "-l $regionFile \\\n";
	}
	if($triosamples ne ""){	#analyse as trio
		my ($child,$father,$mother) = split(' ',$triosamples);
		open(TRIO,">$outdir/trio.txt");
		print TRIO "$child\n$father\n$mother";
		close TRIO;
		$command .= "-s $outdir/trio.txt -T trioauto \\\n";			# TODO: Change trioauto to trioxd/trioxs for ChrX regions and male/female child
	}
	$command .= " - 2> $outdir/".$filePrefix.$region."mpileup.log ";#| grep -v \"#\" >> $outdir/$region.ontarget.varfilter.vcf "; 
}else{
	if($resetFilter){
		$command = " awk '{OFS=\"\\t\";if(\$1 !~ /^#/){\$7=\".\"}; print}' $infiles";
	}else{
		$command = "cat $infiles";
	}
	
}


if($outputAll || $noVarfilter){
	$command .= " > $outFile";
	$logger->info("Running mpileup only...");
}else{
	$command.="| $vcfutils varFilter -Q$vfQ -d$vfd -D$vfD -a$vfa -w$vfw -W$vfW -1$vfp1 -2$vfp2 -3$vfp3 -4$vfp4 -p > $notFiltered 2> $filtered";
	$logger->info("Running mpileup | vcfutils varFilter...");
}

#print "Command: $command";
$logger->info("Command: $command");
system($command);


if(!($outputAll || $noVarfilter)){
	#put the two outputfiles back together, but mark one as filtered
	my $filterString = "";
	$filterString .= qq{##FILTER=<ID=VARQ,Description=\\"Filtered out by vcfutils.pl varFilter, because of to low MQ (< $vfQ)\\">\\n};
	$filterString .= qq{##FILTER=<ID=VARd,Description=\\"Filtered out by vcfutils.pl varFilter, because of to low read depth (< $vfd)\\">\\n};
	$filterString .= qq{##FILTER=<ID=VARD,Description=\\"Filtered out by vcfutils.pl varFilter, because of to high read depth (> $vfD)\\">\\n};
	$filterString .= qq{##FILTER=<ID=VARa,Description=\\"Filtered out by vcfutils.pl varFilter, because of to low read depth of variant bases (< $vfa)\\">\\n};
	$filterString .= qq{##FILTER=<ID=VARG,Description=\\"Filtered out by vcfutils.pl varFilter, because of ??? (probably something with surrounding variants)\\">\\n};
	$filterString .= qq{##FILTER=<ID=VARg,Description=\\"Filtered out by vcfutils.pl varFilter, because of ??? (probably something with surrounding variants)\\">\\n};
	$filterString .= qq{##FILTER=<ID=VARP,Description=\\"Filtered out by vcfutils.pl varFilter, because of bias probability (see PV4)\\">\\n};
	$filterString .= qq{##FILTER=<ID=VARM,Description=\\"Filtered out by vcfutils.pl varFilter, because of ???\\">\\n};
	$filterString .= qq{##FILTER=<ID=VARS,Description=\\"Filtered out by vcfutils.pl varFilter, because of Hardy-Weinberg equilibrium\\">};
	
	
	
	
	system("export PERL5LIB=".$vcftoolslib);
	$command = qq{(awk -F"\\t" '{OFS="\\t";if(\$1~/CHROM/){print "$filterString"};print}' $notFiltered; awk -F"\\t" '\$8="VAR"\$1{OFS="\\t";for(i=2; i<=NF;i++) printf("%s%s", \$i,(i==NF) ? "\\n" : OFS)}' $filtered) | $vcftoolssort > $outFile};
	$logger->info("Joining filtered/non_filtered files with awk...");
	$logger->debug("Command: $command");
	system($command);
	
	system("rm $filtered");
	system("rm $notFiltered");
}






