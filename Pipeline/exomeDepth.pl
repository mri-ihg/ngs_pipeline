#!/usr/bin/perl

use strict;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);
use DBI;


my $prog_path = dirname(abs_path($0));
require $prog_path."/Utilities.pm";

my $outdir     = "";
my $help       = 0;
my $params     = Utilities::getParams();
my $logfile    = "pipeline.log";
my $loglevel   = "INFO";
my $settings   = "";
my $infile     = "";
my $assay      = "";
my $sample     = "";


my $helptext      = 
"
This script runs exomeDepth to discover CNVs.

-i	<infile.bam>; default outdir/merged.rmdup.bam; if a directory is given: take all *.sort.bam files from the directory as input
-o	<outdir>; required
-se	<settings>; required
-s	<samplename>; required
-a	<assay>; required
-lf	log file; default: pipeline.log
-ll	log level: ERROR,INFO,DEBUG; default: INFO
-h	this help\n";


GetOptions(
"i=s"  => \$infile,
"o=s"  => \$outdir, 
"se=s" => \$settings,
"s=s"  => \$sample,
"a=s"  => \$assay,
"lf=s" => \$logfile,
"ll=s" => \$loglevel,
"h"    => \$help);

Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();

if ($help == 1 || $outdir eq "" || $settings eq "" || $sample eq "" || $assay eq "") {
	print $helptext;exit(1);
}


if($infile eq ""){
	$infile = $outdir."/merged.rmdup.bam";
}

my $dbh = Utilities::connectExomeDB($settings);		#normalize the values from the CSV file with the median of all exons that have at least one clearly heterozygous variant
my $snvTable       = $params->{settings}->{$settings}->{exomedb}->{snvtable};
my $snvsampleTable = $params->{settings}->{$settings}->{exomedb}->{snvsampletable};
my $coredb      = $params->{coredb}->{database};
my $sampleTable = $params->{coredb}->{sampletable};

my $gender;
#get gender from sample table
my $sql = "SELECT sex FROM $coredb.$sampleTable WHERE NAME='$sample';";
my $sth = $dbh->prepare($sql) || $logger->error("Can't prepare statement: $DBI::errstr");
$sth->execute() || $logger->error($DBI::errstr);
$gender = $sth->fetchrow_array();
$logger->debug("Sex: $gender");

$gender = "" unless $gender;
$gender = "" if $gender eq "unknown";	

#create a link to the input bam (R script requires this special name)
#my $linkname = $outdir."/Exome$sample";
#system("ln -s $infile $linkname 2> /dev/null");

if(-d $infile){
	my @infiles = glob("$infile/*.sort.bam");
	$infile = join (",",@infiles);
}


#get BAM counts for further processing
my $ref          = $params->{settings}->{$settings}->{reference};
my $getBamCounts = $params->{programs}->{exomeDepth}->{getBamCounts};
my $targetFile   = $params->{settings}->{$settings}->{targets}->{$assay}->{exomedepth}->{target};
my $bamOutfile   = "$outdir/exomeDepth.TestCount";
my $command      = "R --no-restore --no-save --args $infile $targetFile $ref $bamOutfile <$prog_path/$getBamCounts > $outdir/exomeDepth.getBamCounts.log 2>&1";
$logger->debug($command);
$logger->info("Getting counts from BAM file...");
system($command);


#calculate CNVs for autosomes
my $cnv          = $params->{programs}->{exomeDepth}->{cnv};
my $exomCounts   = $params->{settings}->{$settings}->{targets}->{$assay}->{exomedepth}->{counts};
$command		 = "R --no-restore --no-save --args $exomCounts $bamOutfile $outdir/exomeDepth. 0 <$prog_path/$cnv > $outdir/exomeDepth.cnv.log 2>&1";
$logger->debug($command);
$logger->info("Calculating CNVs for autosomes...");
system($command);

#convert CSV to BED
&R2Bed($sample,"$outdir/exomeDepth.all.csv","$outdir/exomeDepth.all.bed");

#convert CNVs to IGV format
&CNV2IGV($sample,"$outdir/exomeDepth.cnv.csv","$outdir/ExomeDepthSingle.seg");

#convert CNVs to VCF format
&BED2VCF($sample,"$outdir/exomeDepth.all.bed","$outdir/exomeDepth.all.vcf");

#convert CNVs to PDF
my $cnv2pdf      = $params->{programs}->{exomeDepth}->{cnv2pdf};
$command		 = "R --no-restore --no-save --args $sample $outdir/exomeDepth.cnv.csv $outdir/exomeDepth.cnv.pdf <$prog_path/$cnv2pdf > $outdir/exomeDepth.cnv2pdf.log 2>&1";
$logger->debug($command);
$logger->info("Converting CNVs to PDF...");
system($command);


#filter BED file
my $genomicSuperDups = $params->{programs}->{exomeDepth}->{genomicSuperDups};
my $dgv				 = $params->{programs}->{exomeDepth}->{dgv};
my $bedtools  		 = $params->{programs}->{bedtools}->{path};
system("$bedtools/intersectBed -a $outdir/exomeDepth.all.bed -b $genomicSuperDups -wa -v -f 0.50 > $outdir/exomeDepth.intersect.bed");
system("$bedtools/intersectBed -a $outdir/exomeDepth.intersect.bed -b $dgv -wa -v -f 0.50> $outdir/exomeDepth.intersect2.bed");


#convert to HTML
&BED2HTML("$outdir/exomeDepth.all.bed","$outdir/exomeDepth.all.html");
&BED2HTML("$outdir/exomeDepth.intersect2.bed","$outdir/exomeDepth.intersect2.html");





#calculate CNVs for chrX
if($gender ne "" && $params->{settings}->{$settings}->{targets}->{$assay}->{exomedepth}->{$gender}){
	$command		 = "R --no-restore --no-save --args $params->{settings}->{$settings}->{targets}->{$assay}->{exomedepth}->{$gender}  $bamOutfile $outdir/exomeDepth.chrX. 1 <$prog_path/$cnv > $outdir/exomeDepth.cnv.log 2>&1";
	$logger->debug($command);
	$logger->info("Calculating CNVs for chrX...");
	system($command);
	
	#convert CSV to BED
	&R2Bed($sample,"$outdir/exomeDepth.chrX.all.csv","$outdir/exomeDepth.chrX.all.bed");
	
	#convert CNVs to IGV format
	&CNV2IGV($sample,"$outdir/exomeDepth.chrX.cnv.csv","$outdir/ExomeDepthSingle.chrX.seg");
	
	#convert CNVs to VCF format
	&BED2VCF($sample,"$outdir/exomeDepth.chrX.all.bed","$outdir/exomeDepth.chrX.all.vcf");
	
	#convert CNVs to PDF
	$command		 = "R --no-restore --no-save --args $sample $outdir/exomeDepth.chrX.cnv.csv $outdir/exomeDepth.chrX.cnv.pdf <$prog_path/$cnv2pdf > $outdir/exomeDepth.cnv2pdf.log 2>&1";
	$logger->debug($command);
	$logger->info("Converting CNVs to PDF...");
	system($command);
	
	
	#filter BED file
	system("$bedtools/intersectBed -a $outdir/exomeDepth.chrX.all.bed -b $genomicSuperDups -wa -v -f 0.50 > $outdir/exomeDepth.chrX.intersect.bed");
	system("$bedtools/intersectBed -a $outdir/exomeDepth.chrX.intersect.bed -b $dgv -wa -v -f 0.50> $outdir/exomeDepth.chrX.intersect2.bed");
	
	
	#convert to HTML
	&BED2HTML("$outdir/exomeDepth.chrX.all.bed","$outdir/exomeDepth.chrX.all.html");
}





#####################################################################
 sub R2Bed {
	my $samplename = shift;
	my $file       = shift;
	my $outfile    = shift;
	# liest den all.csv file
	# wandelt ihn zu bed files um, um fuer intersectBed vorzubereiten
	
	
	my ($startp,$endp,$type,$nexons,$start,$end,$chr,$id,$BF,$readsexpected,$readsobserved,$readsratio);
	
	open(IN, "$file");
	open(OUT, ">$outfile");
	while (<IN>) {
		if (/chromosome/) {next;} # skip header line
		chomp;
		s/\"//g;
		($startp,$endp,$type,$nexons,$start,$end,$chr,$id,$BF,$readsexpected,$readsobserved,$readsratio)=split(/\,/);
		#print OUT "chr$chr\t$start\t$end\t$samplename\t$type\t$nexons\t$readsratio\n";
		print OUT "chr$chr\t$start\t$end\t$samplename\t$type\t".($endp-$startp+1)."\t$readsratio\n";	#TW 16.03.2015 --> nexons does not give the correct number of affected exons if more than one CNVs start at the same exon
	}
	close IN;
	close OUT;
}



#####################################################################
 sub CNV2IGV {
 	my $samplename = shift;
 	my $file       = shift;
 	my $outfile    = shift;
 	
 	
 	my $dbh = Utilities::connectExomeDB($settings);		#normalize the values from the CSV file with the median of all exons that have at least one clearly heterozygous variant
	my $snvTable       = $params->{settings}->{$settings}->{exomedb}->{snvtable};
	my $snvsampleTable = $params->{settings}->{$settings}->{exomedb}->{snvsampletable};
	my $coredb      = $params->{coredb}->{database};
	my $sampleTable = $params->{coredb}->{sampletable};
 	my $idsample;
 	
 	my $sql = "SELECT idsample FROM $coredb.$sampleTable WHERE NAME='$sample';";
    my $sth = $dbh->prepare($sql) || $logger->error("Can't prepare statement: $DBI::errstr");
	$sth->execute() || $logger->error($DBI::errstr);
	exit $logger->error("patient name $sample not found in $coredb.$sampleTable - exiting") unless $idsample = $sth->fetchrow_array();			#get idsample
 	
 	
 	#check if there are any variants for the sample in the database --> otherwise normalization won't be possible
 	$sql = "SELECT count(idsnvsample) FROM $snvsampleTable WHERE idsample=$idsample and coverage>=100 and percentvar<=55 and percentvar>=45;";
 	$sth = $dbh->prepare($sql) || $logger->error("Can't prepare statement: $DBI::errstr");
	$sth->execute() || $logger->error($DBI::errstr);
	
	my $normalize = 1;
	if($sth->fetchrow_array() > 0){
		open (IN, "$file");
		my $i = 0;
		my @normValues;
		while (<IN>) {
			if ($i == 0) {$i++; next;}
			chomp;
			s/\"//g;
			my ($chrom,$start,$end,$width,$gc,$ratio)=split(/\,/);
			if ($ratio ne 'NA') {
				$sql = "select count(ss.idsnvsample) from $snvTable v inner join $snvsampleTable ss on ss.idsnv=v.idsnv where ss.idsample=$idsample and v.chrom='chr".$chrom."' and v.start>=$start and v.start<=$end and ss.coverage>=100 and ss.percentvar<=55 and ss.percentvar>=45;";
				$sth = $dbh->prepare($sql) || $logger->error("Can't prepare statement: $DBI::errstr");
				$sth->execute() || $logger->error($DBI::errstr);
				if($sth->fetchrow_array() > 0){
					push(@normValues,$ratio);
				}
			}
		}
 		close IN;
 		$normalize = &median(\@normValues);
	}
 	
 	
 	open (IN, "$file");
	open (OUT, ">$outfile");
	
	print OUT "#type=GENE_EXPRESSION
#track  graphtype=points color=0,0,0 altColor=0,0,0 viewLimits=-1.2:5 maxHeightPixels=120:120:120
feature	chrom	start	end	value\n";
	
	my $i = 0;
	while (<IN>) {
		if ($i == 0) {$i++; next;}
		chomp;
		s/\"//g;
		my ($chrom,$start,$end,$width,$gc,$ratio)=split(/\,/);
		if ($ratio ne 'NA') {
			$ratio = $ratio - $normalize;
		}
		print OUT "$samplename	$chrom	$start	$end	$ratio\n";
	}
	
	close IN,
	close OUT;
	 	
 }
 
 
 
#####################################################################
sub BED2HTML { 
	my $infile  = shift;
	my $outfile = shift;
	my @row     = ();
	my $tmp     = "";
	my $i       = 0;
	
	open (IN, "$infile");
	open (OUT, ">$outfile");
	
	print OUT qq(
<html>
<head>
<style type="text/css">
body               { background-color:#fffff3; }
td.person          { background-color:#e88a6a; }
td.default         { background-color:#fffff3; }
td.blue            { background-color:#c9e5e3; }
td.dna             { background-color:#aedecf; }
td.dnaReport       { background-color:#d7e2be; }
td.dnaReportMethod { background-color:#c9e5e3; }
td.cytoMaterial    { background-color:#c9e5e3; }
td.cytoTest        { background-color:#dbcbd3; }
td.cytoFaerbung    { background-color:#d7e2be; }
td.dominant        { background-color:#ffeddc; }
td.cytoFish        { background-color:#e5cab5; }
td.cytoReport      { background-color:#bac7db; }
td.barcode1        { background-color:#e2ecef; }
td.barcode2        { background-color:#c2d6db; }
td.barcode3        { background-color:#a6d7e3; }
td.barcode4        { background-color:#7ebdcd; }
td.barcode5        { background-color:#55a6ba; }
td.barcode6        { background-color:#3d7f8f; }
td.barcode7        { background-color:#11667b; }
td.barcode8        { background-color:#11586a; }
input.readonly     { background-color:#CCCCCC; }
table              { border-color:#CCCCCC;border-style:solid;border-width:0px 1px 1px 0px;}
table.outer        { border-color:#fffff3;border-style:solid;border-width:0px 0px 0px 0px;}
td                 { border-color:#CCCCCC;border-style:solid;border-width:1px 0px 0px 1px;}
td.outer           { border-color:#fffff3;border-style:solid;border-width:0px 0px 0px 0px;}
td.n               { border-color:#fffff3;border-style:solid;border-width:0px 0px 0px 0px;}
th                 { border-color:#CCCCCC;border-style:solid;border-width:1px 0px 0px 1px;}
*.big              { font-size: 18px; font-family: Arial, Verdana, Helvetica, sans-serif;}
*                  { font-size: 12px; font-family: Arial, Verdana,  Helvetica, sans-serif;}
a                  { text-decoration:none;color:#2e6385;}
a:hover            { text-decoration:underline;}
a.menu             { font-size:14px;color:#FFFFFF; font-family:Sans-Serif,Arial,Helvetica,Verdana;text-decoration:none;}   
a.menu:hover       { text-decoration:none;color:#99a4b5;}
a.menuactive       { font-size:14px;color:#99a4b5; font-family:Sans-Serif,Arial,Helvetica,Verdana;text-decoration:none;}   
*.header           { border-color:#fffff3;border-style:solid;border-width:0px 0px 0px 0px;font-size:12px;background-color:#53739c; color:white;font-family:Sans-Serif,Arial,Helvetica,Verdana;text-decoration:none;}   
td.header          { border-color:#fffff3;border-style:solid;border-width:0px 2px 0px 0px;}
</style>
</head>
<body>
) ;

	print OUT "CNVs from ExomeDepth<br><br>";
	print OUT "<table border=1>\n";
	print OUT "<tr>";
	print OUT "<td>";
	print OUT "chrom";
	print OUT "</td>\n";
	print OUT "<td>";
	print OUT "ID";
	print OUT "</td>\n";
	print OUT "<td>";
	print OUT "del dup";
	print OUT "</td>\n";
	print OUT "<td>";
	print OUT "n exons";
	print OUT "</td>\n";
	print OUT "<td>";
	print OUT "ratio";
	print OUT "</td>\n";
	print OUT "</tr>\n";
	
	while (<IN>) {
		chomp;
		print OUT "<tr>";
		@row=split(/\t/);
		$i = 0;
		foreach $tmp (@row) {
			if (($i==0) or ($i==1) or ($i==2)) {
			if ($i==0) {
			print OUT "<td>";
			print OUT "<a href='http://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg19&position=$row[0]:$row[1]-$row[2]'>$row[0]:$row[1]-$row[2]</a>";
			print OUT "</td>\n";
			}
			}
			else {
				print OUT "<td>";
				print OUT "$tmp";
				print OUT "</td>\n";
			}
			$i++;
		}
		print OUT "</tr>\n";
	}
	
	print OUT "</table>\n";
	print OUT "</body>\n";
	print OUT "</html>\n";
	
	close IN;
	close OUT;

}


#####################################################################
sub BED2VCF { 
	my $samplename = shift;
	my $infile  = shift;
	my $outfile = shift;
	#my @row     = ();
	#my $tmp     = "";
	#my $i       = 0;
	
	my $dbh      = Utilities::connectVariationDB($settings);
	my $database = $params->{settings}->{$settings}->{variationdb}->{database};
	my $table    = $params->{settings}->{$settings}->{variationdb}->{genetable};
	
	my $filter = "PASS";
	if($params->{settings}->{$settings}->{targets}->{$assay}->{exomedepth}->{maxrsd} && -e $outdir."/exomeDepth.stats"){
		open EDSTATS, $outdir."/exomeDepth.stats" || exit $logger->error("Cannot open $outdir/exomeDepth.stats!");
		my $line = <EDSTATS>;
		$line    = <EDSTATS>;
		my ($tmp1,$tmp2,$exomedepthmad, $meanRsd) = split(' ',$line);
		$filter = "TooManyCNVs" if $meanRsd > $params->{settings}->{$settings}->{targets}->{$assay}->{exomedepth}->{maxrsd};
	}
	
	open (IN, "$infile");
	open (OUT, ">$outfile");

	#print header
	
	print OUT "##fileformat=VCFv4.1\n";
	print OUT "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
	print OUT "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n";
	print OUT "##INFO=<ID=SVLEN,Number=-1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">";
	print OUT "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the CNV.\">\n";
	print OUT "##INFO=<ID=CLASS,Number=1,Type=String,Description=\"Class of the CNV. duplication or deletion\">\n";
	print OUT "##INFO=<ID=EXONS,Number=1,Type=Integer,Description=\"Number of exons affected by this CNV\">\n";
	print OUT "##INFO=<ID=RATIO,Number=1,Type=Float,Description=\"Ratio\">\n";
	print OUT "##INFO=<ID=CNVTRANS,Number=.,Type=String,Description=\"Transcript(s) this CNV overlaps with\">\n";
	print OUT "##INFO=<ID=CNVGENE,Number=.,Type=String,Description=\"Gene(s) this CNV overlaps with\">\n";
	print OUT "##INFO=<ID=CNVSTRAND,Number=.,Type=String,Description=\"Strand(s) of the Transcript(s)/Gene(s) this CNV overlaps with\">\n";	
	print OUT "##INFO=<ID=FUNC,Number=.,Type=String,Description=\"Putative function of this variant, according to all isoforms.\">\n";
	print OUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$samplename\n";
	
	#write variants
	while(<IN>){
		chomp;
		my ($chr,$start,$end,$samplename,$type,$nexons,$readsratio) = split();
		
		#get information of overlapping genes from database
		my $query = "Select distinct name,strand,geneSymbol FROM $table
					 WHERE chrom='$chr' AND (($end>=txStart AND $end<=txEnd) OR (txEnd>=$start AND txEnd<=$end))";
		$logger->debug( "query: " . $query );
		my $out = $dbh->prepare($query) || exit $logger->error($DBI::errstr);
		$out->execute || exit $logger->error($DBI::errstr);
		
		my $trans   = "";
		my $gene    = "";
		my $strands = "";
		
		while (my ($name,$strand,$geneSymbol) = $out->fetchrow_array)
		{
			$trans   .= $name.",";
			$strands .= $strand.",";
			$gene    .= $geneSymbol.",";
		}
		$trans   =~ s/,$//;
		$strands =~ s/,$//;
		$gene    =~ s/,$//;
		$gene    =~ s/\s//g;
		
		my $svtype = "DEL";
		$svtype    = "DUP" if $type eq "duplication";
		my $svlen  = ($end-$start);
		$svlen    *= -1 if $svtype eq "DEL";
		
		print OUT "$chr\t$start\t.\tN\t<$svtype>\t999\t$filter\tEND=$end;SVTYPE=$svtype;SVLEN=".($end-$start).";CLASS=$type;EXONS=$nexons;RATIO=$readsratio;FUNC=unknown";
		print OUT ";CNVTRANS=$trans;CNVGENE=$gene;CNVSTRAND=$strands" if $trans ne "";
		print OUT "\tGT\t0/1\n";
		
	}
	close IN;
	close OUT;


}


##############################################################################################
sub median {
	my $ref = shift;
	my $median;
	my @array = ();
	my $n = 0;
	my $t = 0;
	
	$n = (@{$ref});
	@{$ref} = sort  {$a<=>$b}  @{$ref};
	
	if ($n % 2 != 0) {
		$median=${$ref}[int($n/2)];
	}
	else {
		$median=  (( ${$ref}[int($n/2-1)] + ${$ref}[int($n/2)] ) / 2);
	}

	return($median);
}