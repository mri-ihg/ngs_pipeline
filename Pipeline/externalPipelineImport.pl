#!/usr/bin/perl

use strict;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);


my $prog_path = dirname(abs_path($0));
require $prog_path."/Utilities.pm";

# Options:
my $sample		= "";
my $vcf 		= "";
my $exomedepth 		= "";
my $settings		= "";

my $caller		= "gatk";

my $bam			= "";
my $flagstat		= "";
my $ontarget		= "";
my $calcbasecov		= "";

my $verifybamid		= "";
my $picardmarkdup	= "";
my $metricsalignment	= "";
my $metricsbasedistrib	= "";
my $metricsGC		= "";
my $metricslibcompl	= "";
my $metricsmeanqual     = "";

my $tmp			= "/tmp";

my $logfile  		= "SCREEN";
my $loglevel		= "INFO";

my $help		= 0;


my $helptext      = 
"Data files:
 -sample		<sample name> 			required
 -vcf 			<vcf file to be imported>	required
 -exomedepth 		<exomedepth results>	      	optional

Settings
 -se			<settings>			required

QC Arguments
 -bam			<sample bam file>		required
 -flagstat		<samtools flagstat results>	optional

 -ontarget		<ontarget script results>	optional
 -calcbasecov		<calcbasecov script results>	optional

 -verifybamid		<verifybamid output>		optional
 -picardmarkdup	<mark duplicates output>	optional

 Picard/GATK Metrics output files ( all optional )
 -metricsalignment	
 -metricsbasedistrib
 -metricsGC
 -metricslibcompl
 -metricsmeanqual

Other options
 -tmp 			<where tmp files are stored> 	optional Default: /tmp
 -lf			log file; default: SCREEN
 -ll 			log level: ERROR,INFO,DEBUG; default: INFO
 -h 			this help\n";



GetOptions(
"sample=s"		=> \$sample,
"vcf=s"			=> \$vcf, 
"exomedepth=s"		=> \$exomedepth, 
"se=s"    		=> \$settings,
"bam=s"			=> \$bam,
"flagstat=s"		=> \$flagstat,
"ontarget=s"		=> \$ontarget,
"calcbasecov=s"		=> \$calcbasecov,
"verifybamid=s"		=> \$verifybamid,
"picardmarkdup=s"	=> \$picardmarkdup,
"metricsalignment=s"	=> \$metricsalignment,
"metricsbasedistrib=s"	=> \$metricsbasedistrib,
"metricsGC=s"		=> \$metricsGC,
"metricslibcompl=s"	=> \$metricslibcompl,
"metricsmeanqual=s"	=> \$metricsmeanqual,
"tmp=s"			=> \$tmp,
"lf=s"			=> \$logfile,
"ll=s"			=> \$loglevel,
"h" => \$help);

Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();

if ($help == 1 || $sample eq "" || $vcf eq "" || $settings eq "" ) {
	print $helptext;exit(1);
}

if ( ! -e $vcf ){ print "VCF file $vcf does not exist\n";exit(1);    }
if ( ! -d $tmp ){ print "TMP folder $tmp does not exist\n"; exit(1); }



# Essential steps for VCF Import:

# Cleanup TMP folder 
system ("rm $tmp/$sample.tmp.*");

$logger->info("Annotation");
system("
	perl $prog_path/annotateVCF.pl -i $vcf -o $tmp/$sample.tmp.plus.vcf -se $settings -w 20
");


$logger->info("Check: filtering SNP Quality");
if ( $bam eq "" )
{
	$logger->info("Skipping: filtering SNP Quality - Supply BAM file to activate function");
	system ("ln -s $tmp/$sample.tmp.plus.vcf $tmp/$sample.tmp.plus.checked.vcf");
}
else
{
	system("
		perl $prog_path/filterSNPqual.pl -b $bam -i $tmp/$sample.tmp.plus.vcf -o $tmp/$sample.tmp.plus.checked.vcf -se $settings
	");
}

$logger->info("Insert variants");
system("
	perl $prog_path/snvdbExomeInsert_vcf.pl -i $tmp/$sample.tmp.plus.checked.vcf -c $caller -se $settings
");

$logger->info("Run snpEff");
system("
	perl $prog_path/runsnpEff.pl -i $tmp/$sample.tmp.plus.checked.vcf -se $settings
");

$logger->info("Insert snpEff");
system("
	perl $prog_path/insertsnpEff.pl -i $tmp/$sample.tmp.plus.checked.snpEff.vcf -se $settings
");

$logger->info("Update variants with snpEff");
system("
	perl $prog_path/updateVariants_vcf.pl -i $tmp/$sample.tmp.plus.checked.snpEff.vcf -se $settings
");

$logger->info("lncRNA Annotation");
system("
	perl $prog_path/annotateVCF.pl -i $tmp/$sample.tmp.plus.checked.snpEff.vcf -o $tmp/$sample.tmp.plus.checked.snpEff.lncRNA.vcf -se $settings -l -f lincrna -w 5
");

$logger->info("Update variants with lncRNA Annotation");
system("
	perl $prog_path/updateVariants_vcf.pl -i $tmp/$sample.tmp.plus.checked.snpEff.lncRNA.vcf -c lincrna -f -se $settings
");

$logger->info("miRNA Annotation");
system("
	perl $prog_path/annotateVCF.pl -i $tmp/$sample.tmp.plus.checked.snpEff.vcf -o $tmp/$sample.tmp.plus.checked.snpEff.miRNA.vcf -se $settings -mir -f lincrna -w 0
");

$logger->info("Update variants with lncRNA Annotation");
system("
	perl $prog_path/updateVariants_vcf.pl -i $tmp/$sample.tmp.plus.checked.snpEff.lncRNA.vcf -c mirna -f -se $settings
");






##################################################################################



#The steps

# Necessary

	#annotate
	#filterSNP Qual with BAM 
	#insert
	#snpeff
	#insert snpeff
	#update snpeff
	#annotate lncrna
	#update lncrna
	#annotate mirna
	#update mirna	

# Optional: Exomedepth
	# insert exomedepth

# Optional: Quality control
	# parseStats
	# statsdb



 


