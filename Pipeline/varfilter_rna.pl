#!/usr/bin/perl

use strict;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);
use Pod::Usage;
umask(002);


my $prog_path = dirname(abs_path($0));
require $prog_path."/Utilities.pm";

my $outfile    = "";
my $help       = 0;
my $params     = Utilities::getParams();
my $logfile    = "pipeline.log";
my $loglevel   = "INFO";
my $settings   = "";
my $infiles    = "";
my $outputAll  = 0;
my $regionFile = "";
my $nct			= 1;
my $maxRam		= "4g";
my $man			= 0;

my $forceActive = 0;



GetOptions(
"i=s"  => \$infiles,
"l=s"  => \$regionFile,
"a"	   => \$outputAll,
"fa"   => \$forceActive,
"o=s"  => \$outfile, 
"se=s" => \$settings,
"m=s"  => \$maxRam,
"nct=s"=> \$nct,
"lf=s" => \$logfile,
"ll=s" => \$loglevel,
"man"  => \$man,
"h"    => \$help);


pod2usage( {-exitval => 0  ,-verbose => 1} ) if $help;
pod2usage( {-exitval => 0  ,-verbose => 2} ) if $man;
pod2usage( {-exitval => 1  ,-verbose => 1} ) if $outfile eq "" || $settings eq "";

Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();

my $ref  = $params->{settings}->{$settings}->{reference};
my $gatk = $params->{programs}->{gatk}->{path};
my $gatk4= $params->{programs}->{gatk4}->{path};
my $java = $params->{programs}->{java}->{path};
my $tmp  = $params->{programs}->{gatk}->{tmpdir};
my $sam  = $params->{programs}->{samtools}->{path};

# Set our version of Java in Path with priority so to overcome any current installation (required by GATK4)
$ENV{'PATH'}=dirname(abs_path($java)).":".$ENV{'PATH'};

# Getting infile
$infiles = dirname($outfile)."/merged.rmdup.bam";

my $javalog = $outfile . ".javalog";

#command to SplitNCigarReads:
#my $command ="
#$java -XX:ParallelGCThreads=1 -Xmx$maxRam -jar $gatk \\
#-T SplitNCigarReads \\
#-R $ref \\
#-I $infiles \\
#";
#$infiles =~ s/bam$/split.bam/;
#$command .= " -o $infiles \\
#-rf ReassignOneMappingQuality \\
#-RMQF 255 \\
#-RMQT 60 \\
#-U ALLOW_N_CIGAR_READS \\
#";

#GATK4
# why arguments removed
# -rf 255->60 max mapping quality not needed anymore
# -U ALLOW_N_CIGAR_READS -> allow all reads is applied by default 
my $command ="
$gatk4 --java-options \"-Xmx$maxRam -Xms$maxRam\" \\
SplitNCigarReads \\
-R $ref \\
-I $infiles \\
";
$infiles =~ s/bam$/split.bam/;
$command .= " -O $infiles \\
";
$command .= " >> $javalog 2>&1";
&Utilities::executeCommand($command, "SplitNCigar for variant calling on RNA-seq sample", $logger);



#run HaplotypeCaller

#$command = "
#$java -XX:ParallelGCThreads=1 -Xmx$maxRam -jar $gatk \\
#-T HaplotypeCaller \\
#-R $ref \\
#-I $infiles \\
#-dontUseSoftClippedBases \\
#-stand_call_conf 20.0 \\
#-stand_emit_conf 20.0 \\
#";
#$infiles =~ s/bam$/vcf/;
#$command .= " -o $infiles >> $javalog 2>&1";

#GATK4
#-stand_emit_conf 20.0 \\ not existing anymore
$command = "
$gatk4 --java-options \"-Xmx$maxRam -Xms$maxRam\" \\
HaplotypeCaller \\
-R $ref \\
-I $infiles \\
--dont-use-soft-clipped-bases \\
-stand-call-conf 20.0 \\
";
$infiles =~ s/bam$/vcf/;
$command .= " -O $infiles >> $javalog 2>&1";
&Utilities::executeCommand($command, "Calling variants", $logger);



#filter variants
#$command = "
#$java -XX:ParallelGCThreads=1 -Xmx$maxRam -jar $gatk \\
#-T VariantFiltration \\
#-R $ref \\
#-V $infiles \\
#-window 35 \\
#-cluster 3 \\
#-filterName FS \\
#-filter \"FS > 30.0\" \\
#-filterName QD \\
#-filter \"QD < 2.0\" \\
#";
#$infiles =~ s/vcf$/filtered.vcf/;
#$command .= " -o $infiles >> $javalog 2>&1";

#GATK4
$command = "
$gatk4 --java-options \"-Xmx$maxRam -Xms$maxRam\" \\
VariantFiltration \\
-R $ref \\
-V $infiles \\
-window 35 \\
-cluster 3 \\
--filter-name FS \\
--filter-expression \"FS > 30.0\" \\
--filter-name QD \\
--filter-expression \"QD < 2.0\" \\
";
$infiles =~ s/vcf$/filtered.vcf/;
$command .= " -O $infiles >> $javalog 2>&1";
&Utilities::executeCommand($command, "Filter variants", $logger);

#make link for downstream stuff
$command = "ln -s $infiles " . dirname($infiles) . "/ontarget.varfilter.vcf";
&Utilities::executeCommand($command, "Create link", $logger);



=head1 NAME

varfilter_rna.pl

=head1 SYNOPSIS

 varfilter_rna.pl -i merged.rmdup.bam -se hg19_test

=head1 DESCRIPTION

This is a wrapper script for GATK VariantCalling. By default GATK HaplotypeCaller is used and a outfile.gvcf file is
generated as preliminary result which is then translated into vcf.

=head1 OPTIONS

 -i	<list of infiles> or text file that includes bam files (must have .list as filetype); if empty outdir/merged.rmdup.bam will be taken as infile
 	if a directory is given: use all *.sort.bam files in this directory
 -o	<outfile>; required
 -se	<settings>; required
 -a	output VCF variant entries for all positions in <region> even if no alternative genotype is called; useful for comparing samples; only works for UnifiedGenotyper
 -nct	number of GATK CPU threads; default: 1; NOT SUPPORTED IF -aj IS CHOSEN
 -m	max. memory for each Java GATK job; default 4g
 -lf	<log file>; default: pipeline.log
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page

=head1 AUTHOR

Thomas Wieland

=cut
 