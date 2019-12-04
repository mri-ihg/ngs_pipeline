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

# Options:
my $gvcf="";
my $outfile    = "";
my $settings   = "";
my $maxalt     = 0;
my $logfile    = "pipeline.log";
my $loglevel   = "INFO";
my $maxRam		= "8g";
my $man			= 0;
my $help       = 0;

#my $isArrayJob	= 0;
#my $isBEDArrayJob= 0;

# 
my $params     = Utilities::getParams();


GetOptions(
"i=s"     => \$gvcf,
"gvcfs=s" => \$gvcf,
"o=s"  => \$outfile, 
"maxalt=s" => \$maxalt,
"se=s" => \$settings,
"m=s"  => \$maxRam,
"lf=s" => \$logfile,
"ll=s" => \$loglevel,
"man"  => \$man,
"h"    => \$help);


pod2usage( {-exitval => 0  ,-verbose => 1} ) if $help;
pod2usage( {-exitval => 0  ,-verbose => 2} ) if $man;
pod2usage( {-exitval => 1  ,-verbose => 1} ) if ( $outfile eq "" || $settings eq "" );

Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();


my $ref  = $params->{settings}->{$settings}->{reference};
my $java = $params->{programs}->{java}->{path};
my $gatk4 = $params->{programs}->{gatk4}->{path};

# Set our version of Java in Path with priority so to overcome any current installation (required by GATK4)
$ENV{'PATH'}=dirname(abs_path($java)).":".$ENV{'PATH'};



#my $javalog = $outfile;
#$javalog    =~ s/vcf$/log/;

my $command="";

	
my $isMultisample = 0;

$command = "
	$gatk4 --java-options  \"-Xmx$maxRam\" GenotypeGVCFs\\
	-R $ref \\
	--variant $gvcf \\
	-O $outfile";
	
	$command .=	" --max-alternate-alleles $maxalt " if $maxalt != 0; 
	
    $logger->debug($command);
	if (&Utilities::executeCommand($command, "Running GenotypeGVCF", $logger)) 
	{
   		$logger->error("Error genotyping GVCF $gvcf");
   		exit -1,
	}
	


=head1 NAME

genotypeGVCF.pl

=head1 SYNOPSIS

 genotypeGVCF.pl -i input.gvcf -o output.vcf -se hg19_wholegenome

=head1 DESCRIPTION

This is a wrapper script for GATK GenotypeGVCFs. It sets the required javapath and calls GATK4
I'm here just because it's too annoying to do it manually

=head1 OPTIONS

 -i	<input.gvcf> GVCF File. Required 
 -o	<output.vcf> Output VCF file, genotyped. Required
 -maxalt maximum number of alternate alleles, default from GATK is 6
 -se	<settings>; required
 -m	max. memory for each Java GATK job; default 8g
 -lf	<log file>; default: pipeline.log
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page

=head1 AUTHOR

Riccardo Berutti

=cut
 
