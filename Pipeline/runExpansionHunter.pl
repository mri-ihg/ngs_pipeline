#!/usr/bin/perl

use strict;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;
use File::Copy;
use Cwd qw(abs_path);
use Data::Dumper;
use IO::File;
use Pod::Usage;

my $prog_path = dirname(abs_path($0));
require $prog_path."/Utilities.pm";

umask(002);


# Required Options
my $bamfile = ""; 
my $outdir = "";
my $inheritoutdir = 0; #When an output directory is not specified, this option enables to take the output directory from the bamfile
my $settings = ""; #In the current.config.xml, the <settings> section for the specified analysis, e.g. "hg19_wholegenome"

# additional options:
my $sex = "female"; #will be female , the default of expansionHunter
my $threads = -1;
my $minlocuscoverage = 10; #default of expansionHunter
my $regionextensionlength = 1000; #default of expansionHunter
my $analysismode = "streaming"; #will correspond to non-default mode "streaming" of expansionHunter, but does not require an indexed bam file nor a small variant catalog (compare mode "seeking")
my $aligner = "dag-aligner"; #will correspond to expansionHunter's default aligner "dag-aligner", otherwise switch to "path-aligner"  

# Question: Some analysis settings (like aligner, analysis mode, coverage, extension length) are specified in the command line, some in the current.config.xml (e.g. references). What is the rational behind that?
my $logfile = "pipeline.log";
my $loglevel = "INFO";

# Help/POD
my $help = 0;
my $man = 0;

#update variables with options from the command line..."?=s" means a character string follows the option flag
GetOptions(
"b=s" => \$bamfile,
"o=s" => \$outdir,
"i" => \$inheritoutdir,
"se=s" => \$settings,
"sx=s" => \$sex,
"t=s" => \$threads,
"c=s" => \$minlocuscoverage,
"l=s" => \$regionextensionlength,
"sm=s" => \$analysismode,
"a=s" => \$aligner,
"lf=s" => \$logfile,
"ll=s" => \$loglevel,
"man" => \$man,
"h" => \$help
);


#Show help to the script and exit (when needed)--------------------------------------------------------------------------------------------------------
pod2usage( {-exitval => 0  ,-verbose => 1} ) if $help;
pod2usage( {-exitval => 0  ,-verbose => 2} ) if $man;

#Check whether required options were specified---------------------------------------------------------------------------------------------------------
pod2usage( {-message => "Please specify a bamfile, output directory, and a key to the settings in the current.config.xml.", -exitval => 1  ,-verbose => 1} ) if $bamfile eq "" || ! -e $bamfile || ( $outdir eq "" && not($inheritoutdir) ) || $settings eq "" ;

$bamfile = abs_path($bamfile); 

#Set output directory to the directory of the bam file. But when an output directory in explicitly defined, than that will be taken as output directory.
#Thus, the option inheritoutdir has lower level priority.
$outdir = abs_path(dirname($bamfile)) if $outdir eq "" && $inheritoutdir;

# Initialize the logger---------------------------------------------------------------------------------------------------------------------------------
Utilities::initLogger($logfile, $loglevel);
my $logger = Utilities::getLogger();

# Get settings
my $params = Utilities::getParams();
#get path to the executable expansionHunter
my $expansionHunter = $params->{programs}->{expansionHunter}->{path};
#get settings for a specific analysis with expansionHunter
my $ref = $params->{settings}->{$settings}->{reference};

# Expansion Hunter Variant Catalog
my $variantCatalog = $params->{settings}->{$settings}->{expansionHunter}->{variantCatalog};


# Get threads-------------------------------------------------------------------------------------------------------------------------------------------
if($threads == -1){
        $threads = 1;
        if($ENV{NSLOTS}){               #get number of slots given by SGE
                $threads = $ENV{NSLOTS};
        }
}



#execute the program-----------------------------------------------------------------------------------------------------------------------------------

# Defaults to femnale if sex is not either male or female.
$sex = "female" if ( $sex ne "male" );


#Note: The min-locus-coverage option ($minlocuscoverage) is currently not supported in expansionHunter
my $command = "$expansionHunter \\
--reads $bamfile \\
--reference $ref \\
--variant-catalog $variantCatalog \\
--output-prefix $outdir/expansionhunter \\
--sex $sex \\
--threads $threads \\
--region-extension-length $regionextensionlength \\
--analysis-mode $analysismode \\
--aligner $aligner \\
--log-level info \\
2> $outdir/expansionHunter.err";


#run expansionHunter
if (&Utilities::executeCommand($command,"Running expansionHunter with bam file $bamfile", $logger)) {
  #in case the execution of the command gave an exit status unequal 0, an error message has to be invoked and an exit status of 100, which tells the pipeline to stop subsequent dependant scripts 
  $logger->error("Running expansionHunter with bam file $bamfile failed");
  exit(100)
} else {
  $logger->debug("Running expansionHunter with bam file $bamfile OK");
};


#PODs--------------------------------------------------------------------------------------------------------------------------------------------------



=head1 NAME

runExpansionHunter.pl

=head1 SYNOPSIS

runExpansionHunter.pl -b path/to/merged.rmdup.bam -o expansionhunterdir -se settings [ -sx sex ]

=head1 DESCRIPTION

This script (alpha version) is a wrapper script to run ExpansionHunter.

Caution: The default analysis mode of this script is "streaming", which is not the default of expansionhunter itself. However, "streaming" does not require an indexed bam nor a short variant catalog and is therefore more convinient. But "streaming" needs more memory.

Putative bugs: In principal, expansionhunter can take URLs as an input BAM file. But, as an option of this script, the output directory can be the directory of the BAM file, which could give an error when an URL is supplied.

Note: In priciple, expansionhunter should have an option for minimal locus coverage (according to the documentation). In the current software version, this option does not exist.


=head1 OPTIONS

 -b     <BAM file> to run expansionhunter on; required
 -o     <output directory> to write the results to; required
 -i     sets output directory to the directory of the BAM file (only in effect when no output directory specified with option -o)
 -se    <expression> of the settings as pre-defined by the configuration file (e.g. current.config.xml); required
 -sx    sex: male, female; default: female (only affects analysis of sex chromosomes)
 -t     <threads>, default: take number of slots from SGE environment variable NSLOTS or 1 if no SGE job
 -c     Currently deprecated (integer specifying the minimal locus coverage; default: 10)
 -l     integer specifying the region extension length; default: 1000
 -sm    analysis mode: streaming, seeking; default: streaming (Caution: Not the default of expansionhunter.)
 -a     aligner in expansionhunter: dag-aligner, path-aligner; default: dag-aligner
 -lf    <log file>; default: pipeline.log
 -ll    log level: ERROR,INFO,DEBUG; default: INFO
 -h     print this helptext
 -man   show manual page

=head1 AUTHOR

Erik Tilch

=cut

