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



my $infile 		= "";
my $outfilepfx	= "";
my $regionFile 	= "";
my $settings   	= "";
my $sample      = "";
my $benchmark	= ""; # sample to benchmark to
my $logfile    	= "pipeline.log";
my $loglevel  	= "INFO";
my $threads		= 1;
my $man			= 0;
my $help       	= 0;

my $params     = Utilities::getParams();

GetOptions(
"i=s"  => \$infile,
"o=s"  => \$outfilepfx,
"l=s"  => \$regionFile,
"se=s" => \$settings,
"sample=s" => \$sample,
"benchmark=s" => \$benchmark,
"lf=s" => \$logfile,
"ll=s" => \$loglevel,
"threads=s" => \$threads,
"man"  => \$man,
"h"    => \$help);


pod2usage( {-exitval => 0  ,-verbose => 1} ) if $help;
pod2usage( {-exitval => 0  ,-verbose => 2} ) if $man;
pod2usage( {-exitval => 1  ,-verbose => 1} ) if (($infile eq "" || $outfilepfx eq "") && ( $sample eq "" )) || $settings eq "";


Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();

# Connect to DB
my $dbh = Utilities::connectCoreDB();

my $database       = $params->{coredb}->{database};
my $statTable      = $params->{coredb}->{stattable};
my $sampleTable    = $params->{coredb}->{sampletable};
my $percentcovTable= $params->{coredb}->{percentcoveragetable};
my $solexadb	   = $params->{solexadb}->{database};
my $sample2librarytable = $params->{solexadb}->{sample2librarytable};
my $librarytable   = $params->{solexadb}->{librarytable};
my $libtypetable   = $params->{solexadb}->{libtypetable};
my $libpairtable   = $params->{solexadb}->{libpairtable};
my $assaytable     = $params->{solexadb}->{assaytable};

my $defaultvcfname = "gatk.ontarget.haplotypecaller.filtered.dbSNP.plus.checked.vcf";

# Docker, Ref

	#Docker
	my $docker			= $params->{programs}->{docker}->{path};
	
	#BedTools
	my $bedtools      = $params->{programs}->{bedtools}->{path}."/bedtools";

	#Reference
	my $ref         	= $params->{settings}->{$settings}->{reference};
		$ref			= (abs_path($ref));
		my $refpath 		= dirname(abs_path($ref));


# Infile / Outdir

	my $libtype = "";


	if ( $sample ne "" )
	{
		# Fill in: $infile / $outfilepfx / regionFile
		my $idsample=0;
		$idsample = getIdsample($sample);

		if ( (! defined($idsample)) || ($idsample==0 ) )
	      	{
			print "Sample does not exist\n"; 
			exit(-1);
		}	
		
		my $basepath = getSamplePath($sample, $settings);
	        $infile = $basepath."/".$defaultvcfname if ($infile eq ""); 	
		$outfilepfx = $basepath."/happy/happy" if ($outfilepfx eq "");
		
		# Get Library Type
		$libtype = getLibType($sample);
		
		$regionFile = dirname(getAssayPath($idsample, $settings))."/targets.original.bed" if ( $libtype ne "genomic" );
		
		my $foreignid = getForeignId($sample);
		my $pedigree  = getPedigree($sample);
		
		#If not specified take foreignid as benchmark;
		$benchmark = $foreignid if $benchmark eq "";

	}

	$infile 		= (abs_path($infile));
		my $infilepath  = dirname($infile);
		
	if ( ! -d dirname($outfilepfx) ){
		system('mkdir '.dirname($outfilepfx).';');
	}
	
	$outfilepfx 		= (abs_path($outfilepfx));
		my $outfilepfxpath = dirname($outfilepfx);


	# Benchmark dataset
	
	if (! defined ( $params->{settings}->{$settings}->{benchmark}->{$benchmark} ) )
	{
		print "Specify benchmark filename";
		exit -1;
	}
	
	my $benchmarkfile	= $params->{settings}->{$settings}->{benchmark}->{$benchmark}->{vcf};
		$benchmarkfile	= (abs_path($benchmarkfile));
		my $benchmarkfilepath = dirname(abs_path($benchmarkfile));
	
	my $hiconffile		= $params->{settings}->{$settings}->{benchmark}->{$benchmark}->{hiconf};
		$hiconffile		= (abs_path($hiconffile));
		my $hiconffilepath = dirname(abs_path($hiconffile));
		


# Default region file is NORMAL Chromosomes (Genomes)
	$regionFile 		= $params->{settings}->{$settings}->{normalchromosomes}		if ( $regionFile eq "" || $libtype eq "genomic" );
	
	$regionFile = (abs_path($regionFile));
		my $regionFilePath = dirname(abs_path($regionFile));

	my $regionfilepathentry = "";
	my $regionfileentry = "";

	

	if ( ! ($regionFile eq "" ))
	{
		$regionFile = abs_path($regionFile);
		
		my $cleanRegionFile = "$outfilepfxpath/cleanRegionFile.bed";
		my $cleanRegionFileCMD = "$bedtools sort -i $regionFile | $bedtools merge -i - | $bedtools sort -i - > $cleanRegionFile";
		
		system($cleanRegionFileCMD);
		
		$regionFile=abs_path($cleanRegionFile);
		
		$regionfileentry 		= "-R $regionFile"; 
		$regionfilepathentry 	= "-v $regionFilePath:$regionFilePath";
	}
	else
	{
		$logger->error("No region file? Something wrong here");
		exit(-100);
	}


#print "$infile $outfilepfx $regionFile\n";



# Starting:
$logger->info("Benchmarking $benchmark: starting");


# Intersection:
# bedtools intersect -a HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz -b HG001_GRCh37_1_22_v4.2.1_benchmark.bed -header | awk '{if($1 !~ /^\#/){print "chr"$0}else{print $0}}'  | bedtools intersect -a - -b /data/isilon/users/datasets/target_files/hg19/Twist_Comprehensive_Exome/modified.targets_50b_margin.bed  | wc -l

my $cleanBenchmark 		= "$outfilepfxpath/benchmark.vcf";
my $cleanBenchmarkCMD 	=	"$bedtools intersect -a $benchmarkfile -b $hiconffile -header | awk '{if( \$0 !~ /^\#/){print \"chr\"\$0}else{print \$0}}' |  sed s/contig=\\<ID=/contig=\\<ID=chr/g | sed s/chrMT/chrM/g | grep -v chrGL | grep -v chrNC | grep -v chrhs37 | bedtools intersect -a - -b $regionFile -header | uniq > $cleanBenchmark ";
my $cleanInput	   		= "$outfilepfxpath/input.vcf";
my $cleanInputCMD     	= "cat $hiconffile | awk '{print \"chr\"\$0}' |  $bedtools intersect -a $infile        -b - -header | bedtools intersect -a - -b $regionFile -header | uniq > $cleanInput";


print "input:     $infile\n";
print "outpfx:    $outfilepfx\n";
print "libtype:   $libtype\n";
print "assay:     $regionFile\n";
print "refsample: $benchmark\n";

$logger->info("Filtering Benchmark file");
system($cleanBenchmarkCMD) if (! -e $cleanBenchmark );
$logger->info("Filtering Input file");
system($cleanInputCMD)     if (! -e $cleanInput     );


# Run 

	my $command="
			 docker run -it 								\\
			 -v $outfilepfxpath:$outfilepfxpath 			\\
			 -v $refpath:$refpath 							\\
			 -v $hiconffilepath:$hiconffilepath				\\
			 $regionfilepathentry							\\
			 pkrusche/hap.py 								\\
			 /opt/hap.py/bin/hap.py 						\\
			 	$cleanBenchmark	 							\\
			 	$cleanInput									\\
				-o $outfilepfx 								\\
				-r $ref 									\\
				$regionfileentry							\\
				--threads $threads							\\
				-V
	";
	
	#				--false-positives $hiconffile				\\
	
	print $command."\n\n";


		
    $logger->info("Benchmarking: started");
  	$logger->debug($command);
	system($command); 
    $logger->info("Benchmarking: completed");



#############
#subroutines#
#############

sub getIdsample {

	my $patId = shift;
	my $sql = "";
	my $sth = "";
	my $id = "";
	
	
	#check if patient id is valid
	$sql = qq{select idsample from $database.$sampleTable where name = '$patId' };
	$sth = $dbh->prepare( $sql ) || $logger->error("Can't prepare statement: $DBI::errstr");
	$sth->execute();
	($id)=$sth->fetchrow_array();
	if(!defined($id))
	{
		print "patient name $patId not found in $database.$sampleTable - exiting\n"; 
		exit(1);
	}
	
	return $id;

}

sub getForeignId {
	my $patId = shift;
	my $sql = "";
	my $sth = "";
	my $foreignid = "";
	
	
	#check if patient id is valid
	$sql = qq{select foreignid from $database.$sampleTable where name = '$patId' };
	$sth = $dbh->prepare( $sql ) || $logger->error("Can't prepare statement: $DBI::errstr");
	$sth->execute();
	($foreignid)=$sth->fetchrow_array();
	if(!defined($foreignid))
	{
		print "patient name $patId not found in $database.$sampleTable - exiting\n"; 
		exit(1);
	}
	
	return $foreignid;
}

sub getLibType {
	my $patId = shift;
	my $sql = "";
	my $sth = "";
	my $libtype = "";
	
	
	#check if patient id is valid
	$sql = qq{select LT.ltlibtype from $database.$sampleTable S inner join $solexadb.$sample2librarytable SL on SL.idsample=S.idsample inner join $solexadb.$librarytable L on L.lid=SL.lid inner join $solexadb.$libtypetable LT on LT.ltid=L.libtype  where S.name = '$patId' };
	$sth = $dbh->prepare( $sql ) || $logger->error("Can't prepare statement: $DBI::errstr");
	$sth->execute();
	($libtype)=$sth->fetchrow_array();
	if(!defined($libtype))
	{
		print "patient name $patId not found in $database.$sampleTable - exiting\n"; 
		exit(1);
	}
	
	return $libtype;
}

sub getPedigree {
	my $patId = shift;
	my $sql = "";
	my $sth = "";
	my $pedigree = "";
	
	
	#check if patient id is valid
	$sql = qq{select pedigree from $database.$sampleTable where name = '$patId' };
	$sth = $dbh->prepare( $sql ) || $logger->error("Can't prepare statement: $DBI::errstr");
	$sth->execute();
	($pedigree)=$sth->fetchrow_array();
	if(!defined($pedigree))
	{
		print "patient name $patId not found in $database.$sampleTable - exiting\n"; 
		exit(1);
	}
	
	return $pedigree;
}

sub getAssayPath {
	my $sampleId = shift;
	my $settings = shift;
	my $sql = "";
	my $sth = "";
	my $assay="";
	my $assayPath="";

	$sql="select ASSAY.name from $database.$sampleTable S 
		inner join $solexadb.sample2library S2L on S2L.idsample=S.idsample
		inner join $solexadb.library L on L.lid=S2L.lid
		inner join $solexadb.$assaytable ASSAY on ASSAY.idassay=L.idassay 
		where S.idsample = $sampleId";
	
	$sth = $dbh->prepare( $sql ) || $logger->error("Can't prepare statement: $DBI::errstr");
        $sth->execute();
        ($assay )=$sth->fetchrow_array();
	
	if ($assay ne "" )
	{
		 $assayPath=$params->{settings}->{$settings}->{targets}->{$assay}->{bed};
	}

	return $assayPath;
}

sub getSamplePath {
	my $sample=shift;
	my $settings=shift;

	my $samplePath = "";
	
	$samplePath=dirname(Utilities::getGVCFPath4($sample,$settings));

	return $samplePath;
}


=head1 NAME

happy.pl

=head1 SYNOPSIS

 happy.pl -i sample.vcf -o happybenchmark.tsv -l region.bed -se hg19_plus 

=head1 DESCRIPTION

Runs hap.py Benchmarking tool for NA12878/NA24385 benchmark samples.

=head1 OPTIONS

 -i	     <infile>      if empty outdir/merged.rmdup.bam will be taken as infile. Where outdir is the path to the outfile
 -o	     <outfile_pfx> output benchmark report prefix
 -l	     <region.bed>  BED file with regions for which variants should be called
 -se     <settings>    required
 -sample <sampleid>    optional
 -benchmark NA12878 or NA24385 optional
 -threads <1>		   threads
 -lf     <log file>    default: pipeline.log
 -ll     log level     ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page

=head1 AUTHOR

Riccardo Berutti

=cut
 
