#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;
use Getopt::Long qw(:config no_ignore_case);
use DBI;
use XML::Simple;
use Cwd qw(abs_path);
use File::Basename;
use File::Path qw(make_path);
use DateTime;
use Pod::Usage;
umask(002);

my $prog_path = dirname( abs_path($0) );

my $flowcell          = "";
my $runfolder         = "/data/runs/Runs/";
my $projectfolder     = "";
my $configfile        = $prog_path . "/conf.initAnalysis.xml";
my $skiplanes         = "";
my $alignmentOnly     = "";
my $parseConfig       = 0;
my $sgequeue          = "";
my $help              = 0;
my $removes           = 0;
my $mergeOnly         = "";
my $version           = "";
my $sleep             = 0;
my $skipSleep         = 0;
my $jobCount          = 0;
my $addedSleepTime    = 0;
my $firstTime         = "";
my $debug             = 0;
my $ignoreFailed      = 0;
my $libtype           = "";
my $noAlignment       = 0;
my $noCount			  = 0;
my $clearOut          = 0;
my $convert2BAM       = "";
my $useFailedReads    = 0;
my $parallelPipeline  = 1;
my $isFlowcell        = 0;
my $runFile           = "";
my $libpair           = "paired-end";
my $globalSettings    = "";
my $ignoreRundate     = 0;
my $maxInsertSize     = -1;
my $defaultInsertSize = 1000;
my $noMerge           = 0;
my $varOnly           = 0;
my $dontRemoveFailed  = 0;
my $logfile           = "";
my $loglevel		  = "INFO";
my $printSQL          = 0;
my $aligner           = "";
my $globAligner		  = "";
my $doTrim            = "";
my $mergeFast         = "";
my $man				  = 0;
my $dontRmdup		  = 0;
my $noVariantCalling  = 0;
my $noDB              = 0;
my $clearFCfiles      = 0;
my $parallelEnv		  = "pipeline";
my $priority		  = "";
my $callRNAVariants   = 0;
my $dependsOnSample   = "";
my $justCreateBAM     = 0;

my $tmpArgument = "";

my $removeAdapters = 0;
my $removePhix = 0;

my $isStrandedRNA = 0;


	GetOptions(
		"fs=s"  => \$flowcell,
		"f"     => \$isFlowcell,
		"r=s"   => \$runfolder,
		"p=s"   => \$projectfolder,
		"c=s"   => \$configfile,
		"s=s"   => \$skiplanes,
		"a=s"   => \$alignmentOnly,
		"m=s"   => \$mergeOnly,
		"pc"    => \$parseConfig,
		"pri=s" => \$priority,
		"sge=s" => \$sgequeue,
		"pe=s"  => \$parallelEnv,
		"h"     => \$help,
		"re"    => \$removes,
		"v=s"   => \$version,
		"sl=s"  => \$sleep,
		"n=s"   => \$skipSleep,
		"st=s"  => \$firstTime,
		"d"     => \$debug,
		"sql"   => \$printSQL,
		"i"     => \$ignoreFailed,
		"lt=s"  => \$libtype,
		"lp=s"  => \$libpair,
		"noal"  => \$noAlignment,
		"nom"   => \$noMerge,
		"nocnt" => \$noCount,
		"normdup" => \$dontRmdup,
		"novar" => \$noVariantCalling,
		"novarDB" => \$noDB,
		"var"   => \$varOnly,
		"co"    => \$clearOut,
		"cof"   => \$clearFCfiles,
		"drf"   => \$dontRemoveFailed,
		"b=s"   => \$convert2BAM,
		"u"     => \$useFailedReads,
		"rd"    => \$ignoreRundate,
		"rf=s"      => \$runFile,
		"se=s"      => \$globalSettings,
		"mi=s"      => \$maxInsertSize,
		"lf=s"      => \$logfile,
		"ll=s"      => \$loglevel,
		"al=s"      => \$globAligner,
		"T"         => \$doTrim,
		"fastmerge" => \$mergeFast,
		"crv"		=> \$callRNAVariants,
		"man"       => \$man,
		"sr"		=> \$isStrandedRNA,
		"ta=s" => \$tmpArgument,
		"removeAdapters" => \$removeAdapters,
		"removePhix" => \$removePhix
	);

pod2usage( {-exitval => 0  ,-verbose => 1} ) if $help;
pod2usage( {-exitval => 0  ,-verbose => 2} ) if $man;
pod2usage( {-exitval => 1  ,-verbose => 1} ) if $flowcell eq "" || ($version ne "vcf" && $version ne "gatk") ;
pod2usage( {-exitval => 2  ,-verbose => 1} ) if ($version eq "vcf" && $removePhix ) ;


if ($varOnly) {
	$noAlignment = 1;
	$noMerge     = 1;
	$noCount	 = 1;
}

my %lanesToSkip;

#calculate skiplanes
if ( $skiplanes ne "" ) {
	my @columns = split( ",", $skiplanes );
	foreach (@columns) {
		if ( $_ =~ m/-/ ) {
			my ( $lower, $upper ) = split( "-", $_ );
			for ( my $i = $lower ; $i <= $upper ; $i++ ) {
				$lanesToSkip{$i} = 1;
			}
		}
		else {
			$lanesToSkip{$_} = 1;
		}
	}
}
my %onlyAlign;

#calculate alignment only lanes
if ( $alignmentOnly ne "" ) {
	my @columns = split( ",", $alignmentOnly );
	foreach (@columns) {
		if ( $_ =~ m/-/ ) {
			my ( $lower, $upper ) = split( "-", $_ );
			for ( my $i = $lower ; $i <= $upper ; $i++ ) {
				$onlyAlign{$i} = 1;
			}
		}
		else {
			$onlyAlign{$_} = 1;
		}
	}
}

my %onlyMerge;

#calculate merge only lanes
if ( $mergeOnly ne "" ) {
	my @columns = split( ",", $mergeOnly );
	foreach (@columns) {
		if ( $_ =~ m/-/ ) {
			my ( $lower, $upper ) = split( "-", $_ );
			for ( my $i = $lower ; $i <= $upper ; $i++ ) {
				$onlyMerge{$i} = 1;
			}
		}
		else {
			$onlyMerge{$_} = 1;
		}
	}
}

#read config file
my $xml    = new XML::Simple;
my $params = $xml->XMLin($configfile);

#connect to database
my $db     = $params->{solexadb}->{database};
my $host   = $params->{solexadb}->{host};
my $port   = $params->{solexadb}->{port};
my $user   = $params->{solexadb}->{user};
my $passwd = $params->{solexadb}->{password};

my $dbh = "";
unless ( $dbh = DBI->connect( "DBI:mysql:database=$db;host=$host;port=$port", $user, $passwd ) ) {
	$dbh = DBI->connect( "DBI:mysql:database=$db;host=$host;port=$port",$user, $passwd )  || die print "$DBI::errstr\n";
}

my $align = 0;
my $merge = 0;

if ( $flowcell eq "STDIN" ) {
	while (<STDIN>) {
		chomp;
		unless ($_ =~ m/^#/) {
			&init($_);
		}
	}
}
elsif ( -e $flowcell && -f $flowcell ) {

	open FCFILE, $flowcell or die "Can't open $flowcell!\n";
	while (<FCFILE>) {
		chomp;
		unless ($_ =~ m/^#/) {
			&init($_);
		}
	}
}
else {
	&init($flowcell);
}

#####################################################################################################################
sub init {
	my $flowcell = shift;

	#if ( $flowcell =~ m/^\d+$/ || $flowcell =~ m/^.\d+$/ )
	unless ( ($flowcell =~ m/XX$/ && ( length $flowcell ) == 9) || $isFlowcell || ($flowcell =~ m/XY$/ && ( length $flowcell ) == 9) || ( $flowcell =~ m/000000000\-/ && ( length $flowcell ) == 15)  )
	{   
		
		my $query = "select count(ltid) from libtype where ltlibtype='$libtype';";		#check if given libtype is OK
		
		print "query: $query" if $printSQL;
		my $out = $dbh->prepare($query) || die print "$DBI::errstr\n";
		$out->execute || die print "$DBI::errstr\n";
		my $next = 0;
		unless ( (my $tmp = $out->fetchrow_array) == 1 ){
			print STDERR "\nERROR - Can't find given libtype \"$libtype\" in solexa.libtype!\n ";
			exit(1);
		}
		
		
		
		 #input value is a sample name
		if ( $onlyAlign{1} ) {
			$align = 1;
		}
		if ( $onlyMerge{1} ) {
			$merge = 0;
		}
		&initSample( $flowcell, $align, $merge, $libtype, $libpair );
	}
	else {    #input value is a flowcell name

		$isFlowcell = 1;

		my $query = "
	select exomehg19.sample.name, group_concat(lane.alane),ltlibtype, libpair.lplibpair
	from run 
	left join lane on lane.rid=run.rid 
	left join pool on lane.idpool=pool.idpool 
	left join library2pool on library2pool.idpool=pool.idpool
	left join library on library2pool.lid=library.lid
	left join sample2library on sample2library.lid=library.lid
	left join exomehg19.sample on exomehg19.sample.idsample=sample2library.idsample
	left join libtype on libtype.ltid=library.libtype
	left join libpair on libpair.lpid=library.libpair
	where run.rname='$flowcell' AND exomehg19.sample.name IS NOT NULL
	group by exomehg19.sample.name,ltlibtype, libpair.lplibpair;
	";
		print "query: $query" if $printSQL;
		my $out = $dbh->prepare($query) || die print "$DBI::errstr\n";
		$out->execute || die print "$DBI::errstr\n";
		my $next = 0;
		while ( my ( $sample, $lane, $libtype, $libpair ) =
			$out->fetchrow_array )
		{
			$align = 0;
			$merge = 0;
			$next  = 0;
			my @lanes = split( ",", $lane );
			foreach (@lanes) {
				if ( $onlyAlign{$_} ) {
					$align = 1;
				}
				if ( $onlyMerge{$_} ) {
					$merge = 1;
				}
				if ( $lanesToSkip{$_} ) {
					$next = 1;
				}
			}

			next if $next;
			&initSample( $sample, $align, $merge, $libtype, $libpair );
		}
		if ( $sgequeue ne "" && $convert2BAM ne "" && !$debug )
		{    #convert sequence files into BAM files
			if ( $jobCount >= $skipSleep )
			{    #start the first "skipSleep" jobs right now
				$addedSleepTime += $sleep
				  ;    #add "sleep" minutes to the starttime of all other jobs
			}

			my $starttime    = &getStartTime($addedSleepTime);
			my $formatedTime =
			    substr( $starttime, 6, 2 ) . "."
			  . substr( $starttime, 4,  2 ) . "."
			  . substr( $starttime, 0,  4 ) . ","
			  . substr( $starttime, 8,  2 ) . ":"
			  . substr( $starttime, 10, 2 );

			my $command =
"perl $prog_path/$params->{scripts}->{sequence2BAMpl} -b -f $flowcell -o $convert2BAM -r $runfolder -sge $sgequeue -st $formatedTime -sl $sleep ";
			if ( $skipSleep != 0 ) {
				$command .= "-n $skipSleep ";
			}
			system($command);



		}
	}
}

#############################################################################################
sub initSample {
	my $sample  = shift;
	my $align   = shift;
	my $merge   = shift;
	my $libtype = shift;
	my $libpair = shift;

	my @dirs;
	my @prefixes;
	my @lanes;
	my @flowcells;
	my @libraries;
	my @filetypes;
	my $project      = "";
	my $pdescr		 = "";
	my $exome        = 0;
	my $usedkit      = "";
	my $settings     = "";
	my $count        = 0;
	my $globOrganism = "";
	
	print STDERR "Sample: $sample\n";


	# EXTERNAL LIBRARIES
	#TW 04.09.2013: first, look for external libraries
	my $query =
		"select library.lname, exomehg19.organism.orname, exomehg19.project.pname, exomehg19.project.pdescription, exomehg19.sample.analysis, ldescription, assay.name,lextfilepath
		from library
			inner join sample2library on sample2library.lid=library.lid
			inner join exomehg19.sample on exomehg19.sample.idsample=sample2library.idsample
			inner join exomehg19.project on exomehg19.sample.idproject=exomehg19.project.idproject
			inner join exomehg19.organism on exomehg19.sample.idorganism=exomehg19.organism.idorganism
			inner join libtype on libtype.ltid=library.libtype
			inner join libpair on libpair.lpid=library.libpair
			left join assay on assay.idassay=library.idassay
			where exomehg19.sample.name='$sample' AND ltlibtype='$libtype' AND libpair.lplibpair='$libpair' AND library.lfailed=0 AND library.lstatus='external' AND library.lextfilepath IS NOT NULL
		";

	print "query: $query" if $printSQL;

	my $out = $dbh->prepare($query) || die print "$DBI::errstr\n";
	$out->execute || die print "$DBI::errstr\n";

	while (
		my (
			$lname,        $organism, $pname, $projectdiscr, $analysis,
			$ldescription, $lassay,   $lextfilepath
		)
		= $out->fetchrow_array
	  )    #read info from database for every flowcell of this sample
	{
		$count++;
		if ($debug) {
			$lname        = "" unless $lname;
			$organism     = "" unless $organism;
			$pname        = "" unless $pname;
			$projectdiscr = "" unless $projectdiscr;
			$analysis     = "" unless $analysis;
			$ldescription = "" unless $ldescription;
			$lassay       = "" unless $lassay;
			$lextfilepath = "" unless $lextfilepath;

			print "$lname, $organism, $projectdiscr, $analysis, $ldescription, $lassay, $lextfilepath\n";
		}
		else {
			# Debug info
			$lname        = "" unless $lname;
			$organism     = "" unless $organism;
			$pname        = "" unless $pname;
			$projectdiscr = "" unless $projectdiscr;
			$analysis     = "" unless $analysis;
			$ldescription = "" unless $ldescription;
			$lassay       = "" unless $lassay;
			$lextfilepath = "" unless $lextfilepath;
			print STDERR "\tExternal library found: $lname, $organism, $projectdiscr, $analysis, $ldescription, $lassay, $lextfilepath\n";
			
			my @files = split( ",", $lextfilepath );
			if(@files>1 && $libpair eq "paired-end"){		#if there are more than 1 files and it is a paired end library --> join R1 & R2 fastq files
				my @tmp;
				for(my $counter = 0; $counter < @files; $counter += 2){
					while(  !( $files[$counter] =~ /fastq\.gz$/ || $files[$counter] =~ /fastq$/ ) &&  $counter < @files){		#if there are files in there that are not fastq files (i.e. BAM files) don't join two of them
						push(@tmp,$files[$counter]);
						$counter++;
					}
					last if ($counter+1) >= @files;
					#--> join the two fastq files
					push(@tmp,$files[$counter].",".$files[$counter+1]);
				}
				@files=@tmp;
			}
			
			my $counter = 0;

			foreach my $file (@files) {
				$counter++;
				push( @prefixes, "$lname\_external\_$counter" );
				if ( $lextfilepath =~ /\.bam$/ ) {
					push( @filetypes, "BAMX" );		# External BAMs can have unexpected formats, BAMX filetype tells to the downstream chain to keep it into account
				}
				else {
					push( @filetypes, "FASTQ" );
				}
				push( @dirs,      $file );
				push( @lanes,     "" );
				push( @libraries, $lname );
			}
			
			$projectdiscr =~ s/\//_/g;
			$projectdiscr =~ s/\s//g;
			$pdescr       = $projectdiscr;
			$project	  = $pname if $project eq "";


			if ( $libtype eq "exomic" || $libtype eq "amplicon") {
				$exome   = 1;
				$usedkit = $lassay;
			}
			if ($libtype eq "RNA") {
				$usedkit = $lassay;
			}
			if ( $globOrganism eq "" ) {
				$globOrganism = lc $organism;
			}
			if ( $globalSettings eq "" ) {
				if($params->{settings}->{ lc $organism }->{$libtype}){
					$settings = $params->{settings}->{ lc $organism }->{$libtype}->{$version};		#TW 11.06.2014: there are now special entries in the config file for certain libtypes, e.g. amplicon
				}else{
					$settings = $params->{settings}->{ lc $organism }->{$version};
				}
			}
			else {
				$settings = $globalSettings;
			}
		}
	}


	# LOOK FOR RUNS CONTAINING SAMPLE:
	$query = "
		select rname, rdaterun, library.lname, exomehg19.organism.orname, exomehg19.project.pname , exomehg19.project.pdescription, kit.cdescription,exomehg19.sample.analysis, group_concat(DISTINCT lane.alane SEPARATOR ' '), group_concat(DISTINCT rfailed), ldescription, assay.name
		from run 
			left join lane on lane.rid=run.rid 
			left join pool on lane.idpool=pool.idpool 
			left join library2pool on library2pool.idpool=pool.idpool
			left join library on library2pool.lid=library.lid
			left join sample2library on sample2library.lid=library.lid
			left join exomehg19.sample on exomehg19.sample.idsample=sample2library.idsample
			inner join exomehg19.project on exomehg19.sample.idproject=exomehg19.project.idproject
			left join stock on library.lkit=stock.sid
			left join kit on stock.cid=kit.cid
			left join exomehg19.organism on exomehg19.sample.idorganism=exomehg19.organism.idorganism
			left join libtype on libtype.ltid=library.libtype
			left join libpair on libpair.lpid=library.libpair
			left join assay on assay.idassay=library.idassay
			where exomehg19.sample.name='$sample' AND ltlibtype='$libtype' AND libpair.lplibpair='$libpair' AND library.lfailed=0
		";
		
	if ($ignoreRundate) {
		$query .= "AND rdaterun<>'0000-00-00'\n";
	}
	if ( !$useFailedReads ) {
		$query .= "AND aread1failed='F' AND aread2failed='F'\n";
	}
	if ($ignoreFailed) {
		$query .= "AND (rfailed IS NULL OR rfailed='')\n";
	}

	$query .=
		"group by rname, library.lname, organism.orname,project.pdescription, kit.cdescription,exomehg19.sample.analysis, ldescription;";

	print "query: $query" if $printSQL;

	$out = $dbh->prepare($query) || die print "$DBI::errstr\n";
	$out->execute || die print "$DBI::errstr\n";

	while (
		my (
			$fc,           $rundate, $lname,    $organism, $pname,
			$projectdiscr, $kit,     $analysis, $lane,
			$rfailed,      $ldescription, $assay
		)
		= $out->fetchrow_array
	  )    #read info from database for every flowcell of this sample
	{
		$count++;
		if ($debug) {
			$fc           = "" unless $fc;
			$rundate      = "" unless $rundate;
			$lname        = "" unless $lname;
			$organism     = "" unless $organism;
			$pname		  = "" unless $pname;
			$projectdiscr = "" unless $projectdiscr;
			$kit          = "" unless $kit;
			$analysis     = "" unless $analysis;
			$lane         = "" unless $lane;
			$rfailed      = "" unless $rfailed;
			$assay		  = "" unless $assay;
			print
"$fc, $rundate, $lname, $organism, $projectdiscr, $kit, $analysis, $lane, $rfailed, $assay\n";
		}
		else {
			$fc           = "" unless $fc;
			$rundate      = "" unless $rundate;
			$lname        = "" unless $lname;
			$organism     = "" unless $organism;
			$pname		  = "" unless $pname;
			$projectdiscr = "" unless $projectdiscr;
			$kit          = "" unless $kit;
			$analysis     = "" unless $analysis;
			$lane         = "" unless $lane;
			$rfailed      = "" unless $rfailed;
			$assay		  = "" unless $assay;
			print STDERR "	Found data: $fc, $rundate, $lname, $organism, $projectdiscr, $kit, $analysis, $lane, $rfailed, $assay\n";
			
			push(@flowcells,$fc);
			my $fcpath = qx/ls -d $runfolder\/*$fc\//;

			chomp $fcpath;
			my $gerald = "";

			#############################################
			#29.03.2012: include the possibility that $runfolder points to a directory that contains bam converted flowcells
			#
			#
			$gerald = qx/ls -d $fcpath\/Sample_$sample  2> \/dev\/null/;

			#$gerald = $fcpath if $gerald eq "";			#TW 14.03.2014: doesn't make sense??
			chomp $gerald;
			

			my @inbams = glob( $gerald . "/*.bam" );

			if (@inbams) {    #Demultiplexed, sample structure
				              #get prefix
				if ( ( my $prefix = &parseBamHeader( $inbams[0] ) ) ne "" ) {
					if ($prefix =~ m/XX_.*/g) {
						$prefix =~ s/XX_.*/XX/;
						push( @prefixes,  $prefix );
					} elsif ($prefix =~ m/XY_.*/g) {
						$prefix =~ s/XY_.*/XY/;
						push( @prefixes,  $prefix );
					}
					push( @filetypes, "BAM" );
				}

			}
			else {
				if ( -e "$fcpath/SamplesDirectories.csv" )
				{             #Demultiplexed, bin structure
					open( IN, "$fcpath/SamplesDirectories.csv" );
					while (<IN>) {
						chomp;
						my @columns = split(",");

						if ( $columns[2] eq $sample ) {
							$gerald = "$fcpath/$columns[-1]";

							my @inbams = glob( $gerald . "/*.bam" );
							if (@inbams) {
								if (
									(
										my $prefix =
										&parseBamHeader( $inbams[0] )
									) ne ""
								  )
								{
									if ($prefix =~ m/XX_.*/g) {
										$prefix =~ s/XX_.*/XX_$columns[-1]/;
										push( @prefixes,  $prefix );
									} elsif (($prefix =~ m/XY_.*/g)) {
										$prefix =~ s/XY_.*/XY_$columns[-1]/;
										push( @prefixes,  $prefix );
									}
									push( @filetypes, "BAM" );
								}
							}
							else {
								$gerald = "";
							}
							last;
						}
					}
					close IN;
				}
				else {    #Not demultiplexed
					$gerald = $fcpath;
					my @inbams = glob( $gerald . "/*.bam" );
					if (@inbams) {
						if ( ( my $prefix = &parseBamHeader( $inbams[0] ) ) ne
							"" )
						{
							if ($prefix =~ m/XX_.*/g) {
								$prefix =~ s/XX_.*/XX/;
								push( @prefixes,  $prefix );
							} elsif ($prefix =~ m/XY_.*/g) {
								$prefix =~ s/XY_.*/XY/;
								push( @prefixes,  $prefix );
							}
							push( @filetypes, "BAM" );
						}
						else {
							$gerald = "";
						}
					}
					else {
						$gerald = "";
					}
				}

			}

			##############################################
			# If the given $runfolder doesn't include .bam files, search it for the usual positions of Illumina sequence.txt or .fastq.gz files
			#
			# starting with Casava 1.8 there exists a folder for each sample instead of a folder for all samples
			if ( $gerald eq "" ) {
				$gerald =
qx/ls -d $fcpath\/Demultiplexed\/Project_all\/Sample_$sample*  2> \/dev\/null/;               #TODO: TS add "_" for 10x runs, i.e. .../Sample_$sample*\/ ???
				if ( $gerald =~ /^\/.*\/(.+$fc)\// ) {
					push( @prefixes,  $1 );
					push( @filetypes, "FASTQ" );
				}
			}
			#07.03.2014 TW: adding support for MiSeq --> Demultiplexed files can be found under Data/Intensities/BaseCalls
			if ($gerald eq ""){
				my $sampleFASTQ = qx/ls  $fcpath\/Data\/Intensities\/BaseCalls\/$sample*_R1_001.fastq.gz  2> \/dev\/null/;
				chomp $sampleFASTQ;
				
				if ( ! -e $sampleFASTQ ){		# MiSeq demultiplex substitutes in the fastq name the "_" with "-"
					my $sampledeunderlined=$sample;
					$sampledeunderlined =~ s/\_/\-/g;
					$sampleFASTQ = qx/ls  $fcpath\/Data\/Intensities\/BaseCalls\/$sampledeunderlined*_R1_001.fastq.gz  2> \/dev\/null/;
					chomp $sampleFASTQ;
				}
					
				if(-e $sampleFASTQ){
					$gerald = dirname($sampleFASTQ);
					if ( $gerald =~ /^\/.*\/(.+$fc)\// ) {
						push( @prefixes,  $1 );
						push( @filetypes, "FASTQ" );
					}
				}
			}

			if ( $gerald eq ""
				&& -e "$fcpath/Data/Intensities/BaseCalls/Demultiplexed/SamplesDirectories.csv"
			  )
			{    #if this file exists the flowcell was indexed.
				open( IN,
"$fcpath/Data/Intensities/BaseCalls/Demultiplexed/SamplesDirectories.csv"
				);
				while (<IN>) {
					chomp;
					my @columns = split(",");

					# At the moment: 4 cases where qseq files can lie

					if ( $columns[2] eq $sample ) {
						$gerald =
qx/ls -d $fcpath\/Data\/Intensities\/BaseCalls\/Demultiplexed\/\/$columns[-1]\/GERALD*  2> \/dev\/null/
						  ;    #Case 1: Demultiplexed
						if ( $gerald =~ /^\/.*\/(.+$fc)\// ) {
							push( @prefixes,  $1 . "_$columns[-1]" );
							push( @filetypes, "ILLUMINA" );
							last;
						}
					}
				}

			}
			if ( $gerald eq "" ) {
				$gerald =
qx/ls -d $fcpath\/Data\/Intensities\/Bustard*\/GERALD*  2> \/dev\/null/
				  ;    #Case 2: RTA went wrong
				if ( $gerald =~ /^\/.*\/(.+$fc)\// ) {
					push( @prefixes,  $1 );
					push( @filetypes, "ILLUMINA" );
				}
			}

			if ( $gerald eq "" ) {
				$gerald =
qx/ls -d $fcpath\/Data\/*Firecrest*\/Bustard*\/GERALD*  2> \/dev\/null/
				  ;    #Case 3: RTA went wrong
				if ( $gerald =~ /^\/.*\/(.+$fc)\// ) {
					push( @prefixes,  $1 );
					push( @filetypes, "ILLUMINA" );
				}
			}

			if ( $gerald eq "" ) {    #Case 4: no demultiplexing, RTA ok
				$gerald =
qx/ls -d $fcpath\/Data\/Intensities\/BaseCalls\/GERALD* 2> \/dev\/null/;

#next if $gerald eq "";
# --> no GERALD directory found --> Illumina pipeline has not yet run or sequencing is in progress

				if ( $gerald =~ /^\/.*\/(.+$fc)\// ) {
					push( @prefixes,  $1 );
					push( @filetypes, "ILLUMINA" );
				}
			}

			chomp $gerald;
			if ( $gerald ne "" ) {
				push( @dirs,      $gerald );
				push( @lanes,     $lane );
				push( @libraries, $lname );
			}

			#print "test: $fc $fcpath $gerald @prefixes\n";
			#print "test: $project; $projectdiscr\n";
			$projectdiscr =~ s/\//_/g;
			$projectdiscr =~ s/\s//g;
			$pdescr       = $projectdiscr;
			$project	  = $pname if $project eq "";
#			if ( $project eq "" ) {
#				$projectdiscr =~ s/\//_/g;
#				$projectdiscr =~ s/\s//g;
#				$project = $projectdiscr;
#			}

			#print "test2: $project; $projectdiscr\n";

			if ( $libtype eq "exomic" || $libtype eq "amplicon" || $libtype eq "MIP" ) {
				$exome = 1;
				
				if($assay && $assay ne ""){			#TW 07.03.2014 --> if assay is defined in library table --> take it directly
					$usedkit = $assay;
				}else{
					unless ($kit) {
						$kit =
						  $ldescription
						  ; # for older samples in the library the kit is not correctly defined but is just in the description field
					}
	
					#print "test: $kit\n";
	
					if ( $kit && $usedkit eq "" ) {
						if ( $kit =~ m/v4/i ) {
							$usedkit = "SureSelect50Mbv4";
						}
						elsif ( $kit =~ m/v5/i ) {
							$usedkit = "SureSelect50Mbv5";
						}
						elsif ( $kit =~ m/v6/i ) {
							$usedkit = "SureSelect60Mbv6";
						}						
						elsif ( $kit =~ m/50\s*Mb/i ) {
							$usedkit = "SureSelect50Mb";
						}
						elsif ( $kit =~ m/Mouse/i ) {
							$usedkit = "SureSelectMouse50Mb";
						}
						elsif ( $kit =~ m/38\s*Mb/i ) {
							$usedkit = "SureSelect38Mb";
						}
	
						else {
	
							print
	"\nNo kit defined in DB for sample $sample ($projectdiscr,$fc,$rundate)\n(0) SureSelect60Mbv6\n(1) SureSelect50Mbv5\n(2) SureSelect50Mbv4\n(3) SureSelect50Mb\n(4) SureSelectMouse50Mb\n(5) SureSelect38Mb\nChoose kit manually [5]: ";
							my $selection = <>;
							chomp $selection;
							$selection = 5 if $selection eq "";
	
							if ( $selection == 0 ) {
								$usedkit = "SureSelect60Mbv6";
							} 
							elsif ( $selection == 1 ) {
								$usedkit = "SureSelect50Mbv5";
							}
							elsif ( $selection == 2 ) {
								$usedkit = "SureSelect50Mbv4";
							}
							elsif ( $selection == 3 ) {
								$usedkit = "SureSelect50Mb";
							}
							elsif ( $selection == 4 ) {
								$usedkit = "SureSelectMouse50Mb";
							}
							else {
								$usedkit = "SureSelect38Mb";
							}
						}
					}
				}
			}
			if ($libtype eq "RNA") {
				if($assay && $assay ne "") {
					$usedkit = $assay;
				}
			}
				
			$organism =~ s/\s//g;

			if ( $globOrganism eq "" ) {
				$globOrganism = lc $organism;
			}
			if ( $globalSettings eq "" ) {
				if(defined $params->{settings}->{ lc $organism }->{$libtype}){
					if(defined $params->{settings}->{ lc $organism }->{$libtype}->{$version}){
						$settings = $params->{settings}->{ lc $organism }->{$libtype}->{$version};
						print STDERR "\tUsing settings: $settings\n";
					}else{
						print STDERR "No settings defined for organism: $organism - libtype: $libtype - version: $version: Nothing done for sample $sample\n"; 	#TW 24.02.2016: if the current combination of organism/libtype/version is not defined --> don't do anything for this sample!
						#return;
						$justCreateBAM = 1;
						print STDERR "Just create BAM from FASTQ files in project folder\n";
					}
							
				}elsif(defined $params->{settings}->{ lc $organism }->{$version}){
					$settings = $params->{settings}->{ lc $organism }->{$version};
					print STDERR "\tUsing settings: $settings\n";
				}else{
					print STDERR "No settings defined for organism: $organism - libtype: $libtype - version: $version: Nothing done for sample $sample\n"; 	#TW 24.02.2016: if the current combination of organism/libtype/version is not defined --> don't do anything for this sample!
					#return;
					$justCreateBAM = 1;
					print STDERR "Just create BAM from FASTQ files in project folder\n";
				}
				
			}
			else {
				$settings = $globalSettings;
				print STDERR "\tUsing settings: $globalSettings\n";
			}

		}

	}

	unless ( $debug || $count == 0 ) {

		if ( $maxInsertSize == -1 ) {
			if ( $libpair =~ /\s(\d+)\skb$/ ) {
				$maxInsertSize = $1 * 2000;
			}
			else {
				$maxInsertSize = $defaultInsertSize;
			}
		}

		my $failedlanes = "";
		unless ($dontRemoveFailed) {

		  #get failed lanes/libraries to remove the files from the output folder
			my $query = "
				select rname, lane.alane
				from run 
				left join lane on lane.rid=run.rid 
				left join pool on lane.idpool=pool.idpool 
				left join library2pool on library2pool.idpool=pool.idpool
				left join library on library2pool.lid=library.lid
				left join sample2library on sample2library.lid=library.lid
				left join exomehg19.sample on exomehg19.sample.idsample=sample2library.idsample
				inner join exomehg19.project on exomehg19.sample.idproject=exomehg19.project.idproject
				left join stock on library.lkit=stock.sid
				left join kit on stock.cid=kit.cid
				left join exomehg19.organism on exomehg19.sample.idorganism=exomehg19.organism.idorganism
				left join libtype on libtype.ltid=library.libtype
				left join libpair on libpair.lpid=library.libpair
				where exomehg19.sample.name='$sample' AND ltlibtype='$libtype' AND libpair.lplibpair='$libpair' AND (library.lfailed=1
			";
			if ( !$useFailedReads ) {
				$query .= "OR (aread1failed='T' AND aread2failed='T')\n";
			}
			$query .= ") group by rname, lane.alane;";
			print "query: $query" if $printSQL;
			$out = $dbh->prepare($query) || die print "$DBI::errstr\n";
			$out->execute || die print "$DBI::errstr\n";

			while ( my ( $fc, $lane ) = $out->fetchrow_array ) {
				$failedlanes .= $fc . "_s_" . $lane . ",";
			}
			$failedlanes =~ s/,$//;
		}

		$libpair =~ s/\s//g;
		
		
		#TW 24.02.2015: if $projectfolder is not given --> set based on $version
		if($projectfolder eq ""){
			if($settings eq "hg38") {
				$projectfolder = $params->{dirs}->{hg38projectfolder};
			} elsif($version eq "gatk"){
				$projectfolder = $params->{dirs}->{projectfolder} . "plus";
			}else{
				$projectfolder = $params->{dirs}->{projectfolder};      #default
			}
		}
		
		make_path( "$projectfolder/$project/$sample", { mode => 0775 } );
		symlink("$project","$projectfolder/$pdescr"); 

		open( OUT,
			">$projectfolder/$project/$sample/$libtype.$libpair.pipeline.ini" )
		  or die(
			"Error opening $projectfolder/$project/$sample/$libtype.$libpair.pipeline.ini\n"
		  );

		#print "outfolder: $projectfolder/$project/$sample/$libtype.$libpair.pipeline.ini\n";
		my $folder    = "$libtype$params->{dirs}->{out}";
		my $subfolder = "$libpair$params->{dirs}->{out}";

		if ($clearOut) {
			print "Deleting all files in outputfolder: rm \"$projectfolder/$project/$sample/$folder/$subfolder/*\"\n";    #clear out folder
			system("find $projectfolder/$project/$sample/$folder/$subfolder/ -type f -delete");
			system("find $projectfolder/$project/$sample/$folder/$subfolder/ -type l -delete");
		}elsif ($clearFCfiles){
			foreach my $currFC (@flowcells){
				my $command;
				if(-l "$projectfolder/$project/$sample/$folder/$subfolder/merged.bam"){		#if merged.bam is a symlink --> only one .sort.bam file should exist where it points do --> don't delete this
					$command = "find $projectfolder/$project/$sample/$folder/$subfolder/ -name \"*$currFC*\" | grep -v \".sort.bam\" | xargs -i rm {} ";			
				}else {
					$command = "rm $projectfolder/$project/$sample/$folder/$subfolder/*$currFC*";
				}
				print "Deleting files for flowcell $currFC: $command\n";
				system($command);
			}
		}
		
		for ( my $i = 0 ; $i < @dirs ; $i++ ) {
			print OUT "indir      :  $dirs[$i]\n";
			print OUT "lane       :  $lanes[$i]\n";
			print OUT "prefix     :  $prefixes[$i]\n";
			print OUT "library    :  $libraries[$i]\n";
			print OUT "#type of the inputfiles; possible values \"ILLUMINA\",\"FASTQ\" or \"BAM\"\n";
			print OUT "infiletype :  $filetypes[$i]\n";
		}
		
		if ( $failedlanes ne "" ) {
			print OUT "#if there are failed lanes/libraries in the database put them here so failed files can be deleted\n";
			print OUT "failedlanes:  $failedlanes\n";
		}
		
		print OUT "outdir     :  $projectfolder\n";

		print OUT "project    :  $project\n";
		print OUT "sample     :  $sample\n";

		print OUT "folder     :  $folder\n";
		print OUT "subfolder  :  $subfolder\n";
		print OUT "settings   :  $settings\n";
		print OUT "organism   :  $globOrganism\n";
		print OUT "#max insert size for bwa sampe\n";
		print OUT "maxinssize :  $maxInsertSize\n";
		print OUT "#version of the pipeline; currently \"vcf\" and \"gatk\" are valid options\n";
		print OUT "version    :  $version\n";
		print OUT "strandedRNA	:	$isStrandedRNA\n";
		print OUT "tmpArg	  :  $tmpArgument\n";
		print OUT "removeAdapters : ".($removeAdapters ? "TRUE" : "FALSE" )."\n";
		print OUT "removePhix : ".($removePhix ? "TRUE" : "FALSE")."\n";
		print OUT "normdup    :  ";
		if($dontRmdup){
			print OUT "TRUE\n";
		}else{
			print OUT "FALSE\n";
		}
		print OUT "dotrim     :  ";

		if($globAligner eq ""){
			if($libtype eq "RNA"){
				$aligner = "star";
			}else{
				if($version eq "vcf"){
					$aligner = "bwa";
				}elsif($version eq "gatk"){
					$aligner = "bwamem";
				}
				
			}
		}else{
			$aligner = $globAligner;
		}

		if ( $doTrim && $aligner eq "gem" ) {
			print OUT "TRUE\n";
		}
		else {
			print OUT "FALSE\n";
		}
		print OUT "fastmerge	: ";
		if ( $mergeFast && $aligner eq "gem" ) {
			print OUT "TRUE\n";
		}
		else {
			print OUT "FALSE\n";
		}

		if ( $exome == 1 || $libtype eq "MIP" ) {
			if ( $usedkit ne "" ) {
				if( ref $params->{targets}->{$usedkit} eq ref {} && $params->{targets}->{$usedkit}->{$settings}){
					print OUT "target     :  $params->{targets}->{$usedkit}->{$settings}\n";
				}else{
					print OUT "target     :  $params->{targets}->{$usedkit}\n";
				}
				print OUT "#window that should be added around each target region for variant calling\n";
				print OUT "targetWin  :  $params->{windows}->{$usedkit}\n";
				print OUT "assay      :  $usedkit\n";
			}
		}elsif($libtype eq "genomic" && $globOrganism eq "human" && ($settings eq "hg19_plus")){
			print OUT "target     :  $params->{targets}->{genomic}\n";
			print OUT  "#window that should be added around each target region for variant calling\n";
			print OUT "targetWin  :  $params->{windows}->{genomic}\n";
			#print OUT #"#region to filter pindel calls (different in genomic samples)\n";
			#print OUT "pindelreg  :  $params->{targets}->{genomicpindel}\n";
		}elsif($libtype eq "mtDNA") {
			print OUT "target     :  $params->{targets}->{mtDNA}\n";	
		}elsif($libtype eq "genomic" && $globOrganism eq "mouse" && ($settings eq "mm10")){
			print OUT "target     :  $params->{targets}->{genomicmouse}\n";
			print OUT "#window that should be added around each target region for variant calling\n";
			print OUT "targetWin  :  $params->{windows}->{genomicmouse}\n";
		}
		if ($libtype eq "RNA" && $usedkit ne "") {
			print OUT "assay      :  $usedkit\n";
		}

		print OUT "aligner    :  $aligner\n";
		
		print OUT "\n## Parts of the pipeline that should be run\n";
		my $arrayRef;

		if ($justCreateBAM && $version eq "vcf") {
			# if no settings are found for a sample, just convert the FASTQ- to a BAM-file in the according project folder 
			# so that the automated file download for collaborators can be standardized
			make_path( "$projectfolder/$project/$sample/$folder/$subfolder", { mode => 0775 } );

			my $command = "perl $prog_path/$params->{scripts}->{sequence2BAMpl} -b -f $flowcell -o $projectfolder/$project/$sample/$folder/$subfolder/ -r $runfolder -sge $sgequeue -css -ssn $sample -lt $libtype";
			system($command);
			if ($libtype ne "scATAC-Seq") {
				$command = "perl $prog_path/bwa_merge_gatk.pl -o $projectfolder/$project/$sample/$folder/$subfolder/ -a fastq2bam -s $sample -d";
				
				my $jobName = $sample . "_FASTQ2BAM_bwamerge";
				#generate qsub string
				my $qsubstring = "qsub -q $sgequeue -N $jobName -hold_jid \"BAM" . $sample ."*\" ";
				#create temp script to submit
				open( OUT, ">$projectfolder/$project/$sample/$folder/$subfolder/sgeSubmTmp.sh" )|| exit print("Cannot open $projectfolder/$project/$sample/$folder/$subfolder/sgeSubmTmp.sh");
				print OUT "#\$  -S /bin/bash\n";
				print OUT $command;
				print OUT "\nif [ \$? -ne 0 ]
then
  $command
else";
				close OUT;
				$qsubstring .= "$projectfolder/$project/$sample/$folder/$subfolder/sgeSubmTmp.sh";
				print "Submitting: $qsubstring\n";
				system($qsubstring);
				#system("rm $projectfolder/$project/$sample/$folder/$subfolder/sgeSubmTmp.sh");		#TODO: uncomment this line!!!
			}
		} else {
			if ( $runFile eq "" ) {
				$arrayRef =
				  $params->{runs}->{alignmentOnly}
				  ->{run};    #insert alignmentOnly runs
				if ($noAlignment) {
					print OUT "#run  :  $arrayRef\n";
				}
				else {
					print OUT "run  :  $arrayRef\n";
				}
	
				$arrayRef =
				  $params->{runs}->{standard}->{run};    #insert merge & stats
				my @standardRuns = @$arrayRef;
				if ( $align == 1 || $noMerge ) {
					foreach (@standardRuns) {
						print OUT "#run  :  $_\n";
					}
				}
				else {
					foreach (@standardRuns) {
						print OUT "run  :  $_\n";
					}
				}
				
				$arrayRef =
				  $params->{runs}->{genome}->{run};    #insert genomic only runs
				my @genomicRuns = @$arrayRef;
				if(($libtype eq "genomic" ) && (!$noVariantCalling)){
					foreach (@genomicRuns) {
						print OUT "run  :  $_\n";
					}				
				}else{
					foreach (@genomicRuns) {
						print OUT "#run  :  $_\n";
					}
				}
				
				$arrayRef =
				  $params->{runs}->{genomedb}->{run};    #insert genomic only runs - DB insertion
				my @genomicDBRuns = @$arrayRef;
				if(($libtype eq "genomic" ) && ( (!$noVariantCalling) && (!$noDB) )){
					foreach (@genomicDBRuns) {
						print OUT "run  :  $_\n";
					}				
				}else{
					foreach (@genomicDBRuns) {
						print OUT "#run  :  $_\n";
					}
				}
	
				$arrayRef =
				  $params->{runs}->{exome}->{run};    #insert exome runs, if needed
				my @exomeRuns = @$arrayRef;
				if ( (($exome == 1 || ($libtype eq "genomic" && ($globOrganism eq "human" || $globOrganism eq "mouse" )) || $libtype eq "MIP" || ($libtype eq "RNA" && $callRNAVariants == 1 && $globOrganism eq "human")) && $merge == 0 && $align == 0 ) && (!$noVariantCalling) ){
					foreach (@exomeRuns) {
						print OUT "run  :  $_\n";
					}
				}
				else {
					foreach (@exomeRuns) {
						print OUT "#run  :  $_\n";
					}
				}
				
				$arrayRef =
				  $params->{runs}->{exomedb}->{run};    #insert exome runs, if needed - DB insertion
				my @exomeDBRuns = @$arrayRef;
				if ( (($exome == 1 || ($libtype eq "genomic" && ($globOrganism eq "human" || $globOrganism eq "mouse" )) || $libtype eq "MIP" || ($libtype eq "RNA" && $callRNAVariants == 1 && $globOrganism eq "human")) && $merge == 0 && $align == 0 ) && ( (!$noVariantCalling) && (!$noDB) ) ){
					foreach (@exomeDBRuns) {
						print OUT "run  :  $_\n";
					}
				}
				else {
					foreach (@exomeDBRuns) {
						print OUT "#run  :  $_\n";
					}
				}
				
				
				#ChIP-Seq
				$arrayRef =
				  $params->{runs}->{chipseq}->{run};
				
				if ( ( $libtype eq "ChIP-Seq" || $libtype eq "xC-Seq" || $libtype eq "ATAC-Seq") && (!$noVariantCalling) ) {
					print OUT "run  :  $arrayRef\n";
				}
	
				
				#RNA-seq stuff
				$arrayRef =
				  $params->{runs}->{rnaseq}->{run};
				my @rnaseqRuns = @$arrayRef;
				if ( ( $exome == 0 && $noCount == 0 && $libtype eq "RNA" ) && (!$noVariantCalling) ) {
					foreach (@rnaseqRuns) {
						print OUT "run  :  $_\n";
					}
				}
				else {
					foreach (@rnaseqRuns) {
						print OUT "#run  :  $_\n";
					}
				}
				
		
				#mtDNA stuff
					$arrayRef =
				  $params->{runs}->{mtDNA}->{run};    #insert exome runs, if needed
				my @mtDNARuns = @$arrayRef;
				if ( ( $libtype eq "mtDNA" ) && (!$noVariantCalling) ) {
					foreach (@mtDNARuns) {
						print OUT "run  :  $_\n";
					}
				}
	
				
				
				
			}
			else {
				if ( $runFile eq "STDIN" )
				{    #don't take scripts to run from config file but from <STDIN>
					while (<STDIN>) {
						print OUT $_;
					}
				}
				elsif ($runFile ne "EMPTY") {
					open( RUNS, $runFile )
					  or die "Can't open runFile: $runFile!\n"
					  ; #don't take scripts to run from config file but from a given file
					while (<RUNS>) {
						print OUT $_;
					}
					close RUNS;
				}
			}
		}

		close OUT;
		
		if ($justCreateBAM) {
			return;
		}
		
		#if ChIP-seq sample, check if control sample is available:
		my $query = "SELECT s.chipseqcontrol 
						from exomehg19.sample s 
						inner join solexa.sample2library sl on sl.idsample=s.idsample
						inner join solexa.library l on sl.lid=l.lid
						left join solexa.libtype lt on lt.ltid=l.libtype
						where s.name = '$sample'
						and lt.ltlibtype='$libtype';";
		$out = $dbh->prepare($query) || die print("$DBI::errstr");
		$out->execute || die print("$DBI::errstr");
		if ($out->rows == 1) {
			$dependsOnSample = $out->fetchrow_array;
			unless (defined $dependsOnSample) {
				$dependsOnSample = "";
			}
		}

		if ( ($parseConfig == 1 || $sgequeue ne "") && $runFile ne "EMPTY" ) {
			if ($debug) {
				print
"perl $prog_path/$params->{scripts}->{parseConfig} $projectfolder/$project/$sample/$libtype.$libpair.pipeline.ini\n";
			}
			else {
				system(
"perl $prog_path/$params->{scripts}->{parseConfig} $projectfolder/$project/$sample/$libtype.$libpair.pipeline.ini"
				);
			}

		}

		if ( $jobCount >= $skipSleep )
		{    #start the first "skipSleep" jobs right now
			$addedSleepTime +=
			  $sleep;    #add "sleep" minutes to the starttime of all other jobs
		}
		my $starttime = &getStartTime($addedSleepTime);

		if ( $sgequeue ne "" && $runFile ne "EMPTY") {
			my $currentlog;
			if ( $logfile eq "" ) {
				$currentlog = "$projectfolder/$project/$sample/pipeline.log";
			}
			else {
				$currentlog = $logfile;
			}
			
			my $command =
"perl $prog_path/parallelpipeline.pl -i $projectfolder/$project/$sample/$libtype.$libpair.pipeline.cfg -sge $sgequeue -lf $currentlog -ll $loglevel -sl $sleep -pe $parallelEnv ";
			$command .= "-pri \"$priority\" " if $priority ne "";
			$command .= "-do $dependsOnSample " if ($dependsOnSample ne "");

			if ($isFlowcell) {

#if this script was called for a flowcell
#calculation of start times gets somehow complicated
#lets say -n 3 -sl 10
#this means that the first three samples start right now and the others start in 10,20,... minutes
#but in these three samples only 1 job starts at the beginning and the next ones also in 10,20,... minutes
				if ( $skipSleep != 0 ) {

					#201107211137 --> 21.07.2011,11:37
					#reformat starttime
					my $formatedTime =
					    substr( $starttime, 6, 2 ) . "."
					  . substr( $starttime, 4,  2 ) . "."
					  . substr( $starttime, 0,  4 ) . ","
					  . substr( $starttime, 8,  2 ) . ":"
					  . substr( $starttime, 10, 2 );
					$command .= "-st $formatedTime -n 1 ";
				}

			}
			else {
				if ( $firstTime ne "" ) {
					$command .= "-st $firstTime ";
				}
				if ( $skipSleep != 0 ) {
					$command .= "-n $skipSleep ";
				}
			}

			if ($debug) {
				$command .= "-d ";
			}

			#print "test: $command\n";
			system($command);

		

		$jobCount++;
		}
	}
}

#################################################################################################
sub parseBamHeader {
	my $bamfile = shift;
	open HEAD, "samtools view -H $bamfile |"
	  or die
	  "Found BAM file $bamfile, but can't read header!\n"; #fixed samtools path!
	while (<HEAD>) {
		chomp;
		if ( $_ =~ /^\@RG/ ) {
			if ( $_ =~ /ID:(.*?)\s/ ) {

				#my $prefix = $1;
				#$prefix =~ s/XX_.*/XX/;
				return $1;
			}
		}
	}
	return "";
}

#################################################################################################
sub getStartTime {
	my $minOffset = shift;

	my $dt;
	my $tz = DateTime::TimeZone->new( name => 'local' );
	if ( $firstTime eq "" ) {
		$dt = DateTime->now;
		$dt->set_time_zone($tz);
	}
	else {
		my ( $date, $time ) = split( ",", $firstTime );
		my ( $day, $month, $year ) = split( /\./, $date );
		my ( $hour, $minute ) = split( ":", $time );
		$dt = DateTime->new(
			year      => $year,
			month     => $month,
			day       => $day,
			hour      => $hour,
			minute    => $minute,
			time_zone => $tz
		);
	}

	$dt->add( minutes => $minOffset );

	my $hour   = sprintf( "%02d", $dt->hour );
	my $minute = sprintf( "%02d", $dt->minute );
	my $date   = $dt->ymd("");

	return "$date$hour$minute";
}

=head1 NAME

initAnalysis.pl

=head1 SYNOPSIS

 initAnalysis.pl -fs <FLOWCELLNAME> -sge custom.q -b /data/runs/Runs/SequenceBams/

=head1 DESCRIPTION

This script is the central script of the pipeline. It retrieves all the information that is needed to
start a new analysis from the database and is controlled by the parameters shown below and the file
"conf.initAnalysis.xml" in the same folder as this script.
The standard way to start the analysis of new flowcells is to retrieve their names with the script
"checkForNewRuns" and forward them to initAnalysis. 

The current status of running analysis can then be viewed by typing: checkForNewRuns -s

=head1 OPTIONS

 Required:
 -fs	<flowcellname> or <samplename>; give a filename containing a list of samples/flowcells or choose "STDIN" to take a list of samples/flowcells from <STDIN>
 -v	<version> of the pipeline to run; can be "vcf" or "gatk"; decides the following options:
     -) used variant caller (this is hardcodeded in "pipelincfg_gatk.pl")
     -) together with organism & library type: used settings (can be overwriten by specifying -se)
     -) together with library type: aligner --> bwa for genomic DNA using "vcf" and bwa mem using "gatk"; star for RNA using "vcf" (not specified for "gatk")
        (can be overwritten by specifying -al)
     -) the default project folder: /.../isilon/..../exomehg19 for "vcf" & /.../isilon/..../exomehg19plus for "gatk"
        (can be overwritten by specifying -p)
 -lt	<libtype>; only if -fs contains a sample id, possible entries: exomic, genomic, RNA, ChIP-Seq, 'ChIP-Seq RNA', hMeDIP, mtDNA, amplicon, MIP, ATAC-Seq
 
 Optional:
 -p	<projectfolder> folder where project folders will be created,
 -f	entry from "-fs" is a flowcell name; otherwise the scripts tries to determine if its a sample or a flowcell name automatically by structure and length
 -lp	<libpair>; only if -fs contains a sample id, possible entries: single-end, paired-end, 'mate pair 3 kb', 'mate pair 10 kb', 'paired-end nextera', 'mate pair 1 kb'; default: paired-end
 -se	<settings> normally settings are taken from the config file according to the species of each sample and the pipeline version (-v). This option overwrites these settings. e.g. to use a testdatabase
 -r	<runfolder>; default: /data/runs/Runs/
 -c	<configfile>; default: conf.initAnalysis.xml in current folder
 -mi	<max_insert_size>; taken from library pair settings; e.g. if 'mate pair 3 kb' --> max insert size 6000; otherwise default: 500
 -i	ignore failed description; i.e. only use lanes where the failed description is empty
 -u	use failed reads; i.e. use also lanes where the "failed" flag in the database has been set
 -rd	ignore flowcells with a rundate '0000-00-00'
 -s	<lanesToSkip> only valid for flowcells!; format: 1,8,3-5
 -a	<lanesAlignmentOnly> perform only alignment; if -fs is only a sample, just give "1" to choose this option; format: 1,8,3-5
 -m	<lanesMergeOnly> perform only alignment and merge; if -fs is only a sample, just give "1" to choose this option; format: 1,8,3-5
 -rf	<Run File> file that contains scripts to run (instead of scripts from the config file); choose "STDIN" to take the list from <STDIN> (only if -fs is not STDIN); 
 	choose "EMPTY" if you don't want to run anything (useful to run only the BAM conversion for a flowcell) 
 -noal	don't align files; i.e. start with merging
 -nom	don't perform merging and stat computation; i.e. start with variant calling
 -novar just alignment and stats but no variant calling
 -novarDB do everything but do not insert variants into the database
 -nocnt for RNA samples only; don't perform htseq count and FPKM calculations
 -normdup	don't remove duplicates after merging (e.g. for amplicons)
 -removeAdapters	perform adapter clipping
 -removePhix	remove read pairs in which one or both align on Phix (use only with gatk)
 -var	run variant calling parts only (includes -noal, -nom and -nosdb); useful to run vcf pipeline directly after normal pipeline
 -co	remove all files from OUT folder; useful if pipeline was run before and files should be cleaned
 -cof	remove all files including the flowcell name(s); i.e. all files generated at initial alignment. Useful to free diskspace. 
 -drf	don't remove failed files. Usually the pipeline removes aligned bam files that originate from failed lanes/libraries to exclude them from further analysis.
 -pc	run parse config script for each pipeline.ini file
 -sge	<sge_queue> submit pipeline job to specified SGE queue
 -pe	<SGE parallel environment> parallel environment that is used to submit multihtreaded jobs (number of threads is defined in conf.initanalysis.xml); default: pipeline
 -pri	priority of the new jobs. Overwrites standard values (-10,-100). Only negative values are possible. NOTE: negative values must be in apostrophes.
 -st	time at which the job(s) should start; format: DD.MM.YYYY,hh:mm; default: now
 -sl	<minutes_to_sleep> time the script waits until it submits the next SGE job, default: 0
 -n	number of jobs that should be started at the beginning without sleeping (i.e. there is no sleeping time between them)
 -re	include standard removes at the end of the pipeline
 -b	<bam_dir> if this dir is given the sequence.txt files from the given flowcell get converted into .bam files via a SGE job; only works if -fs
	contains a flowcellname and -sge is given
 -lf	<logfile.log> specify a central log file; otherwise use a pipeline.log file in the output directory of each sample
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -d	don't change or create anything, just print info
 -sql	print SQL queries
 -al	aligner; possible values "bwa", "bwamem" and "gem"; default: determined by libtype: DNA --> bwa; RNA --> star
 -sr	specify if input are stranded RNA-seq data (default: no)
 -T	just for gem alignment; if specified, alignment will be performed: 1st with full length reads, 2nd with 20bp trimmed reads and 3rd with 5bp and 20bp trimmed reads; resulting bam files will be merged at the end
 -fastmerge just for gem alignment and "-T" option, performs the merging of the trimmed results faster
 -crv	perform variant calling for RNA-seq sample (it is recommended that alignment is performed with STAR); only possible if sample libtype is RNA; default: no
 -h	print this helptext
 -man	show man page

=head1 AUTHOR

Thomas Wieland
Riccardo Berutti

=cut
