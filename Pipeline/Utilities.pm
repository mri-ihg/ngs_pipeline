package Utilities;

use strict;
use XML::Simple;
use File::Basename;
use Sys::Hostname;
use DBI;
use Log::Log4perl qw(get_logger :levels);
use Cwd qw(abs_path);
use Bio::DB::Fasta;

my $prog_path = dirname(abs_path(__FILE__));

my $logger;
my $fasta;


# XML::Parser instead of XML::Simple needed
# to parse correctly:
# current.config.xml
$XML::Simple::PREFERRED_PARSER='XML::Parser';


#######################################################################################
#initLogger: Initialise the log4perl environment
#######################################################################################
sub initLogger{
	my $path  = shift;
	my $loglevel = shift;
	
	
	my $level = $INFO;
	if($loglevel eq "DEBUG"){
		$level = $DEBUG;
	}elsif($loglevel eq "ERROR"){
		$level = $ERROR;
	}
	
	
	$logger = get_logger();
    $logger->level($level); 
     
    # Appenders
    my $appender;
    
    if($path eq "SCREEN"){
    	#$appender = Log::Log4perl::Appender->new("Log::Dispatch");
		$appender = Log::Log4perl::Appender->new("Log::Log4perl::Appender::ScreenColoredLevels");
    }else{
    	$appender = Log::Log4perl::Appender->new(
	    	"Log::Dispatch::File",
	    	filename => $path,
	    	mode => "append"
	    );
    }
    
   
    $logger->add_appender($appender);
     
    # Layouts
    my $layout = 
    	Log::Log4perl::Layout::PatternLayout->new(
    	"%d %p (%H)> %F{1}:%L %M - %m%n");
    $appender->layout($layout);
	
}

#######################################################################################
#getLogger: Return log4perl object
#######################################################################################
sub getLogger{
	return $logger;
}


#######################################################################################
#getParams: Read in  a config file
#######################################################################################
sub getParams {
	my $file = "";

	if(@_){
		$file = shift;
	}
	
	
	if($file eq ""){
			
		if(-e $prog_path."/current.config.xml"){
		
		
			$file = $prog_path."/current.config.xml";
		}else{
			print STDERR "\nNo config file found! Either define one in the call of Utilities::getParams() or create $prog_path/current.config.xml!\n\n";
			exit(1);
		}
	}
	

	my $xml = new XML::Simple;

	my $params = $xml->XMLin($file);

	#TODO:
	# Introduce alias mechanism for
	#    exome kits
	#    settings that need to override just a small number of parameters

	return $params;
}

#######################################################################################
#connectExomeDB: Connect to the ExomeDB
#######################################################################################
sub connectExomeDB {
	my $setting = shift;
	
	my $params = Utilities::getParams();
	my $db     = $params->{settings}->{$setting}->{exomedb}->{database};
	my $host   = $params->{settings}->{$setting}->{exomedb}->{host};
	my $port   = $params->{settings}->{$setting}->{exomedb}->{port};
	my $user   = $params->{settings}->{$setting}->{exomedb}->{user};
	my $passwd = $params->{settings}->{$setting}->{exomedb}->{password};
	my $dbh;
	
	$dbh = DBI->connect("DBI:mysql:database=$db;host=$host;port=$port;mysql_local_infile=1", $user, $passwd,{PrintError => 0, RaiseError => 0});
	unless($dbh)
	 {
		$dbh = DBI->connect("DBI:mysql:database=$db;host=localhost;port=$port;mysql_local_infile=1", $user, $passwd,{PrintError => 0, RaiseError => 0}) or		#connection over IP address doesn't work locally --> try localhost connection
		 
		die print "$DBI::errstr\n"};
	
	return $dbh;
}

#######################################################################################
#connectSolexaDB: Connect to the SolexaDB
#######################################################################################
sub connectSolexaDB {
	my $params = Utilities::getParams();
	my $db     = $params->{solexadb}->{database};
	my $host   = $params->{solexadb}->{host};
	my $port   = $params->{solexadb}->{port};
	my $user   = $params->{solexadb}->{user};
	my $passwd = $params->{solexadb}->{password};
	my $dbh;
	
	$dbh = DBI->connect("DBI:mysql:database=$db;host=$host;port=$port", $user, $passwd,{PrintError => 0, RaiseError => 0});
	unless($dbh)
	 {
		$dbh = DBI->connect("DBI:mysql:database=$db;host=localhost;port=$port", $user, $passwd,{PrintError => 0, RaiseError => 0}) or		#connection over IP address doesn't work locally --> try localhost connection
		 
		die print "$DBI::errstr\n"};
	
	return $dbh;
}

#######################################################################################
#connectCoreDB: Connect to the coreDB (most of the thime exome db)
#######################################################################################
sub connectCoreDB {
	my $params = Utilities::getParams();
	my $db     = $params->{coredb}->{database};
	my $host   = $params->{coredb}->{host};
	my $port   = $params->{coredb}->{port};
	my $user   = $params->{coredb}->{user};
	my $passwd = $params->{coredb}->{password};
	my $dbh;
	
	$dbh = DBI->connect("DBI:mysql:database=$db;host=$host;port=$port", $user, $passwd,{PrintError => 0, RaiseError => 0});
	unless($dbh)
	 {
		$dbh = DBI->connect("DBI:mysql:database=$db;host=localhost;port=$port", $user, $passwd,{PrintError => 0, RaiseError => 0}) or		#connection over IP address doesn't work locally --> try localhost connection
		 
		die print "$DBI::errstr\n"};
	
	return $dbh;
}

#######################################################################################
#connectVariationDB: Connect to a variation DB
#######################################################################################
sub connectVariationDB {
	my $setting = shift;
	
	my $params = Utilities::getParams();
	my $db     = $params->{settings}->{$setting}->{variationdb}->{database};
	my $host   = $params->{settings}->{$setting}->{variationdb}->{host};
	my $port   = $params->{settings}->{$setting}->{variationdb}->{port};
	my $user   = $params->{settings}->{$setting}->{variationdb}->{user};
	my $passwd = $params->{settings}->{$setting}->{variationdb}->{password};
	my $dbh;
	
	$dbh = DBI->connect("DBI:mysql:database=$db;host=$host;port=$port", $user, $passwd,{PrintError => 0, RaiseError => 0	});
	unless($dbh) {
		$dbh = DBI->connect("DBI:mysql:database=$db;host=localhost;port=$port", $user, $passwd,{PrintError => 0, RaiseError => 0	}) 		#connection over IP address doesn't work locally --> try localhost connection
		|| 
		die print "$DBI::errstr\n"};
	
	return $dbh;
}

#######################################################################################
#connectRnaDB: Connect to the RnaDB
#######################################################################################
sub connectRnaDB {
	my $setting = shift;
	
	my $params = Utilities::getParams();
	my $db     = $params->{settings}->{$setting}->{rnadb}->{database};
	my $host   = $params->{settings}->{$setting}->{rnadb}->{host};
	my $port   = $params->{settings}->{$setting}->{rnadb}->{port};
	my $user   = $params->{settings}->{$setting}->{rnadb}->{user};
	my $passwd = $params->{settings}->{$setting}->{rnadb}->{password};
	my $dbh;
	
	$dbh = DBI->connect("DBI:mysql:database=$db;host=$host;port=$port", $user, $passwd,{PrintError => 0, RaiseError => 0});
	unless($dbh)
	 {
		$dbh = DBI->connect("DBI:mysql:database=$db;host=localhost;port=$port", $user, $passwd,{PrintError => 0, RaiseError => 0}) or		#connection over IP address doesn't work locally --> try localhost connection
		 
		die print "$DBI::errstr\n"};
	
	return $dbh;
}

###############################
##prepare output prefix
##############################
sub prepareOutputPrefix() {
	my $prefix = shift;
	my $dirName = shift;
	my $fileName = shift;
	
	my $op = $dirName;
	
	if ( $fileName =~ /(s_\d)/ ) {    #old format
		$op .= $prefix . "_" . $1;
	}
	elsif ( $fileName =~ /_L(\d+)_/ ) {    #new format
		my $lane = sprintf( "%d", $1 );
		$op .= $prefix . "_s_" . $lane;
		if ( $fileName =~ /_(\d+).fastq/ ) {		#add bin
			$op .= "_" . $1;
		}elsif( $fileName =~ /_(\d+).bam/ ) {
			$op .= "_" . $1;
		}
	}else{
		$op .= $prefix;			#other prefix is given
	}
	return $op;
}


###############################
##for trimming in GEM alignment
##parse the "x, y" parameter
###############################
sub parseTrimmingParameter() {
	my $trim = shift;
	my %tParam = ();
	
	$trim =~ /^([0-9]+),\s*([0-9]+)$/;
	$tParam{"x"} = $1;
	$tParam{"y"} = $2;
	
	return %tParam;
}
##################################################################
##for trimming in GEM alignment
##prepare the name for the trimmed sequence files (fastq, fq, ...)
##################################################################
sub prepareTrimmingFileName() {
	my $infile = shift;
	my $t = shift;
	
	my $ofn = "";
	my %trimParam = &parseTrimmingParameter($t);
	
	if ($infile =~ /\.trim$/) {
		$ofn = substr($infile, 0, -5);
		$ofn .= "_trimmed".$trimParam{"x"}."x".$trimParam{"y"};
		return $ofn;
	}
	else {
		#prepare outputfilename
		if ($infile =~ /\.gz$/) {
			if ($infile =~ /\.fastq.gz$/) {
				$ofn = substr($infile, 0, -9);
			} else {
				$ofn = substr($infile, 0, -6)
			}
		} else {
			if ($infile =~ /\.fastq$/) {
				$ofn = substr($infile, 0, -6);
			}
			else {
				$ofn = substr($infile, 0, -3);
			}
		}
		
		$ofn .= "_trimmed".$trimParam{"x"}."x".$trimParam{"y"}.".fastq";
		return $ofn;
	}
}


###############################################
##Convert an interval in seconds to an interval
##in dd:hh:mm:ss
###############################################
sub seconds_to_ddhhmmss {
	    my $seconds = shift;
	
	    my ($d, $h, $m, $s)  = 0;
	    $d = int($seconds/3600/24);
	    $h = int(($seconds/3600)%24);
	    $m = int(($seconds/60)%60);
	    $s = int($seconds%60);
	
	    return sprintf "%d:%02d:%02d:%02d", $d,$h,$m,$s;
}


###############################################
##Execute a command on the command line and
##check whether it returns successful or with
##errors
###############################################
sub executeCommand {
	my $command = shift;
	my $commandInfo = shift;
	my $l = shift;
	
	if ($l) {
		$l->info("###" . $commandInfo . "###");
		$l->info("CMD: $command");
	} else {
		print STDERR $commandInfo . "\n";
		print STDERR "CMD: $command\n";
	}
	my $start_time = time();
	my $ret = system($command);
	my $end_time = time();
	if ($ret) {
		if ($l) {
			$l->error("CMD $command died with return value: $ret");
		} else {
			print STDERR "CMD $command died with return value: $ret\n";
		}
	} else {
		if ($l) {
			$l->info("Finished successfully in " . &seconds_to_ddhhmmss($end_time - $start_time) . " (ddd:hh:mm:ss)");
		} else {
			print STDERR "Finished in " . &seconds_to_ddhhmmss($end_time - $start_time) . " (ddd:hh:mm:ss)\n";
		}
	}
	#if ($ret) { ERROR } else { SUCCESS }
	return $ret;
}


##################################################################################################################
#referenceHasChr : checks if reference has chr Prefix
##################################################################################################################
# Check if the reference genome uses "chr" in front of chromosome names unless provided from command line argument
# Useful to install the online DB where you don't need the reference file

sub referenceHasChr {
	my $setting = shift;
	my $params = Utilities::getParams();
	
	my $hasChr = 0;
	# 	Autotest if the reference needs chr prefix
	open FAI, $params->{settings}->{$setting}->{reference}.".fai" or die("Can't open $params->{settings}->{$setting}->{reference}.fai!");
	my $line = <FAI>;
	$hasChr  = 1 if $line =~ /^chr/;
	close FAI;
	
	return $hasChr;
}


##################################################################################################################
#twoBit: get reference sequence from fasta file. Now uses the BioPerl API instead of the twoBit program.
##################################################################################################################
sub twoBit {
	my $chrom   = shift;
	my $start   = shift;
	my $end     = shift;
	my $strand  = shift;
	my $setting = shift;
	my $cds     = "";
	
	my $params = Utilities::getParams();
	
	unless($fasta){
		$fasta = Bio::DB::Fasta->new($params->{settings}->{$setting}->{reference});
	}
	$cds = uc $fasta->seq("$chrom:$start,$end");

	#reverse complement
	if($strand eq "-")
	{
		
		
		$cds = reverse($cds);					#TODO: check if this actually works!!!!!!!!!!!!!!!
		my @nucs = split("",$cds);
		my $revComp = "";
		foreach my $nuc (@nucs)
		{
			$nuc = comp_rev($nuc);
			$revComp .= $nuc;
		}
		$cds = $revComp;
	}
	
	return $cds;
}

################################################
#iub: translate iub nomenclature to nucleotides#
################################################
sub iub {
my $iub = shift;
my %iub = ();
$iub{A}="A";
$iub{C}="C";
$iub{T}="T";
$iub{G}="G";
$iub{M}="AC";
$iub{R}="AG";
$iub{W}="AT";
$iub{S}="CG";
$iub{Y}="CT";
$iub{K}="GT";
$iub{V}="ACG";
$iub{H}="ACT";
$iub{D}="AGT";
$iub{B}="CGT";

$iub=$iub{$iub};

return $iub;
}
###################################################
#comp_rev: return reverse complement of nucleotide#
###################################################
sub comp_rev {
my $nuc = shift;
if ($nuc eq "A") {
	$nuc = "T";
}
elsif ($nuc eq "T") {
	$nuc = "A";
}
elsif ($nuc eq "G") {
	$nuc = "C";
}
elsif ($nuc eq "C") {
	$nuc = "G";
}

return ($nuc);
}


############################################
#init_code: translate codons to amino acids#
############################################
sub init_code {
my $code = shift;
$$code{TTT}="Phe";
$$code{TTC}="Phe";
$$code{TTA}="Leu";
$$code{TTG}="Leu";
$$code{TCT}="Ser";
$$code{TCC}="Ser";
$$code{TCA}="Ser";
$$code{TCG}="Ser";
$$code{TAT}="Tyr";
$$code{TAC}="Tyr";
$$code{TAA}="Stop";
$$code{TAG}="Stop";
$$code{TGT}="Cys";
$$code{TGC}="Cys";
$$code{TGA}="Stop";
$$code{TGG}="Trp";
$$code{CTT}="Leu";
$$code{CTC}="Leu";
$$code{CTA}="Leu";
$$code{CTG}="Leu";
$$code{CCT}="Pro";
$$code{CCC}="Pro";
$$code{CCA}="Pro";
$$code{CCG}="Pro";
$$code{CAT}="His";
$$code{CAC}="His";
$$code{CAA}="Gln";
$$code{CAG}="Gln";
$$code{CGT}="Arg";
$$code{CGC}="Arg";
$$code{CGA}="Arg";
$$code{CGG}="Arg";
$$code{ATT}="Ile";
$$code{ATC}="Ile";
$$code{ATA}="Ile";
$$code{ATG}="Met";
$$code{ACT}="Thr";
$$code{ACC}="Thr";
$$code{ACA}="Thr";
$$code{ACG}="Thr";
$$code{AAT}="Asn";
$$code{AAC}="Asn";
$$code{AAA}="Lys";
$$code{AAG}="Lys";
$$code{AGT}="Ser";
$$code{AGC}="Ser";
$$code{AGA}="Arg";
$$code{AGG}="Arg";
$$code{GTT}="Val";
$$code{GTC}="Val";
$$code{GTA}="Val";
$$code{GTG}="Val";
$$code{GCT}="Ala";
$$code{GCC}="Ala";
$$code{GCA}="Ala";
$$code{GCG}="Ala";
$$code{GAT}="Asp";
$$code{GAC}="Asp";
$$code{GAA}="Glu";
$$code{GAG}="Glu";
$$code{GGT}="Gly";
$$code{GGC}="Gly";
$$code{GGA}="Gly";
$$code{GGG}="Gly";
#return $code;
}

############################################
#init_code_mtdna: translate codons to amino acids#
############################################
sub init_code_mtdna {
my $code = shift;
$$code{TTT}="Phe";
$$code{TTC}="Phe";
$$code{TTA}="Leu";
$$code{TTG}="Leu";
$$code{TCT}="Ser";
$$code{TCC}="Ser";
$$code{TCA}="Ser";
$$code{TCG}="Ser";
$$code{TAT}="Tyr";
$$code{TAC}="Tyr";
$$code{TAA}="Stop";
$$code{TAG}="Stop";
$$code{TGT}="Cys";
$$code{TGC}="Cys";
$$code{TGA}="Trp";
$$code{TGG}="Trp";
$$code{CTT}="Leu";
$$code{CTC}="Leu";
$$code{CTA}="Leu";
$$code{CTG}="Leu";
$$code{CCT}="Pro";
$$code{CCC}="Pro";
$$code{CCA}="Pro";
$$code{CCG}="Pro";
$$code{CAT}="His";
$$code{CAC}="His";
$$code{CAA}="Gln";
$$code{CAG}="Gln";
$$code{CGT}="Arg";
$$code{CGC}="Arg";
$$code{CGA}="Arg";
$$code{CGG}="Arg";
$$code{ATT}="Ile";
$$code{ATC}="Ile";
$$code{ATA}="Met";
$$code{ATG}="Met";
$$code{ACT}="Thr";
$$code{ACC}="Thr";
$$code{ACA}="Thr";
$$code{ACG}="Thr";
$$code{AAT}="Asn";
$$code{AAC}="Asn";
$$code{AAA}="Lys";
$$code{AAG}="Lys";
$$code{AGT}="Ser";
$$code{AGC}="Ser";
$$code{AGA}="Stop";
$$code{AGG}="Stop";
$$code{GTT}="Val";
$$code{GTC}="Val";
$$code{GTA}="Val";
$$code{GTG}="Val";
$$code{GCT}="Ala";
$$code{GCC}="Ala";
$$code{GCA}="Ala";
$$code{GCG}="Ala";
$$code{GAT}="Asp";
$$code{GAC}="Asp";
$$code{GAA}="Glu";
$$code{GAG}="Glu";
$$code{GGT}="Gly";
$$code{GGC}="Gly";
$$code{GGA}="Gly";
$$code{GGG}="Gly";
#return $code;
}


########################################################
#indel2VCF: convert standard indel format to VCF
########################################################
sub indel2VCF {
	my $chrom   = shift;
	my $start   = shift;
	my $strand  = shift;
	my $indel   = shift;
	my $setting = shift;
	my $bam		= shift;
	
		
		
	my $params = Utilities::getParams();
	
	unless($fasta){
		$fasta = Bio::DB::Fasta->new($params->{settings}->{$setting}->{reference});
	}
	
	my $type = 0;		#insertion
	
	if($indel =~ /^-/){
		$type = 1;		#deletion
	}
	$indel =~ s/\+//;
	$indel =~ s/-//;
	$indel = uc $indel;
	
	my ($window, $windowStart, $windowEnd, $ref, $alt, $startPos);
	$window = length $indel;
	
	#print "orig start: $start\n";
	
	if($strand eq "-"){		#means that the indel is on - strand AND that it has to be reversed AND will only be shifted to the right, since 
							#the format says that it should be on the left most position
		$indel = reverse($indel);
		$indel =~ tr/AGCT/TCGA/; 
	}elsif($strand eq ""){			#if no strand is given assume that the sequence is in forward notation, but not necessarily on the right most position --> shift it to the right, if possible
									#!!!!!!!!!!!!! NOTE: this conversion is only tested with coordinates from our old snv DB table. 
		
		if($type == 1){
			$windowStart = $start;
			$windowEnd   = $windowStart + $window;
			#my $shiftCounter = 0;
			while(1){		
				my $firstNuc = uc $fasta->seq("$chrom:$windowStart,$windowStart");
				my $lastNuc  = uc $fasta->seq("$chrom:$windowEnd,$windowEnd");
				if($firstNuc ne $lastNuc){	#no more shifting possible
					last;
				}
				$indel = substr($indel,1,(length($indel)-1)).$lastNuc;
				#$shiftCounter++;
				$windowStart++;
				$windowEnd++;
				#if($shiftCounter % $window == 0){	# choose the last complete occurence as start position
				#	$start = $windowStart;
				#}
			}
			$start = $windowStart;
			
		}else{
			$windowStart = $start;
			$windowEnd   = $windowStart + ($window -1);
			while(1){
				my $nuc = uc $fasta->seq("$chrom:$windowStart,$windowEnd");
				if($nuc ne $indel){
					last;
				}

				$windowStart = ($windowEnd + 1);
				$windowEnd   = $windowStart + ($window -1);
				
			}
			$start = $windowStart -1;
		}	
	}
	#print "converted start: $start\n";
	
	
	
	
	
	if($type == 1){			#DELETION
		
			
		if($strand eq "-"){	#deletion is on backward strand --> shift it to the right as far as possible
			
			$startPos = $start - 1;						#since the position given is already the leftmost possible position, the VCF start position is one nucleotide in front of it
			$alt = uc $fasta->seq("$chrom:$startPos,$startPos");	#we can already retrieve this position from the reference genome since it will be the first letter of the reference and the alternative allele
			$ref = $alt.$indel;								#the given indel is defintely the first possible part that is deleted
			
			$windowStart = $start;
			$windowEnd	 = $windowStart + $window;
			
			
			while(1){		
				my $firstNuc = uc $fasta->seq("$chrom:$windowStart,$windowStart");
				my $lastNuc  = uc $fasta->seq("$chrom:$windowEnd,$windowEnd");
				if($firstNuc ne $lastNuc){	#no more shifting possible
					return ($startPos,$ref,$alt);
				}
				$ref .= $firstNuc;	#if the two nucleotides are the same a shift to the right is possible 
				$alt .= $firstNuc;
				$windowStart++;
				$windowEnd++;
			}
		}else{		#deletion is on the forward strand --> it is on the rightmost position, so shift it to the left
			$ref = $indel;	#the given indel must be part of the deleted string
			$alt = "";
			$windowStart = $start -1;
			$windowEnd   = $windowStart + $window;
			while(1){		
				my $firstNuc = uc $fasta->seq("$chrom:$windowStart,$windowStart");
				my $lastNuc  = uc $fasta->seq("$chrom:$windowEnd,$windowEnd");
				$ref = $firstNuc.$ref;	#add firstNuc to the left of both ref and alt
				$alt = $firstNuc.$alt;
				if($firstNuc ne $lastNuc){	#no more shifting possible --> current windowStart ist the position in front of the possible deletions
					return ($windowStart,$ref,$alt);
				}
				
				$windowStart--;
				$windowEnd--;
			}
		}
	}else{ #INSERTION
		if($strand eq "-"){	#insertion is on backward strand --> shift it to the right as far as possible
			$startPos = $start;								#for an insertion $start is already at the nucleotide in front of the first possible insertion 
			$ref = uc $fasta->seq("$chrom:$startPos,$startPos");	#we can already retrieve this position from the reference genome since it will be the first letter of the reference and the alternative allele
			$alt = $ref.$indel;								#the given indel is defintely the first possible part that is inserted
			
			$windowStart = $startPos +1;
			$windowEnd   = $windowStart + ($window -1);
			while(1){
				my $nuc = uc $fasta->seq("$chrom:$windowStart,$windowEnd");
				if($nuc ne $indel){
					return ($startPos,$ref,$alt);
				}
				$ref .= $nuc;	#if the window is the same as the insertion a shifting by the size of the window is possible
				$alt .= $nuc;
				
				$windowStart = ($windowEnd + 1);
				$windowEnd   = $windowStart + ($window -1);
				
			}
		}else{	#insertion is on forward strand
			$alt = $indel;
			$ref = "";
			$windowEnd   = $start;
			$windowStart = ($windowEnd - $window) + 1;
			while(1){
				my $nuc = uc $fasta->seq("$chrom:$windowStart,$windowEnd");
				if($nuc ne $indel){
					$nuc = substr($nuc,($window -1), 1);		#the last nucleotide of the current (non-matching) window is the first base in front of the insertion
					$ref = $nuc.$ref;
					$alt = $nuc.$alt;
					return ($windowEnd,$ref,$alt);
				}
				$ref = $nuc.$ref;
				$alt = $nuc.$alt;
				$windowEnd   = $windowStart -1;
				$windowStart = ($windowEnd - $window) + 1;
			}
			
		}
	}
}


########################################################
########################################################
#getMinimalIndel: converts samtools representation of indels (including repeats) to version without repeats
########################################################
sub getMinimalIndel {
	my $ref   = shift;
	my $alt   = shift;
	$ref = uc $ref;
	my @tmp = split(",",$alt);			#WARNING: returns only most likely allele when minimizing!
	$alt = $tmp[0];
	$alt = uc $alt;
	

	return ($ref,$alt) if (length($ref)==1 || length($alt)==1) || !($ref =~/[acgt]+/i && $alt =~/[acgt]+/i);	#if indel is already minimal -> return it

	if(length($ref)>length($alt)){	# deletion
		$ref = substr($ref,0,length($ref)-length($alt)+1);
		$alt = substr($alt,0,1);
	}else{				# insertion
		$alt = substr($alt,0,length($alt)-length($ref)+1);
		$ref = substr($ref,0,1);
	}
	return ($ref,$alt);
}

	
############################################################
#getVCF: get a VCF file from a list of idsnvs and idsamples
############################################################
sub getVCF {
	my $dbh      = shift;
	my $sampledb = shift;
	my $tmp      = shift;
	my @idsnvs   = @$tmp;
	$tmp         = shift;
	my @idsamples= @$tmp;
	
	
	my $ret = qq{##fileformat=VCFv4.1
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
##INFO=<ID=AF,Number=1,Type=Float,Description="Allele frequency of this variant in the Munich Exome DB">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="# high-quality bases">
##FORMAT=<ID=MQ,Number=1,Type=Integer,Description="Root-mean-square mapping quality of covering reads">
##FORMAT=<ID=SQ,Number=1,Type=Integer,Description="Variant quality (samtools)">
##FORMAT=<ID=PF,Number=1,Type=Integer,Description="Percent of reads on the forward strand">
##FORMAT=<ID=PV,Number=1,Type=Integer,Description="Percent of reads showing the variant">
##FILTER=<ID=VARQ,Description="Filtered out by vcfutils.pl varFilter, because of to low MQ (< 25)">
##FILTER=<ID=VARd,Description="Filtered out by vcfutils.pl varFilter, because of to low read depth (< 3)">
##FILTER=<ID=VARD,Description="Filtered out by vcfutils.pl varFilter, because of to high read depth (> 9999)">
##FILTER=<ID=VARa,Description="Filtered out by vcfutils.pl varFilter, because of to low read depth of variant bases (< 2)">
##FILTER=<ID=VARG,Description="Filtered out by vcfutils.pl varFilter, because of ??? (probably something with surrounding variants)">
##FILTER=<ID=VARg,Description="Filtered out by vcfutils.pl varFilter, because of ??? (probably something with surrounding variants)">
##FILTER=<ID=VARP,Description="Filtered out by vcfutils.pl varFilter, because of bias probability (see PV4)">
##FILTER=<ID=VARM,Description="Filtered out by vcfutils.pl varFilter, because of ???">
##FILTER=<ID=VARS,Description="Filtered out by vcfutils.pl varFilter, because of Hardy-Weinberg equilibrium">
##FILTER=<ID=Q20,Description="Filtered out by filterSNPqual.pl because median quality < 20">
##FILTER=<ID=Q15,Description="Filtered out by filterSNPqual.pl because median quality < 15">
##FILTER=<ID=Q10,Description="Filtered out by filterSNPqual.pl because median quality < 10">
##FILTER=<ID=Q5,Description="Filtered out by filterSNPqual.pl because median quality < 5">
##FILTER=<ID=Q3,Description="Filtered out by filterSNPqual.pl because median quality < 3">
##SamplesInMunichExomeDB=};


	#get number of samples in database
	my $query = "SELECT count(DISTINCT idsample) FROM snvsample;";
	my $out = $dbh->prepare($query) || die print "$DBI::errstr";
	$out->execute() || die print "$DBI::errstr";
	my ($sampleCount) =  $out->fetchrow_array;
	$ret .= "$sampleCount\n#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT";
	
	#print STDERR "getNumber: ".(time - $^T)."\n";
	
	#get names of samples
	foreach my $currSample(@idsamples) {
		$query = "SELECT name FROM $sampledb.sample WHERE idsample=?;";
		$out = $dbh->prepare($query) || die print "$DBI::errstr";
		$out->execute(($currSample)) || die print "$DBI::errstr";
		my ($currName) =  $out->fetchrow_array;
		$ret .= "\t$currName";
	}
	#print STDERR "getSamples: ".(time - $^T)."\n";
	
	#get entries for each SNV
	foreach my $currSNV (@idsnvs){
		
		#get entries
		$query = "SELECT chrom,start,if(rs='','.',rs),refallele,allele,freq FROM snv WHERE idsnv=?"; #
		$out = $dbh->prepare($query) || die print "$DBI::errstr";
		$out->execute(($currSNV)) || die print "$DBI::errstr";
		my @columns =  $out->fetchrow_array;
		$ret .= "\n".join("\t",@columns[0..4]);
		my $info = "AF=".sprintf("%.3f",($columns[5]/($sampleCount*2)));
		#print STDERR "getSNV: ".(time - $^T)."\n";
		
		
		my $format = "GT:GQ:DP:MQ:SQ:PF:PV";
		my $dp   = 0;
		my $qual = 0;
		my $foundSamples = 0;
		my %filter;
		#get info from snvsample
		foreach my $currSample(@idsamples) {
			$query = "SELECT alleles,percentvar,percentfor,snvqual,mapqual,coverage,filter,gtqual FROM snvsample WHERE idsnv = ? AND idsample = ?";
			$out = $dbh->prepare($query) || die print "$DBI::errstr";
			$out->execute(($currSNV,$currSample)) || die print "$DBI::errstr";
			if(my ($alleles,$percentvar,$percentfor,$snvqual,$mapqual,$coverage,$filter,$gtqual) =  $out->fetchrow_array){
				$foundSamples++;
				
				if($alleles == 0){
					$format .= "\t0/0";
				}elsif($alleles == 1){
					$format .= "\t0/1";
				}else{
					$format .= "\t1/1";
				}
				
				$format .= ":$gtqual:$coverage:$mapqual:$snvqual:$percentfor:$percentvar";
				
				my @columns = split(",",$filter);
				foreach(@columns){
					$filter{$_} = 1;		#set filters
				}
				
				$dp   += $coverage;
				$qual += $snvqual;
				
				
			}else{
				$format .= "\t.";
			}
			#print STDERR "getSnvsample: ".(time - $^T)."\n";
			
			
		}
		
		#mean of qual
		if($foundSamples>0){
			$qual = sprintf("%d",($qual/$foundSamples));
		}
		$ret .= "\t$qual";
		
		#remove PASS if more than one filter and PASS was set
		if( (keys %filter)>1 && $filter{"PASS"}){
			delete $filter{"PASS"};
		}
		$ret .= "\t".join(";",keys %filter);
		
	
		$info .= ";DP=$dp";
		$ret .= "\t$info\t$format";
	}

	
	return $ret;
}
		
		
############################################################
#getReferenceLen
############################################################
sub getReferenceLen		
{
	my $settings = shift;
	my $params 	= 	Utilities::getParams();
	my $ref =		$params->{settings}->{$settings}->{reference};
	my $fai =		$ref.".fai";
		
	
	open (FASTAINDEX, "$fai") or die "Not indexed reference";
	
	my $length = 0;
	
	while (<FASTAINDEX>) {
		chomp($_);
		my @data = split("\t", $_);
		if ( defined $data[1] )	{
			$length += $data[1];
		}
	}
	
	return $length;
}

############################################################
#validateChrPos				chr[XYM]:[1-9][0-9]*-[1-9][0-9]*
############################################################
sub validateChrPos
{
	my $pos = shift; 

	return ( $pos =~ /chr([1-9]|1[0-9]|2[1-2]|M|X|Y)(\:[1-9][0-9]*\-[1-9][0-9]*$)/ );
}

############################################################
#validateChrPosFull			chr[XYM] and optionally :[1-9][0-9]*-[1-9][0-9]*
############################################################
sub validateChrPosFull
{
	my $pos = shift; 

	return ( $pos =~ /chr([1-9]|1[0-9]|2[1-2]|M|X|Y)($|\:[1-9][0-9]*\-[1-9][0-9]*$)/ );
}

############################################################
#validateChrPosGeneric			whatever:[1-9][0-9]*-[1-9][0-9]*
############################################################
sub validateChrPosGeneric
{
	my $pos = shift; 
	return ( $pos =~ /[a-zA-z0-9\_\.\|]*(\:[1-9][0-9]*\-[1-9][0-9]*$)/ );
}

############################################################
#getMinimalRepresentation
# converts a vcf entry to its minimal representation:
# 1000 AAAT AAAA  -> 1003 T A
#
############################################################
sub getMinimalRepresentation
{
	my ($start, $ref, $alt) = (shift,shift,shift);

    if ( length($ref)==1 && length($alt)==1 )
    {
		# do nothing
	}
	else
	{
		# strip off identical suffixes
		while (substr($alt,length($alt)-1) eq substr($ref,length($ref)-1) &&
				min(length($alt),length($ref)) > 1 ) 
		{
			$alt = substr($alt,0,length($alt)-1);
			$ref = substr($ref,0,length($ref)-1);
		}
		
        # strip off identical prefixes and increment position
        while (substr($alt,0) eq substr($ref,0)  &&
				min(length($alt),length($ref)) > 1)
		{
			$alt = substr($alt, 1);
			$ref = substr($ref, 1);
			$start = $start + 1;
		}

	}
        return ($start, $ref, $alt);
}
sub min {
        my ($x,$y) = (shift,shift);
        my $min = ($x, $y)[$x > $y];
}


sub getMinimalRepresentation_RB
{
	my ($start, $ref, $alt) = (shift,shift,shift);

	if ( length($ref)==1 && length($alt)==1 )
	{
		#Well formatted
	}
	elsif ( length($ref) == length($alt) )
	{
		#Multiple substitution

		#CLIP Tail
		# In identical length identical tails are to be clipped!
		# ACCTA->CTCTA  becomes AC->CT
		while ( length($ref)>1 && ( substr($ref, -1, 1) eq substr($alt, -1, 1) ))
		{
			$ref=substr($ref, 0, -1);
			$alt=substr($alt, 0, -1);
		}

		#CLIP Head
		# AT->AG becomes T->G
		while ( length($ref)>1 && ( substr($ref, 0, 1) eq substr($alt, 0, 1) ))
		{
			$ref=substr($ref, 1);
			$alt=substr($alt, 1);
			$start++;
		}
	}
	elsif ( length($alt) > length($ref) )
	{
		#Insertion
		#CLIP Head
		# AAC->AACGT becomes C->CGT
		while ( ( length($ref)>1 && length($alt)>1 )&& ( substr($ref, 0, 1) eq substr($alt, 0, 1) ))
		{
			$ref=substr($ref, 1);
			$alt=substr($alt, 1);
			$start++;
		}
		
	}
	elsif ( length($alt) < length($ref) )
	{
		#Deletion
		#CLIP Head
		# AACTT->AA becomes ACTT->A
		while ( ( length($ref)>1 && length($alt)>1 )&& ( substr($ref, 0, 1) eq substr($alt, 0, 1) ))
		{
			$ref=substr($ref, 1);
			$alt=substr($alt, 1);
			$start++;
		}
		

	}
	else
	{
		#Impossible error		
	}

	return ($start, $ref, $alt);
}


############################################################
#getGVCFPath
# gets the path for a GVCF, given settings and sample name
# if bedARRAYJob settings are given
# 
# Params
# - sample name
# - settings
# - bedarray 0 if not, SGE_TASK_ID if yes
# - libtype if ambiguous
# - libpair if ambiguous
#
############################################################
sub getGVCFPath
{
	# Name
	my $sample_name = shift;
		$sample_name = "" if ! defined $sample_name;
	# Settings ( to get path etc )
	my $settings = shift;
		$settings = "" if ! defined $settings;
	# Is bedarray, if yes number if not, 0
	my $bedarray = shift;
		$bedarray = "" if ! defined $bedarray;
	# Libtype (in case of ambiguous libs (many types))
	my $libtypeconstraint = shift;
		$libtypeconstraint = "" if ! defined $libtypeconstraint;
	# Libpair status (in case of ambiguous libs (different paired statuses))
	my $libpairconstraint = shift;
		$libpairconstraint = "" if ! defined $libpairconstraint;
	
	# Logger
	my $logger = shift;	
	
	# output
	my $path = "";
	
	#Cleanup
	return "" if ( $settings eq ""  || $sample_name eq "" );
	my $gatk4 = ( $bedarray eq "4" ? 1 : 0 ); 
	$bedarray = 0 if ( $bedarray eq "" || $bedarray eq "4");
	
	
	# Db and params and table names
	my $dbh = Utilities::connectCoreDB();
	my $params = Utilities::getParams;
	
	#Databases
	my $coredb   = $params->{coredb}->{database};
	my $solexadb = $params->{solexadb}->{database};
	#Tables
	my $sampletable 	= $coredb.".".	$params->{coredb}->{sampletable};
	my $projecttable 	= $coredb.".".	$params->{coredb}->{projecttable};
	my $organismtable   = $coredb.".".	$params->{coredb}->{organismtable};
	my $sample2librarytable = $solexadb.".".$params->{solexadb}->{sample2librarytable};
	my $librarytable	= $solexadb.".".$params->{solexadb}->{librarytable};
	my $libtypetable	= $solexadb.".".$params->{solexadb}->{libtypetable};
	my $libpairtable	= $solexadb.".".$params->{solexadb}->{libpairtable};
	
	
	# Get Analysis folder
	my $analysis_folder="";
	if ( $params->{settings}->{$settings}->{analysis}->{folder} )
	{
		$analysis_folder = $params->{settings}->{$settings}->{analysis}->{folder};
	}
	else
	{
		$logger->error("Analysis folder not defined for settings $settings") if defined $logger;
		return "";
	}
	
	# Find sample
	
	# Constraint conditions for libtype (genomic, exomic, etc) and for libpair (paired-end, single-end)
		# necessary when a sample has multiple libs with different type or paired states.
		my $libtypecondition= ( $libtypeconstraint ne "") ? "and LT.ltlibtype=\"$libtypeconstraint\"" : "";
		my $libpaircondition= ( $libpairconstraint ne "") ? "and LP.lplibpair=\"$libpairconstraint\"" : "";
		
		# Find sample
		my $sql = qq{
				select 
				    S.name, P.pname, O.orname, group_concat(distinct(LT.ltlibtype)), group_concat(distinct(LP.lplibpair))
				from 
				    $sampletable S 
				    inner join $projecttable P on P.idproject = S.idproject 
				    inner join $organismtable O on O.idorganism=S.idorganism
				    inner join $sample2librarytable S2L on S2L.idsample=S.idsample 
				    inner join $librarytable L on S2L.lid=L.lid  
				    inner join $libtypetable LT on L.libtype=LT.ltid 
				    inner join $libpairtable LP on LP.lpid=L.libpair  
				where 
				    S.name=\"$sample_name\"
				    and S.sbam<>\"\" 
				    $libtypecondition
				    $libpaircondition
				    and S.nottoseq=0 
				    and L.lfailed=0 
				group by 
				    S.name 
		};
		
		my $sth = $dbh->prepare($sql);
		$sth->execute();
		
		#Get libtype; get libpair; get organism
		#  libtype and liborganism serve to check that the batch is made of compatible samples
		#  libtype and libpair can be multiple (group conact distinct, in case they are the script crashes and asks for specification)
		my ($sample, $projectname, $organism, $libtype, $libpair)=$sth->fetchrow_array();

				
		# Sample not found
		if ( !defined($sample) )
		{
			$logger->error("Sample $sample not found") if defined $logger;
			return "";
		}
		
	# GVCF name
	my $gvcf_subpath="";
	
	# Defaults for the path inside sample folders:
 		# BedArray
		my $gvcfSubFolder_BEDArray	="HaplotypeCaller";
		my $gvcfSuffix_BEDArray		=".haplotypecaller.gvcf.gz";
		# Normal
		my $gvcfSubFolder_Exome		=".";
		my $gvcfSuffix_Exome		="gatk.ontarget.haplotypecaller.gvcf.gz";
		   $gvcfSuffix_Exome		="gatk4.ontarget.haplotypecaller.gvcf.gz" if ($gatk4);
	
	if($bedarray && $params->{settings}->{$settings}->{genome_splits} )
	{
		# Build gvcf subpath
		$gvcf_subpath=$gvcfSubFolder_BEDArray."/".$bedarray.$gvcfSuffix_BEDArray;
	}
	elsif($bedarray && !($params->{settings}->{$settings}->{genome_splits}) )
	{
		$logger->error("Bed array job needs genome splits defined") if ( defined $logger );
		return "";
	}
	else
	{
		# Build gvcf subpath
		$gvcf_subpath=$gvcfSubFolder_Exome."/".$gvcfSuffix_Exome;
	}

	# Check if there are commas in libtype and libpair
		#But first apply constraint:
		$libtype = ( $libtypeconstraint ne "" ? $libtypeconstraint : $libtype );
		$libpair = ( $libpairconstraint ne "" ? $libpairconstraint : $libpair );
	
		my @libtypes = split(",", $libtype);
		my @libpairs = split(",", $libpair);
		
		if ((@libtypes != 1) || (@libpairs != 1))
		{
			$logger->error("Sample $sample has multiple libtypes ($libtype) or libpairs ($libpair), please specify libtype and libpair as arguments of the analysis") if ( defined $logger );
		} 
	
	# Build Gvcf path - Kept into account BEDARRAY JOB ($gvcf_subpath pre calculated)
	$path = $analysis_folder."/".$projectname."/".$sample_name."/".$libtype."out"."/".$libpair."out"."/".$gvcf_subpath;
	
	return $path;
} 


# Call GetGVCFPath set bedarray=4 to flag GATK4.
sub getGVCFPath4
{
	# Name
	my $sample_name = shift;
		$sample_name = "" if ! defined $sample_name;
	# Settings ( to get path etc )
	my $settings = shift;
		$settings = "" if ! defined $settings;
	#4
	# Libtype (in case of ambiguous libs (many types))
	my $libtypeconstraint = shift;
		$libtypeconstraint = "" if ! defined $libtypeconstraint;
	# Libpair status (in case of ambiguous libs (different paired statuses))
	my $libpairconstraint = shift;
		$libpairconstraint = "" if ! defined $libpairconstraint;
		
		
	getGVCFPath($sample_name, $settings, 4 ) if $libtypeconstraint eq "";
	getGVCFPath($sample_name, $settings, 4, $libtypeconstraint ) if $libpairconstraint eq "";
	getGVCFPath($sample_name, $settings, 4, $libtypeconstraint, $libpairconstraint ) if $libpairconstraint ne "";
	
}




1;
