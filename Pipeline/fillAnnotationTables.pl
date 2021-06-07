#!/usr/bin/perl

use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use DBI;
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);
use List::MoreUtils qw/ uniq /;
use Pod::Usage;
use Tabix;

my $prog_path = dirname( abs_path($0) );
require $prog_path . "/Utilities.pm";

# database
my $dbh          = "";
my $sql          = "";
my $sth          = "";
my $delete       = 0;

my $settings     = "";
my $logfile      = "SCREEN";
my $loglevel     = "INFO";
my $help         = 0;
my $man		     = 0;

my $chrprefix = 0;
my $nochrprefix = 0;

my $caddFile     = "";
my $caddRegion   = "";

my $clinVarFile  = "";

my $hifile	     = "";

my $dbNSFP       = "";
my $doPPH	     = 0;
my $doSIFT       = 0;

my $exacFile     = "";
my $gnomadFolder = "";
my $gnomadExomesFile= "";
my $gnomadConstraintsFile = "";

my $tgenomesFile = "";

my $kaviarFile   = "";

GetOptions(
	"d"    => \$delete,
	"se=s" => \$settings,
	"chrprefix" => \$chrprefix,
	"nochrprefix" => \$nochrprefix,
	"c=s"  => \$caddFile,
	"cr=s" => \$caddRegion,
	"cv=s" => \$clinVarFile,
	"hi=s" => \$hifile,
	"db=s" => \$dbNSFP,
	"p"    => \$doPPH,
	"s"    => \$doSIFT,
	"e=s"  => \$exacFile,
	"g=s"  => \$gnomadFolder,
	"gex=s"=> \$gnomadExomesFile,
	"gc=s" => \$gnomadConstraintsFile,
	"t=s"  => \$tgenomesFile,
	"kav=s"=> \$kaviarFile,
	"lf=s" => \$logfile,
	"ll=s" => \$loglevel,
	"h"    => \$help,
	"man"  => \$man
);

pod2usage( { -exitval => 0, -verbose => 1 } ) if $help;
pod2usage( { -exitval => 0, -verbose => 2 } ) if $man;
pod2usage( { -exitval => 1, -verbose => 1 } ) if ($caddFile eq "" && $clinVarFile eq "" && $hifile eq ""  && $dbNSFP eq ""  && $exacFile eq "" && $tgenomesFile eq "" && $kaviarFile eq "" && $gnomadFolder eq "" && $gnomadExomesFile eq "" && $gnomadConstraintsFile eq "") || (($doPPH || $doSIFT) && ! -d $dbNSFP ) || ($settings eq "" );
pod2usage( { -exitval => 1, -verbose => 1 } ) if ( $gnomadFolder ne "" && $gnomadExomesFile ne "");

my $params         = Utilities::getParams();

my $database				= $params->{coredb}->{database};
my $host					= $params->{coredb}->{host};
my $user					= $params->{coredb}->{user};
my $password				= $params->{coredb}->{password};

# VariantDB is the REFERENCE DB
my $variantdb			    = $params->{settings}->{$settings}->{variationdb}->{database};

# All subtables into Reference DB:
my $caddtable               = $variantdb . "." . $params->{settings}->{$settings}->{variationdb}->{caddtable};
my $pphtable                = $variantdb . "." . $params->{settings}->{$settings}->{variationdb}->{pphtable};
my $sifttable               = $variantdb . "." . $params->{settings}->{$settings}->{variationdb}->{sifttable};
my $clinvartable            = $variantdb . "." . $params->{settings}->{$settings}->{variationdb}->{clinvartable};
my $haploinsufficiencytable = $variantdb . "." . $params->{settings}->{$settings}->{variationdb}->{haploinsufficiencytable};
my $exactable				= $variantdb . "." . $params->{settings}->{$settings}->{variationdb}->{exactable};
my $gnomadtable				= $variantdb . "." . $params->{settings}->{$settings}->{variationdb}->{exactable};	
my $gnomadconstraintstable  = $variantdb . "." . $params->{settings}->{$settings}->{variationdb}->{gnomadconstraints};

#TODO gnomadtable=exactable; #Should be renamed
#   $gnomadtable				= $variantdb . "." . "gnomad";	
#   $exactable = $gnomadtable; 
#### TODO: COMMENT THE UP ONE BACK
my $kaviartable				= $variantdb . "." . $params->{settings}->{$settings}->{variationdb}->{kaviartable};
my $tgenomestable		    = $variantdb . "." . $params->{settings}->{$settings}->{variationdb}->{tgenomestable};


# TODO: Transform it into an option
#### Temporarly change gnomadTable for import (when updating)
#$gnomadtable = "hg19.evs_tmp";
#$gnomadtable = "hg19.evs2";
####

my $exomeprefix="";

Utilities::initLogger( $logfile, $loglevel );
my $logger = Utilities::getLogger();

$dbh = Utilities::connectCoreDB();


# Check if the reference genome uses "chr" in front of chromosome names unless provided from command line argument
# Useful to install the online DB where you don't need the reference file
my $useChr = 0;
my $chr = "";

if ( $chrprefix || $nochrprefix ) 
{
	if ($chrprefix)
	{
		$useChr = 1;
		$chr="chr";
	}
}
else
{
	# Autotest if the reference needs chr prefix
	open FAI, $params->{settings}->{$settings}->{reference}.".fai" or exit $logger->error("Can't open $params->{settings}->{$settings}->{reference}.fai!");
	my $line = <FAI>;
	$useChr  = 1 if $line =~ /^chr/;
	close FAI;

	$chr    = "chr" if $useChr;
}

#delete required entries
if($delete){
	
	executeQuery("delete from $caddtable") if $caddFile ne "";
	executeQuery("delete from $clinvartable") if $clinVarFile ne "";
	executeQuery("delete from $haploinsufficiencytable") if $hifile ne "";
	executeQuery("delete from $pphtable") if $doPPH;
	executeQuery("delete from $sifttable") if $doSIFT;
	executeQuery("delete from $exactable") if $exacFile ne "";
	executeQuery("delete from $gnomadtable") if $gnomadFolder ne "";
	executeQuery("delete from $kaviartable") if $kaviarFile ne "";
}


# If GNOMAD Exomes, enable the ExAC routine to import the gnomad.exomes file
if ( $gnomadFolder ne "" || $gnomadExomesFile ne "")
{
	if($gnomadFolder)
	{
		foreach(glob("$gnomadFolder/gnomad.exomes.*vcf*gz"))
		{
			#It's only one anyway
			$exacFile=$_;
			$gnomadExomesFile=$exacFile;
		}
	}
	else
	{
		$exacFile=$gnomadExomesFile;			
	}
	
	$exactable=$gnomadtable;
	#No extra columns
	#$exomeprefix="exome_";
	
	
}


# CADD
if($caddFile ne ""){ #get data from CADD
	
	open MYSQL, "| mysql -h$host -u$user -p$password $database --local-infile=1 -e 'LOAD DATA LOCAL INFILE \"/dev/stdin\" INTO TABLE $caddtable;'" or exit $logger->error("Can't open pipe to MySQL: | mysql -h$host -u$user -p$password $database --local-infile=1 -e 'LOAD DATA LOCAL INFILE \"/dev/stdin\" INTO TABLE $caddtable;'");;
	
	if($caddRegion ne ""){
		my $tabix = Tabix->new( '-data' => $caddFile ) or exit $logger->error("Can't open $caddFile!") ;
		open CR, $caddRegion or exit $logger->error("Can't open $caddRegion!");
		while(<CR>){
			next if $_ =~ /^#/; #skip header 
			chomp;
			my @region = split("\t");
			my $iter = $tabix->query($region[0],$region[1],$region[2]);
			$logger->info("CADD - Processing region $region[0]:$region[1]-$region[2]...");
			while ( ( my $line = $tabix->read($iter) ) ) {
				print MYSQL $chr.$line;
			}
			
		}
		close CR;
	}else {
		if($caddFile =~ /^http/ || $caddFile =~ /^ftp/){		# --> if URL use curl
			$caddFile = "curl -s ".$caddFile." | gunzip -c |"
		}else{
			$caddFile = "gunzip -c ".$caddFile." |" if $caddFile =~ /\.gz$/;
		}
		open CADD, $caddFile or exit $logger->error("Can't open $caddFile!");
			
		my $counter = 0;
		while ( <CADD>) {
			next if $_ =~ /^#/; #skip header
			print MYSQL $chr.$_;
			$counter++;
			$logger->info("CADD - Inserted $counter entries - last entry: $chr$_") if $counter % 1000000 == 0;	
		}
		
		close CADD;
				
	}	
	close MYSQL;
} 


# POLYPHEN2 / SIFT fron dbNSFP3 
if($doPPH || $doSIFT) { #insert Polyphen2 and/or SIFT scores
	foreach(glob("$dbNSFP/dbNSFP3.*_variant.chr*")){
		open DBNSFP, $_ or exit $logger->error("Can't open $_!");
		$logger->info("dbNSFP - Processing file $_...");
		while(my $line = <DBNSFP>){
			next if $line  =~ /^#/; #skip header 
			chomp $line;
			my @columns = split("\t",$line);
			if($doPPH){
				if($columns[34] eq "B"){
					$columns[34] = "benign";
				}elsif($columns[34] eq "D"){
					$columns[34] = "probably damaging";
				}elsif($columns[34] eq "P"){
					$columns[34] = "possibly damaging";
				}else{
					$columns[34] = "unknown";
				}
				my ($score,@tmp) = split(";",$columns[32]);
				executeQuery("insert into $pphtable (chrom,start,ref,alt,hvar_prediction,hvar_prob) values ('$chr$columns[7]',$columns[8],'$columns[2]','$columns[3]','$columns[34]',$score)") if $score ne ".";
			}
			if($doSIFT){
				my ($score,@tmp) = split(";",$columns[23]);
				executeQuery("insert into $sifttable (chrom,start,ref,alt,score) values ('$chr$columns[7]',$columns[8],'$columns[2]','$columns[3]',$score)") if $score ne ".";
			}
			
		}
		
		close DBNSFP;
	}
}


# ClinVar
if($clinVarFile ne ""){	#insert ClinVar
	my %clnsigtrans = (
		0 , "Uncertain significance", 
		1 , "not provided", 
		2 , "Benign", 
		3 , "Likely benign", 
		4 , "Likely pathogenic", 
		5 , "Pathogenic",
		6 , "drug response", 
		7 , "histocompatibility", 
		255 , "other"
	);
	
	
	$clinVarFile = "gunzip -c ".$clinVarFile." |" if $clinVarFile =~ /\.gz$/;
	open CV, $clinVarFile or exit $logger->error("Can't open $clinVarFile!");
	
	while(<CV>){
		next if $_ =~ /^#/; #skip header 
		chomp;
		my @columns = split("\t");
		
		my @info = split(";",$columns[7]);
		my @clnsig;
		my @clnacc;
		foreach my $entry (@info){
			my ($key,$value) = split("=",$entry);
			if($key eq "CLNSIG"){
				@clnsig = split(",",$value);
			}elsif($key eq "CLNACC"){
				@clnacc = split(",",$value);
			}
		}
	
		
		my @altalleles = split(",",$columns[4]);
		
		for(my $i = 0; $i<@altalleles;$i++){
			my @subclnacc;
			my @subclnsig;
			if(defined $clnacc[$i]){
				@subclnacc = split(/\|/,$clnacc[$i]);
			}else{
				@subclnacc = split(/\|/,$clnacc[0]);
			}
			if(defined $clnacc[$i]){
				@subclnsig = split(/\|/,$clnsig[$i]);
			}else{
				@subclnsig = split(/\|/,$clnsig[0]);
			}
			
			for(my $j = 0; $j < @subclnacc; $j++){
				executeQuery("insert into $clinvartable (chrom,start,ref,alt,rcv,path) values ('$chr$columns[0]',$columns[1],'$columns[3]','$altalleles[$i]','$subclnacc[$j]','".$clnsigtrans{$subclnsig[$j]}."')");
			}
		}
	}
	
	close CV;
}

#Haploinsufficiency
if($hifile ne ""){ #insert haploinsufficiency
	open HI, $hifile or exit $logger->error("Can't open $hifile!");
	
	while(<HI>){
		next if $_ =~ /^track/; #skip header 
		chomp;
		my @columns = split("\t");
		$columns[0] =~ s/^chr// unless $useChr;
		my ($gene,$hiscore,$hipercent) = split(/\|/,$columns[3]);
		$hipercent  =~ s/\%$//;
		
		executeQuery("insert into $haploinsufficiencytable (chrom,start,end,genesymbol,hiscore,hipercent) values ('$columns[0]',$columns[1],$columns[2],'$gene',$hiscore,$hipercent)");		
		
		
	}
	
	close HI;
}

#Haploinsufficiency gene related from GnomAD
# GnomAD 
if ($gnomadConstraintsFile ne "" )
{
	if ( ! -e $gnomadConstraintsFile )
	{
		$logger->error("Gnomad Constraints table does not exist!");
		exit (-1);
	}

	# Straight import the file (tail -n +2 skips header)
	my $command="zcat $gnomadConstraintsFile | tail -n +2 | mysql -h$host -u$user -p$password $database --local-infile=1 -e 'LOAD DATA LOCAL INFILE \"/dev/stdin\" INTO TABLE $gnomadconstraintstable;'";
	
	if (&Utilities::executeCommand($command, "Importing GnomAD constraints table", $logger)) {
			exit;
	} 

}


#GnomAD
#	foreach(glob("$dbNSFP/dbNSFP3.*_variant.chr*")){
#		open DBNSFP, $_ or exit $logger->error("Can't open $_!");
#		$logger->info("dbNSFP - Processing file $_...");
#		while(my $line = <DBNSFP>){
if ($gnomadFolder ne "")
{
	# LOAD DATA INFILE HANDLER
	open MYSQL, "| mysql -h$host -u$user -p$password $database --local-infile=1 -e 'LOAD DATA LOCAL INFILE \"/dev/stdin\" INTO TABLE $gnomadtable;'" or exit $logger->error("Can't open pipe to MySQL: | mysql -h$host -u$user -p******* $database --local-infile=1 -e 'LOAD DATA LOCAL INFILE \"/dev/stdin\" INTO TABLE $gnomadtable;'");
	#open MYSQL, " | cat";
	
	foreach(glob("$gnomadFolder/gnomad.genomes.*vcf*")) 
	{	
		my $gnomadFile=$_;
		
		if($gnomadFile ne ""){

			# If gzipped gunzip on the fly, otherwise keep its name
			$gnomadFile = "gunzip -c ".$gnomadFile." |" if ( $gnomadFile =~ /\.gz$/ || $gnomadFile =~ /\.bgz$/ ) ;
			
			$logger->debug("Opening now file $gnomadFile - to load into GnomAD-table $gnomadtable");
			open GN, $gnomadFile or exit $logger->error("Can't open $gnomadFile!");
			
			while(<GN>){
				
				next if $_ =~ /^#/; #skip header 
				chomp;
				my @columns = split("\t");
				
				my @alt = split(",",$columns[4]);		
				my $n = @alt;		#get number of alternative alleles

                                my $filter = $columns[6];
                                $filter = "VQSR" unless $filter eq "PASS";	
			
				my @info = split(";",$columns[7]);
	
				# All datasets
				my @ac;
				my @hom;
				my @het;
				my $an;		
				
				# Africans
				my @ac_afr;
				my @hom_afr;
				my @het_afr;
				my $an_afr;
				
				# Non-finnish europeans
				my @ac_nfe;
				my @hom_nfe;
				my @het_nfe;
				my $an_nfe;
				
				#my @ac_popmax;
				#my @hom_popmax;
				#my @het_popmax;
				#my $an_popmax;
				my @af_popmax;
				
				foreach my $entry (@info){
					my ($key,$value) = split("=",$entry);
					if($key eq "AC"){
                                                @ac = split(",",$value);
                                        }elsif($key eq "nhomalt"){
                                                @hom = split(",",$value);
                                        }elsif($key eq "AN"){
                                                $an = $value;
					}elsif($key eq "AC_afr"){
						@ac_afr = split(",",$value);
					}elsif($key eq "nhomalt_afr"){
						@hom_afr = split(",",$value);
					}elsif($key eq "AN_afr"){
						$an_afr = $value;
					}elsif($key eq "AC_nfe"){
						@ac_nfe = split(",",$value);
					}elsif($key eq "nhomalt_nfe"){
						@hom_nfe = split(",",$value);
					}elsif($key eq "AN_nfe"){
						$an_nfe  = $value;
					}elsif($key eq "AF_popmax"){
						@af_popmax = split(",",$value);
					}
					
					#}elsif($key eq "AC_popmax"){
					#	@ac_popmax = split(",",$value);
					#}elsif($key eq "nhomalt_popmax"){
					#	@hom_popmax = split(",",$value);
					#}elsif($key eq "AN_popmax"){
					#	$an_popmax  = $value;
					
					# New GNOMAD format to mark segdup and lcr in INFO field instead of filter
					# it's not a key,value pair - as it's parsed here
					#TODO: Why don't you use a VCF parser instead?
					if ( $entry eq "segdup" || $entry eq "lcr" )
					{
						$filter = "VQSR";
					} 
				}
				
				
				my $sum = 0;
				my $afr_sum = 0;
				my $nfe_sum = 0;
				#my $popmax_sum = 0;
				
				for(my $i = 0; $i < $n; $i++){

                                        $het[$i]     = ($ac[$i]    - ($hom[$i]*2));                 #create numbers per allele
                                        $sum        += $hom[$i]     + $het[$i];

					$het_afr[$i] = ($ac_afr[$i]- ($hom_afr[$i]*2));
					$afr_sum    += $hom_afr[$i] + $het_afr[$i];
					
					$het_nfe[$i] = ($ac_nfe[$i]- ($hom_nfe[$i]*2));
					$nfe_sum    += $hom_nfe[$i] + $het_nfe[$i];
					
					#$het_popmax[$i] = ($ac_popmax[$i]- ($hom_popmax[$i]*2)) 	if( defined ($ac_popmax[$i]) && defined ($hom_popmax[$i]) );
					#$popmax_sum    += $hom_popmax[$i] + $het_popmax[$i]			if( defined ($hom_popmax[$i]) && defined ($het_popmax[$i]));
				}
				
				my $homref     = (($an     - ($sum*2) )/2);
				my $afr_homref = (($an_afr - ($afr_sum*2) )/2);
				my $nfe_homref = (($an_nfe - ($nfe_sum*2) )/2);
				#my $popmax_homref = 0;
				#	$popmax_homref = (($an_popmax - ($popmax_sum*2) )/2)			if( defined($an_popmax));
				
				for(my $i = 0; $i < $n; $i++){
					
					my $ref = $columns[3];
					my $alt = $alt[$i];
					if(length($ref) > 1 &&  length($ref) == length($alt)){		#special case of biallelic site with substitution and deletion at same site e.g. CGTT    C,TGTT 
						$ref = substr($ref,0,1);
						$alt = substr($alt,0,1);
					}
					
					#RB 20160719: Minimal representation for alleles
					my $entry_start=$columns[1];
					my ($newentry_start, $newref, $newalt)=Utilities::getMinimalRepresentation($entry_start, $ref,$alt);
					
					if ( $newentry_start != $entry_start || $newref ne $ref || $newalt ne $alt )
					{
						print "$entry_start -> $newentry_start | $ref -> $alt || $newref -> $newalt\n";
					}
					
					#$het_popmax[$i]=0 if (! defined $het_popmax[$i]);
					#$hom_popmax[$i]=0 if (! defined $hom_popmax[$i]);
					
					#print $af_popmax[$i] if (  defined $af_popmax[$i]);
					$af_popmax[$i]=0 if ( ! defined $af_popmax[$i]);
					
					#executeQuery("insert ignore into $exactable (chrom,start,refallele,allele,filter,ea_homref,ea_het,ea_homalt,aa_homref,aa_het,aa_homalt) values ('$chr".$columns[0]."',".$newentry_start.",'$newref','$newalt','".$filter."',$nfe_homref,".$het_nfe[$i].",".$hom_nfe[$i].",".$afr_homref.",".$het_afr[$i].",".$hom_afr[$i].");" );
					# Wants an ID to be set as autoincrement
					# TODO: ADD new zeroes for exomes of GNOMAD
					
					#print MYSQL "0"."\t".$chr.$columns[0]."\t".$newentry_start."\t".$newref."\t".$newalt."\t".$filter."\t".$nfe_homref."\t".$het_nfe[$i]."\t".$hom_nfe[$i]."\t".$afr_homref."\t".$het_afr[$i]."\t".$hom_afr[$i]."\t".$popmax_homref."\t".$het_popmax[$i]."\t".$hom_popmax[$i]."\n";
					print MYSQL "0"."\t".$chr.$columns[0]."\t".$newentry_start."\t".$newref."\t".$newalt."\t".$filter."\t".$nfe_homref."\t".$het_nfe[$i]."\t".$hom_nfe[$i]."\t".$afr_homref."\t".$het_afr[$i]."\t".$hom_afr[$i]."\t".$af_popmax[$i]."\t".$homref."\t".$het[$i]."\t".$hom[$i]."\n";
				}
			}
			close GN;
		}
	}
}


#ExAC - Moved after GnomAD - so I can update the table with exomic values
if($exacFile ne ""){
	if($exacFile =~ /^http/ || $exacFile =~ /^ftp/){		# --> if URL use curl
		$exacFile = "curl -s ".$exacFile." | gunzip -c |"
	}else{
		$exacFile = "gunzip -c ".$exacFile." |" if (  $exacFile =~ /\.gz$/ || $exacFile =~ /\.bgz$/  );
	}
	open EX, $exacFile or exit $logger->error("Can't open $exacFile!");
	
	while(<EX>){
		
		next if $_ =~ /^#/; #skip header 
		chomp;
		my @columns = split("\t");
		
		my $filter = $columns[6];
		$filter = "VQSR" unless $filter eq "PASS";

		my @alt = split(",",$columns[4]);		
		my $n = @alt;		#get number of alternative alleles
		
		my @info = split(";",$columns[7]);
	

				my @ac;
                                my @hom;
                                my @het;
                                my $an;	
		
				my @ac_afr;
				my @hom_afr;
				my @het_afr;
				my $an_afr;
				
				my @ac_nfe;
				my @hom_nfe;
				my @het_nfe;
				my $an_nfe;
				
				#my @ac_popmax;
				#my @hom_popmax;
				#my @het_popmax;
				#my $an_popmax;			
				my @af_popmax;
				
				foreach my $entry (@info){
					my ($key,$value) = split("=",$entry);
					if($key eq "AC"){
                                                @ac = split(",",$value);
                                        }elsif($key eq "nhomalt"){
                                                @hom = split(",",$value);
                                        }elsif($key eq "AN"){
                                                $an = $value;
					}elsif($key eq "AC_afr"){
						@ac_afr = split(",",$value);
					}elsif($key eq "nhomalt_afr"){
						@hom_afr = split(",",$value);
					}elsif($key eq "AN_afr"){
						$an_afr = $value;
					}elsif($key eq "AC_nfe"){
						@ac_nfe = split(",",$value);
					}elsif($key eq "nhomalt_nfe"){
						@hom_nfe = split(",",$value);
					}elsif($key eq "AN_nfe"){
						$an_nfe  = $value;
					}elsif($key eq "AF_popmax"){
						@af_popmax = split(",",$value);
					}
					
					#}elsif($key eq "AC_popmax"){
					#	@ac_popmax = split(",",$value);
					#}elsif($key eq "nhomalt_popmax"){
					#	@hom_popmax = split(",",$value);
					#}elsif($key eq "AN_popmax"){
					#	$an_popmax  = $value;
				}
		
		
		my $sum = 0;
		my $afr_sum = 0;
		my $nfe_sum = 0;
		#my $popmax_sum = 0;
		
		for(my $i = 0; $i < $n; $i++){
                        $het[$i]     = ($ac[$i]    - ($hom[$i]*2));                 #create numbers per allele
                        $sum        += $hom[$i]     + $het[$i];

			$het_afr[$i] = ($ac_afr[$i]- ($hom_afr[$i]*2));
			$afr_sum    += $hom_afr[$i] + $het_afr[$i];
			
			$het_nfe[$i] = ($ac_nfe[$i]- ($hom_nfe[$i]*2));
			$nfe_sum    += $hom_nfe[$i] + $het_nfe[$i];
						
			#$het_popmax[$i] = ($ac_popmax[$i]- ($hom_popmax[$i]*2)) 	if( defined ($ac_popmax[$i]) && defined ($hom_popmax[$i]) );
			#$popmax_sum    += $hom_popmax[$i] + $het_popmax[$i]			if( defined ($hom_popmax[$i]) && defined ($het_popmax[$i]));
		}
		my $homref     = (($an     - ($sum*2) )/2);
		my $afr_homref = (($an_afr - ($afr_sum*2) )/2);
		my $nfe_homref = (($an_nfe - ($nfe_sum*2) )/2);
		#my $popmax_homref = 0; #(($an_popmax - ($popmax_sum*2) )/2);
		#$popmax_homref = (($an_popmax - ($popmax_sum*2) )/2)			if( defined($an_popmax));
		
		for(my $i = 0; $i < $n; $i++){
			
			my $ref = $columns[3];
			my $alt = $alt[$i];
			if(length($ref) > 1 &&  length($ref) == length($alt)){		#special case of biallelic site with substitution and deletion at same site e.g. CGTT    C,TGTT 
				$ref = substr($ref,0,1);
				$alt = substr($alt,0,1);
			}
			
			#RB 20160719: Minimal representation for alleles
			my $entry_start=$columns[1];
			my ($newentry_start, $newref, $newalt)=Utilities::getMinimalRepresentation($entry_start, $ref,$alt);
			
			if ( $newentry_start != $entry_start || $newref ne $ref || $newalt ne $alt )
			{
				print "$entry_start -> $newentry_start | $ref -> $alt || $newref -> $newalt\n";
			}
			
			#$het_popmax[$i]=0 if (! defined $het_popmax[$i]);
			#$hom_popmax[$i]=0 if (! defined $hom_popmax[$i]);
			$af_popmax[$i]=0 if ( ! defined $af_popmax[$i]);
			
			if ( $gnomadExomesFile eq "" )
			{
				#It's exac
				executeQuery("insert ignore into $exactable (chrom,start,refallele,allele,".$exomeprefix."filter,".$exomeprefix."ea_homref,".$exomeprefix."ea_het,".$exomeprefix."ea_homalt,".$exomeprefix."aa_homref,".$exomeprefix."aa_het,".$exomeprefix."aa_homalt,".$exomeprefix."popmax_af,".$exomeprefix."homref,".$exomeprefix."het,".$exomeprefix."homalt) values ('$chr".$columns[0]."',".$newentry_start.",'$newref','$newalt','".$filter."',$nfe_homref,".$het_nfe[$i].",".$hom_nfe[$i].",".$afr_homref.",".$het_afr[$i].",".$hom_afr[$i].",".$af_popmax[$i].",".$homref.",".$het[$i].",".$hom[$i].") ON DUPLICATE KEY UPDATE ".$exomeprefix."filter='".$filter."', ".$exomeprefix."ea_homref=". $nfe_homref .", ".$exomeprefix."ea_het=".$het_nfe[$i].", ".$exomeprefix."ea_homalt=".$hom_nfe[$i].", ".$exomeprefix."aa_homref=".$afr_homref.", ".$exomeprefix."aa_het=".$het_afr[$i].", ".$exomeprefix."aa_homalt=".$hom_afr[$i].", ".$exomeprefix."popmax_af=".$af_popmax[$i].", ".$exomeprefix."homref=".$homref.", ".$exomeprefix."het=".$het[$i].", ".$exomeprefix."homalt=".$hom[$i].";" );
				#executeQuery("insert ignore into $exactable (chrom,start,refallele,allele,".$exomeprefix."filter,".$exomeprefix."ea_homref,".$exomeprefix."ea_het,".$exomeprefix."ea_homalt,".$exomeprefix."aa_homref,".$exomeprefix."aa_het,".$exomeprefix."aa_homalt,".$exomeprefix."popmax_homref,".$exomeprefix."popmax_het,".$exomeprefix."popmax_homalt ) values ('$chr".$columns[0]."',".$newentry_start.",'$newref','$newalt','".$filter."',$nfe_homref,".$het_nfe[$i].",".$hom_nfe[$i].",".$afr_homref.",".$het_afr[$i].",".$hom_afr[$i].",".$popmax_homref.",".$het_popmax[$i].",".$hom_popmax[$i].") ON DUPLICATE KEY UPDATE ".$exomeprefix."filter='".$filter."', ".$exomeprefix."ea_homref=". $nfe_homref .", ".$exomeprefix."ea_het=".$het_nfe[$i].", ".$exomeprefix."ea_homalt=".$hom_nfe[$i].", ".$exomeprefix."aa_homref=".$afr_homref.", ".$exomeprefix."aa_het=".$het_afr[$i].", ".$exomeprefix."aa_homalt=".$hom_afr[$i].", ".$exomeprefix."popmax_homref=".$popmax_homref.", ".$exomeprefix."popmax_het=".$het_popmax[$i].", ".$exomeprefix."popmax_homalt=".$hom_popmax[$i].";" );
			}
			else
			{
				#Sum on already present values
				executeQuery("insert ignore into $exactable (chrom,start,refallele,allele,".$exomeprefix."filter,".$exomeprefix."ea_homref,".$exomeprefix."ea_het,".$exomeprefix."ea_homalt,".$exomeprefix."aa_homref,".$exomeprefix."aa_het,".$exomeprefix."aa_homalt,".$exomeprefix."popmax_af,".$exomeprefix."homref,".$exomeprefix."het,".$exomeprefix."homalt) values ('$chr".$columns[0]."',".$newentry_start.",'$newref','$newalt','".$filter."',$nfe_homref,".$het_nfe[$i].",".$hom_nfe[$i].",".$afr_homref.",".$het_afr[$i].",".$hom_afr[$i].",".$af_popmax[$i].",".$homref.",".$het[$i].",".$hom[$i].") ON DUPLICATE KEY UPDATE ".$exomeprefix."filter='".$filter."', ".$exomeprefix."ea_homref=".$exomeprefix."ea_homref + ".$nfe_homref .", ".$exomeprefix."ea_het=".$exomeprefix."ea_het + ".$het_nfe[$i].", ".$exomeprefix."ea_homalt=".$exomeprefix."ea_homalt + ".$hom_nfe[$i].", ".$exomeprefix."aa_homref=".$exomeprefix."aa_homref + ".$afr_homref.", ".$exomeprefix."aa_het=".$exomeprefix."aa_het + ".$het_afr[$i].", ".$exomeprefix."aa_homalt=".$exomeprefix."aa_homalt + ".$hom_afr[$i].", ".$exomeprefix."popmax_af=GREATEST(".$exomeprefix."popmax_af, ".$af_popmax[$i]."), ".$exomeprefix."homref=".$exomeprefix."homref + ".$homref.", ".$exomeprefix."het=".$exomeprefix."het + ".$het[$i].", ".$exomeprefix."homalt=".$exomeprefix."homalt + ".$hom[$i].";" );
				#executeQuery("insert ignore into $exactable (chrom,start,refallele,allele,".$exomeprefix."filter,".$exomeprefix."ea_homref,".$exomeprefix."ea_het,".$exomeprefix."ea_homalt,".$exomeprefix."aa_homref,".$exomeprefix."aa_het,".$exomeprefix."aa_homalt,".$exomeprefix."popmax_homref,".$exomeprefix."popmax_het,".$exomeprefix."popmax_homalt ) values ('$chr".$columns[0]."',".$newentry_start.",'$newref','$newalt','".$filter."',$nfe_homref,".$het_nfe[$i].",".$hom_nfe[$i].",".$afr_homref.",".$het_afr[$i].",".$hom_afr[$i].",".$popmax_homref.",".$het_popmax[$i].",".$hom_popmax[$i].") ON DUPLICATE KEY UPDATE ".$exomeprefix."filter='".$filter."', ".$exomeprefix."ea_homref=".$exomeprefix."ea_homref + ".$nfe_homref .", ".$exomeprefix."ea_het=".$exomeprefix."ea_het + ".$het_nfe[$i].", ".$exomeprefix."ea_homalt=".$exomeprefix."ea_homalt + ".$hom_nfe[$i].", ".$exomeprefix."aa_homref=".$exomeprefix."aa_homref + ".$afr_homref.", ".$exomeprefix."aa_het=".$exomeprefix."aa_het + ".$het_afr[$i].", ".$exomeprefix."aa_homalt=".$exomeprefix."aa_homalt + ".$hom_afr[$i].", ".$exomeprefix."popmax_homref=".$exomeprefix."popmax_homref + ".$popmax_homref.", ".$exomeprefix."popmax_het=".$exomeprefix."popmax_het + ".$het_popmax[$i].", ".$exomeprefix."popmax_homalt=".$exomeprefix."popmax_homalt + ".$hom_popmax[$i].";" );
			}
		}
	}
	close EX;
}


#1000genomes
if($tgenomesFile ne ""){
	if($tgenomesFile =~ /^http/ || $tgenomesFile =~ /^ftp/){		# --> if URL use curl
		$tgenomesFile = "curl -s ".$tgenomesFile." | gunzip -c |"
	}else{
		$tgenomesFile = "gunzip -c ".$tgenomesFile." |" if $tgenomesFile =~ /\.gz$/;
	}
	open TG, $tgenomesFile or exit $logger->error("Can't open $tgenomesFile!");
	
	while(<TG>){
		next if $_ =~ /^#/; #skip header 
		chomp;
		my @columns = split("\t");
		
		my @info = split(";",$columns[7]);
		my @af;
		my $dp = 0;
		foreach my $entry (@info){
			my ($key,$value) = split("=",$entry);
			if($key eq "VT" && !($value eq "SNP" || $value eq "INDEL")){
				next;
			}elsif($key eq "AF"){
				@af = split(",",$value);
			}elsif($key eq "DP"){
				$dp = $value;
			}
		}
		
		my @alleles = split(",",$columns[4]);
	
		for(my $i = 0; $i < @alleles; $i++){
			executeQuery("insert into $tgenomestable (chrom,pos,rssnp,ref,alt,quality,pass,dp,af) values ('$chr$columns[0]',$columns[1],'$columns[2]','$columns[3]','$alleles[$i]',$columns[5],'$columns[6]',$dp,$af[$i])");
		}
	}
	
	close TG;
}

#Kaviar
if ($kaviarFile ne "")
{
	#
	# Data source:
	#
	# 	/data/mirror/goldenpath/hg19/database/Kaviar/Kaviar-160204-Public/vcfs/Kaviar-160204-Public-hg19.vcf.gz
	#
	#	#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
	#	1	10001	.	T	C	.	.	AF=0.0000379;AC=1;AN=26378;DS=SS6004475
	#	1	10002	.	A	C,T	.	.	AF=0.0001137,0.0000379;AC=3,1;AN=26378;DS=HGDP00927|HGDP00998|HGDP01284,HGDP01029
	#
	# 	INFO Fields to parse Kaviar
	#
	# 	AF - MULTIPLE - Allele frequency(-ies)
	# 	AC - MULTIPLE - RAW ALLELE COUNT
	# 	AN - SINGLE, TO DUPLICATE - TOTAL NUMBER OF SAMPLES
	# 	DS - MULTIPLE - Reference sources (cooperations - PIPE "|" separated)
	#
	# 	END - VARIABLE NUMBER OF VALUES - IMPOSSIBLE TO GET WHICH ALLELE IT IS REFERRED TO (TO TRASH!)
	#
	
	if($kaviarFile =~ /^http/ || $kaviarFile =~ /^ftp/){		# --> if URL use curl
		$kaviarFile = "curl -s ".$kaviarFile." | gunzip -c |"
	}else{
		$kaviarFile = "gunzip -c ".$kaviarFile." |" if $kaviarFile =~ /\.gz$/;
	}
	open KAVIAR, $kaviarFile or exit $logger->error("Can't open $kaviarFile!");
	
	# Cycle through entries 
	while(<KAVIAR>){
		
		#skip header 
		next if $_ =~ /^#/;
		chomp;
		my @columns = split("\t");
		

		my $chrom = $columns[0];
		my $start = $columns[1];
		my $ref   = $columns[3];
		my @alt = split(",",$columns[4]);		
		my $n_alt = @alt;		#get number of alternative alleles		
		
		my @info = split(";",$columns[7]);
		
		
		# Process info entries

		# Data store - renewed every cycle
		my @af;	#allele_frequency
		my @ac;	#raw allele count
		my $an; #sample count
		my @ds; #references
		
		foreach my $entry (@info)
		{
			# Split INFO (key value)
			my ($key,$value) = split("=",$entry);

			if($key eq "AF"){
				@af = split(",",$value);
			}elsif($key eq "AC"){
				@ac = split(",",$value);
			}elsif($key eq "AN"){
				$an = $value;
			}elsif($key eq "DS"){
				@ds = split(",",$value);
			}
		}
		
		
		# Cycling through alt alleles and load them into DB
		for(my $alt_loop = 0; $alt_loop < $n_alt; $alt_loop++)
		{	

			my $entry_start=$start;
			my $alt = $alt[$alt_loop];
			# Entry format must be reduced to minimal representation
			($entry_start, $ref, $alt)=Utilities::getMinimalRepresentation($entry_start, $ref,$alt);
			
			if ( defined $ds[$alt_loop] )
			{ 
				$ds[$alt_loop] =~ s/\|/\,/g;
				$ds[$alt_loop] =~ s/\'/\''/g;
				$ds[$alt_loop] = "'".$ds[$alt_loop]."'";
			}
			else
			{
				$ds[$alt_loop]="NULL";
			}
				
			
			executeQuery("insert ignore into $kaviartable (chrom,start,refallele,allele,af,ac,an,ds) values ('$chr".$chrom."',".$entry_start.",'$ref','$alt',".$af[$alt_loop].",".$ac[$alt_loop].",".$an.",".$ds[$alt_loop].");");
		}
	}
	close KAVIAR;

}




################################### subs #########################################
sub executeQuery {
	my $sql = shift;
	#print $sql."\n";
	$logger->debug($sql);
    $sth = $dbh->prepare($sql) || $logger->error("Can't prepare statement: $DBI::errstr");
    $sth->execute() || $logger->error($DBI::errstr.":::::::::::::::: Query: ".$sql);
}


=head1 NAME

fillAnnotationTables.pl

=head1 SYNOPSIS

fillAnnotationTables.pl -hi pgen.1001154.s002.tx -cv variant_summary.txt -d

=head1 DESCRIPTION

This script can be used to fill tables in the database containing additional external annotations shown in the
results tables of the Web-Interface. Table names must be defined in the current.config.xml file and data sources
must be defined in the respective comments. 

CADD:
   Since the precomputed table of CADD is very large (~ 80GB) this script accepts directly the HTTP link to the CADD
   server. You can either give a BED file with regions you want to insert (-cr) into the database (e.g. your exome target region)
   or insert all positions. Find the current file here: http://cadd.gs.washington.edu/download
   NOTE: If you are importing the whole file it takes several hours and the table uses ~260 GB of disk space!
   
Polyphen2 & SIFT:
   For simplicity Polyphen2 & SIFT predictions are taken from dbNSFP, a precomputed collection of many prediction and
   conservation scores. dbNSFP can be dowloaded here: https://sites.google.com/site/jpopgen/dbNSFP
   NOTE: This script is written for version 3 of dbNSFP. In other versions the numbering of the columns may be different
   		 and this script might not work.

ClinVar:
   The script requires the ClinVar VCF file. The current version of the file can be found here:
   ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz
   
Haploinsufficiency:
   The haploinsufficiency for each gene can be downoladed from the Supporting Information of Huang et al. PLoS genetics 2010;6;10;e1001154:
   http://europepmc.org/articles/PMC2954820/bin/pgen.1001154.s002.txt
   NOTE: The file is for hg19/GRCh37 only!

ExAC:
   The script requires the ExAC VCF file which can be found here:
   ftp://ftp.broadinstitute.org/pub/ExAC_release/
   
GnomAD:
   Same format as ExAC, but for genomes. Plus exome tracks. Download all the VCFs into a folder.

1000 Genomes:
   The script requires the 1000 Genomes VCF file (only the file with summed allele frequencies; no individual genotypes are required):
   ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz
   NOTE : Only SNPs and Indels are imported (VT=SNP or VT=INDEL)
   NOTE2: If you only want to import variants from a certain region (e.g. exome target region) you can use tabix to download only these variants:
          tabix -B ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz targets.bed > regions_of_interest.vcf
          
Kaviar:
   The script requires the Kaviar VCF file Kaviar-160204-Public/vcfs/Kaviar-160204-Public-hg19.vcf.gz contained into the KAVIAR bundle:
   http://s3-us-west-2.amazonaws.com/kaviar-160204-public/Kaviar-160204-Public-hg19.vcf.tar
   Please note that the version with data sources is required.
   	          

=head1 OPTIONS

 General:
 -se	name of the settings in the current.config.xml file that holds path to reference genome, 
 	to the annotation file and to possible additional annotation files; use default settings if nothing is given
 -chrprefix  reference contigs have "chr" as prefix (overrides automatic check)
 -nochrprefix reference contigs don't have "chr" as prefix (overrides automatic check)
 -d     delete entries of the choosen data sets before inserting the new entries
 -lf	log file; default: print to screen
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page
 
 CADD:
  -c	Path to CADD file. Either local or HTTP (from http://cadd.gs.washington.edu/download)
  	Must be zipped with bgzip and indexed with tabix (original files are by default).
  -cr	<regions.bed> regions from the CADD file that should be inserted. default: insert everything
  	NOTE: Chromosomes in the original CADD file don't start with "chr"
 
 Polyphen2 & SIFT:
  -db	</path/to/dbNSFP/> folder containing the unzipped dbNSFP tables (one per chromosome)
  -p	insert polyphen scores
  -s	insert sift scores
 
 ClinVar:
  -cv	<clinvar.vcf.gz> Tab seperated ClinVar release file 
	
 Haploinsufficiency:
  -hi	<pgen.1001154.s002.txt> haploinsufficiency per gene file

 ExAC:
  -e	<ExAC.rX.X.sites.vep.vcf.gz> ExAC VCF file (can be streamed directly from the web)
  
 GnomAD:
  -g    </path/to/GnomAD> where GnomAD VCF files have been downloaded. Replaces ExAC
  -gex	</path/to/GnomAD/gnomad.exome.XXX.vcf> direct path to gnomad exome file. Not necessary if -g enabled it finds it automatically
  -gc	</path/to/GnomAD/constraints_gene_file.bgz> direct path to gnomad gene constraints file.
  
 1000 Genomes:
  -t	<ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz> 1000 Genomes VCF file (can be streamed directly from the web)

 Kaviar:
  -kav  <Kaviar-160204-Public-hg19.vcf.gz> Kaviar VCF file

=head1 AUTHOR

Riccardo Berutti, Thomas Wieland

=cut
