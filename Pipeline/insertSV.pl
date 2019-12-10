#!/usr/bin/perl 

#Shut up Tabix and VCFTools garbage message
BEGIN { $SIG{__WARN__} = sub {warn $_[0] unless( $_[0] =~ m/Vcf.pm/ || $_[0] =~ m/Tabix.pm/  ) }; };


use strict;
use Getopt::Long;
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);
use Pod::Usage;
use Tabix;
use Data::Dumper;
use List::Util qw(min max);

my $prog_path = dirname( abs_path($0) );
require $prog_path . "/Utilities.pm";
require $prog_path . "/BEDRecord.pm";
require $prog_path . "/SequenceDictionary.pm";


# Debug timers
my %elapsed_time;
my %counter;

my $settings = "default";


my $infile		 = "";
my $sample 		 = "";
my $delete       = 0;
my $minInterSampleOverlap = 0.6; # TODO: IF WE CHANGE minOverlap!   WARNING: This function works only if minOverlap>0.5!!!!
my $tgoverlap    = 0.6;
my $dgvOverlap   = 0.5;

my $help         = 0;
my $man          = 0;
my $logfile      = "SCREEN";
my $loglevel     = "INFO";



GetOptions(
	"i=s"  => \$infile,
	"v=s"  => \$minInterSampleOverlap,
	"e=s"  => \$tgoverlap,
	"r=s"  => \$dgvOverlap,
	"se=s" => \$settings,
	"s=s"  => \$sample,
	"d"    => \$delete,
	"lf=s" => \$logfile,
	"ll=s" => \$loglevel,
	"h"    => \$help,
	"man"  => \$man
);

pod2usage( { -exitval => 0, -verbose => 1 } ) if $help;
pod2usage( { -exitval => 0, -verbose => 2 } ) if $man;
pod2usage( { -exitval => 1, -verbose => 1 } ) if $sample eq "" || $infile eq "";


my $params = Utilities::getParams();
Utilities::initLogger( $logfile, $loglevel );		#initialize logging system
my $logger = Utilities::getLogger();

require $params->{programs}->{vcftools}->{pm}
  ;    #VCFTools library to parse VCF file
  
my $dbh = Utilities::connectExomeDB($settings);

my $exomedb        = $params->{settings}->{$settings}->{exomedb}->{database};
my $svTable        = $params->{settings}->{$settings}->{exomedb}->{svtable};
my $svSampleTable  = $params->{settings}->{$settings}->{exomedb}->{svsampletable};
			
my $geneTable      = $params->{settings}->{$settings}->{exomedb}->{genetable};
my $svgeneTable    = $params->{settings}->{$settings}->{exomedb}->{svgenetable};

my $coredb      = $params->{coredb}->{database};
my $sampleTable = $params->{coredb}->{sampletable};


#get idsample
my $sql =  qq{select idsample from $coredb.$sampleTable where name = '$sample' };
$logger->debug($sql);
my $sth = $dbh->prepare($sql)
	  || $logger->error("Can't prepare statement: $DBI::errstr");
$sth->execute() || $logger->error($DBI::errstr);
my $idsample = $sth->fetchrow_array();
if ( !defined( $idsample ) ) {
	$logger->error(
		"sample name $sample not found in $coredb.$sampleTable - exiting");
	exit(1);
}



my $variationdb   = $params->{settings}->{$settings}->{variationdb}->{database};
my $knowngene     = $params->{settings}->{$settings}->{variationdb}->{genetable};




$logger->info("Opening tabix files with annotations...");
my $giaBtabix      = Tabix->new('-data' => $params->{settings}->{$settings}->{svannotation}->{giab});
my $tGenomestabix  = Tabix->new('-data' => $params->{settings}->{$settings}->{svannotation}->{tgenomessv});
my $lowcompltabix  = Tabix->new('-data' => $params->{settings}->{$settings}->{svannotation}->{lowcompl});
my $dgvtabix       = Tabix->new('-data' => $params->{settings}->{$settings}->{svannotation}->{dgv});
my $gaptabix       = Tabix->new('-data' => $params->{settings}->{$settings}->{svannotation}->{gap});
my $gensupduptabix = Tabix->new('-data' => $params->{settings}->{$settings}->{svannotation}->{gensupdup});
my $regtabix       = Tabix->new('-data' => $params->{settings}->{$settings}->{svannotation}->{regulation});
my $dbsnptabix     = Tabix->new('-data' => $params->{settings}->{$settings}->{svannotation}->{dbsnp});
#my $codingtabix    = Tabix->new('-data' => $params->{settings}->{$settings}->{svannotation}->{codingsequence});
my $codingtabix    = Tabix->new('-data' => $params->{settings}->{$settings}->{svannotation}->{codingexons});

#
# OriginalSVs
# Keeping track of the original SVs for the sample - if I re-run a sample and there are no more entries in SVSample which support old entries in SVs
# these old entries are recalculated in the end so that if they lost a supporting svsample they are properly recalculated or deleted
#
my %originalSVs;

#delete entries from svsample table
if($delete){

	#get idsvs to change the merging if necessary
	$sql = "select idsv from $svSampleTable where idsample=$idsample;";
	$logger->debug($sql);
	$sth = $dbh->prepare($sql) || $logger->error("Can't prepare statement: $DBI::errstr");
	$sth->execute() || $logger->error($DBI::errstr);
	while(my ($idsv) = $sth->fetchrow_array()){
		# Remember which are the original SVs which are supported by this sample
		$originalSVs{$idsv} = 1;
	}

	$sql = "delete from $svSampleTable where idsample=$idsample;";
	$logger->debug($sql);
	$sth = $dbh->prepare($sql) || $logger->error("Can't prepare statement: $DBI::errstr");
	$sth->execute() || $logger->error($DBI::errstr);
	
	
	# Delete orphan SVs
	$sql = "delete from $svTable where (SELECT count(distinct SVS.idsv) from $svSampleTable SVS where SVS.idsv=$svTable.idsv)=0";
	#$sql = "delete from $svTable where idsv not in (SELECT distinct SVS.idsv from $svSampleTable SVS)";
	$logger->debug($sql);
	$sth = $dbh->prepare($sql) || $logger->error("Can't prepare statement: $DBI::errstr");
	$sth->execute() || $logger->error($DBI::errstr);


	# Delete orphan SVGene
	#TODO: CHANGE
	$sql = "delete from $svgeneTable where (SELECT count(distinct SV.idsv) from $svTable SV where SV.idsv=$svgeneTable.idsv)=0";
	$logger->debug($sql);
	$sth = $dbh->prepare($sql) || $logger->error("Can't prepare statement: $DBI::errstr");
	$sth->execute() || $logger->error($DBI::errstr);
	
}

#my $idsvsample     = "";
my $svstype        = "";
my $svschrom       = "";
my $idsvsample     = "";
my $newsvsstart    = "";
my $newsvsend      = "";
my $newsvslength   = "";
my $newIdSV        = 0;
my @row            = ();
my $overlap        = "";
my $oldidsv        = "";
my $oldidsvsample  = "";
my $oldsvsstart    = "";
my $oldsvsend      = "";
my $oldsvslength   = "";



#prepare query
#SELECT svs.idsv,svs.idsvsample,svs.start,svs.end
#FROM svsample svs 
#WHERE svtype  = '$svstype'
#AND   chrom   = '$svschrom'
#AND   start  >=  ($newsvsstart - $newsvslength)
#AND   start  <=  ($newsvsend)
# TODO: IF WE CHANGE minOverlap!   WARNING: This function works only if minOverlap>0.5!!!!
$sql = "
SELECT svs.idsv,svs.idsvsample,svs.start,svs.end, svs.svlen
FROM $svSampleTable svs 
WHERE svtype  = ?
AND   chrom   = ?
AND   start  >=  (? - ?)
AND   start  <=  ?
AND   abs(svlen)  <=  ? 
AND   abs(svlen)  >=  ? 
";
my $sth2 = $dbh->prepare($sql) || $logger->error("Can't prepare statement: $DBI::errstr");


my $vcf = Vcf->new( file => $infile );
$vcf->parse_header();
my $counter = 0;
$logger->info("Inserting SVs...");
while ( my $vcfline = $vcf->next_data_hash() ) {
	
	$newIdSV = 0;
	($idsvsample,$svschrom,$newsvsstart,$newsvsend,$svstype, $newsvslength) = &insertSVSample($vcfline); #function in same script
	
	# Need absolute value for comparisons - DELs are <0
	$newsvslength=abs($newsvslength);
	$counter++;
	$logger->info("$counter variants inserted...") if $counter%5000 == 0;

	# get all overlapping SVs in svsample
	#IF
	if ( $svstype eq "INS")
	{
		# Insert overlap is defined as whatever crosses INS position and the region surrounding of size min(500,newsvslen)
		$sth2->execute( $svstype , $svschrom , $newsvsstart , 0.5*min(500,$newsvslength) , $newsvsend+0.5*min(500,$newsvslength), $newsvslength*(1/$minInterSampleOverlap), $newsvslength*($minInterSampleOverlap) ) || $logger->error($DBI::errstr);
	}
	else
	{
		$sth2->execute($svstype,$svschrom,$newsvsstart,$newsvslength,$newsvsend, $newsvslength*(1/$minInterSampleOverlap), $newsvslength*($minInterSampleOverlap)  ) || $logger->error($DBI::errstr);
	}

	while (@row = $sth2->fetchrow_array) {
		# select all overlapping SVs in svsample >= minoverlap
		my ($oldidsv,$oldidsvsample,$oldsvsstart,$oldsvsend,$oldsvslength)=@row;
		
		#TODO:
		# The search finds always the newly inserted svsample - (with no oldidsv ) - skip it
		next if ( $oldidsvsample == $idsvsample );
		#next if ( ! defined $oldidsv );
		#next if ( $oldidsv eq "" || $oldidsv eq "NULL" );
		
		# Need absolute value for comparisons - DELs are <0
		$oldsvslength = abs($oldsvslength);
		
		# Differentiate INS and the rest
		if ( $svstype ne "INS")
		{
			$overlap = min($newsvsend,$oldsvsend) - max($newsvsstart,$oldsvsstart);
			next if ($overlap <= 0);
		}
		
		next if $newsvslength == 0;
		next if $oldsvslength == 0;
		
		# Treat INS differently ()
		if (
				( $svstype eq "INS" ) ?
					(	(min($newsvslength,$oldsvslength)/max($newsvslength,$oldsvslength)) >= $minInterSampleOverlap )
					:
					(	($overlap/$newsvslength >= $minInterSampleOverlap) and ($overlap/$oldsvslength >= $minInterSampleOverlap)	)	
		) 
		{
			if($newIdSV == 0){		#first overlapping SV
				$newIdSV = $oldidsv;
			}else{					#more than one overlapping SV --> delete them
				if ( $newIdSV != $oldidsv )
				{
					
					# Update oldidsv to newidsv
					$sql = "update $svSampleTable set idsv=$newIdSV where idsv=$oldidsv";
					$logger->debug($sql);
					my $sth = $dbh->prepare($sql)
						|| $logger->error("Can't prepare statement: $DBI::errstr");
					$sth->execute() || $logger->error($DBI::errstr." - $sql");
					
					# Delete SV oldidsv
					&deleteSV($oldidsv);
				}
			}
		}
	}



	#while(my ($idsv) = $sth->fetchrow_array()){
		
	#	if($newIdSV == 0){		#first overlapping SV
	#		$newIdSV = $idsv;
	#	}else{					#more than one overlapping SV --> delete them
	#		#before deleting them --> point all svsample entries pointing on it to $newIdSV
	#		$sql = "update $svSampleTable set idsv=$newIdSV where idsv=$idsv";
	#		$logger->debug($sql);
	#		my $sth = $dbh->prepare($sql)
	#			  || $logger->error("Can't prepare statement: $DBI::errstr");
	#		$sth->execute() || $logger->error($DBI::errstr);
			
	#		$originalSVs{$idsv} = 0;
	#		#delete old entry
	#		&deleteSV($idsv);
	#	}
	#}
	
	
	#if no overlapping SV has been found --> enter a new entry
	# this is fixing insertions -> svlen is calculated after 
	if($newIdSV == 0){
		# svslen < 0 if DEL otherwise > 0
		$newsvslength = ( $svstype eq "DEL" ? -1 * abs($newsvslength) : abs($newsvslength)); 
		
		#$sql = "insert into $svTable (chrom,start,end,svtype) values ('$svschrom',$newsvsstart,$newsvsend,'$svstype')
		$sql = "insert into $svTable (chrom,start,end,svtype,svlen) values ('$svschrom',$newsvsstart,$newsvsend,'$svstype',$newsvslength)
				on duplicate key update idsv=LAST_INSERT_ID(idsv)";
		$logger->debug($sql);
		$dbh->do($sql) ||  $logger->error("Can't insert statement $sql: $DBI::errstr");
		$newIdSV = $dbh->last_insert_id(undef, undef, qw($svTable idsv)) or die $DBI::errstr;
		# Deleting svs not necessary here - it will be done by updateOrDeleteSV
	}
	
	#add link of current svsample to sv
	$sql = "update $svSampleTable set idsv=$newIdSV where idsvsample=$idsvsample";
	$logger->debug($sql);
	$sth = $dbh->prepare($sql)
		  || $logger->error("Can't prepare statement: $DBI::errstr");
	$sth->execute() || $logger->error($DBI::errstr);
	
	#update annotations
	#$logger->info("newIDSv= $newIdSV");
	&updateOrDeleteSV($newIdSV);
	# This SV must not be processed again
	$originalSVs{$newIdSV} = 0;
	
	
	

}
#
# Update all originalSVs{x}=1 
# the old SVs of the sample that were not processed are updated
# maybe there's no supporting svsample anymore for them
#
# > update all SVs that might have lost a link
$counter = 0;
foreach (keys %originalSVs){
	if($originalSVs{$_} && $_ ne ""){
		&updateOrDeleteSV($_);
		$counter++;
	}
}
$logger->info("$counter not used SVs updated or deleted!");


foreach my $key ( keys %elapsed_time )
{
	$logger->info( "$key - ". Utilities::seconds_to_ddhhmmss($elapsed_time{$key}) );
}



############## insert SV sample ###################
sub insertSVSample {
	my $vcfline = shift;
	
	my ($cndosage,$cnscore1,$cnscore2,$cnscore3,$cnscore4,$cnunique,$pialleles,$pidp,$pipercentvar,$bddp,$bdor1,$bdor2,$lppe,$lpsr,$mtgt,$mtgq) = ("NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL","NULL");
	
	
	$cndosage     = $vcfline->{INFO}->{CNDOSAGE} if defined $vcfline->{INFO}->{CNDOSAGE};
	$cnscore1     = $vcfline->{INFO}->{CNSCORE1} if defined $vcfline->{INFO}->{CNSCORE1};
	$cnscore2     = $vcfline->{INFO}->{CNSCORE2} if defined $vcfline->{INFO}->{CNSCORE2};
	$cnscore3     = $vcfline->{INFO}->{CNSCORE3} if defined $vcfline->{INFO}->{CNSCORE3};
	$cnscore4     = $vcfline->{INFO}->{CNSCORE4} if defined $vcfline->{INFO}->{CNSCORE4};
	$cnunique     = $vcfline->{INFO}->{CNUNIQUE} if defined $vcfline->{INFO}->{CNUNIQUE};
	$pialleles    = $vcfline->{INFO}->{PIALLELES} if $vcfline->{INFO}->{PIALLELES};
	$pidp         = $vcfline->{INFO}->{PIDP} if $vcfline->{INFO}->{PIDP};
	$pipercentvar = $vcfline->{INFO}->{PIPERCENTVAR} if defined $vcfline->{INFO}->{PIPERCENTVAR};
	$bddp         = $vcfline->{INFO}->{BDDP} if defined $vcfline->{INFO}->{BDDP};
	$bdor1        = "'".$vcfline->{INFO}->{BDOR1}."'" if defined $vcfline->{INFO}->{BDOR1};
	$bdor2	      = "'".$vcfline->{INFO}->{BDOR2}."'" if defined $vcfline->{INFO}->{BDOR2};
	$lppe         = $vcfline->{INFO}->{LPPE} if defined $vcfline->{INFO}->{LPPE};
	$lpsr		  = $vcfline->{INFO}->{LPSR} if defined $vcfline->{INFO}->{LPSR};
	$mtgt		  = $vcfline->{INFO}->{MTGT} if defined $vcfline->{INFO}->{MTGT};
	$mtgq		  = $vcfline->{INFO}->{MTGQ} if defined $vcfline->{INFO}->{MTGQ};
	
	
	# Recalculate or cure SVLEN	
	if ( $vcfline->{INFO}->{SVTYPE} eq "INS" )
	{
		# ins - keep reported svlen
		$vcfline->{INFO}->{SVLEN} = abs($vcfline->{INFO}->{SVLEN});
	}
	elsif ( $vcfline->{INFO}->{SVTYPE} eq "DEL" )
	{
		# del -1 * (end - start + 1)
		if ( defined $vcfline->{INFO}->{END} ) # otherwise keep svlen
		{
			$vcfline->{INFO}->{SVLEN} = -1 * abs( $vcfline->{INFO}->{END} - $vcfline->{POS} + 1 );
		}
	}
	else
	{
		# else end-start +1
		if ( defined $vcfline->{INFO}->{END} ) # otherwise keep svlen
		{
			$vcfline->{INFO}->{SVLEN} = abs( $vcfline->{INFO}->{END} - $vcfline->{POS} + 1 );
		}
	}
	
	
	#insert variant
	my $sql = "insert into $svSampleTable (idsample,chrom,start,end,svtype,svlen,caller,cndosage,cnscore1,cnscore2,cnscore3,cnscore4,cnunique,pialleles,pidp,pipercentvar,bddp,bdor1,bdor2,lppe,lpsr,mtalleles,mtgq)
	values                 ($idsample,'$vcfline->{CHROM}',$vcfline->{POS},$vcfline->{INFO}->{END},'$vcfline->{INFO}->{SVTYPE}',$vcfline->{INFO}->{SVLEN},'$vcfline->{INFO}->{CALLER}', $cndosage,$cnscore1,$cnscore2,$cnscore3,$cnscore4,$cnunique,$pialleles,$pidp,$pipercentvar,$bddp,$bdor1,$bdor2,$lppe,$lpsr,$mtgt,$mtgq)";
	$logger->debug($sql);
	my $sth = $dbh->prepare($sql)
		  || $logger->error("Can't prepare statement: $DBI::errstr");
	$sth->execute() || $logger->error($DBI::errstr);

    my $id = $dbh->last_insert_id(undef, undef, qw($svSampleTable idsvsample)) or exit $logger->error("Can't retrieve last inserted id!");
    return ($id,$vcfline->{CHROM},$vcfline->{POS},$vcfline->{INFO}->{END},$vcfline->{INFO}->{SVTYPE},$vcfline->{INFO}->{SVLEN});	
}






################### update sv table
sub updateOrDeleteSV{
	my $idsv = shift;
	

	
	#check if there are still variants pointing to this SV
	my $sql = "select count(distinct(idsample)) from $svSampleTable where idsv=$idsv;";
	$logger->debug($sql);
	my $sth = $dbh->prepare($sql) || $logger->error("Can't prepare statement: $DBI::errstr");
	$sth->execute() || $logger->error($DBI::errstr);
	my $frequency = $sth->fetchrow_array();
	
	#if no variants are pointing at SV --> delete it
	if($frequency == 0){
		$elapsed_time{"deleteSV"} -= time();
		&deleteSV($idsv);
		$elapsed_time{"deleteSV"} += time();
		$counter{"deleteSV"}++;
	}else{

		# I have to update SV
		# Updating with the new coordinates can create an entry that violates a unique key
		# If so, I must find it, and delete it - after deleting recover svsamples belonging
		# to the old SV (should not be the case)
		  
		my $checked=0;
		while($checked eq 0){
			$sql = "select idsv from $svTable v 
					where 
							v.chrom = ( select chrom from $svTable where idsv=$idsv ) and 
							v.start = (select round(avg(ss.start)) from $svSampleTable ss where ss.idsv=$idsv) and
							v.end   = (select round(avg(ss.end))   from $svSampleTable ss where ss.idsv=$idsv) and 
							v.svlen = (select round(avg(ss.svlen)) from $svSampleTable ss where ss.idsv=$idsv) and 
							v.idsv != $idsv and
							v.svtype = ( select svtype from $svTable where idsv=$idsv )
							";
			$logger->debug($sql);
			$sth = $dbh->prepare($sql) || $logger->error("Can't prepare statement: $DBI::errstr");
			#$elapsed_time{"getNewSvCoords"} -= time();
			$sth->execute() || $logger->error($DBI::errstr);
			#$elapsed_time{"getNewSvCoords"} += time();
			#$counter{"getNewSvCoords"}++;
			
			$checked=1;
			
			while ( my $existingIdSv = $sth->fetchrow_array() )
			{
				$logger->info("CUcu");
				# Delete
				$sql = "delete from $svTable where idsv=$existingIdSv";
				$logger->debug($sql);
				my $sth = $dbh->prepare($sql)
					|| $logger->error("Can't prepare statement: $DBI::errstr");
				$sth->execute() || $logger->error($DBI::errstr." - $sql");
			
				# Point old existing on svsmaple to idsv
				$sql = "update $svSampleTable set idsv=$idsv where idsv=$existingIdSv";
				$logger->debug($sql);
				my $sth = $dbh->prepare($sql)
					|| $logger->error("Can't prepare statement: $DBI::errstr");
				$sth->execute() || $logger->error($DBI::errstr." - $sql");
				$checked=0;
			}
		}			
		
		#otherwise --> update the entry
		#first with new coordinates
		$sql = "update $svTable v set
		v.start=(select round(avg(ss.start)) from $svSampleTable ss where ss.idsv=$idsv),
		v.end  =(select round(avg(ss.end))   from $svSampleTable ss where ss.idsv=$idsv),
		v.svlen=(select round(avg(ss.svlen)) from $svSampleTable ss where ss.idsv=$idsv),
		v.freq = $frequency
		where v.idsv=$idsv";
		
#		
#		$sql = "insert into $svTable (idsv) values ($idsv ) 
#		on duplicate key update
#		start=(select avg(ss.start) from $svSampleTable ss where ss.idsv=$idsv),
#		end  =(select avg(ss.end) from $svSampleTable ss where ss.idsv=$idsv),
#		svlen=(select avg(ss.svlen) from $svSampleTable ss where ss.idsv=$idsv),
#		freq = $frequency,
#		idsv=LAST_INSERT_ID(idsv)";
#		
		
		
		$logger->debug($sql);
		$elapsed_time{"calcSVsize"} -= time();
		$dbh->do($sql) || $logger->error($DBI::errstr);
#		my $updatedIdSV = $dbh->last_insert_id(undef, undef, qw($svTable idsv)) or die $DBI::errstr;
		$elapsed_time{"calcSVsize"} += time();
		$counter{"calcSVsize"}++;
		
		
#		# Merged to an already existing SV then reassing idsv to svsamples and idsv is the new (old) one
#		if ( $updatedIdSV ne $idsv )
#		{
#			$logger->info("Duplicate idsv on updateordelete: $sql");
#			$sql = "update $svSampleTable set idsv=$updatedIdSV where idsv=$idsv";
#			$logger->debug($sql);
#			$sth = $dbh->prepare($sql) || $logger->error("Can't prepare statement: $DBI::errstr");
#			$logger->info("Duplicate idsv on updateordelete SQL: $sql");
#			$sth->execute() || $logger->error($DBI::errstr);	
#			$idsv=$updatedIdSV;
#		}
		
		
		#delete old svgene entries
		$sql = "delete from $svgeneTable where idsv=$idsv;";
		$logger->debug($sql);
		$sth = $dbh->prepare($sql) || $logger->error("Can't prepare statement: $DBI::errstr");
		$elapsed_time{"deleteSVgene"} -= time();
		$sth->execute() || $logger->error($DBI::errstr);
		$elapsed_time{"deleteSVgene"} += time();
		$counter{"deleteSVgene"}++;
		
		#add new snvgene
		$sql = "INSERT IGNORE INTO $svgeneTable (idgene,idsv) SELECT distinct g.idgene,v.idsv
				FROM $svTable v
				INNER JOIN $variationdb.$knowngene kg ON (kg.chrom=v.chrom AND    v.end>=kg.txStart and v.start<=kg.txEnd )
				INNER JOIN $geneTable g on g.geneSymbol=kg.geneSymbol
				WHERE v.idsv=$idsv;";
		$logger->debug($sql);
		$sth = $dbh->prepare($sql) || $logger->error("Can't prepare statement: $DBI::errstr");
		$elapsed_time{"addSVgene"} -= time();
		$sth->execute() || $logger->error($DBI::errstr);
		$elapsed_time{"addSVgene"} += time();
		$counter{"addSVgene"}++;
		
		
		#get coordinates of SV and calculate annotations based on BED hashes
		$sql = "select v.chrom,v.start,v.end
		from $svTable v
		where v.idsv=$idsv;";
		$logger->debug($sql);
		$sth = $dbh->prepare($sql) || $logger->error("Can't prepare statement: $DBI::errstr");
		$elapsed_time{"getNewSvCoords"} -= time();
		$sth->execute() || $logger->error($DBI::errstr);
		$elapsed_time{"getNewSvCoords"} += time();
		$counter{"getNewSvCoords"}++;
		
		
		my ($chrom,$start,$end) = $sth->fetchrow_array();
		
		#$start--;	#subtract 1 from start for annotation using tabix
		#$end++;		#add 1 to end for annotation using tabix

		my $currSV = BEDRecord->new($chrom,$start,$end);
		
		#calculate 1KG overlap
		$elapsed_time{"tabix1KG"} -= time();
		my @overlap = $currSV->getOverlappingSequencesTabix($tGenomestabix,($tgoverlap*100) );		#get list of reciprocal overlapping indices
		my $af1kg   = 0;
		my $num1kg  = 0;
		my $type1kg = "";
		foreach(@overlap){
			$af1kg   += $_->name();
			$num1kg++;
			$type1kg .= ${$_->rest()}.",";
		}
		if($type1kg eq ""){
			$type1kg = "NULL";
		}else{
			$type1kg =~ s/,$//;
			$type1kg = "'".$type1kg."'";
		}
		$elapsed_time{"tabix1KG"} += time();
		$counter{"tabix1KG"}++;
		
		
		#calculate genome in a bottle overlap
		$elapsed_time{"tabixGIAB"} -= time();
		$currSV->getOverlappingSequencesTabix($giaBtabix);
		my $giaboverlap = $currSV->overlap()/100;
		$elapsed_time{"tabixGIAB"} += time();
		$counter{"tabixGIAB"}++;
		
		#calculate low complexity overlap
		$elapsed_time{"tabixLowComplexity"} -= time();
		$currSV->getOverlappingSequencesTabix($lowcompltabix);
		my $lowcomplexityoverlap = ($currSV->overlap()/100);
		$elapsed_time{"tabixLowComplexity"} += time();
		$counter{"tabixLowComplexity"}++;
		
		#calculate dgv overlap
		$elapsed_time{"tabixDGV"} -= time();
		my $dgv     = 0;
		my @dgvover = $currSV->getOverlappingSequencesTabix($dgvtabix);
		foreach(@dgvover){
			$dgv++ if ($_->overlap()/100) > $dgvOverlap;
		}
		$elapsed_time{"tabixDGV"} += time();
		$counter{"tabixDGV"}++;
		
		#calculate gap overlap
		$elapsed_time{"tabixGAP"} -= time();
		$currSV->getOverlappingSequencesTabix($gaptabix);
		my $gapoverlap = ($currSV->overlap()/100);
		$elapsed_time{"tabixGAP"} += time();
		$counter{"tabixGAP"}++;
		
		
		#calculate genomic superduplication overlap
		$elapsed_time{"tabixSuperdup"} -= time();
		my $genomicsuperduplicationoverlap     = 0;
		@dgvover = $currSV->getOverlappingSequencesTabix($gensupduptabix);
		foreach(@dgvover){
			$genomicsuperduplicationoverlap++ if ($_->overlap()/100) > $dgvOverlap;
		}
		$elapsed_time{"tabixSuperdup"} += time();
		$counter{"tabixSuperdup"}++;
		
		#calculate overlaps with coding sequence and promotors
		my $overlaps = "";
		$elapsed_time{"tabixProm"} -= time();
		$currSV->getOverlappingSequencesTabix($regtabix);
		$overlaps .= "Promoter," if $currSV->overlap() > 0;
		$elapsed_time{"tabixProm"} += time();
		$counter{"tabixProm"}++;

		$elapsed_time{"tabixCoding"} -= time();
		$currSV->getOverlappingSequencesTabix($codingtabix);
		$overlaps .= "Coding_Region," if $currSV->overlap() > 0;
		$overlaps  =~ s/,$//;
		if($overlaps eq ""){
			$overlaps  = "NULL";
		}else{
			$overlaps  = "'".$overlaps."'";
		}
		$elapsed_time{"tabixCoding"} += time();
		$counter{"tabixCoding"}++;
		
		
		
		#calculate DBSNP overlap
		$elapsed_time{"tabixDBSNP"} -= time();
		@overlap = $currSV->getOverlappingSequencesTabix($dbsnptabix,($tgoverlap*100));			#get list of reciprocal overlapping indices
		my $afdbsnp   = 0;
		foreach(@overlap){
			$afdbsnp  += $_->name();
		}
		$afdbsnp = 999 if $afdbsnp > 999;
		$elapsed_time{"tabixDBSNP"} += time();
		$counter{"tabixDBSNP"}++;
		
		
		
		#update all other annotations
		$sql ="UPDATE $svTable v SET
			   v.af1kg             = $af1kg,
     		   v.num1KG            = $num1kg,
     		   v.type1KG           = $type1kg,
     		   v.giaboverlap       = $giaboverlap,
     		   v.lowcomploverlap   = $lowcomplexityoverlap,
     		   v.dgvoverlap        = $dgv,
     		   v.gapoverlap        = $gapoverlap,
     		   v.gensupdupsoverlap = $genomicsuperduplicationoverlap,
     		   v.overlaps          = $overlaps,
     		   v.afdbsnp		   = $afdbsnp
     		   WHERE v.idsv=$idsv
			   ;";
		$logger->debug($sql);
		#print STDERR "$sql\n\n";
		$sth = $dbh->prepare($sql) || $logger->error("Can't prepare statement: $DBI::errstr");
		$elapsed_time{"dbUpdateTabix"} -= time();
		$sth->execute() || $logger->error($DBI::errstr);
		$elapsed_time{"dbUpdateTabix"} += time();
		$counter{"dbUpdateTabix"}++;
	}
}

################### delete entry in SV table
sub deleteSV{
	my $idsv = shift;
	
	#first, delete svgene entries
	my $sql = "delete from $svgeneTable where idsv=$idsv;";
	$logger->debug($sql);
	my $sth = $dbh->prepare($sql) || $logger->error("Can't prepare statement: $DBI::errstr");
	$sth->execute() || $logger->error($DBI::errstr);
	
	
	#delete sv entry
	$sql = "delete from $svTable where idsv=$idsv;";
	$logger->debug($sql);
	$sth = $dbh->prepare($sql) || $logger->error("Can't prepare statement: $DBI::errstr");
	$sth->execute() || $logger->error($DBI::errstr);
	
}






=head1 NAME

insertSV.pl

=head1 SYNOPSIS

 insertSV.pl -i merged.sv.vcf -s [SAMPLE] -se settings 

=head1 DESCRIPTION

This script inserts structural variants (SVs) that from the merged VCF file produced by "mergeSV.pl".
It calculates the overlap to SVs from other samples already in the database and all kind of overlaps 
and statistics for the variant

=head1 OPTIONS

 -s	sample name; REQUIRED
 -i	<merged.sv.vcf> VCF file containing SVs, REQUIRED
 -se	name of the settings in the current.config.xml file that holds path to reference genome, 
 	to the annotation file and to possible additional annotation files; use default settings if nothing is given
 -v	required reciprocal overlap to join two variants of two DIFFERENT samples; default: 0.9
 -e required reciprocal overlap to annotate 1000 Genome overlap; default: 0.6
 -r required reciprocal overlap to annotate GAP, DGV and GenomicSuperdup; default: 0.5
 -d	delete entries for given sample before inserting. NOTE: also deletes entries from svgene table
 -lf	log file; default: print to screen
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page

=head1 AUTHOR

Thomas Wieland, Tim Strom, Riccardo Berutti

=cut
