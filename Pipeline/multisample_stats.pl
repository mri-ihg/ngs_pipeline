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
use DateTime;
use Tie::IxHash;
use Term::ANSIColor;
use Scalar::Util qw(looks_like_number);

my $prog_path = dirname( abs_path($0) );
require $prog_path . "/Utilities.pm";
#require $prog_path . "/BEDRecord.pm";


# database
my $dbh           = "";
my $sql           = "";
my $sth           = "";
my $logfile       = "SCREEN";
my $loglevel      = "INFO";

# Options
my $sampleselect_sql  = "";
my $list          = "";

my $name          = "";
my $description   = "";

my $grouping      = 0;
my $group_affected = 0;
my $group_diseasename = "NO_DISEASE";
my $group_external = 0;
my $group_cooperation = 0;
my $group_project = 0;
my $group_assay   = 0;

# Boxplots
my $boxplots = 0;

# Landscape
my $landscape = 0;

my $suppresswarnings = 0;

my $help          = 0;
my $man           = 0;

# Extravars
my $sqlquery      = "";

# Names
my $case_name	  = "Case";
my $ctrl_name	  = "Control";

# Opts
GetOptions(
	"name=s"		     => \$name,
	"description=s"      => \$description,
    "sql=s"              => \$sampleselect_sql,
    "group_disease=s"    => \$group_diseasename,
    "group_affected"     => \$group_affected,
    "group_external"     => \$group_external,
    "group_cooperation"	 => \$group_cooperation,
    "group_project"      => \$group_project,
    "group_assay"		 => \$group_assay,
    "case_name=s"		 => \$case_name,
    "control_name=s"	 => \$ctrl_name,
    "boxplots"		=> \$boxplots,
    "landscape"		=> \$landscape,
    "silent"		=> \$suppresswarnings,
	"h"                  => \$help,
	"man"                => \$man
);

my $debug = 0;


# Grouping? 
$grouping = 1 if ( $group_diseasename ne "NO_DISEASE" || $group_affected || $group_external || $group_cooperation || $group_project );

# Boxplots and Grouping not possible together
if ( $grouping && $boxplots )
{
	die "Grouping and Boxplot are not possible together";
	exit (-1);
}

# Orientation PDF
my $orientation = "";
	$orientation = "--landscape" if $landscape; 

# Quality Parameters
my %QualParameters;

tie %QualParameters, "Tie::IxHash"; # This thingy keeps the hash sorted as I want it to be

# Stats:
%QualParameters = (
	"ST.seq" => "Sequenced (Gb)",
	"ST.reads" => "Read count",
	"ST.mapped" => "Mapped read count",
	"ST.percentm" => "Mapped read (%)",
	"ST.properlyp" => "Properly paired reads",
	"ST.duplicates" => "Duplicate reads (%)",
	"ST.opticalduplicates" => "Optical duplicate reads (%)",
	"ST.onbait" => "Reads on bait (%)",
	"ST.avgcov" => "Average coverage (X)",
	"ST.mediancov" => "Median coverage (X)",
	"ST.uncovered" => "Not covered (%)",
	"ST.cov1x" => "Percent covered >1X",
	"ST.cov4x" => "Percent covered >4X",
	"ST.cov8x" => "Percent covered >8X",
	"ST.cov20x" => "Percent covered >20X",
	"ST.tstv" => "Calculated TS/TV",
	"ST.mix" => "Contamination",
	"ST.exomedepthrsd" => "Exomedepth noise",
	"ST.mismatchrate" => "Mismatch rate",
	"ST.avgqual" => "Average base quality",
	"ST.libcomplexity" => "Estimated library complexity",
	"ST.q30fraction" => "Bases Q>30 (%)"
); 

# Service vars

	my @QualParametersList;				# Array list, to build SQL string and to get them back in order
	my %QualParametersAccumulator;		# To accumulate ALL values for histogram
	my %QualParametersAverage;			# To calculate the average
	my %QualParametersTotal;			# To know how many items are in the plot
	
	foreach my $item ( keys %QualParameters )
	{
		push @QualParametersList,  $item;
	
			$QualParametersAccumulator{$item}=[];
			$QualParametersAverage{$item}=0;
			$QualParametersTotal{$item}=0;
	}
	
	my $QualParametersListSQL=join(",", @QualParametersList);


# Stop conditions
pod2usage( { -exitval => 0, -verbose => 1 } ) if $help;
pod2usage( { -exitval => 0, -verbose => 2 } ) if $man;
pod2usage( { -exitval => 1, -verbose => 1 } ) if (( $sampleselect_sql eq "" && $list eq "" ) || ( $sampleselect_sql ne "" && $list ne "" ));


# Escape quotes
$sampleselect_sql =~ s/\"/\\\"/g;
	

# Getting parameters:
my $params = Utilities::getParams();

	# Tools    
	# Get the CoreDB
	my $coredb                    = $params->{coredb}->{database};
	my $sampletable               = $params->{coredb}->{sampletable};

	# Multisample db and tables
	my $multisampledb			  = $params->{multisampledb}->{database};
	my $batchtable                = $params->{multisampledb}->{batchtable};
    my $datafreezetable           = $params->{multisampledb}->{datafreezetable};
    my $sample2batchtable         = $params->{multisampledb}->{sample2batchtable};
    my $settingstable             = $params->{multisampledb}->{settingstable};

	# Get the multisample settings if they exist
	my $coredbh = Utilities::connectCoreDB();

	#FUNCTIONAL
	###$sql = "SELECT $QualParametersListSQL from $coredb.$sampletable S inner join $coredb.exomestat ST on S.idsample=ST.idsample where S.idsample in ( SELECT idsample from ( $sampleselect_sql ) SQ ) "; #TODO inner join variantstat (requires settings)
	$sql = "SELECT 
						$QualParametersListSQL, 
						IF(S.saffected=1, '$case_name', '$ctrl_name') as AFFECTED, 
						IF(D.name='$group_diseasename', '$case_name', '$ctrl_name') as DISEASE, 
						IF( L.lextfilepath <> '', 'External', 'Internal')  as EXTERNAL,
						C.name as COOPERATION,
						P.pname as PROJECT,
						A.name as ASSAY
			from 
						$coredb.$sampletable S 
						inner join $coredb.exomestat ST on S.idsample=ST.idsample 
						inner join exomehg19.disease2sample D2S on D2S.idsample=S.idsample
				        inner join exomehg19.disease D on D.iddisease=D2S.iddisease
        				inner join exomehg19.diseasegroup DG on DG.iddiseasegroup=D.iddiseasegroup
        				inner join exomehg19.organism O on O.idorganism=S.idorganism
        				inner join exomehg19.project P on P.idproject=S.idproject 
        				inner join exomehg19.cooperation C on C.idcooperation=S.idcooperation 
        				inner join solexa.sample2library S2L on S2L.idsample=S.idsample
        				inner join solexa.library L on L.lid=S2L.lid
        				left  join solexa.assay A on A.idassay=L.idassay 	
			where 
						S.idsample in ( SELECT idsample from ( $sampleselect_sql ) SQ ) 
						group by S.idsample "; #TODO inner join variantstat (requires settings) 		
		
	$sth = $coredbh->prepare($sql) || die "Error in preparing $sql";
	$sth->execute() || die "Error in executing $sql";
	
		
    while( my @data=$sth->fetchrow_array() )
    {
    	my $id=0;
    	
    	# Group Info
    	my $aff=$data[(scalar @QualParametersList)];
    	my $dis=$data[(scalar @QualParametersList)+1];
    	my $ext=$data[(scalar @QualParametersList)+2];
    	my $cop=$data[(scalar @QualParametersList)+3];
    	my $prj=$data[(scalar @QualParametersList)+4];
    	my $ass=$data[(scalar @QualParametersList)+5];
    	
    	print "Count:".((scalar @QualParametersList)-3)."\n";
    	print "$aff - $dis - $ext //";
    	my @group;
    		push @group, $aff if $group_affected;
    		push @group, $dis if $group_diseasename ne "NO_DISEASE";
    		push @group, $ext if $group_external;
    		push @group, $cop if $group_cooperation;
    		push @group, $prj if $group_project;
    		push @group, $ass if $group_assay;
    		
    		my $grouping_line = "";
    			$grouping_line = " ".( join "_", @group ) if $grouping;
    	
    	#print $#data."\n";
    	foreach my $item (@data)
    	{
    		my $itemname=$QualParametersList[$id];
    		
    		if ( looks_like_number($item) )
    		{
    			push @{$QualParametersAccumulator{$itemname}}, $item.$grouping_line;
    		}
  
  			if ( $debug )
  			{
				print $itemname." [".(scalar @{$QualParametersAccumulator{$itemname}} )."] ".$item."\n";  				
  			}
  			
  
    		$id++;
    		last if $id == (scalar @QualParametersList); 
    	}
    }



# Creating histograms
my $left=1;
my $html="<html><title>Stats for Dataset $name</title><body><table>";
my $warnings="";

print "=========================================================================================================\n";

foreach my $parameter (@QualParametersList)
{
	# Go:
	
	my $dataset = "";
	
	foreach my $value ( @{$QualParametersAccumulator{$parameter}} )
	{
		$dataset.=$value."\n" if defined $value;
	}
	
	open(my $parFileOut, ">", "$parameter.dat");
	print $parFileOut $dataset;
	close $parFileOut;
	
	my $parameterdescription=$QualParameters{$parameter};
	
	print "Plotting parameter: ".color('bold')."$parameterdescription \n".color('reset');


	my $silencer = "";
		$silencer = " 2>/dev/null " if $suppresswarnings; 

	if ( $dataset ne "")
	{
		system ("cat $parameter.dat | Rscript /data/isilon/users/berutti/snippets/R/plotHistToFileFlex.R \"$parameterdescription\" \"$parameterdescription\" $parameter.png         $silencer");
		system ("cat $parameter.dat | Rscript /data/isilon/users/berutti/snippets/R/plotHistToFileFlexCairo.R \"$parameterdescription\" \"$parameterdescription\" $parameter.svg    $silencer");
	
		system ("cat $parameter.dat | Rscript /data/isilon/users/berutti/snippets/R/plotHistToFileGrouping.R \"$parameterdescription\" \"$parameterdescription\" bygroup_$parameter.png  $silencer ");
	
		if ( $boxplots )
		{
			system ("cat $parameter.dat | awk -vPARAMETER=\"$parameter\" '{print PARAMETER\" \"\$1}' | Rscript /data/isilon/users/berutti/snippets/R/boxplothorizontal.R \"\" \"$parameterdescription\" $parameter.png     $silencer");
		}
		
		# Grouping for PDF
		if ( $grouping ) 
		{
			$html.="<tr>";
			$html.="<td><img src=\"$parameter.png\" width=\"50%\"></td>";
			$html.="<td><img src=\"bygroup_$parameter.png\" width=\"50%\"></td>";
			$html.="</tr>";
		}
		elsif ( $boxplots )
		{
			$html.="<tr>";
                        $html.="<td><img src=\"$parameter.png\" width=\"90%\"></td>";
			$html.="</tr>";
		}
		else
		{
			$html.="<tr>" if $left eq 1;
			$html.="<td><img src=\"$parameter.png\" width=\"50%\"></td>";
			$html.="</tr>" if $left eq 0;
			
			$left = ( $left eq 0 ) ? 1 : 0 ;	
		}
		#$html.="<td><img src=\"$parameter.svg\" width=\"50%\"></td>";

	}
	else
	{
		print "[WARNING] NO DATA FOR $parameterdescription \n";
		$warnings .="[WARNING] NO DATA FOR $parameterdescription<br>";
	}
	
	print "=========================================================================================================\n";
}

$html.="</table>";
$html.="<b>Warnings:</b><br>".$warnings if $warnings ne "" && $suppresswarnings == 0;
$html.="</body></html>";
system "echo -e \"$html\" > dataset.html";
#system "/usr/local/packages/htmldoc/bin/htmldoc dataset.html ....
system "/data/isilon/users/software/utils/htmldoc-1.9.12/htmldoc/htmldoc dataset.html --quiet --webpage --no-title --no-toc --outfile dataset.pdf $orientation 2>/dev/null /dev/null";


=head1 NAME

multisample_stats.pl

=head1 SYNOPSIS

multisample_stats.pl -sql "whatever query returns a table with an idsample field defined"
	to launch me use: 
	multisample_settings.pl -name "mydataset" -makestats
  

=head1 DESCRIPTION

Creates quality plots

=head1 OPTIONS

 -sql 		  makes stats in place

 -boxplots        draw boxplots instead of histograms ( more compact - no grouping option possible for now )

 -landscape       pdf in landscape mode
 -silent          do not print error messages  
 -h	          print this helptext
 -man	          show man page

=head1 AUTHOR

Riccardo Berutti

=cut

