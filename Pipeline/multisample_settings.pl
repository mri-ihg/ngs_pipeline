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
use Term::ANSIColor;

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
my $isnew         = 0;
my $isupdatesql   = 0;
my $isdelete      = 0;
my $isList		  = 0;
my $force         = 0;
my $name          = "";
my $settings      = "";
my $description   = "";
my $sql_in        = "";
my $isBEDArray    = 0;

my $isMakePed     = 0;
	my $phenostring ="";
	my $disease     ="";
my $isMakeStats   = 0;
	my $group_diseasename = "NO_DISEASE";
	my $group_affected = 0;
	my $group_external = 0;
	my $group_cooperation = 0;
	my $group_project     = 0;
	my $group_assay       = 0;

my $help          = 0;
my $man           = 0;

# Extravars
my $sqlquery      = "";

# Opts
GetOptions(
    "new"            => \$isnew,
    "updatesql"		 => \$isupdatesql,
    "delete"         => \$isdelete,
    "list"			 => \$isList,
    "f"              => \$force,
    "name=s"         => \$name,
    "settings=s"     => \$settings,
    "se=s"           => \$settings,
    "description=s"  => \$description,
    "sql=s"          => \$sql_in,
    "bedArray"		 => \$isBEDArray,
    "getped"		 => \$isMakePed,
    "makeped"		 => \$isMakePed,
    "pheno=s"		 => \$phenostring,
    "disease=s"		 => \$disease,
    "makestats"      => \$isMakeStats,
    "cov_cooperation"=> \$group_cooperation,
    "group_disease=s"=> \$group_diseasename,
    "group_affected" => \$group_affected,
    "group_external" => \$group_external,
    "group_cooperation"	 => \$group_cooperation,
    "group_project"  => \$group_project,
    "group_assay"    => \$group_assay,
	"h"              => \$help,
	"man"            => \$man
);

pod2usage( { -exitval => 0, -verbose => 1 } ) if $help;
pod2usage( { -exitval => 0, -verbose => 2 } ) if $man;
pod2usage( { -exitval => 0, -verbose => 1 } ) if (( $isnew eq 0 && $isupdatesql eq 0 && $isdelete eq 0 && $isMakePed eq 0 && $isMakeStats eq 0 && $isList eq 0 ) || ( $isnew eq 1 && $isdelete eq 1 ) ) || ( $name eq "" && $isList eq 0 );


# Getting parameters:
my $params = Utilities::getParams();

	# Tools    
	# Get the CoreDB
	my $coredb                    = $params->{coredb}->{database};
	my $sampletable               = $params->{coredb}->{sampletable};
	
	my $diseasetable              = $params->{coredb}->{diseasetable};
	my $disease2sampletable       = $params->{coredb}->{disease2sampletable};
	
	my $cooperationtable          = $params->{coredb}->{cooperationtable};

	# Multisample db and tables
	my $multisampledb			  = $params->{multisampledb}->{database};
	my $batchtable                = $params->{multisampledb}->{batchtable};
    my $datafreezetable           = $params->{multisampledb}->{datafreezetable};
    my $sample2batchtable         = $params->{multisampledb}->{sample2batchtable};
    my $settingstable             = $params->{multisampledb}->{settingstable};

	# Get the multisample settings if they exist
	my $coredbh = Utilities::connectCoreDB();

	$sql = "SELECT idmssettings from $multisampledb.mssettings where name=\"$name\"";
	$sth = $coredbh->prepare($sql) || die "Error in preparing $sql";
	$sth->execute() || die "Error in executing $sql";
	
	my ($idmssettings) = $sth->fetchrow_array();
	
	if ( defined $idmssettings && $isnew )
	{
		print "New settings\nERROR: Already existing settings with name: $name.\n";
		exit(-1);
	}

	if ( ! defined $idmssettings && $isdelete )
	{
		print "Delete settings\nERROR: No settings with name: $name.\n";
		exit(-1);
	}


# Insert or DELETE
if ( $isnew )
{	
	# NEW SETTINGS
	
	if ( ( $settings eq "" ) || ( $description eq "" ) || ( $sql eq "" ))
	{
		print "To create settings, arguments are required\n";
		exit -1;
	} 
		
	#Determine if SQL is a file or not:
	if ( $sql_in =~ m/^[Ss][Ee][Ll][Ee][Cc][Tt]\ / )
	{
		$sqlquery=$sql_in;
	}
	else
	{
		if ( ! -f $sql_in )
		{
			print "Error: $sql_in not a file nor an SQL query";
			exit(-1);
		}
		
		open( SQLIN, "<$sql_in" );
		
		while(<SQLIN>)
		{
			my $line=$_;
			chomp $line;
			
			$sqlquery .= "$line ";
		}	
	}
	
	# Escape quotes
	$sqlquery =~ s/\"/\\\"/g;
	
	$sql = "
			INSERT INTO $multisampledb.$settingstable
				( name, selectquery, settings, bedarray, description )
			VALUES
				( \"$name\", \"$sqlquery\", \"$settings\", $isBEDArray, \"$description\" )
	";

	$sth = $coredbh->prepare($sql) || die ("Can't prepare statement: $DBI::errstr");
	$sth->execute() || die "Cannot create multisettings $name. SQL Error. $sql";
	
	print "Multisample settings $name created.\n";
	
	exit(0);
}
elsif ( $isupdatesql )
{
	# UPDATE SQL STRING
	
	if ( $sql eq "" )
	{
		print "To update settings, -sql \"SELECT...\" or -sql query.sql argument is required\n";
		exit -1;
	} 
		
	#Determine if SQL is a file or not:
	if ( $sql_in =~ m/^[Ss][Ee][Ll][Ee][Cc][Tt]\ / )
	{
		$sqlquery=$sql_in;
	}
	else
	{
		if ( ! -f $sql_in )
		{
			print "Error: $sql_in not a file nor an SQL query";
			exit(-1);
		}
		
		open( SQLIN, "<$sql_in" );
		
		while(<SQLIN>)
		{
			my $line=$_;
			chomp $line;
			
			$sqlquery .= "$line ";
		}
	}
	
	# Escape quotes
	$sqlquery =~ s/\"/\\\"/g;
	
	$sql = "
			UPDATE $multisampledb.$settingstable
				set selectquery=\"$sqlquery\"
				where name=\"$name\"
			";

	$sth = $coredbh->prepare($sql) || die ("Can't prepare statement: $DBI::errstr");
	$sth->execute() || die "Cannot update multisettings $name. SQL Error. $sql";
	
	print "Multisample settings $name SQL Query updated.\n";
	
	exit(0)	
}
elsif ( $isdelete )
{
	# DELETE SETTINGS
		
	# If ! force
	if (! $force)
	{
		$sql = "SELECT count(B.idmsbatch) from $multisampledb.$batchtable B 
						INNER JOIN $multisampledb.$settingstable MS on MS.idmssettings = B.idmssettings 
						WHERE
							MS.name=\"$name\"
				";
			$sth = $coredbh->prepare($sql) || die ("Can't prepare statement: $DBI::errstr");
			$sth->execute() || die ("Can't execute statement: $DBI::errstr");		
		
		my ($count_batches) = $sth->fetchrow_array();
					
		if ( $count_batches > 0 )
		{
			print "Cannot delete multisettings $name. Detected $count_batches batches.\nUse -force option to delete, batches will be deleted as well";
			exit(0);
		}
	}	
	
	# Delete datafreezes
	$sql = "DELETE from $multisampledb.$datafreezetable where idmssettings=$idmssettings";
		$sth = $coredbh->prepare($sql) || die ("Can't prepare statement: $DBI::errstr");
		$sth->execute() || die ("Can't execute statement: $DBI::errstr");		

	# Delete sample2batches
	$sql = "DELETE from $multisampledb.$sample2batchtable where idmsbatch in ( SELECT idmsbatch from $multisampledb.$batchtable where idmssettings=$idmssettings)";
		$sth = $coredbh->prepare($sql) || die ("Can't prepare statement: $DBI::errstr");
		$sth->execute() || die ("Can't execute statement: $DBI::errstr");
	
	# Delete batches
	$sql = "DELETE from $multisampledb.$batchtable where idmssettings=$idmssettings";
		$sth = $coredbh->prepare($sql) || die ("Can't prepare statement: $DBI::errstr");
		$sth->execute() || die ("Can't execute statement: $DBI::errstr");
	
	# Delete settings
	$sql = "DELETE from $multisampledb.$settingstable where idmssettings=$idmssettings";
		$sth = $coredbh->prepare($sql) || die ("Can't prepare statement: $DBI::errstr");
		$sth->execute() || die ("Can't execute statement: $DBI::errstr");


	print "Multisample settings $name deleted.\n";
	
	exit(0);
}
elsif ( $isMakePed )
{
	my $settings_sql = qq{select idmssettings, selectquery, settings, bedarray from $multisampledb.$settingstable where name=?};
	my $settings_sth = $coredbh->prepare($settings_sql);
	$settings_sth->execute($name) || die ("Error in getting settings: $settings_sql");
	
	my ($idmultisettings, $sampleselect_sql, $settings, $perBEDArray)=$settings_sth->fetchrow_array();
	
	# No multisample settings found
	if ( ! defined $idmultisettings )
	{
		print "No multisample-calling settings named $name. Create settings first\n";
		exit 1;
	}
	
	$sql = "SELECT XS.name, XS.pedigree, (SELECT XF.name from $coredb.$sampletable XF where XF.idsample=XS.father) fath, (SELECT XM.name from $coredb.$sampletable XM where XM.idsample=XS.mother) moth, XS.sex, XS.saffected, XD.symbol, XD.name, XC.name from $coredb.$sampletable XS left join $coredb.$disease2sampletable XD2S on XD2S.idsample=XS.idsample left join $coredb.$diseasetable XD on XD.iddisease=XD2S.iddisease inner join $coredb.$cooperationtable XC on XC.idcooperation=XS.idcooperation where XS.idsample in ( SELECT idsample from ( $sampleselect_sql ) XSQ ) order by XS.pedigree asc, fath asc, moth asc";
	#$sql = "SELECT S.name, S.pedigree, (SELECT F.name from $coredb.$sampletable F where F.idsample=S.father) fath, (SELECT M.name from $coredb.$sampletable M where M.idsample=S.mother) moth, sex, saffected, D.symbol, D.name from $coredb.$sampletable S left join $coredb.$disease2sampletable D2S on D2S.idsample=S.idsample left join $coredb.$diseasetable D on D.iddisease=D2S.iddisease where S.idsample in ( SELECT idsample from ( $sampleselect_sql ) SQ ) order by pedigree asc, fath asc, moth asc";
	
	$sth = $coredbh->prepare($sql) || die ("Can't prepare statement: $DBI::errstr");
	$sth->execute() || die ("Can't execute statement: $DBI::errstr");
	
	

	# Print PED Header
	my $header = "#FAM_ID\tIND_ID\tFAT_ID\tMOT_ID\tSEX\t".( $phenostring ne "" ? $phenostring : "PHENO" );
	$header .="\tCOOP" if $group_cooperation;
	$header .="\n";
	
	print $header;
	
	while ( my ($sample_name, $sample_pedigree, $sample_father, $sample_mother, $sample_sex, $sample_affected, $diseasesymbol, $diseasename, $coopname ) = $sth->fetchrow_array()  ) 
	{
		my $sample_affected_out="-9";
		
		$sample_father = 0 if ! defined $sample_father;
		$sample_mother = 0 if ! defined $sample_mother;
		
		if ( defined $sample_sex )
		{
			$sample_sex=( $sample_sex eq "male" ? 1 : $sample_sex eq "female" ? 2 : 0 );
		}
		else
		{
			$sample_sex=0;
		}

		# If disease is specified then set as affected all samples that have the disease AND are affected, if disease is not specified use saffected status
		# For cancer samples tumor/healthy tissue this should hold.
		if ( $disease ne "" )
		{
			# If Disease matches AND affected, affected, if not affected or disease do not match then unaffected
			$sample_affected_out = ( 
										( ($diseasesymbol eq $disease) || ($diseasename eq $disease) )
										?
										( $sample_affected eq "1" ? "2" : "1"  )
										:
										"1"  
									);
		}
		else
		{
			if ( defined $sample_affected )
			{
					$sample_affected_out = ( $sample_affected eq "0" ? 1 : $sample_affected eq "1" ? 2 : -9 );
			}
			else
			{
				$sample_affected_out=-9;
			}
		}		
		# Header is #FAM_ID\tIND_ID\tFAT_ID\tMOT_ID\tSEX\tPHENO(name)
		my $row = $sample_pedigree."\t".$sample_name."\t".$sample_father."\t".$sample_mother."\t".$sample_sex."\t".$sample_affected_out;
		$row .= "\t".$coopname if $group_cooperation;
		$row .="\n";
		
		print $row;
	}
	
}
elsif ( $isMakeStats ) 
{
	my $settings_sql = qq{select idmssettings, selectquery, settings, bedarray from $multisampledb.$settingstable where name=?};
	my $settings_sth = $coredbh->prepare($settings_sql);
	$settings_sth->execute($name) || die ("Error in getting settings: $settings_sql");
	
	my ($idmultisettings, $sampleselect_sql, $settings, $perBEDArray)=$settings_sth->fetchrow_array();
	
	# No multisample settings found
	if ( ! defined $idmultisettings )
	{
		print "No multisample-calling settings named $name. Create settings first\n";
		exit 1;
	}
	

	# Send Makestats
	print "Creating stat plots for settings ".color('bold').$name.color('reset')."\nSQL QUERY:\n$sampleselect_sql\n\n";
	$sampleselect_sql =~ s/\"/\'/g;
	
	# Extra opts:
	#"group_disease=s"=> \$group_diseasename,
    #"group_affected" => \$group_affected,
    #"group_external" => \$group_external,
    my $extraopts="";
    	$extraopts.=" -group_affected " if $group_affected;
    	$extraopts.=" -group_external " if $group_external;
    	$extraopts.=" -group_diseasename =  \"$group_diseasename\"" if $group_diseasename ne "NO_DISEASE";
    	$extraopts.=" -group_cooperation " if $group_cooperation;
    	$extraopts.=" -group_project "     if $group_project;
    	$extraopts.=" -group_assay "       if $group_assay;
    	
	
	system("$prog_path/multisample_stats.pl -name $name -sql \"$sampleselect_sql\" $extraopts ");
	print "[DONE]";
	
}
elsif ( $isList )
{
	print "\n";
	print "List of the available multisample settings: \n";
	
	my $settings_sql = qq{select idmssettings, name, selectquery, settings, bedarray, description from $multisampledb.$settingstable order by idmssettings desc};
	my $settings_sth = $coredbh->prepare($settings_sql);
	$settings_sth->execute() || die ("Error in getting settings: $settings_sql");

			printf "\n";
			print color("bold");
			printf "-------------------------------------------------------------------------------------------------------------------\n";
			printf "\| %-35s \| %-35s \| %-35s \|\n", "Name", "Pipeline settings", "Description";
			printf "-------------------------------------------------------------------------------------------------------------------\n";
            print color("reset");                                                                                                                           
	
	while ( my ($idmultisettings, $msname, $sampleselect_sql, $settings, $perBEDArray, $description)=$settings_sth->fetchrow_array() )
	{
			printf "\| %-35s \| %-35s \| %-35s \|\n", "$msname", "$settings", "$description";
			printf "-------------------------------------------------------------------------------------------------------------------\n";
	}
	
	print "\n";
	
	
}
else
{
	print "Error\n";
	exit (-1);
}


 
# #FUTURE COMMANDS
# -organism    
# -cov.min
# -cov.max
# -mix.max
# -project
# -coop
# -name.like
# -libtype
 


=head1 NAME

multisample_settings.pl

=head1 SYNOPSIS

multisample_settings.pl -new -name "Genomes" -settings hg19_wholegenome -description "Whole genome sequencing" -sql /path/to/select.sql || -sql "SELECT idsample from..."
multisample_settings.pl -delete -name "Genomes"
multisample_settings.pl -getPed -name "Genomes"  

=head1 DESCRIPTION

Adds, removes, edits multisample settings

=head1 OPTIONS

 -new         add new multisample settings (requires -name, -settings, -description, -sql)
 -updatesql   updates sql string (requires -name, -sql) after quality controls, for example
 -delete      deletes multisample settings (requires -name, -f to delete if multisample batches and freezes exist)
 -f           force delete
 -list		  list available multisample settings
 -name        multisample settings name
 -settings    pipeline settings ( hg19_plus, hg19_wholegenome, etc.. )
 -description add text description. \" Quotes \" are required
 -sql         path to text file containing sql selection for sample ids, or sql syntax
 -bedArray    select if genomesplits are available and sample GVCFs were produced per genome split
 
 Special Functions: 
 -makeped	  creates ped file for multisample calling post-processing
 -getped	  (syn makePed)
 -makestats   creates plots with quality parameters from exomestat.
 
 -h	          print this helptext
 -man	      show man page

=head1 AUTHOR

Riccardo Berutti

=cut

