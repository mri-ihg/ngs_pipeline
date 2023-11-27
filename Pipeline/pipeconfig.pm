package pipeconfig;

###########################################################
#
# pipeconfig.pm
#
# Pipeline making class, with class
#
# Author: 	Riccardo Berutti
# Date:		2017-07-07 
#
#
#	USAGE:
#
#   1) DECLARE OBJECT
#		my $pipelineconfig = pipeconfig->new("pipename", "/my/path/to/config.cfg", ["my comment to pipeline :)"] , [0 or 1: if not running slot must be printed]);
#
#	2) DECLARE HEADER (AMENDABLE WITH MORE PARAMETERS AT ANY TIME BEFORE WRITE)
#		$pipelineconfig->header( 
#					{ "outfolder" => "path/to/outfolder" },
#					{ "pipetype"  => "multisample" }
#		);
#
#	3) DECLARE PROGRAM
#		$pipelineconfig->program(
#					"program_name.pl",
#					"comment, this program makes fun",
# 					1, 												(if this program runs, otherwise, 0, or whatever)
#					"waitforme.pl,waitalsoforme.pl,andforme.pl", 	(no spaces please)
#					1,												(number of slots, or "perBedArray", or whatever supported by parallelpipeline.pl)
#					{ "outfolder" => "path/to/outfolder" },
#					{ "pipetype"  => "multisample" }
#		);
#	4) WRITE CONFIG FILE
#		$pipelineconfig->write();
#
############################################################
use Tie::IxHash;

# Variables
tie my %pipecfg, 'Tie::IxHash';
my @pipecfgheaderlines;
my @pipecfglines;
my $pipename = "";
my $pipecomment = "";
my $pipecfgfile="";
my $HANDLE;
my $header;
my $printnonrunningslots;

############################################################
# New : init it
############################################################
sub new
{
	
	my $class = shift;
	my $inpipename = shift;
	my $outfile = shift;
	my $comment = shift;
	my $nonrunningprint = shift;
	
	my $self  = {};
	
	if ( $inpipename eq "" )
	{
		die "Error: specify pipeline name";
	}
	
	if ( $outfile eq "" || (! defined($outfile)) ) 
	{
		die "Error specify config file."; 
	}
	
	$nonrunningprint=1 if ( ! defined $nonrunningprint ); 
	
	# Nonrunning slots are ignored by default
	$printnonrunningslots = ( $nonrunningprint eq "1" ? 1 : 0);
	
	# I'm a class now
	bless $self, $class;
	
	# Set pipeline name
	$pipename = $inpipename;
	
	# Set config file 
	$pipecfgfile=$outfile;
		
	# Header not yet inserted
	$header=0;
			
	return $self;
}

############################################################
# Pipeline header : declare pipeline header (amendable)
############################################################
sub header
{
	
	my $class = shift;

	# If I need some argument later to amend header I can re-invoke header without writing the title or the header again
	if ( $header eq 0 )
	{
		# Title and comment lines
		if ( $pipecomment eq "" )
		{
			addheaderline(title($pipename));
		}
		else
		{
			addheaderline(title($pipename,$comment));
		}
		
		# Pipeline name
		addheaderline(param("pipename", $pipename));
	}
	
	# Header key value pairs
	while (my %inargs=%{shift @_} )
    {
		foreach my $key (keys %inargs)
        {
        	$pipecfg{$key}=$inargs{$key};
        	
        	addheaderline( 
        		param($key, $pipecfg{$key})
        	);
        }
    }

	# Header has been written    
    $header=1;
	
	return 0;
}

############################################################
# Pipeline Program : declare program
############################################################
sub program
{
	my $class = shift;
	my $program = shift;
	my $comment = shift;
	my $runs = shift;
	my $waitfor = shift;
	my $slots = shift;
	
	# Do not bother to write it
	if ( $runs ne 1 && $printnonrunningslots eq 0 )
	{
		return;
	}
	
	# Write comment header
	addline(title($program,$comment));
	
	# Write program
	addline(param("pgr", $program));
	
	# Write arguments
	while (my %inargs=%{shift @_} )
    {
		foreach my $key (keys %inargs)
        {
        	next if $key eq "#";
        	$pipecfg{$key}=$inargs{$key};
        	addline(
        		param("param", $key, $inargs{$key})
        	);
        }
    }
    
    # Write dependson
    addline(param( "dependson", $waitfor ));
    
    # Write slots
    addline(param( "slots", $slots ));
    
    # Write run on config
    addline(param( $runs eq 1 ? "run" : "#run" ));
    
    # Program ends
    addline(programend());
    
    return 0;
}


############################################################
# Write : write config file - last step
############################################################
sub write
{
	my $class = shift;
	
	if ( $header == 0 )
	{
		die "Error: Insert header";
	}
	
	# Write to file:
	open (HANDLE, ">$pipecfgfile");
		
	# Print header
	foreach my $headerline ( @pipecfgheaderlines )
	{
		print HANDLE $headerline."\n";
	}
	
	# Print Programs
	foreach my $programline ( @pipecfglines )
	{
		print HANDLE $programline."\n";
	}	
	
	# Close the config gracefully
	print HANDLE spacer().title($pipename.": end of pipeline");
	
	# Close file
	close (HANDLE);
	
	#return 0;
}


############################################################
# Internal functions
############################################################

#
# ARRAYS are used to store separately header and program list
# header can be added parameters after declared
# header and programs are written in the order they are declared  
#

# Addheaderline (to header array )
sub addheaderline()
{
	my $line=shift;
	push (@pipecfgheaderlines, $line);
}

# Addline ( to program array )
sub addline()
{
	my $line=shift;
	push (@pipecfglines, $line);
}

#Parameter
sub param
{
	my $arg1=shift;
	my $arg2=shift;
	my $arg3=shift;
		$arg1="" if ! defined $arg1;
		$arg2="" if ! defined $arg2;
		$arg3="" if ! defined $arg3;
	
	# Runs or not run/#run
	return "run" if ( $arg1 eq "run" || $arg1 eq "#run" );
	# Comment if "# blah"
	return comment($arg2) if ( $arg1 eq "#" && $arg2 ne "");
	# Key : val if "key val"
	return $arg1."\t:\t".$arg2 if ($arg1 ne "" && $arg2 ne "" && $arg3 eq "");
	# Key : key : val if "key key val"
	return $arg1."\t:\t".$arg2."\t:\t".$arg3 if ($arg1 ne "" && $arg2 ne "" && $arg3 ne "");
	# Nothing if dunnowhat
	return "";
}

#Comment
sub comment
{
	my $arg = shift;
	return "# $arg\n";
}

#Spacer
sub spacer
{
	return "############################################################\n";
}

#Program end
sub programend
{
	return spacer()."\n\n";
}

#Title
sub title
{
	my $title=shift;
	
	my $out=spacer().
			comment($title);
	
	while (my $line=shift)
	{
		$out.=comment($line);	
	}
	
	$out.=spacer();
	
	return $out; 
}




1;