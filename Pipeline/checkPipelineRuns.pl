#!/usr/bin/perl

use strict;
use Getopt::Long;
use File::Find;
use DBI;
use XML::Simple;
use Pod::Usage;
use File::Basename;
#use Mail::Send;
use MIME::Lite;
use Cwd qw(abs_path);
use DateTime;


my $prog_path = dirname( abs_path($0) );
require $prog_path."/Utilities.pm";

my $help        = 0;
my $man         = 0;
my $mailto      = "";
my $noStdout    = 0;
my $hideEntries = 0;
my $weeks		= -1;

GetOptions(
"m=s"  => \$mailto,
"w=s"  => \$weeks,
"n"   => \$noStdout,
"i"   => \$hideEntries,
"h"   => \$help,
"man" => \$man
);

#pod2usage( {-exitval => 1  ,-verbose => 1} ) if $infile eq "" || $settings eq "";
pod2usage( {-exitval => 0  ,-verbose => 1} ) if $help;
pod2usage( {-exitval => 0  ,-verbose => 2} ) if $man;

my @mailTo = split(",",$mailto);

#get running jobs
my %runningJobs;
open JOBS, "perl $prog_path/checkForNewRuns.pl -os|" or die "Can't open stream from checkForNewRuns.pl!\n";
while(<JOBS>){
	chomp;
	my ($type, $sample, $project, $current_stage, $running) = split("\t");
	$runningJobs{$type."_".$sample}->{stage}   = $current_stage ;
	$runningJobs{$type."_".$sample}->{running} = $running ;
	#print "test: $type\_$sample\n";
}
close JOBS;


#get not hidden jobs from database
my $params      = Utilities::getParams();
my $solexadb    = $params->{solexadb}->{database};
my $sampleTable = $params->{coredb}->{sampletable};

#if only the data from the last N weeks should be shown, get the date
my $dateoffset = "";
if($weeks > -1){
	
	my $tz = DateTime::TimeZone->new( name => 'local' );
	my $dt = DateTime->now;
	$dt->subtract( weeks => $weeks );
	$dateoffset = "and u.rdate>='".$dt->year."-".$dt->month."-".$dt->day."'";
	
}


my $sql = "select s.name,t.ltlibtype,a.lplibpair,s.sex,r.pdescription,p.status,p.starttime,p.endtime,p.numofvars,p.errors,p.logfile,p.seq,p.duplicates,p.cov20x,p.exomedepthrsd,p.sry,p.mix,p.currentsettings,p.idpipeline,s.idsample,p.idlibtype,p.idlibpair,p.numofsvs,p.numofcnvs
from $sampleTable s
inner join pipeline p on s.idsample=p.idsample
inner join project r on s.idproject=r.idproject
inner join $solexadb.libtype t on t.ltid=p.idlibtype
inner join $solexadb.libpair a on a.lpid=p.idlibpair
where p.hide=0;";

my $dbh = Utilities::connectCoreDB();
my $out = $dbh->prepare($sql) || die print "$DBI::errstr\n";
$out->execute || die print "$DBI::errstr\n";

my (@finishedSamples,@runningSamples,@failedSamples);

while ( my @columns = $out->fetchrow_array ) {
	for(my $i = 0; $i < @columns; $i++){
		$columns[$i] = "" unless defined $columns[$i];
	}
	
	#get run information
	$sql = "select rdate, rname 
from solexa.sample2library sl  
inner join solexa.library i on i.lid=sl.lid
inner join solexa.library2pool lp on lp.lid=i.lid 
inner join solexa.lane l on lp.idpool=l.idpool 
inner join solexa.run u on u.rid=l.rid 
where sl.idsample=$columns[19] and i.libtype=$columns[20] and i.libpair=$columns[21] $dateoffset
order by u.rdate 
desc limit 1;";
	
	#print STDERR "$sql \n";
	
	my $out2 = $dbh->prepare($sql) || die print "$DBI::errstr\n";
	$out2->execute || die print "$DBI::errstr\n";
	my ($rdate,$rname);
	unless(($rdate,$rname) = $out2->fetchrow_array ){
		next if $weeks > 0; 							#skip sample if no new run can be found
		$rdate = "";
		$rname = "";
	}
	push(@columns,($rdate,$rname));

	if($columns[5] eq "finished"){	#pipeline finished correctly
		push(@finishedSamples,\@columns);
		
	}else{		#samples marked as running
		#check if sample is really running
		if($runningJobs{uc(substr($columns[1],0,1).substr($columns[2],0,1))."_$columns[0]"}){
			push(@runningSamples,\@columns);
		}else{
			my $sql2 = "update pipeline set status='failed' where idpipeline=$columns[18]";
			my $out2 = $dbh->prepare($sql2) || die print "$DBI::errstr\n";
			$out2->execute || die print "$DBI::errstr\n";
			push(@failedSamples,\@columns);
		}
		
	}
}

exit(0) unless @finishedSamples || @runningSamples || @failedSamples;	#exit if the pipeline hasn't been run

&printToSTDOUT() unless $noStdout;		#print to stdout if chosen
&sendMail() if $mailto ne "";			#send e-mail if chosen


#hide entries if chosen
if($hideEntries){
	$sql = "update pipeline set hide=1 where status='finished' or status='failed'";
	$out = $dbh->prepare($sql) || die print "$DBI::errstr\n";
	$out->execute || die print "$DBI::errstr\n";
}



#########################################################################
sub printToSTDOUT {
	#failed samples
	if(@failedSamples){
		print "\nFailed Samples\n\nSample\tLibtype\tLibpair\tProject\tRundate\tRunname\tStarttime\tLogfile\tSettings\tidpipeline\n--------------------------------------------------------------------------------------\n";
		foreach my $tmp(@failedSamples){
			my @columns = @$tmp;
			print "$columns[0]\t$columns[1]\t$columns[2]\t$columns[4]\t$columns[24]\t$columns[25]\t$columns[6]\t$columns[10]\t$columns[17]\t$columns[18]\n"
		}
	}
	if(@runningSamples){
		print "\nRunning Samples\n\nSample\tLibtype\tLibpair\tProject\tRundate\tRunname\tStarttime\tCurrent Program\tJob Running\tidpipeline\n-------------------------------------------------------------------------------------------------------------\n";
		foreach my $tmp(@runningSamples){
			my @columns = @$tmp;
			print "$columns[0]\t$columns[1]\t$columns[2]\t$columns[4]\t$columns[24]\t$columns[25]\t$columns[6]\t".$runningJobs{uc(substr($columns[1],0,1).substr($columns[2],0,1))."_$columns[0]"}->{stage}."\t".$runningJobs{uc(substr($columns[1],0,1).substr($columns[2],0,1))."_$columns[0]"}->{running}."\t$columns[18]\n"
		}
	}
	if(@finishedSamples){
		print "\nFinished Samples\n\nSample\tLibtype\tLibpair\tProject\tRundate\tRunname\tStarttime\tEndtime\tErrors\tLogfile\tSeq\tDuplicates\tCov20x\tRsd\tSex\tSRY\tMix\tVariants\tidpipeline\n----------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n";
		foreach my $tmp(@finishedSamples){
			my @columns = @$tmp;
			print "$columns[0]\t$columns[1]\t$columns[2]\t$columns[4]\t$columns[24]\t$columns[25]\t$columns[6]\t$columns[7]\t$columns[9]\t$columns[10]\t$columns[11]\t$columns[12]\t$columns[13]\t$columns[14]\t$columns[3]\t$columns[15]\t$columns[16]\t$columns[8],$columns[22],$columns[23]\t$columns[18]\t\n"
		}
	}
}

#########################################################################
sub sendMail {

	my $msg = MIME::Lite->new(
         To      =>$mailto,
         Subject =>'Pipeline Update',
         Type    =>'multipart/related'
    );
    my $body = "";
	
	$body .= "\n". &getDocHeader();
	if(@failedSamples){
		$body .= "\n". qq(
		<span class="bigbold">Failed Samples</span>
		
		<table style="border:1" cellspacing="0" cellpadding="2"> 
	<tr><th align="center">Sample</th><th align="center">Libtype</th><th align="center">Libpair</th><th align="center">Project</th><th align="center">Rundate</th><th align="center">Runname</th><th align="center">Starttime</th><th align="center">Logfile</th><th align="center">Settings</th><th align="center">idpipeline</th></tr> );
		
		foreach my $tmp(@failedSamples){
			my @columns = @$tmp;
			$body .= "\n". qq(<tr> <td >$columns[0]</td><td >$columns[1]</td><td >$columns[2]</td><td >$columns[4]</td><td >$columns[24]</td><td >$columns[25]</td><td >$columns[6]</td><td ><a href="https://ihgseq3.helmholtz-muenchen.de$columns[10]">LOG</a></td><td >$columns[17]</td><td >$columns[18]</td> </tr> 
			);
		}
		
		$body .= "\n". qq(</table> <br><br><br>
		);
	}
	if(@runningSamples){
		$body .= "\n". qq(
		<span class="bigbold">Running Samples</span>
		
		<table style="border:1" cellspacing="0" cellpadding="2"> 
	<tr><th align="center">Sample</th><th align="center">Libtype</th><th align="center">Libpair</th><th align="center">Project</th><th align="center">Rundate</th><th align="center">Runname</th><th align="center">Starttime</th><th align="center">Current Program</th><th align="center">Running</th><th align="center">idpipeline</th></tr> );
		foreach my $tmp(@runningSamples){
			my @columns = @$tmp;
			my $stage   = $runningJobs{uc(substr($columns[1],0,1).substr($columns[2],0,1))."_$columns[0]"}->{stage};
			my $running = $runningJobs{uc(substr($columns[1],0,1).substr($columns[2],0,1))."_$columns[0]"}->{running};
			my $runningClass = "class=\"bad\"";
			$runningClass    = "class=\"good\"" if $running eq "Yes";
			$body .= "\n". qq(<tr> <td >$columns[0]</td><td >$columns[1]</td><td >$columns[2]</td><td >$columns[4]</td><td >$columns[24]</td><td >$columns[25]</td><td >$columns[6]</td><td >$stage</td><td $runningClass>$running</td><td >$columns[18]</td> </tr> 
			);
		}
		$body .= "\n". qq(</table><br><br><br>
		);
	}
	if(@finishedSamples){
		$body .= "\n". qq(
		<span class="bigbold">Finished Samples</span>
		
		<table style="border:1" cellspacing="0" cellpadding="2"> 
	<tr><th align="center">Sample</th><th align="center">Libtype</th><th align="center">Libpair</th><th align="center">Project</th><th align="center">Rundate</th><th align="center">Runname</th><th align="center">Starttime</th><th align="center">Endtime</th><th align="center">Errors</th><th align="center">Logfile</th><th align="center">Seq</th><th align="center">Duplicates</th><th align="center">Cov20x</th><th align="center">Rsd</th><th align="center">Sex</th><th align="center">SRY</th><th align="center">Mix</th><th align="center">Variants(SNVs,pindel,exomedepth)</th><th align="center">Settings</th><th align="center">idpipeline</th></tr> 
	);
		
		
		foreach my $tmp(@finishedSamples){
			my @columns = @$tmp;
			my $errorClass = "";
			$errorClass    = "class=\"bad\"" if $columns[9]>0;
			my $seqClass   = "";
			$seqClass      = "class=\"bad\"" if $columns[11] < 8;
			my $dupClass   = "";
			$dupClass      = "class=\"bad\"" if $columns[12] > 0.3;
			my $cov20Class = "";
			$cov20Class    = "class=\"bad\"" if $columns[13] < 93;
			my $rsdClass   = "";
			$rsdClass      = "class=\"bad\"" if $columns[14] > 2.5;
			my $sryClass   = "";
			$sryClass      = "class=\"bad\"" if ($columns[15] > 10 && $columns[3] eq "female") || ($columns[15] < 200 && $columns[3] eq "male");
			my $mixClass   = "";
			$mixClass      = "class=\"bad\"" if $columns[16] > 0.02;
			my $varClass   = "";
			$varClass      = "class=\"bad\"" if  ($columns[8] < 70000 || $columns[8] > 80000 ) || ($columns[22]<200 || $columns[22]>1000) || ($columns[23]<10 || $columns[23]>800);
			
			$body .= "\n". qq(<tr> <td >$columns[0]</td><td >$columns[1]</td><td >$columns[2]</td><td >$columns[4]</td><td >$columns[24]</td><td >$columns[25]</td><td >$columns[6]</td><td >$columns[7]</td><td $errorClass>$columns[9]</td><td ><a href="https://ihgseq3.helmholtz-muenchen.de$columns[10]">LOG</a></td><td $seqClass>$columns[11]</td><td $dupClass>$columns[12]</td><td $cov20Class>$columns[13]</td><td $rsdClass>$columns[14]</td><td >$columns[3]</td><td $sryClass>$columns[15]</td><td $mixClass>$columns[16]</td><td $varClass>$columns[8],$columns[22],$columns[23]</td><td >$columns[17]</td><td >$columns[18]</td>);
		}
		$body .= "\n". qq(</table><br><br><br>
		);
	}
	
	$body .= "\n". &getDocFooter();
	
	$msg->attach(
        Type => 'text/html',
        Data => $body
    );

    $msg->send();
}




#########################################################################
sub getDocHeader {
	return qq(<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN" 
	"http://www.w3.org/TR/html4/loose.dtd">

<html>
<head>

	
        

<style type="text/css">
body               { background-color:#fffff3; }
td.person          { background-color:#e88a6a; }
td.default         { background-color:#fffff3; }
td.blue            { background-color:#c9e5e3; }
td.dna             { background-color:#aedecf; }
td.dnaReport       { background-color:#d7e2be; }
td.dnaReportMethod { background-color:#c9e5e3; }
td.cytoMaterial    { background-color:#c9e5e3; }
td.cytoTest        { background-color:#dbcbd3; }
td.cytoFaerbung    { background-color:#d7e2be; }
td.cytoFish        { background-color:#e5cab5; }
td.cytoReport      { background-color:#bac7db; }
td.barcode1        { background-color:#e2ecef; }
td.barcode2        { background-color:#c2d6db; }
td.barcode3        { background-color:#a6d7e3; }
td.barcode4        { background-color:#7ebdcd; }
td.barcode5        { background-color:#55a6ba; }
td.barcode6        { background-color:#3d7f8f; }
td.barcode7        { background-color:#11667b; }
td.barcode8        { background-color:#11586a; }
td.nodata		   { background-color:#FFFFFF; width: 30px }
td.bad   		   { background-color:#FF0000; width: 30px }
td.ok   		   { background-color:#FFFF00; width: 30px }
td.good   		   { background-color:#00940A; width: 30px }
td.excellent	   { background-color:#7DFF86; width: 30px }
table.fixed		   { table-layout: fixed }
input.readonly     { background-color:#CCCCCC; }
table              { border-color:#CCCCCC;border-style:solid;border-width:0px 1px 1px 0px;}
table.outer        { border-color:#fffff3;border-style:solid;border-width:0px 0px 0px 0px;}
td                 { border-color:#CCCCCC;border-style:solid;border-width:1px 0px 0px 1px;}
td.outer           { border-color:#fffff3;border-style:solid;border-width:0px 0px 0px 0px;}
th                 { border-color:#CCCCCC;border-style:solid;border-width:1px 0px 0px 1px;}
*.big              { font-size: 18px; font-family: Arial, Verdana, Helvetica, sans-serif;}
*.bigbold          { font-weight: bold; font-size: 18px; font-family: Arial, Verdana, Helvetica, sans-serif;}
*                  { font-size: 12px; font-family: Arial, Verdana,  Helvetica, sans-serif;}
a                  { text-decoration:none;color:#2e6385;}
a:hover            { text-decoration:underline;}
a.menu             { font-size:14px;color:#FFFFFF; font-family:Sans-Serif,Arial,Helvetica,Verdana;text-decoration:none;}   
a.menu:hover       { text-decoration:none;color:#99a4b5;}
a.menuactive       { font-size:14px;color:#99a4b5; font-family:Sans-Serif,Arial,Helvetica,Verdana;text-decoration:none;}   
*.header           { border-color:#fffff3;border-style:solid;border-width:0px 0px 0px 0px;font-size:12px;background-color:#53739c; color:white;font-family:Sans-Serif,Arial,Helvetica,Verdana;text-decoration:none;}   
td.header          { border-color:#fffff3;border-style:solid;border-width:0px 2px 0px 0px;}
</style>
</head>
<body> 

	
);
}


#########################################################################
sub getDocFooter {
	return qq(</body> </html> );
}

=head1 NAME

 checkPipelineRuns.pl

=head1 SYNOPSIS

 checkPipelineRuns.pl

=head1 DESCRIPTION

This script checks the table "pipeline" for all entries that have not the status "hide". It prints some stats on the commandline
or sends an E-Mail. It's meant to run as a cron job once a day to send summarys on analysed samples.

There are three types of samples:
 1) Samples that are running properly (entry "started" in DB and jobs in SGE)
 2) Samples NOT running properly (entry "started" in DB and NO jobs in SGE)
 3) Samples have finished properly (entry "finished" in DB)

=head1 OPTIONS

 -m	<mailto1@server.com,mailto2@server.com> comma seperated list of mail addresses that should receive changes
 -n	no output on commandline
 -i	hide all "finished" and "failed" entries in the DB, so they won't be displayed in the next call.
 -w	weeks; if this is given only samples from runs that have been started after today - weeks are displayed
 -h	show help text
 -man	show man page

=head1 AUTHOR

Thomas Wieland

=cut