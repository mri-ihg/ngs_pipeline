#!/usr/bin/perl

use strict;
use Getopt::Long qw(:config no_ignore_case);
use File::Basename;
use File::Copy;
use Cwd qw(abs_path);
#use Log::Log4perl qw(get_logger :levels);
use IO::File;
use Pod::Usage;
umask(002);

my $prog_path = dirname( abs_path($0) );
require $prog_path."/Utilities.pm";
require $prog_path."/CollectMetrics.pm";


my $bamfile = "";
my $outdir  = "";
my $help	= 0;
my $man     = 0;
my $settings= "";
my $logfile  	  = "pipeline.log";
my $loglevel 	  = "INFO";
my $bins    = 500;


GetOptions(
"b=s"  => \$bamfile,
"o=s"  => \$outdir,
"se=s" => \$settings,
"i=s"  => \$bins,
"lf=s" => \$logfile,
"ll=s" => \$loglevel,
"man"  => \$man,
"h"	   => \$help
);

pod2usage( {-exitval => 0  ,-verbose => 1} ) if $help;
pod2usage( {-exitval => 0  ,-verbose => 2} ) if $man;
pod2usage( {-exitval => 1  ,-verbose => 1} ) if $bamfile eq "" || $settings eq "";

Utilities::initLogger($logfile,$loglevel);
my $dbh = Utilities::connectExomeDB($settings);

my $logger = Utilities::getLogger();

$outdir = dirname($bamfile) if $outdir eq "";



my $params    = Utilities::getParams();
my $cnvnator  = $params->{programs}->{cnvnator}->{path};
my $rootfile  = $outdir."/".basename($bamfile);
$rootfile     =~ s/bam$/root/;
my $outfile   = $outdir."/".basename($bamfile);
$outfile      =~ s/bam$/cnvnator\.out/;
my $vcffile   = $outfile;
$vcffile      =~ s/out$/vcf/;
my $uniqueout = $outfile;
$uniqueout    =~ s/out$/unique\.out/;
my $refname   = $params->{settings}->{$settings}->{refname};


my $coredb                    = $params->{coredb}->{database};
my $refdb					  = $params->{settings}->{$settings}->{variationdb}->{database};


#set environment variable
$ENV{LD_LIBRARY_PATH} .= ":".$params->{programs}->{cnvnator}->{rootlibs};

#only run cnvnator on normal chromosomes
my $chromString = "";
if($params->{settings}->{$settings}->{normalchromosomes}){
	$chromString = "-chrom";
	open CHRS, $params->{settings}->{$settings}->{normalchromosomes} or exit $logger->error("Can't open  $params->{settings}->{$settings}->{normalchromosomes}!");
	while(<CHRS>){
		chomp;
		my @columns = split(/\t/);
		#$columns[0]   =~ s/^chr//;
		$chromString .= " $columns[0]";
	}
	close CHRS;
}

die $logger->error("No referenceperchrom folder defined for settings $settings!") unless $params->{settings}->{$settings}->{referenceperchrom};
my $refperchr = $params->{settings}->{$settings}->{referenceperchrom};

my $command  = "$cnvnator -root $rootfile $chromString -genome $refname -unique -tree $bamfile";
$logger->info("Running CNVnator tree...");
$logger->debug("Command: ".$command);
system($command);

$command  = "$cnvnator -root $rootfile $chromString -genome $refname -his $bins -d $refperchr";
$logger->info("Running CNVnator his...");
$logger->debug("Command: ".$command);
system($command);

$command  = "$cnvnator -root $rootfile $chromString -genome $refname -stat $bins";
$logger->info("Running CNVnator stat...");
$logger->debug("Command: ".$command);
system($command);

$command  = "$cnvnator -root $rootfile $chromString -genome $refname -partition $bins";
$logger->info("Running CNVnator partition...");
$logger->debug("Command: ".$command);
system($command);

#$command  = "$cnvnator -root $rootfile $chromString -genome $refname -stat $bins";
#$logger->info("Running CNVnator stat...");
#$logger->debug("Command: ".$command);
#system($command);

$command  = "$cnvnator -root $rootfile $chromString -genome $refname -call $bins > $outfile";
$logger->info("Running CNVnator call...");
$logger->debug("Command: ".$command);
system($command);

my $cnvnator2vcf  = $params->{programs}->{cnvnator}->{cnvnator2vcf};
$command  = "$cnvnator2vcf $outfile $refperchr > $vcffile";
$logger->info("Running CNVnator2VCF...");
$logger->debug("Command: ".$command);
system($command);

$command  = "awk -F\"\t\" '((\$9 < 0.1) && (\$9 >=0))' $outfile > $uniqueout";
$logger->info("Running CNVnator awk...");
$logger->debug("Command: ".$command);
system($command);


#convert unique.out file into BED format
my $uniquebed = $uniqueout;
$uniquebed    =~ s/out$/bed/;
open(IN, "$uniqueout");
open(OUT, ">$uniquebed");

while (<IN>) {
	chomp;
	my @line=split(/\t/);
	my ($chrom,$tmp)=split(/\:/,$line[1]);
	my ($start,$end)=split(/\-/,$tmp);
	print OUT "$chrom\t$start\t$end";
	foreach (@line) {
		print OUT "\t$_";
	}
	print OUT "\n";
}


close IN;
close OUT;

#filter variants if files are defined

my $gapOut		 = $uniquebed;
$gapOut			 =~ s/bed$/minusGAP\.bed/;
my $superDupsOut = $gapOut;
$superDupsOut    =~ s/bed$/minusDup\.bed/;
my $dgvOut       = $superDupsOut;
$dgvOut			 =~ s/bed$/minusDGV\.bed/;
my $exonOut      = $dgvOut;
$exonOut		 =~ s/bed$/exons\.bed/;

my $bedtools     = $params->{programs}->{bedtools}->{path};


if($params->{settings}->{$settings}->{gap}){
	$command = "$bedtools/intersectBed -a $uniquebed -b $params->{settings}->{$settings}->{gap} -wa -v -f 0.50 > $gapOut";
	$logger->info("Filtering GAP...");
	$logger->debug("Command: ".$command);
	system($command);
}else{
	symlink($uniquebed,$gapOut);
}
if($params->{settings}->{$settings}->{genomicsuperdups}){
	$command = "$bedtools/intersectBed -a $gapOut -b $params->{settings}->{$settings}->{genomicsuperdups} -c -f 0.50 > $superDupsOut";
	$logger->info("Filtering genomicSuperDups...");
	$logger->debug("Command: ".$command);
	system($command);
}else{
	symlink($gapOut,$superDupsOut);
}
if($params->{settings}->{$settings}->{dgv}){
	$command = "$bedtools/intersectBed -a $superDupsOut -b $params->{settings}->{$settings}->{dgv} -c -f 0.50 > $dgvOut";
	$logger->info("Filtering DGV...");
	$logger->debug("Command: ".$command);
	system($command);
}else{
	symlink($superDupsOut,$dgvOut);
}

if($params->{settings}->{$settings}->{refseqexons}){
	$command = "$bedtools/intersectBed -a $dgvOut -b $params->{settings}->{$settings}->{refseqexons} -loj -f 1E-9 > $exonOut";
	$logger->info("Annotating exons...");
	$logger->debug("Command: ".$command);
	system($command);
}else{
	symlink($dgvOut,$exonOut);
}

&parse_exons($exonOut);

$outfile  = $exonOut;
$outfile  =~ s/bed$/html/;

#print HTML file
#####################################################################
my $item   = "";
my $i      = 0;
my $n      = 1;
my @labels = (
'n',
'Chr',
'Start',
'Stop',
'Type',
'UCSC',
'Lenght',
'Dosage',
'Score',
'Score',
'Score',
'Score',
'Unique',
'n Dup',
'n DGV',
'n Exons',
'Genes',
'Omim'
);

open(IN, "$exonOut");
open(OUT, ">$outfile");

#print OUT "<html>\n<table border=1>\n";
printHeader();
print OUT qq(<table border="0" cellspacing="0" cellpadding="0">
<tr>
<td>Minimum exons:</td><td><input type="text" id="minexons" name="minexons" value="0.9"></td>
<td>Maximum exons:</td><td><input type="text" id="maxexons" name="maxexons" value=""></td>
</tr>
<tr>
<td>Minimum DGV:</td><td><input type="text" id="mindgv" name="mindgv" value=""></td>
<td>Maximum DGV:</td><td><input type="text" id="maxdgv" name="maxdgv" value="5.1"></td>
</tr>
<tr>
<td>Minimum dup:</td><td><input type="text" id="mindup" name="mindup" value=""></td>
<td>Maximum dup:</td><td><input type="text" id="maxdup" name="maxdup" value="0.1"></td>
</tr>
</table>);
&tableheader("");
print OUT "<thead><tr>";
foreach (@labels) {
	print OUT "<th align=\"center\">$_</th>";
}
print OUT "</tr></thead><tbody>\n";

while (<IN>) {
	chomp;
	my @line=split(/\t/);
	print OUT "<tr><td>$n</td>";
	$i = 0;
	foreach my $item (@line) {
		if ($i == 4) {
			print OUT qq(<td><a href='http://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg19\&position=$item' title='UCSC Browser'>$item</a></td>);
		}
		elsif ($i == 16) { #omim
			$item=&omim_link($item);
			print OUT "<td>$item</td>";
		}
		else {
			print OUT "<td>$item</td>";
		}
		$i++;
	}
	print OUT "</tr>\n";
	$n++;
}

#print OUT "<table>\n</html>\n";
print OUT "</tbody></table></div>";
&tablescript("6","");

close IN;
close OUT;
########################################################################
# parse exons
########################################################################
sub parse_exons {
my $exonFile = shift;
my @line  = ();
my @tmp   = ();
my %all   = ();
my %genes = ();
my @genes = ();
my %genestmp = ();
my %exons = ();
my $chrom = "";
my $start = "";
my $end   = "";
my $genes = "";
my $tmp   = "";
my $tmpomim = "";
my $line  = "";
my $omim  = "";
my $i     = 0;

open (IN, "$exonFile");
open (OUT, ">$exonFile.tmp");

# alle Zeilen mit der gleichen strukturellen Variante in eine Zeile vereinen
while (<IN>) {
	chomp;
	@line  = split(/\t/);
	$genes = $line[-1];
	$chrom = $line[0];
	$start = $line[1];
	#for better sorting
	$start = substr("0000000000",0,10-length("$start")) . $start;
	$end   = $line[2];
	$start = substr("0000000000",0,10-length("$end")) . $end;
	if ($genes eq ".") {
		$exons{"$chrom\_$start\_$end"}=0;
	}
	else {
		$exons{"$chrom\_$start\_$end"}++;
	}
	if (!exists($genes{"$chrom\_$start\_$end"})) { 
		$genes{"$chrom\_$start\_$end"}=$genes;
		pop(@line);
		pop(@line);
		pop(@line);
		pop(@line);
		$i    = 0;
		$line = "";
		foreach (@line) {
			if ($i >= 1) {
				$line .= "\t";
			}
			$line .= "$_";
			$i++;
		}
		$all{"$chrom\_$start\_$end"}=$line;
	}
	else {
		$genes{"$chrom\_$start\_$end"}.= ";";
		$genes{"$chrom\_$start\_$end"}.= $genes;
	}
	
}

# beseitige alle doppelten Gennamen, suche Omim-Gennummer
#sortieren
foreach $chrom (sort keys %all) {
	@genes = split(/\;/,$genes{$chrom});
	%genestmp = ();
	$omim = "";
	foreach $genes (@genes) {
		#$genes = &genesymbol($genes);
		$genestmp{$genes}=$genes;
	}
	$i = 0;
	foreach $genes (sort keys %genestmp) {
		if ($i == 0) {
			$tmp = $genes;
			$omim = &omim($genes);
		}
		else {
			$tmp .= " ";
			$tmp .= $genes;
			$tmpomim   = &omim($genes);
			if ($tmpomim ne "") {
				$omim .= " ";
				$omim .= $tmpomim;
			}
		}
		$i++;
	}
	if ($omim eq "") {
		$omim = "-";
	}
	if ($i < 20) { # gib, wenn mehr als 20 Gene, nur die Anzahl an.
		print OUT "$all{$chrom}\t$exons{$chrom}\t$tmp\t$omim\n";
	}
	else {
		print OUT "$all{$chrom}\t$exons{$chrom}\t$i\t$omim\n";
	}
}

close IN;
close OUT;
move("$exonFile.tmp",$exonFile);

}
########################################################################
# genesymbol
########################################################################
sub genesymbol {
my $gene  = shift;

my $query = "SELECT geneSymbol FROM ".$refdb.".kgXref WHERE kgID='$gene'";
my $out = $dbh->prepare($query) || die print "$DBI::errstr";
$out->execute() || die print "$DBI::errstr";
$gene = "";
$gene = $out->fetchrow_array;
print "$gene\n";

return($gene);
}
########################################################################
# omim
########################################################################
sub omim {
my $gene  = shift;
my $tmp   = "";

my $query = "SELECT omim FROM gene WHERE genesymbol='$gene'";
my $out = $dbh->prepare($query) || die print "$DBI::errstr";
$out->execute() || die print "$DBI::errstr";
$gene = "";
while ($tmp = $out->fetchrow_array) {
	$gene .= $tmp;
}
return($gene);
}
########################################################################
# omim_link
########################################################################
sub omim_link {
my $omim        = shift;
my @omim        = ();
my $tmp         = "";
my $query       = "";
my $out         = "";
my $res         = "";

if (($omim ne "") and ($omim != 0)) {
	$_=$omim;
	s/^\s+//g;
	s/\s+$//g;
	$omim=$_;
	@omim = split(/\s+/,$omim);
	$omim = "";
	foreach $tmp (@omim) {
		$query = "SELECT group_concat(DISTINCT omimdisease,' ',inheritance,' ',disease separator '\n')
			  FROM ".$coredb.".omim WHERE omimgene='$tmp'";
		$out = $dbh->prepare($query) || die print "$DBI::errstr";
		$out->execute() || die print "$DBI::errstr";
		$res = $out->fetchrow_array;
		$omim .= "<a href='http://www.ncbi.nlm.nih.gov/omim/$tmp' title='$res'>$tmp</a> ";
	}
}
else {
	$omim="-";
}
return($omim);
}

########################################################################
# tableheader
########################################################################
sub tableheader {
my $width = shift;
if ($width ne "") {
	$width = "style=\"width:$width\"";
}

print OUT qq(
<style type="text/css" title="currentStyle">
\@import "/jquery/DataTables-1.9.4/media/css/demo_page.css";
\@import "/jquery/DataTables-1.9.4/media/css/demo_table.css";
\@import "/jquery/DataTables-1.9.4/extras/TableTools/media/css/TableTools.css";
</style>
<div id="container" $width>
<table  border="1" cellspacing="0" cellpadding="0" class="display" id="example"> 
);

}
########################################################################
# tablescript
########################################################################
sub tablescript {
my $numeric = shift;
my $string  = shift;

print OUT qq(
<script type="text/javascript" charset="utf-8">

\$.fn.dataTableExt.afnFiltering.push(
    function( oSettings, aData, iDataIndex ) {
        var iColumn = 15;
        var iMin = document.getElementById('minexons').value * 1;
        var iMax = document.getElementById('maxexons').value * 1;
        
        var iVersion = aData[iColumn] == "-" ? 0 : aData[iColumn]*1;
        if ( iMin == "" \&\& iMax == "" )
        {
            return true;
        }
	
        else if ( iMin == "" \&\& iVersion < iMax )
        {
            return true;
        }
        else if ( iMin < iVersion \&\& "" == iMax )
        {
            return true;
        }
        else if ( iMin < iVersion \&\& iVersion < iMax )
        {
            return true;
        }
        return false;
    }
);
\$.fn.dataTableExt.afnFiltering.push(
    function( oSettings, aData, iDataIndex ) {
        var iColumn = 14;
        var iMin = document.getElementById('mindgv').value * 1;
        var iMax = document.getElementById('maxdgv').value * 1;
         
        var iVersion = aData[iColumn] == "-" ? 0 : aData[iColumn]*1;
        if ( iMin == "" \&\& iMax == "" )
        {
            return true;
        }
        else if ( iMin == "" \&\& iVersion < iMax )
        {
            return true;
        }
        else if ( iMin < iVersion \&\& "" == iMax )
        {
            return true;
        }
        else if ( iMin < iVersion \&\& iVersion < iMax )
        {
            return true;
        }
        return false;
    }
);
\$.fn.dataTableExt.afnFiltering.push(
    function( oSettings, aData, iDataIndex ) {
        var iColumn = 13;
        var iMin = document.getElementById('mindup').value * 1;
        var iMax = document.getElementById('maxdup').value * 1;
         
        var iVersion = aData[iColumn] == "-" ? 0 : aData[iColumn]*1;
        if ( iMin == "" \&\& iMax == "" )
        {
            return true;
        }
        else if ( iMin == "" \&\& iVersion < iMax )
        {
            return true;
        }
        else if ( iMin < iVersion \&\& "" == iMax )
        {
            return true;
        }
        else if ( iMin < iVersion \&\& iVersion < iMax )
        {
            return true;
        }
        return false;
    }
);


\$(document).ready(function() {
 var oTable = \$('#example').dataTable({
 	"bPaginate":      true,
  	"bLengthChange":  true,
	"bFilter":        true,
 	"bSort":          true,
	"bInfo":          false,
	"bAutoWidth":     true,
	"iDisplayLength": -1,
	"aLengthMenu": [[-1, 100, 50, 25], ["All", 100, 50, 25]],
	"sDom": 'T<"clear">lfrtip',
	"oTableTools": {
		"sRowSelect": "multi",
		"aButtons": [  ]
	},
	"aoColumnDefs": [
		{ "sType": "numeric",
		  "aTargets": [$numeric]
		},
		{ "sType": "string",
		  "aTargets": [$string]
		}
	]
});
//column filter
\$("thead input").keyup( function () {
	oTable.fnFilter( this.value, \$("thead input").index(this) );
} ); 

\$('#minexons').keyup( function() { oTable.fnDraw(); } );
\$('#maxexons').keyup( function() { oTable.fnDraw(); } );
\$('#mindgv').keyup( function() { oTable.fnDraw(); } );
\$('#maxdgv').keyup( function() { oTable.fnDraw(); } );
\$('#mindup').keyup( function() { oTable.fnDraw(); } );
\$('#maxdup').keyup( function() { oTable.fnDraw(); } );



//oTable.fnAdjustColumnSizing(); 
//oTable.width("auto");
});
</script>
);
#		"aButtons": [ "select_all", "select_none" ]

print OUT q(
<style type="text/css">
table.dataTable tr.odd  { background-color: #f9f9f1; }
table.dataTable tr.even { background-color: #fffff3; }
table.dataTable tr.odd  td.sorting_1 { background-color: #efefef; }
table.dataTable tr.even td.sorting_1 { background-color: #f9f9f5; }
table.dataTable td {padding: 2px; }
</style>
);

}
########################################################################
# printHeader
########################################################################

sub printHeader {
my $self        = shift;
my $background  = shift;
#my $user        = shift;
#my $igvport     = shift;
#my $cgi         = new CGI;


print OUT qq(
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN"
    "http://www.w3.org/TR/html4/loose.dtd">
<html>
<head>
<title>Exome</title>
) ;


print OUT qq(
<script type="text/javascript" language="JavaScript" src="/jquery/jquery-1.8.3.min.js"></script>
<script class="jsbin" src="/jquery/DataTables-1.9.4/media/js/jquery.dataTables.js"></script>
<script class="jsbin" src="/jquery/DataTables-1.9.4/extras/TableTools/media/js/TableTools.js"></script>
<script class="jsbin" src="/jquery/DataTables-1.9.4/extras/TableTools/media/js/ZeroClipboard.js"></script>
<script type="text/javascript" language="JavaScript" src="/jquery/DataTables-1.9.4/extras/dataTables.numHtmlSort.js"></script>
<script type="text/javascript" language="JavaScript" src="/jquery/DataTables-1.9.4/extras/dataTables.numHtmlTypeDetect.js"></script>
<script language="JavaScript" src="/cal/calendar_db.js"></script>
<script type="text/javascript" src="/jquery/medialize-jQuery-contextMenu-09dffab/src/jquery.contextMenu.js"></script>
<script type="text/javascript" src="/jquery/medialize-jQuery-contextMenu-09dffab/src/jquery.ui.position.js"></script>
<script type="text/javascript" src="/jquery/medialize-jQuery-contextMenu-09dffab/prettify/prettify.js"></script>
<script type="text/javascript" src="/jquery/medialize-jQuery-contextMenu-09dffab/screen.js"></script>
<link rel="stylesheet" href="/cal/calendar.css">
<link rel="stylesheet" type="text/css" href="/jquery/medialize-jQuery-contextMenu-09dffab/src/jquery.contextMenu.css">
<link rel="stylesheet" type="text/css" href="/jquery/medialize-jQuery-contextMenu-09dffab/prettify/prettify.sunburst.css">
<link rel="stylesheet" type="text/css" href="/jquery/medialize-jQuery-contextMenu-09dffab/sreen.css">
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
td.dominant        { background-color:#ffeddc; }
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
input.readonly     { background-color:#CCCCCC; }
table              { border-color:#CCCCCC;border-style:solid;border-width:0px 1px 1px 0px;}
table.outer        { border-color:#fffff3;border-style:solid;border-width:0px 0px 0px 0px;}
td                 { border-color:#CCCCCC;border-style:solid;border-width:1px 0px 0px 1px;}
td.outer           { border-color:#fffff3;border-style:solid;border-width:0px 0px 0px 0px;}
td.n               { border-color:#fffff3;border-style:solid;border-width:0px 0px 0px 0px;}
th                 { border-color:#CCCCCC;border-style:solid;border-width:1px 0px 0px 1px;}
*.big              { font-size: 18px; font-family: Arial, Verdana, Helvetica, sans-serif;}
*                  { font-size: 12px; font-family: Arial, Verdana,  Helvetica, sans-serif;}
a                  { text-decoration:none;color:#2e6385;}
a:hover            { text-decoration:underline;}
a.menu             { font-size:14px;color:#FFFFFF; font-family:Sans-Serif,Arial,Helvetica,Verdana;text-decoration:none;}   
a.menu:hover       { text-decoration:none;color:#99a4b5;}
a.menuactive       { font-size:14px;color:#99a4b5; font-family:Sans-Serif,Arial,Helvetica,Verdana;text-decoration:none;}   
*.header           { border-color:#fffff3;border-style:solid;border-width:0px 0px 0px 0px;font-size:12px;background-color:#53739c; color:white;font-family:Sans-Serif,Arial,Helvetica,Verdana;text-decoration:none;}   
td.header          { border-color:#fffff3;border-style:solid;border-width:0px 2px 0px 0px;}
.box2 {
    height: 12px;
    overflow: auto;
    padding: 5px;
}
</style>
</head>
) ;

if ($background eq "white") {
	print OUT qq(<body bgcolor=\"#ffffff\" >);
}
else {
	print OUT qq(<body bgcolor=\"#CCCCCC\" >);
}


}

########################################################################
# printFooter
########################################################################

sub printFooter {
my $self        = shift;

print OUT qq(
</body>
</html>
);

}
########################################################################


=head1 NAME

runCNVnator.pl

=head1 SYNOPSIS

 runCNVnator.pl -b merged.rmdup.bam -se hg19

=head1 DESCRIPTION

This script is a wrapper script to run cnvnator (http://sv.gersteinlab.org/) on WG datasets. 

=head1 OPTIONS

 -b	<bamfile> to run cnvnator on; required
 -se	<settings>; required
 -o	outdir; default: directory where bamfile lies in
 -i	bin size; default: 500
 -lf	<log file>; default: pipeline.log
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page

=head1 AUTHOR

Thomas Wieland

=cut