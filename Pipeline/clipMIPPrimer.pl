#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use File::Basename;
use Cwd qw(abs_path);
use Log::Log4perl qw(get_logger :levels);
use Pod::Usage;

my $prog_path = dirname(abs_path($0));
require $prog_path."/Utilities.pm";


my $infile   = "";
my $outfile  = "";
my $help     = 0;
my $man		 = 0;
my $logfile  = "SCREEN";
my $loglevel = "INFO";
my $settings = "default";
my $offset   = 1;
my $oneBased = 0;
my $doHardClip = 0;
my $countFile= "";

my ($sample,$dup,$reads,$mapped,$mappedPer,$seq,$uncov,$x1,$x4,$x8,$x20,$avgCov,$std,$medianCov,$mstd,$autosomalCov,$onBait);
my $idsample = 0;

GetOptions(
"i=s"  => \$infile, 
"o=s"  => \$outfile,
"1"    => \$oneBased,
"hard" => \$doHardClip,
"c=s"  => \$countFile,
"lf=s" => \$logfile,
"se=s" => \$settings,
"ll=s" => \$loglevel, 
"h"    => \$help,
"man"  => \$man);

pod2usage( {-exitval => 0  ,-verbose => 1} ) if $help;
pod2usage( {-exitval => 0  ,-verbose => 2} ) if $man;
pod2usage( {-exitval => 1  ,-verbose => 1} ) if $infile eq "" || $outfile eq "";

Utilities::initLogger($logfile,$loglevel);
my $logger = Utilities::getLogger();

my $params = Utilities::getParams();
my $samtools = $params->{programs}->{samtools}->{path};

$offset = 0 if $oneBased;

my $primers;

unless($primers  = $params->{settings}->{$settings}->{mipprimers}){
	$logger->error("No MIP primer file specified in settings!");
	exit(1);
}


#read primer file and store forward and reverse primers with according length in hash
my %primer;
my %primerNames;
my %primerCounts;

open PRIM, $primers or exit $logger->error("Can't open primer file $primers!");
my $tmp = <PRIM>;
while(<PRIM>){
	chomp;
	my @columns = split();
	$columns[2] += $offset;		#add 0 based offset if needed
	$columns[3] += $offset;
	$columns[4] += $offset;
	$columns[5] += $offset;
	
	$primerCounts{$columns[-1]} = 0;
	
	if($columns[2] < $columns[4]){ #reverse primer comes first
		my @tmp = (($columns[3]-$columns[2]),($columns[5]-$columns[4]) );
		$primer{$columns[1].$columns[2]} = \@tmp;
		$primerNames{$columns[1].$columns[2]} = $columns[-1];
	}else{
		my @tmp = (($columns[5]-$columns[4]),($columns[3]-$columns[2]) );
		$primer{$columns[1].$columns[4]} = \@tmp;
		$primerNames{$columns[1].$columns[4]} = $columns[-1];
	}
}
close PRIM;

#open input and output streams
# READS are first sorted by NAME to get PAIRS
open IN, "$samtools sort -o -n -m 2G $infile $infile.tmp | $samtools view -h -f 2 -F 256 - |" or exit $logger->error("Can't pipe input BAM $infile!");

# On output READS are sorted again by position
$outfile =~ s/\.bam$//;
open OUT, "| $samtools view -hSb - | $samtools sort -m 2G - $outfile" or exit $logger->error("Can't pipe output BAM $outfile!");


#get properly paired read pairs and clip primers by adjusting the cigar string
while(<IN>){
	if($_ =~ /^\@/){
		print OUT $_;
		next;
	}
	chomp;
	my @columns = split("\t");
	
	my @forward;
	my @reverse;
	
	my $newfwlen=0;
	my $newrwlen=0;
	
	# Get forward and reverse reads
	# Since reads are sorted by NAME pairs are always one after the other	
	if($columns[1] & 16){ 		#read is on the reverse strand
		@reverse = @columns;
		my $tmp  = <IN>;
		chomp $tmp;
		@forward = split("\t",$tmp)
	}else{						#read is on the forward strand
		@forward = @columns;
		my $tmp  = <IN>;
		@reverse = split("\t",$tmp)
	}
	
	# Search for primer defined by position
	# Stored as chrX123456
	#
	#TODO bug by TW to correct!!
	#TODO chr21:12345 == chr2:112345
	
	if(defined $primer{$forward[2].$forward[3]}){		#if mapping position of read can be found in primer file 
		$primerCounts{$primerNames{$forward[2].$forward[3]}} += 2;			#add read count for mip
		
		my ($fwLength,$rwLength) = @{$primer{$forward[2].$forward[3]}};	# Get FW and REVERSE primer lengths 
		($forward[5],$offset,$newfwlen) = &softclip($forward[5],$fwLength,1,$doHardClip);		# Alter FW CIGAR (Pass CIGAR, PRIMER LEN, isForward )
		
		$forward[3] += $offset;		#set offset in mapping position if read includes insertion/deletion
		$reverse[7] += $offset;		# change offset of the mapping position of the mate for the reverse primer ## NOTE: the reverse won't by definition change its alignment offset!
		
		($reverse[5],$offset,$newrwlen) = &softclip($reverse[5],$rwLength,0,$doHardClip);		# Alter RW CIGAR (Pass CIGAR, PRIMER LEN , isForward - no -> 0 )
		
		# Change sequences
		
		if ( $doHardClip )
		{
			#$forward[9]  = behead($forward[9],  $fwLength);
			$forward[9]  = behead($forward[9],  $newfwlen);
			#$forward[10] = behead($forward[10], $fwLength);
			$forward[10]  = behead($forward[10],  $newfwlen);
			
			#$reverse[9]  = betail($reverse[9],  $rwLength);
			$reverse[9]  = betail($reverse[9],  $newrwlen);
			#$reverse[10] = betail($reverse[10], $rwLength);
			$reverse[10] = betail($reverse[10], $newrwlen);
		}
	}

	#De'Bug
	#print join("\t",@forward)."\n";
	#print join("\t",@reverse)."\n";	
	print OUT join("\t",@forward)."\n";
	print OUT join("\t",@reverse)."\n";

	
}
close IN;
close OUT;

if($countFile ne ""){		#output read counts
	open OUT, ">$countFile"or exit $logger->error("Can't open $countFile!");
	foreach my $mip(sort keys %primerCounts){
		print OUT $mip."\t".$primerCounts{$mip}."\n";
	}
	
	close OUT;	
}


sub behead {
	my $string = shift;
	my $cutlen = shift;
	return substr($string,$cutlen,length($string)-$cutlen); 
}

sub betail {
	my $string = shift;
	my $cutlen = shift;
	return substr($string,0,length($string)-$cutlen);
}



sub softclip {
	my $cigar      = shift;
	my $length     = shift;
	my $forward    = shift;
	my $hardclip   = shift;
	
	my $ret = "";
	my $offset = 0;
	my $remainingLength;
	$remainingLength = $length;
	
	 
	
	my $totLength = 0;
	
	if($forward){
		while ( $cigar =~ m/(\d+\D)/g ) {
			my $type = substr( $1, ( length $1 ) - 1, 1 );
			my $curr = substr( $1, 0, ( length $1 ) - 1 );
			
			
			if($remainingLength <= 0){			# clipping is done --> just add the rest of the cigar string to the end
				$ret .= $1;
				$totLength += ($curr );
			}else{
				if($type eq "M" || $type eq "S"){		#current part of CIGAR is "match" or "softclip"
					
					if($curr > $remainingLength){					# current stretch of M or S is long enough to harbor all remaining
						if( $type eq "M" ){
							$ret .= ($curr - $remainingLength)."M";
						}elsif( ($type eq "S") and $hardclip ){
							$ret .= ($curr - $remainingLength)."S";
						}else{
							$length += ($curr - $remainingLength);
						}
						
						$totLength += ($curr -$remainingLength);
					}
					$remainingLength -= $curr;
					
				}elsif($type eq "D"){	#current part of CIGAR is "deletion"
					if($curr > $remainingLength){					# current stretch of D is at the boundaries of the primer -> split it
						$ret .= ($curr - $remainingLength)."D";
						$length -= $remainingLength;
						$offset += $remainingLength;
						
					}else{
						$length -= $curr;
						$offset += $curr;
					}
					$remainingLength -= $curr;
				}elsif($type eq "I"){	#current part of CIGAR is "insertion"
						$length += $curr;
						$offset -= $curr;
				}
			}
		}
		#add soft clipping to begin of cigar string
		$ret = $length.($hardclip ? "H" : "S").$ret;
		$totLength += $length;
		$offset += $length; 	# we have to add the length of the softclipped bases because it gets adjusted automatically
		
		
	}else{		#reverse read
		my @revCigar;
		while ( $cigar =~ m/(\d+\D)/g ) {		#for the reverse read --> turn the cigar string around
			push(@revCigar,$1);
		}
		for(my $i = (@revCigar -1);$i>=0; $i--){
			my $type = substr( $revCigar[$i], ( length $revCigar[$i] ) - 1, 1 );
			my $curr = substr( $revCigar[$i], 0, ( length $revCigar[$i] ) - 1 );
			
			
			if($remainingLength <= 0){			# clipping is done --> just add the rest of the cigar string to the end
				$ret = $revCigar[$i].$ret;
				$totLength += ($curr );
			}else{
				if($type eq "M" || $type eq "S"){		#current part of CIGAR is "match" or "softclip"
					
					if($curr > $remainingLength){					# current stretch of M or S is long enough to harbor all remaining
						if($type eq "M"){
							$ret = ($curr - $remainingLength)."M".$ret;
						}elsif( ($type eq "S") and $hardclip ){
							$ret = ($curr - $remainingLength)."S".$ret;
						}else{
							$length += ($curr - $remainingLength);
						}
						
						
						
						$totLength += ($curr -$remainingLength);
					}
					$remainingLength -= $curr;
					
				}elsif($type eq "D"){	#current part of CIGAR is "deletion"
					if($curr > $remainingLength){					# current stretch of D is at the boundaries of the primer -> split it
						$ret .= ($curr - $remainingLength)."D";
						$length -= $remainingLength;
						
					}else{
						$length -= $curr;
					}
					$remainingLength -= $curr;
				}elsif($type eq "I"){	#current part of CIGAR is "insertion"
						$length += $curr;
				}
			}
			
		}
		
		$ret .= $length.($hardclip ? "H" : "S");
	}
	
	return ($ret,$offset,$length);
}


=head1 NAME

clipMIPPrimer.pl

=head1 SYNOPSIS

 clipMIPPrimer.pl -i in.bam -o out.bam -se settings

=head1 DESCRIPTION

This script softclips the primer regions of reads from a MIP experiment. It takes the positions
of the primers from a file specified in the settings.xml.

=head1 OPTIONS
 -i	<infile.bam>; REQUIRED
 -o	<outfile.bam>; REQUIRED
 -hard 	perform hard clip (default soft clip)
 -1	primer file is one based instead of 0 based --> don't add offset of one to coordinates
 -c	<counts.out> output number of reads per primer; default don't print
 -se	name of the settings in the current.config.xml file that holds path to reference genome, 
 	to the annotation file and to possible additional annotation files; use default settings if nothing is given
 -lf	log file; default: print to screen
 -ll	log level: ERROR,INFO,DEBUG; default: INFO
 -h	print this helptext
 -man	show man page

=head1 AUTHOR

Thomas Wieland
Riccardo Berutti

=cut