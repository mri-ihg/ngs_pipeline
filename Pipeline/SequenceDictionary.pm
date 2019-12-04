package SequenceDictionary;
use strict;

####################################################################################
# Thomas Wieland, 12.2011 														   #
# This file represents a class to handle a SAM sequence dictionary. The main       #
# purpose at the moment is to determine the order of chromosomes within a          #
# reference genome.															   	   #
####################################################################################



####################################################################################
# Constructor; requires a file that contains a SAM/BAM header including the        #
# sequence dictionary entries (@SQ). It doesn't have to be an actual file but can  #
# be a piped samtools call, like 'samtools view -H file.bam |'                     #
####################################################################################
sub new {
    my $proto = shift;
	my $class = ref($proto) || $proto;
	my $self  = {};
	bless ($self, $class);
	$self->{FILE} = shift;
	$self->readDictionary();
	
	return $self;
}


####################################################################################
# Getter and setter methods                                                        #
####################################################################################
sub file {
	my $self = shift;
	if (@_) { $self->{FILE} = shift }
	return $self->{FILE};
}

sub readDictionary {
	my $self = shift;
	my $filename = $self->{FILE}; #generate a sequence dictionary either from a .bed or a .bam file
	my $isBAM    = 0;
	if($filename =~ /\.bam$/){
		$filename = "samtools view -H $filename |";
		$isBAM = 1;
	}
	
	open IN, $filename or print "Can't open $filename!";
	my $counter = 0;
	$self->{CHROMS} = {};
	if($isBAM){
		#parse BAM file
		while(my $line = <IN>){
			next if !($line =~ /^\@SQ/);
			if($line =~ /SN:(.*)\s+LN:(\d*)/){
				$counter++;
				$self->{CHROMS}->{$1}->{LENGTH} = $2;		# save length of chromosome
				$self->{CHROMS}->{$1}->{ORDER}  = $counter;
				
			}		
		}
	}else{
		#parse FAI or BED file
		while(my $line = <IN>){
			chomp $line;
			my ($chr,@rest) = split("\t",$line);
			unless($self->{CHROMS}->{$chr}){
				$counter++;
				if($filename =~ /\.fai$/){
					$self->{CHROMS}->{$chr}->{LENGTH} = $rest[1];		# save length of chromosome
				}else{
					$self->{CHROMS}->{$chr}->{LENGTH} = 0;		# save length of chromosome
				}
				
				$self->{CHROMS}->{$chr}->{ORDER}  = $counter;
			}
		}
	}
}


####################################################################################
# compareChrs: compares two chromosomes according to their order in the sequence   #
# dictionary. returns:															   #
# <0 if $chr1 <  $chr2															   #
# 0  if $chr1 == $chr2                                                             #
# >0 if $chr1 >  $chr2															   #
####################################################################################
sub compareChrs {
	my $self = shift;
	my $chr1 = shift;
	my $chr2 = shift;
	return ($self->{CHROMS}->{$chr1}->{ORDER} - $self->{CHROMS}->{$chr2}->{ORDER});
}

1;