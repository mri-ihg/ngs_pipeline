package BEDRecord;
use strict;
use Data::Dumper;
use IO::File;
use Storable qw(dclone);
use Tabix;

####################################################################################
# Thomas Wieland, 12.2011 														   #
# This file represents a class to handle single records, i.e. lines of a BED file  #
# or similar genomic regions. It provides methods to compare and find Records.     #
####################################################################################

####################################################################################
# Global variables   															   #
####################################################################################
my $globSD;


####################################################################################
# Constructor; requires definition of chromosome, start- and endposition           #
# These are also the minimum requirements of an entry in the BED format descr.     #
####################################################################################
sub new {
    my $proto = shift;
	my $class = ref($proto) || $proto;
	my $self  = {};
	bless ($self, $class);
	$self->{CHR}     = shift;
	$self->{STARTPOS}= shift;
	$self->{ENDPOS}  = shift;
	$self->{NAME}    = shift;
	$self->{OVERLAP} = undef;			#needed for overlap calculations
	
	return $self;
}


####################################################################################
# parseBED: parse a line from a BED file                                           #
####################################################################################
sub parseBED {
	my $proto = shift;
	my $line  = shift;
	chomp $line;
	my ($chr,$start,$end,$name,@rest) = split("\t",$line);
	my $self= $proto->new($chr,$start,$end);
	if($name){
		$name =~ s/\s/_/g;
		$self->name($name);
	}
	if(@rest){
		$self->rest(\@rest);
	}
	return $self;
}


####################################################################################
# readOverlappingSequences: returns the overlapping sequences from a file. This    #
# method is particulary useful if you want to get the overlapping sequences from   #
# $file for a SORTED list of regions, since it keeps the memory footprint low by   #
# just "scrolling" through $file without caching whole chromosomes. It also returns#
# the last read entry because otherwise it would get lost.						   #
####################################################################################
sub readOverlappingSequences {
	my $self       = shift;
	my $file       = shift;
	my $sd         = shift;
	my $readBuffer = shift;
	my @ret;
	my $currBed;
	#if the currently cached read is still to the right of $self --> no overlaps
	if($readBuffer){
		for(my $i = 0; $i < @$readBuffer; $i++){
			my $overlap = $self->calcOverlap(@$readBuffer[$i],$sd);
			if($overlap == -1){							#$currEntry is strictly to the right of $self --> leave it in readbuffer, following entries might overlap it
				return ($readBuffer,\@ret) if @ret;
				return ($readBuffer,undef);
			}
			if($overlap > 0){
				push(@ret,@$readBuffer[$i]);	#overlap detected, push it into @ret and leave it in readbuffer, because following entries might also overlap
			}else{
				shift @$readBuffer;		#remove entry from read buffer if its strictly to the left of $self --> will be strictly to the left of all following variants
				$i--;
			}
		}
	}
	
	
	
	
	while(my $line = <$file>){
		$currBed = BEDRecord->parseBED($line);
		
		my $overlap = $self->calcOverlap($currBed,$sd);
		push(@$readBuffer,$currBed) if $overlap != 0;	#add $currBed to the read buffer if it isn't strictly to the left of $self
		last if $overlap == -1;	#if the current read is already to the right of $self --> scrolled by
		if($overlap>0){
			#print "overlap >0 : $overlap\n";
			push(@ret,$currBed);	
		}
	}
	return ($readBuffer,\@ret) if @ret;
	return ($readBuffer,undef);
	
}

####################################################################################
# Getter and setter methods                                                        #
####################################################################################
sub chr {
	my $self = shift;
	if (@_) { $self->{CHR} = shift }
	return $self->{CHR};
}

sub startpos {
	my $self = shift;
	if (@_) { $self->{STARTPOS} = shift }
	return $self->{STARTPOS};
}

sub endpos {
	my $self = shift;
	if (@_) { $self->{ENDPOS} = shift }
	return $self->{ENDPOS};
}

sub name {
	my $self = shift;
	if (@_) { $self->{NAME} = shift }
	return $self->{NAME};
}

sub rest {
	my $self = shift;
	if (@_) { $self->{REST} = shift }
	return $self->{REST};
}

sub overlap {
	my $self = shift;
	if (@_) { $self->{OVERLAP} = shift }
	return $self->{OVERLAP};
}

sub toString {
	my $self = shift;
	my $ret = $self->chr()."\t".$self->startpos()."\t".$self->endpos();
	if($self->name()){
		$ret .= "\t".$self->name();
	}
	if($self->rest()) {
		my @rest = @{$self->rest()};
		if(@rest){
			foreach my $curr(@rest){
				if($curr){
					$ret .= "\t".$curr;
				}		
			}
		}
	}
	return $ret;
}

####################################################################################
# sortBySD: sort an array of BED entries by order of a sequence ditionary          #
####################################################################################
sub sortBySD {
	my $self         = shift;
	my $arrayPointer = shift;
	$globSD          = shift;
	
	my @ret = sort compareToSort @$arrayPointer;
	return \@ret;
}

####################################################################################
# compare function for sort: first chromosome according to sequence dictionary,    #
# then startpos, then endpos													   #
####################################################################################
sub compareToSort {
	return &compareTo($a,$b,$globSD);
}

sub compareTo {
	my $a  = shift;
	my $b  = shift;
	my $sd = shift;
	
	my $chrCompare;
	if($sd){
		$chrCompare = $sd->compareChrs($a->chr(),$b->chr());
	}else{
		$chrCompare = $a->chr() cmp $b->chr();
	}
	 
	return -1 if $chrCompare < 0;
	return  1 if $chrCompare > 0;
	return -1 if $a->startpos()<$b->startpos;
	return  1 if $a->startpos()>$b->startpos;
	return -1 if $a->endpos()<$b->endpos;
	return  1 if $a->endpos()>$b->endpos;
	return  0;
}

####################################################################################
# calcOverlap: calculate overlap of two BEDRecord objects. Needs a sequence 	   #
# dictionary object to calculate positions of chromosomes relative to their order  #
# in the sequence dictionary.													   #
# returns:                                                                         #
# -1 if $self lies to the left of $other                                           #
# 0  if $self lies to the right of $other										   #
# 0 < x <= 100 if $self and $other overlap, such that x gives the percentage of    #
#			   how many bases of $self overlap with $other						   #
####################################################################################
sub calcOverlap{
	my $self  = shift;
	my $other = shift;
	my $sd    = shift;
	
	#check if the regions are on the same chromosome
	my $chrCompare;
	if($sd){
	 	$chrCompare = $sd->compareChrs($self->chr(), $other->chr());
	}else{
		$chrCompare = $self->chr() cmp $other->chr();
	}
	if($chrCompare == 0){		#same chromosome
	
		my $selfStart  = $self->startpos();
		my $selfEnd    = $self->endpos();
		
		my $otherStart = $other->startpos();											
		my $otherEnd   = $other->endpos();
		
		$selfEnd-- if $self->endpos()>$self->startpos();		#-1 is needed because the last coordinate of a BED feature is actually not part of the feature (e.g. a SNPs coordinates are for instance chr1	5	6)
		$otherEnd-- if $other->endpos()>$other->startpos(); 	#EDIT: TW 23.03.2016: some SV callers produce startpos == endpos values for insertions --> don't subtract -1 from end
		
		# CASE 1:	
		# $self :		######
		# $other:				####
		return -1 if $selfEnd < $otherStart;
		
		# CASE 2:	
		# $self :				######
		# $other:		####
		return 0 if $selfStart > $otherEnd;
		
		if($selfEnd <= $otherEnd ){
			
			# CASE 3:	
			# $self :		 ######
			# $other:		########		
			return 100 if $selfStart >= $otherStart;
			
			# CASE 4:	
			# $self :		 ######
			# $other:		   ######		
			return (($selfEnd-$otherStart)/($selfEnd-$selfStart))*100;
		}
		
		# CASE 5:	
		# $self :		#######
		# $other:		 ####
		return  ((($otherEnd-$otherStart)/($selfEnd-$selfStart))*100) if $selfStart <= $otherStart;
		
		# CASE 6:	
		# $self :		  #######
		# $other:		 ####
		return  ((($otherEnd-$selfStart)/($selfEnd-$selfStart))*100);
		
		
		
		
	}elsif($chrCompare < 0){	#chromosome of $self comes first in sequence dictionary
		return -1;
	}
	return 0;
}


####################################################################################
# intersect: calculates the intersection of two BEDRecords. Returns a new          #
#   BEDRecord with the overlapping area if the records are overlapping and         #
#   nothing otherwise.                                                             #
####################################################################################
sub intersect{
	my $self  = shift;
	my $other = shift;
	my $sd    = shift;
	
	#check if the regions are on the same chromosome
	my $chrCompare;
	if($sd){
	 	$chrCompare = $sd->compareChrs($self->chr(), $other->chr());
	}else{
		$chrCompare = $self->chr() cmp $other->chr();
	}
	if($chrCompare == 0){		#same chromosome
	
		my $selfStart  = $self->startpos();
		my $selfEnd    = $self->endpos()-1;		#-1 is needed because the last coordinate of a BED feature is actually not part of the feature (e.g. a SNPs coordinates are for instance chr1	5	6)
		my $otherStart = $other->startpos();
		my $otherEnd   = $other->endpos()-1; 
	
		# CASE 1:	
		# $self :		######
		# $other:				####
		return if $selfEnd < $otherStart;
		
		# CASE 2:	
		# $self :				######
		# $other:		####
		return if $selfStart > $otherEnd;
		
		if($selfEnd <= $otherEnd ){
			
			# CASE 3:	
			# $self :		 ######
			# $other:		########		
			return $self if $selfStart >= $otherStart;
			
			# CASE 4:	
			# $self :		 ######
			# $other:		   ######		
			return BEDRecord->new($self->chr(),$otherStart,($selfEnd+1));
		}
		
		# CASE 5:	
		# $self :		#######
		# $other:		 ####
		return  $other if $selfStart <= $otherStart;
		
		# CASE 6:	
		# $self :		  #######
		# $other:		 ####
		return BEDRecord->new($self->chr(),$selfStart,($otherEnd+1));
	}
}


####################################################################################
# getOverlappingSequences: returns all overlapping sequences and the percentage	   #
# with which $self overlaps with them. BEDRecords to test for must provided in a   #
# SORTED array since this method uses binary searching. Returns an array of        #
# BEDRecords where the overlap of $self with each record is stored in 			   #
# record->{OVERLAP}																   #
####################################################################################
sub getOverlappingSequences {
	my $self         = shift;
	my $hash         = shift;
	my $seqDict      = shift;
	return () unless $hash->{$self->chr()};				#check if chromosome that is currently looked for is in the hash
	
	
	my @array		 = @{$hash->{$self->chr()}};			# get the array for the current chromosome out from the hash
	
	my $lastPos      = @array;
	my $firstPos     = 0;
	my $currPos		 = int($lastPos / 2);
	my $oldCurr      = -1;
	
	#do binary search
	while(1){
		
		my $overlap = $self->calcOverlap($array[$currPos],$seqDict);
		if($overlap<0){
			$lastPos  = $currPos;
		}elsif($overlap == 0){
			$firstPos = $currPos;
		
		}else{			#overlap found
		
			#it is possible that there are more than one overlapping BEDRecords in the array --> search nearby records to find overlaps
			
			while($currPos>0 && $self->calcOverlap($array[$currPos-1],$seqDict) > 0){		#go to the left as long as regions are overlapping	
				$currPos--;
			}
			
			my @ret;
			while($currPos<@array && ($overlap = $self->calcOverlap($array[$currPos],$seqDict))>0 ){
				$array[$currPos]->overlap($overlap);
				push(@ret, $array[$currPos]);
				$currPos++;
			}
			return @ret;
			
		}
		
		$oldCurr = $currPos;
		$currPos = $firstPos+int(($lastPos-$firstPos) / 2);
		return () if $oldCurr==$currPos;			#if the new calculated position is the same as the old --> no overlap found
	}
	
}



#######################################################################################
# getOverlappingSequencesTabix: gets the overlapping sequences from a tabix indexed   #
#   using Heng Li's Tabix module. This method returns an array of BEDRecords where    #
#   each BEDRecord element in the array contains the calculated overlap. If           #
#   if minOverlap is given the array will only include BEDRecords with a RECIPROCAL   #
#   overlap > $minOverlap.															  #
#   $self->overlap() will hold the sum of all overlaps (Can be >1 if entries in the   #
#   tabix file are overlapping!!!! Use sumOverlap if this is a problem!)    		  #
#	This method is especially useful to create annotations							  #  
#																					  #
#   NOTE: tabix can handle multiple file formats, but this methods assumes that the   #
#		  current file is in BED format!											  #
#######################################################################################
sub getOverlappingSequencesTabix {
	my $self         = shift;
	my $tabix        = shift;
	my $minOverlap   = shift;

	$self->overlap(0);
	my @ret;
	my $iter = $tabix->query($self->chr(), $self->startpos()-1,$self->endpos()+1);			# +/-1 because tabix doesn't handle that correctly
	eval {$tabix->read($iter)}; return @ret if $@;
	while (my $line = $tabix->read($iter)){	
		chomp $line;
		my @columns = split("\t",$line);
		my $curr = BEDRecord->new($columns[0],$columns[1],$columns[2]);
		$curr->name($columns[3]) if @columns > 3;
		$curr->rest(\@columns[4..(@columns-1)]) if @columns > 4;
		
		$curr->overlap($self->calcOverlap($curr));
		
		if( !defined($minOverlap) || ($curr->overlap() > $minOverlap && $curr->calcOverlap($self) > $minOverlap) ){		#check overlap if required
			$self->overlap($self->overlap() + $curr->overlap());
			push(@ret,$curr);
		}
	}
	
	return @ret;
}






####################################################################################
# getOverlappingIndices: returns all indices within the given array of BEDRecords  #
# with a RECIPROCAL overlap >= $minOverlap; $minOverlap must be >0 and <1		   #
####################################################################################
sub getOverlappingIndices {
	my $self         = shift;
	my $hash         = shift;
	my $minOverlap   = shift;
	my $seqDict      = shift;
	return () unless $hash->{$self->chr()};				#check if chromosome that is currently looked for is in the hash
	
	
	my @array		 = @{$hash->{$self->chr()}};			# get the array for the current chromosome out from the hash
	
	my $lastPos      = @array;
	my $firstPos     = 0;
	my $currPos		 = int($lastPos / 2);
	my $oldCurr      = -1;
	
	#do binary search
	while(1){
		
		my $overlap = $self->calcOverlap($array[$currPos],$seqDict);
		if($overlap<0){
			$lastPos  = $currPos;
		}elsif($overlap == 0){
			$firstPos = $currPos;
		
		}else{			#overlap found
		
			#so far we found only a single overlapping BEDRecord; it could be that this overlapping BEDRecord is e.g. spanning the whole chromosome
			#we need to find all elements that can overlap $self sufficiently
			#our array is sorted by the start position --> we go only so far to the left until we reach a start position that is so far away from $self
			#that it is impossible to get an reciprocal overlap > $minOverlap
			
			while($currPos>0 &&  (($self->endpos()-$self->startpos())/($self->endpos()-$array[$currPos-1]->startpos()))> $minOverlap){	
				$currPos--;
			}
			
			my @ret;
			while($currPos<@array && $array[$currPos]->startpos()<=$self->endpos() ){		#scroll through the array until the start of the current BEDRecord lies to the right of $self 
				
				
				#check if the reciptrocal overlap is > $minOverlap
				if($self->calcOverlap($array[$currPos],$seqDict) > ($minOverlap*100) && $array[$currPos]->calcOverlap($self,$seqDict) > ($minOverlap*100)){
					push(@ret, $currPos);
				}
				$currPos++;
			}
			return @ret;
			
		}
		
		$oldCurr = $currPos;
		$currPos = $firstPos+int(($lastPos-$firstPos) / 2);
		return () if $oldCurr==$currPos;			#if the new calculated position is the same as the old --> no overlap found
	}
	
}







####################################################################################
# sumOverlap: returns how many percent of $self are overlapped by regions in an    #
# array. The method colapses the single overlapping regions into distinct regions, #
# because otherwise some bases may be summed up more than once.                    #
####################################################################################
sub sumOverlap {
	my $self         = shift;
	my $arrayPointer = shift;
	my $seqDict      = shift;
	
	#print "neuertest: ".$self->toString()."\n";
	my @regions = $self->getOverlappingSequences($arrayPointer,$seqDict);
	
	#print Dumper(@regions);
	return 0 unless @regions;
	
	my @colapsedRegions;
	push(@colapsedRegions, $regions[0]);
	
	for(my $i=1;$i<@regions;$i++){
		if($colapsedRegions[-1]->calcOverlap($regions[$i],$seqDict)>0){ #regions are overlapping -> collapse
			$colapsedRegions[-1]->endpos($regions[$i]->endpos());
		}else{
			push(@colapsedRegions,$regions[$i]);
		}
	}
	
	my $hash;
	$hash->{$self->chr()} = \@colapsedRegions;							#build chr hash
	
	#print Dumper($hash);
	
	my @finalRegions = $self->getOverlappingSequences($hash,$seqDict);		#calculate overlaps with collapsed regions

	#sum up overlap percent
	my $sum = 0;
	
	foreach(@finalRegions){
		#print $self->toString()." to ".$_->toString()." = ".$_->overlap();
		$sum += $_->overlap();
	}
	
	return $sum;
}


####################################################################################
# sumSuperregion: sum the single overlaps of a super region (e.g. overlap for all  #
# the exons of a gene), i.e. all regions of an array. But since the single overlap #
# values are percentages of this region, the overall overlap must be calculated    #
# by calculating the sum of the overlapping bases                                  #
####################################################################################
sub sumSuperregion {
	my $self         = shift;
	my $arrayPointer = shift;
	my @regions = @$arrayPointer;
	
	my $allBases        = 0;
	my $overlappedBases = 0;
	
	foreach(@regions){
		$allBases        += ($_->endpos()-$_->startpos());
		$overlappedBases += ($_->endpos()-$_->startpos())*$_->overlap();
	}
	
	return $overlappedBases/$allBases;

}

####################################################################################
# maxIntersect: calculates the intersection of two BEDRecords. Returns a new       #
#   BEDRecord with the maximum overlapping area if the records are overlapping and #
#   nothing otherwise.                                                             #
####################################################################################
sub maxIntersect{
	my $self  = shift;
	my $other = shift;
	my $sd    = shift;
	
	#check if the regions are on the same chromosome
	my $chrCompare;
	if($sd){
	 	$chrCompare = $sd->compareChrs($self->chr(), $other->chr());
	}else{
		$chrCompare = $self->chr() cmp $other->chr();
	}
	if($chrCompare == 0){		#same chromosome
	
		my $selfStart  = $self->startpos();
		my $selfEnd    = $self->endpos();
		my $otherStart = $other->startpos();
		my $otherEnd   = $other->endpos(); 
	
		# CASE 1:	
		# $self :		######
		# $other:				####
		return if $selfEnd < $otherStart;
		
		# CASE 2:	
		# $self :				######
		# $other:		####
		return if $selfStart > $otherEnd;
		
		# SPECIAL CASE:
		# $self :		########
		# $other:		########
		return if ($selfStart == $otherStart && $selfEnd == $otherEnd);
		
				
		if($selfEnd <= $otherEnd ){
			
			# CASE 3:	
			# $self :		 ######
			# $other:		########		
			return $other if $selfStart >= $otherStart;
			
			# CASE 4:	
			# $self :		 ######
			# $other:		   ######		
			return BEDRecord->new($self->chr(),$selfStart,$otherEnd);
		}
		
		# CASE 5:	
		# $self :		#######
		# $other:		 ####
		return  $self if $selfStart <= $otherStart;
		
		# CASE 6:	
		# $self :		  #######
		# $other:		 ####
		return BEDRecord->new($self->chr(),$otherStart,$selfEnd);
	}
}

#####################################################################################
# collapseOverlappingArrayEntries: takes an array with BEDRecord objects and        #
#   compares all entries against each other and returns the max overlapping regions #
#####################################################################################
sub collapseOverlappingArrayEntries {
	my $arrayRef = shift;
	my @array = @$arrayRef;
	my @merged_array = ();
	my $overlap = 0;
	my $tmpRec;
	
	push(@merged_array, shift @array);
	foreach my $arrEntry (@array) {
		$overlap = 0;
		foreach my $mergedArrEntry (@merged_array) {
			$tmpRec = $arrEntry->maxIntersect($mergedArrEntry);
			if ($tmpRec) {
				$mergedArrEntry = $tmpRec;
				$overlap = 1;
			}
		}
		if (!$overlap) {
			push(@merged_array, $arrEntry);
		}
	}
	return (\@merged_array);
}


####################################################################################
# insertBEDRecord: inserts $self into the hash at the correct position             #
# insert position is returned													   #
####################################################################################
sub insertBEDRecord {
	my $self         = shift;
	my $hash         = shift;
	my $seqDict      = shift;
	if(!(defined $hash->{$self->chr()}) || @{$hash->{$self->chr()}} == 0){
		my @tmp = ($self);
		$hash->{$self->chr()} = \@tmp;
		return 0;
	}				#check if chromosome that is currently inserted is in the hash --> insert it if not
	
	my $array		 = $hash->{$self->chr()};			# get the array for the current chromosome out from the hash
	
	my $lastPos      = @{$array};
	my $firstPos     = 0;
	my $currPos		 = int($lastPos / 2);
	my $oldCurr      = -1;
	
	#do binary search in the array to find insert position
	while(1){
		
		$oldCurr = $currPos;

		my $compare = $self->compareTo($array->[$currPos],$seqDict);
		if($compare<0){
			$lastPos  = $currPos;
			$currPos  = $firstPos+int(($lastPos-$firstPos) / 2);
		}elsif($compare>0){
			$firstPos = $currPos;
			$currPos  = $firstPos+int(($lastPos-$firstPos) / 2);
		}
		
		
		#position to insert variant found
		if($oldCurr == $currPos){
			
			if($compare < 0){		#insert before current entry
				splice(@{$array},$currPos,0,$self);
				return $currPos;
			}else{					#insert after current entry
				splice(@{$array},($currPos+1),0,$self);
				return ($currPos+1);
			}
		}
	}
}

####################################################################################
# removeBEDRecords: removes a list of BEDRecords with given indices from an array  #
####################################################################################
sub removeBEDRecords {
	my $arrayRef = shift;
	my $indices  = shift;
	
	my @ind = @{$indices};
	for(my $i = (@ind-1);$i>=0; $i-- ){
		splice(@{$arrayRef},$ind[$i],1);			#start deleting from the back because obviously indices are different when deleting from the front
	}
}



1;
