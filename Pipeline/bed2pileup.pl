#!/usr/bin/perl
use strict;
use warnings;
use diagnostics;
use Getopt::Long;

while(<STDIN>){
	if($_ =~ /track/){print $_; next;}
	chomp;
	my ($chr,$start,$end,@rest) = split();
	for(my $i=$start; $i <= $end; $i++){
		print "$chr\t$i\n";
	}
}