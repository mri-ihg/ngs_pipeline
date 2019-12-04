#!/usr/bin/perl

#print "track name=exome38Mb description=\"Agilent 38Mb Exome Kit Target Regions\" useScore=1\n";

while(<STDIN>){
	chomp;
	my ($chr,$pos) = split(":");
	my ($pos1,$pos2) = split("-", $pos);
	print "$chr\t$pos1\t$pos2\n";
}