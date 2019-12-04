#!/usr/bin/perl

use Getopt::Long qw(:config no_ignore_case);
use Bio::DB::Sam;

my $ref = "PATHTO/hg19.fa";
GetOptions(
"i=s" => \$bedfile,
"r=s" => \$ref,
"b=s" => \$bamfile
);


my $sam = Bio::DB::Sam->new(									#needs a bam file --> dummy.bam must be in the current path
			-fasta => $ref,	
			-bam   => $bamfile
		);

#open input/output files, if specified
if($bedfile ne ""){
	close STDIN;
	$bedfile =~ s/(.*\.gz)\s*$/gzip -dc < $1|/;
	
	open(STDIN,$bedfile) or die "Can't open $bedfile!\n";
	
	
}


while(<STDIN>){
	chomp;
	my ($chrom,$start,$end,@rest) = split();
	my $seq = uc $sam->seq($chrom,$start,$end);
	print ">$chrom:$start-$end\n$seq\n";
}
