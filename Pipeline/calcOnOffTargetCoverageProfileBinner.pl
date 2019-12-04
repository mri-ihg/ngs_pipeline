#!/usr/bin/perl

# Cumulative histogram / normalised to total

#TODO: Improve it with custom parameters as a general tool

my @hist;

my $total=0;

for ( my $iLoop = 0; $iLoop<=1000; $iLoop++ )
{
	$hist[$iLoop]=0;
}

while(<>)
{
	chomp;
	$hist[ ($_<1000) ? $_ : 1000 ]++;
	$total++;
}

for ( my $iLoop=1000; $iLoop>=1; $iLoop-- )
{
	if ( $iLoop>=1 )
	{
		$hist[$iLoop-1]+=$hist[$iLoop];
	}
}

for ( my $iLoop=0; $iLoop<=200; $iLoop++ )
{
	print $iLoop." ".($hist[$iLoop]/$total)."\n";
}


=head1 NAME

calcOnOffTargetCoverageProfileBinner.pl

=head1 SYNOPSIS

 CoverageValues (topped to 100) | calcOnOffTargetCoverageProfileBinner.pl

=head1 DESCRIPTION

Bins from 0 to 200 for each bin returns the fraction of entries >= bin value.

=head1 AUTHOR

Riccardo Berutti

=cut
