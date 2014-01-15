#!/usr/bin/perl -w

# calculate the difference between mutants and WT
# note only Col0 #2 is used

use strict;

sub trimmean
{
	my $arr_ref = $_[0];
	my $length = scalar (@$arr_ref);
	
	if ($length < 2)
	{
		die "@$arr_ref";	
	}
	
	my $sum = 0;
	for(my $i = 1; $i < $length-1 ; $i++)
	{
		$sum += $arr_ref->[$i];	
	}	
	my $trimmean = $sum/($length - 2);
	return $trimmean;
}

my @a = (4,2,1,6,7,3);
my @b = sort @a;
my $t = trimmean(\@b);

print "$t\n";

print STDERR "\a";
exit;