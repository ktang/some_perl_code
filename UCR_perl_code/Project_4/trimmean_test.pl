#!/usr/bin/perl -w

use strict;

sub trimmean
{
	my $arr_ref = $_[0];
	my $length = scalar (@$arr_ref);
	
	if ($length < 3)
	{
		print STDERR join "\t",@{$arr_ref};
		die "less than 3,@$arr_ref";	
	}
	
	my $sum = 0;
	for(my $i = 1; $i < $length-1 ; $i++)
	{
		$sum += $arr_ref->[$i];	
	}	
	my $trimmean = $sum/($length - 2);
	return $trimmean;
}

my @a = qw(-4 5 9);

print join("\t",@a),"\n";

my $b = trimmean(\@a);

print "$b\n";


exit;

