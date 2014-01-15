#!/usr/bin/perl -w
use strict;

my @nums;
for my $i (0..7){
	for my $j(0..7){
		$nums[$i, $j] = $i * $j ;
	}
}

format 