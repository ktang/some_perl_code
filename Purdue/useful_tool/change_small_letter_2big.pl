#!/usr/bin/perl -w

use strict;
use Bio::SeqIO;


my $usage = "$0 <infile> \n";
die $usage unless (@ARGV == 1);

my $infile = shift or die "input";

open(IN, $infile);
while(<IN>){
	if(/^>/){
		print;
	}
	else{
		print uc;
	}
}
exit;