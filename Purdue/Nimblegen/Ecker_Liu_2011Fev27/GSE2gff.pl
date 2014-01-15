#!/usr/bin/perl -w

use strict;

while(<>){
	my @a = split "\t";
	my $chr = "Chr".$a[0];
	print join("\t",($chr,".",".",$a[2],$a[3],".",".",".",".")) , "\n";
}