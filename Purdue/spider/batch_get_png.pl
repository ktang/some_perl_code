#!/usr/bin/perl -w
use strict;

my $debug = 0;

my $script = "/Users/tang58/try/Learn/spidering_Hacks/get_png_from_efp.pl";
die unless (-e $script);

my $usage = "$0 \n <inlist> \n\n";

die $usage unless (@ARGV == 1);

my $file = shift or die;

die unless (-e $file);

open(IN, $file) or die;

while (<IN>){
	chomp;
	my $cmd = "perl $script $_";
	print STDERR $cmd ,"\n\n";
	if(!$debug){
		`$cmd`;
	}
}
close(IN);