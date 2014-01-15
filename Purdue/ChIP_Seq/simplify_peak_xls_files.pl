#!/usr/bin/perl -w

#simplify_peak_xls_files.pl
#input MACS2 NAME_peaks.xls file
# output a txt file with no heading # and no Pt and Mt

use utf8;#可以吗？
use strict;
use File::Spec;

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $genome_size = "1.2e8";
my $keep_dup = "all";

my $usage = "\n$0 \n\n <input> <output> \n\n";
die $usage unless(@ARGV == 2);

my $input = shift or die;
my $output = shift or die;

die unless(-e $input);
die if(-e $output);

open(IN, $input) or die;
open(OUT, ">>$output" ) or die;

while (<IN>) {
	next if(/^#/ or /^$/);
	chomp;
	my @a = split "\t";
	print OUT join("\t", @a), "\n" if($a[0] =~ /\d/ or $a[0] eq "chr");
}



close IN;
close OUT;
exit;