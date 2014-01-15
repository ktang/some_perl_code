#!/usr/bin/perl -w

use strict;
my $usage = "$0 <input> <output>";

die $usage unless (@ARGV == 2);
my $C24_file = "/Users/tang58/deep_seq_analysis/SNP_pos_in_C24andRos1.txt";

my ($input, $output) = @ARGV[0..1];

die unless (-e $input and !(-e $output));

open (C24, $C24_file) or die "cannot open C24:$!";

open (IN, $input) or die "cannot open $input:$!";

open (OUT, ">$output") or die "cannot open $output:$!";

my %SNPs;

while(<C24>){
	chomp;
	my @a = split "\t";
	my $chr = lc $a[0];
	my $pos = $a[1];
	$SNPs{$chr}->[$pos] = 1;
}

while(<IN>){
	next if /^\#/;
	chomp;
	my @a = split "\t";
	my $chr = lc $a[0];
	my $pos = $a[3];
	print OUT $_,"\n" unless (defined $SNPs{$chr}->[$pos]);
}

close(IN);
close(OUT);

exit;