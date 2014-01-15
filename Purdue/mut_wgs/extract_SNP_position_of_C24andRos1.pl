#!/usr/bin/perl -w

use strict;

my $debug = 0;

my $usage = "$0 <pileup_1> <pileup_2> <output>";
die $usage unless(@ARGV == 3);

my ($input_1, $input_2, $output) = @ARGV[0..2];

open (IN1,$input_1)
	or die "cannot open $input_1:$!";
	
open (IN2,$input_2)
	or die "cannot open $input_2:$!";
	
open (OUT,">$output")
	or die "cannot open $output:$!";
	
my %snps;
	
while(<IN1>){
	chomp;
	my @pts = split "\t";
	$snps{$pts[0]}->{$pts[1]} = 1;
}

close(IN1);

while(<IN2>){
	chomp;
	my @pts = split "\t";
	$snps{$pts[0]}->{$pts[1]} = 1;
}

close(IN2);

foreach my $chr (sort keys %snps){
	foreach my $pos (sort {$a <=> $b} keys %{$snps{$chr}}){
		print OUT join("\t", ($chr,$pos) ), "\n";
	}
}