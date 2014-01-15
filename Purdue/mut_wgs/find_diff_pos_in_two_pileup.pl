#!/usr/bin/perl -w

use strict;

my $debug = 0;

my $usage = "$0 <pileup_1> <pileup_2> <output_1> <output_2> <same>";
die $usage unless(@ARGV == 5);

my ($input_1, $input_2, $output_1, $output_2, $same_out) = @ARGV[0..4];

open (IN1,$input_1)
	or die "cannot open $input_1:$!";
	
open (IN2,$input_2)
	or die "cannot open $input_2:$!";
	
open (OUT1,">$output_1")
	or die "cannot open $output_1:$!";
	
open (OUT2,">$output_2")
	or die "cannot open $output_2:$!";
	
open (SAME,">$same_out")
	or die "cannot open $same_out:$!";
	
my %snps;
	
while(<IN1>){
	chomp;
	my @pts = split "\t";
	$snps{$pts[0]}->{$pts[1]} = [@pts[0..4]];
}

close(IN1);

while(<IN2>){
	chomp;
	my @pts = split "\t";
	if (defined $snps{$pts[0]}->{$pts[1]}){
		my $dep1 = ${$snps{$pts[0]}->{$pts[1]}}[3];
		my $seq1 = ${$snps{$pts[0]}->{$pts[1]}}[4];
		print SAME join("\t", (@pts[0..4], $dep1,$seq1)), "\n";
		delete $snps{$pts[0]}->{$pts[1]};
	}
	else{
		print OUT2 join("\t", @pts[0..4]), "\n";
	}
}

foreach my $chr (sort keys %snps){
	foreach my $pos (sort {$a <=> $b} keys %{$snps{$chr}}){
		print OUT1 join("\t", @{$snps{$chr}->{$pos}}), "\n";
	}
}