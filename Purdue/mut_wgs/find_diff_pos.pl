#!/usr/bin/perl -w
# already has /Users/tang58/deep_seq_analysis/SNP_pos_in_C24andRos1.txt
# filter ros1_mutant get pos diff from ros1.
use strict;

my $debug = 0;

my $usage = "$0 <input> STDOUT";
die $usage unless(@ARGV == 1);

my $pos_file = "/Users/tang58/deep_seq_analysis/SNP_pos_in_C24andRos1.txt";

my $marker_file = "/Users/tang58/data/C24/Col0_SNP_in_C24.txt";

my $input = $ARGV[0];

open (SNP,$pos_file)
	or die "cannot open $pos_file:$!";
	
open (IN,$input)
	or die "cannot open $input:$!";

open(MA, $marker_file)
	or die "cannot open $marker_file:$!";	
	
my %snps;
	
while(<SNP>){
	chomp;
	my @pts = split "\t";
	$snps{$pts[0]}->{$pts[1]} = 1;
}

while(<MA>){
	chomp;
	my @a = split "\t";
	my $chr = "Chr".$a[1];
	$snps{$chr}->{$a[2]} = 1;
}

close(SNP);

while(<IN>){
	chomp;
	my @pts = split "\t";
	if (defined $snps{$pts[0]}->{$pts[1]}){
	}
	else{
		print join("\t", @pts[0..4]), "\n";
	}
}

close (IN);