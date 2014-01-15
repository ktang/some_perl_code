#!/usr/bin/perl -w

use strict;

# purpose:
# substract the potential_TAIR_ref_error_and_heterozygous_sites_v0.1
# from 3998loci_mpileup_in_6MAP20bam_and_4muts_with_head.txt
my %h;

my $usage = "$0 \n <long_list> <heterozygous_sites_list> STDOUT\n\n";

die $usage unless ( @ARGV == 2);
my $long_list = shift or die;
my $heterozygous_sites_list = shift or die;

open (IN, $heterozygous_sites_list) or die;
my $h = <IN>;
while(<IN>){
	chomp;
	my @a = split "\t";
	$h{$a[0]}->{$a[1]} = 1;
}
close IN;

open(IN, $long_list) or die;
$h = <IN>;
print $h;

while (<IN>){
#	chomp;
#	next if (/^#/);
	my @a = split "\t";
	my ($chr, $pos) = @a[0..1];
	unless( defined $h{$chr}->{$pos}){
		print $_;
	}
}

