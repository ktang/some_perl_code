#!/usr/bin/perl -w

use strict;

my %records;

#/Users/tang58/misc/Zhu_Xiaohong/Mar27_filter_MAPQ20/mpileup/4056loci_mpileup_in_6MAP20bam.txt

#/Users/tang58/misc/Zhu_Xiaohong/Mar27_filter_MAPQ20/mpileup/4056loci_mpileup_in_4mutants.txt
my $usage = "$0 \n <input_mpileup_in_cell_culture> <input_mpileup_in_mut> STDOUT \n\n";
print STDERR "files should have head\n\n";

die $usage unless (@ARGV == 2);

my $file_cell_culture = shift or die;
my $file_mut	      = shift or die;

open(CEL, $file_cell_culture ) or die;
open(MUT, $file_mut ) or die;

my $h_cell = <CEL>;
chomp $h_cell;
my $h_mut = <MUT>;
chomp $h_mut;
my @a_h_mut = split "\t", $h_mut;

print join("\t", ($h_cell, @a_h_mut[3..$#a_h_mut])), "\n";
my ($l_cell, $l_mut);

while ($l_cell = <CEL>, $l_mut = <MUT>){
#	chomp;
	chomp $l_cell;
	chomp $l_mut;
	my @a_cell = split "\t", $l_cell;
	my @a_mut  = split "\t", $l_mut;
	if( $a_mut [1] != $a_cell[1] or   $a_mut [2] ne $a_cell[2]  ){
		print STDERR "wrong:\n";
		print STDERR $l_cell, "\n";
		print STDERR $l_mut, "\n\n";
	}
	if($a_mut[2] !~ /[ACTG]/ ){
		print STDERR join("\t", (@a_cell, @a_mut[3..$#a_mut] )), "\n";
	}
	print join("\t", (@a_cell, @a_mut[3..$#a_mut] )), "\n";
}

close CEL;
close MUT;