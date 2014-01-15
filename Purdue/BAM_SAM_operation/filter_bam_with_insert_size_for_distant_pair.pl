#!/usr/bin/perl -w
use strict;
use File::Spec;

my $debug = 0;

#my $min_quality = 20;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <input_bam> <output_sam> <insert_size_cutoff>\n\n";
die $usage unless(@ARGV == 3);

my $input = shift or die;
my $output = shift or die;
#my $min_quality = shift or die;
my $insert_size_cutoff = shift or die;
die unless (-e $input) ;
die if (-e $output);

if($debug){
	print STDERR "input: $input\n";
	print STDERR "output: $output\n";
	#print STDERR "min_quality: $min_quality\n\n";
	print STDERR "insert_size_cutoff: $insert_size_cutoff\n\n";
	exit;
}

open( IN, "samtools view -h $input |") or die;
open( OUT,  "| samtools view -b -S -  > $output") or die;

#HWI-ST531R:161:D12G0ACXX:6:1103:17302:77589     163     chr1    1       42      101M    =       181     281     CCCTAAACCCTAAACCCTAAACCCTAAACCTCTGAATCCTTAATCCCTAAATCCCTAAATCTTTAAATCCTACATCCATGAATCCCTAAATACCTAATTCC   CCCFFFFFHHGHHJJJGGIJCHHGAHI<FEEF<<FF>HEGIJIIIJIDFHIHHIEDGIGCCGHJJJGDGEACGAEEEEHBECDEF@CCC@CACDCCDD>CA   AS:i:0  XN:i:0  XM:i:0  XO:i:0  XG:i:0  NM:i:0  MD:Z:101        YS:i:0  YT:Z:CP
#					0	 1	 2	 3	 4	 5	 6	 7	 8

while(<IN>){
	if(/^@/){
		print OUT $_;
	}else{
		#chomp;
		my @a = split "\t";
		#print OUT $_ if($a[4] >= $min_quality);
		if($a[6] ne "=" or $a[8] == 0 or abs($a[8]) > $insert_size_cutoff){
			print OUT $_;
		}
	}
}
close IN;
close OUT;
exit;
