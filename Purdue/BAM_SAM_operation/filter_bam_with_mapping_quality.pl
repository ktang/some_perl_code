#!/usr/bin/perl -w
use strict;
use File::Spec;

my $debug = 0;

#my $min_quality = 20;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <input_bam> <output> <min_quality>\n\n";
die $usage unless(@ARGV == 3);

my $input = shift or die;
my $output = shift or die;
my $min_quality = shift or die;
die unless (-e $input) ;
die if (-e $output);

if($debug){
	print STDERR "input: $input\n";
	print STDERR "output: $output\n";
	print STDERR "min_quality: $min_quality\n\n";
	exit;
}

open( IN, "samtools view -h $input |") or die;
open( OUT,  "| samtools view -b -S -  > $output") or die;



while(<IN>){
	if(/^@/){
		print OUT $_;
	}else{
		#chomp;
		my @a = split "\t";
		print OUT $_ if($a[4] >= $min_quality);
	}
}



close IN;
close OUT;


exit;
