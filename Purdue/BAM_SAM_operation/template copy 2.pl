#!/usr/bin/perl -w
use strict;
use File::Spec;

my $debug = 1;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <input_bam> <outdir> <outpre>\n\n";
die $usage unless(@ARGV == 3);

my $input = shift or die;
my $outdir = shift or die;
my $outpre = shift or die;

my $output;

die unless (-e $input) ;
die unless (-d $outdir);

open( IN, "samtools view -h $input |") or die;
open( OUT,  "| samtools view -b -S -  > $output") or die;



while(<IN>){
	if(/^@/){
		print OUT $_;
		
	}else{
		chomp;
		my @temp = split "\t";
	}
}



close IN;
close OUT;


exit;
