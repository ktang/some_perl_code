#!/usr/bin/perl -w
use strict;
use File::Spec;
my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <input_bam> <output_bam>\n\n";
die $usage unless(@ARGV == 2);

my $input = shift or die;
my $output = shift or die;

open( IN, "samtools view -h $input |") or die;
open(OUT,  "| samtools view -b -S -  > $output") or die;


while(<IN>){
	if(/^@/){
		print OUT $_;
		#print  $_;
	}else{
		chomp;
		my @a = split "\t";
		if($a[3] == 16539980){
			print OUT $_, "\n";
			#print $_, "\n";
		}
	}
	
}

close IN;
close OUT;