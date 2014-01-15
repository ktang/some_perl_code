#!/usr/bin/perl -w

use strict;
use Bio::SeqIO;

my $usage = '$0 <input> <output>';

die $usage unless(@ARGV ==2);

my $input = $ARGV[0];
my $output = $ARGV[1];
my $out = "";
#sanger # to illumina B
my $in = Bio::SeqIO->new( -format	=> "fastq",
						  -variant	=> 'sanger',
						  -file		=>  $input);
if (-e $output){
	die "$output exist";
}else{						  
	$out = Bio::SeqIO->new( -format	=>'fastq',
						   -variant	=>'illumina',
						   -file	=> ">$output");
}						   
while(my $seq = $in->next_seq){
	$out->write_seq($seq);
}