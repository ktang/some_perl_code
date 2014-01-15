#!/usr/bin/perl -w
use strict;
# chr1    3       .       CTAAACCCTAAACCCTAAACCCTAAACC    CTAAACCCTAAACCCTAAACCCTAAACCCTAAACC     34.2    .       INDEL;DP=3;VDB=0.0060;AF1=1;AC1=2;DP4=0,0,2,0;MQ=40;FQ=-40.5    GT:PL:GQ        1/1:73,6,0:10
# chr1    59214   .       C       A       7.8     .       DP=14;VDB=0.0216;AF1=0.5;AC1=1;DP4=2,2,0,4;MQ=41;FQ=10.4;PV4=0.43,8.8e-05,0,1   GT:PL:GQ        0/1:37,0,113:39
# 0       1       2        3       4       5      6       7

my %records;
# input and output is std
# cat *vcf | perl $0 > output.txt

while (<>){
#	chomp;
	next if (/^#/);
	my @a = split "\t";
	my ($chr, $pos) = @a[0..1];
	$records{$chr}->{$pos} = 1;
}

foreach my $chr (sort keys %records){
	foreach my $pos(sort {$a <=> $b} keys %{$records{$chr}}){
		print join("\t", ($chr, $pos )), "\n";
	}
}