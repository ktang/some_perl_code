#!/usr/bin/perl -w

use strict;
my $usage = "$0 <infq> <out_fq>";

die $usage unless (@ARGV == 2);

my ($input, $output) = @ARGV[0..1];

open (IN, $input)
	or die "cannot open input file $input:$!";

open (OUT,">$output")
	or die "cannot open output file $output:$!";
my $no = 0;

my $line;

while ($line = <IN>){
	$no++;
	if ($no % 4 == 2){
		$line =~ tr/./N/;
		print OUT $line;
	}
	else{
		print OUT $line;
	}	
	
}

close(IN);
close(OUT);
exit;
