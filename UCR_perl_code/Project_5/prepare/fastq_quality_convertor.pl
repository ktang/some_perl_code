#!/usr/bin/perl -w
=head2
Because FASTX use B(66) as 2, but some illumina fastq use #(35) as 2.
chr(ord()+31); 
=cut
use strict;
my $usage = "$0 <infq> <out_fq>";

die $usage unless (@ARGV == 2);

my ($input, $output) = @ARGV[0..1];

open (IN, $input)
	or die "cannot open input file:$!";

open (OUT,">$output")
	or die "cannot open output file:$!";
my $no = 0;

my $line;

while ($line = <IN>){
	$no++;
	if ($no % 4 == 0){
		chomp $line;
		my @pts = split "",$line;
		for (my $i = 0; $i <= $#pts; $i++)
		{$pts[$i] = chr( ord( $pts[$i] ) + 31 );}
		$line = join "",@pts;
		print OUT $line,"\n";
	}
	else{
		print OUT $line;
	}	
	
}

close(IN);
close(OUT);
exit;
