#!/usr/bin/perl
use warnings;
use strict;

my $usage = "<$0> <input>";

die $usage unless(@ARGV == 1);

my $out = $ARGV[0].".gff";

open (IN, "<$ARGV[0]")
	or die "cannot open input file:$!";
	
open (OUT ,">$out")
	or die "cannot open output file:$!";

my $line;
my @pts;

while($line = <IN>)
{
	chomp $line;
	@pts = split "\t", $line;
	my $ch = "chr".$pts[0];
	if ($pts[4] == 0) {$pts[4] = 0}
	if ($pts[5] == 0) {$pts[5] = 0}
	print OUT "$ch\t.\t.\t$pts[1]\t$pts[2]\t$pts[3]\t$pts[4]\t$pts[5]\t$pts[6]\n";	
}


close(IN);
close(OUT);
print STDERR "\a";
exit;