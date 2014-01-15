#!/usr/bin/perl
use warnings;
use strict;

my $usage = "<input> <output>";

die $usage unless(@ARGV == 2);

open (SORTED,"/Users/kaitang/Desktop/Nimblegen/probe_files/version2_Aug12/sorted_position_ID_v2")
	or die "cannot open sorted probe file:$!";


open (IN, "<$ARGV[0]")
	or die "cannot open input file:$!";
	
open (OUT ,">$ARGV[1]")
	or die "cannot open output file:$!";

my $line;
my @pts;
my %vals;


while($line = <IN>)
{
	chomp $line;
	@pts = split "\t", $line;
    $vals{$pts[2]} = $pts[5];
}

while ($line = <SORTED>)
{
	chomp $line;
	@pts = split "\t", $line;
	print OUT "$pts[0]\t.\t.\t$pts[1]\t$pts[2]\t$vals{$pts[3]}\t$pts[3]\t.\t.\n";
}

close(IN);
close(OUT);
print STDERR "\a";
exit;