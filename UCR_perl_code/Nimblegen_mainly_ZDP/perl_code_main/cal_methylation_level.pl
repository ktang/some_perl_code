#!/usr/bin/perl
use warnings;
use strict;

#set a window size, count mC number in 
#each window

my $usage = "<$0> <window_size> <input_file>";

die $usage unless(@ARGV == 2);

$ARGV[1] =~ /(\S+)\.gff/;

my $pre = $1;

my $out = $pre."_wsize_".$ARGV[0];

open (IN, "<$ARGV[1]")
	or die "cannot open input file:$!";
	
open (OUT ,">$out")
	or die "cannot open output file:$!";

my $line;
my @pts;
my %nums;

my $index;
my $wsize = $ARGV[0];

while($line = <IN>)
{
	chomp $line;
	if ($line =~ /^\#/)
	{next;}
	@pts = split "\t", $line;
	$index = int $pts[3]/$wsize;
	
	$nums{$pts[0]}->{$index}++;
}

foreach my $chr (sort keys %nums)
{
	foreach my $pos (sort {$a <=> $b} keys %{$nums{$chr}})
	{
		my $start = $pos * $wsize;
		my $end = ($pos + 1) * $wsize - 1;
		print OUT "$chr\t$start\t$end\t$nums{$chr}->{$pos}\n";
	}
}


close(IN);
close(OUT);
print"\a";
exit;