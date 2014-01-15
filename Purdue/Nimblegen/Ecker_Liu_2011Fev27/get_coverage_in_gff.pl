#!/usr/bin/perl -w

use strict;

#my $usage = "$0 <cutoff> STDOUT"; # gap between to regions if less than
#cutoff bp, consider as the same one;

#die $usage unless(@ARGV);

my $cutoff = 100; #$ARGV[0];

my($last_chr, $last_start, $last_end) = ("Chr0", -1, -1);

my($start, $end) = (0, 0);
my @regions = ();
while(<>){
	chomp;
	my @a = split "\t";
	if ($a[0] ne $last_chr){
		push @regions, [$last_chr, $start, $end];
		$start = $a[3];
		$end = $a[4];
		$last_start = $start;
		$last_end = $end;
		$last_chr = $a[0];
	}else{
		if ($a[3] < $start){
			die "a3 < start\n $a[3], $start\n";
		}
		if($a[3] - $last_end <= $cutoff){
			if($end < $a[4]){
				$end = $a[4];
				$last_end = $end;
			}
		}
		else{
			push @regions, [$a[0], $start, $end];
			$start = $a[3];
			$end = $a[4];
			$last_start = $start;
			$last_end = $end;
		}
	}
}
		
push @regions, [$last_chr, $start, $end];

#my $index = $#regions;
#print STDERR "index = $index";

for(my $i = 1; $i <= $#regions; $i++){
	my ($chr, $s, $e) = @{$regions[$i]};
	print join ("\t", ($chr, $s, $e) ), "\n";
}