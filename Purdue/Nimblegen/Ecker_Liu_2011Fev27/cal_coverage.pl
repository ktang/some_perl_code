#!/usr/bin/perl -w

use strict;

my $usage = "$0 <NCBI_coverage_file> <input_bed_file>";
#format BED
#Chr	Start	End

die $usage unless(@ARGV == 2);

my $ncbi = $ARGV[0];
my $input = $ARGV[1];

open (NCBI, $ncbi)
	or die "cannot open $ncbi file:$!";
	
open (IN, $input)
	or die "cannot open $input file:$!";

my @ncbis = <NCBI>;
my @beds = <IN>;

close (NCBI);
close (IN);	
	
my ($i_n, $i_b) = (0, 0);
my ($total_cover, $part_cover, $no_cover) = (0, 0, 0);

LOOP: while( $i_n <= $#ncbis and $i_b <= $#beds){
	my $thisN = $ncbis[$i_n];
	my $thisB = $beds[$i_b];
	chomp $thisN;
	chomp $thisB;
	my ($chr_n, $start_n, $end_n) = split "\t", $thisN;
	my ($chr_b, $start_b, $end_b) = split "\t", $thisB;
	$chr_n = lc $chr_n;
	$chr_b = lc $chr_b;
	
	if($chr_b eq $chr_n){
		if($start_b >= $start_n  and $end_b <= $end_b){
			$total_cover ++;
			$i_b++;
			goto LOOP;
		}
		if($end_b < $end_n){
			$no_cover++;
			$i_b++;
			goto LOOP;
		}
		if($end_b >= $start_n and $end_b <= $end_n and $start_b < $start_n){
			$part_cover++;
			$i_b++;
			goto LOOP;
		}
		if($start_b >= $start_n and $start_b <= $end_n and $end_b > $end_n){
			$part_cover++;
			$i_b++;
			$i_n++;
			goto LOOP;
		}
		
		if($start_b < $start_n and $end_b > $end_n){
			$part_cover++;
			$i_b++;
			$i_n++;
			goto LOOP;
		}
		
		
	}elsif($chr_b lt $chr_n){
		$no_cover++;
		$i_b++;
	}else{
		$i_n++;
	}
}

my $peak_no = $#beds + 1; 
print "total peak Num: $peak_no\n";
print "no cover: $no_cover\npartly: $part_cover\ncomplete: $total_cover\n";