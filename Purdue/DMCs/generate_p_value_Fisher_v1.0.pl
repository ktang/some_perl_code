#! /usr/bin/perl -w
#		   WT	    mut
#		  word2   ~word2
#mC	 word1    n11      n12 | n1p
#No.T	~word1    n21      n22 | n2p
#          --------------
#		  np1      np2   npp
use strict;
use Text::NSP::Measures::2D::Fisher2::twotailed;

#v1.0 use Dr. Liu's output as input to prove the correctness
#of the code.

#v1.1 use Col-0 high WT

#v1.2 use C24 input

#my ($npp, $n1p,    $np1,      $n11);
#   total  mC_sum   WT_depth  WT_mC;

my $usage = "$0 <input> STDOUT";

die $usage unless (@ARGV == 1);
my $input = $ARGV[0];

#0chr	 1pos	2strand	 3type   4WT	 5WT_per	 6WT_met 7mu1	 8	 9	 10mu2	11	12	13p_value
#chr1    29      +       CHH     7       0.142857        0       2       0       0       3       0       0       1

open (IN, $input) or die "cannot open $input: $!";

while (<IN>){
    my ($npp, $n1p,    $np1,      $n11);
    chomp;
    my @a = split "\t";
    my $wt_depth = $a[4];
    my $mut_depth = $a[7] + $a[10];
    if($wt_depth >= 4 and $mut_depth >= 4){
	my $wt_mC = round($wt_depth * $a[5] );
	my $mut_mC = round($a[7] * $a[8]) + round($a[10] * $a[11]);
	$npp = $wt_depth + $mut_depth;
	$n1p = $wt_mC + $mut_mC;
	$np1 = $wt_depth;
	$n11 = $wt_mC;
	my  $left_value = calculateStatistic( n11=>$n11,
					      n1p=>$n1p,
					      np1=>$np1,
					      npp=>$npp);
	print join("\t", (@a[0..12], $left_value)), "\n";
    }
}

#my  $left_value = calculateStatistic( n11=>$n11,
#                                     n1p=>$n1p,
#                                      np1=>$np1,
#                                      npp=>$npp);

exit;

sub round{
    my($number) = shift;
    return int($number + .5 * ($number <=> 0)); # take care of negative numbers too
}
