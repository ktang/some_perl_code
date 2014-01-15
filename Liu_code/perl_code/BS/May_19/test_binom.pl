#!/usr/bin/perl -w
# test binomial function
use strict;
use Math::CDF qw(:all);

my $n = 10;
my $rate = 0.03;
if(@ARGV > 0){
	$n = $ARGV[0];
}
if(@ARGV > 1){
	$rate = $ARGV[1];
}
print "n=$n, rate=$rate\n";
foreach my $i(0..$n){
	my $p = pbinom($i,$n, $rate);
	my $ip = $p;
	if($i > 0){
		$ip = $p - pbinom($i-1, $n, $rate);
	}
    print "i=$i,accum. p=$p, p=$ip\n";
}
