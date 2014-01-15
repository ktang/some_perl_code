#!/usr/bin/perl -w

=head2
this is version 2

the change is the output is formmated to 10 digits.
############################

input two intergers.
one indicate how many regions to
generate,
the other indicate the length of 
the region.

out put bed format in stdout.

Notice: the beginning is 0;

average length = 673
=cut

use strict;

my $usage = "$0 <num_regions> <region_len>";

die $usage unless (@ARGV == 2);
my $num_regions = $ARGV[0];
my $region_len  = $ARGV[1];

my ($len_chr1,$len_chr2,$len_chr3,$len_chr4,$len_chr5)=
(30427671,19698289,23459830,18585056,26975502);

my $total_chr = $len_chr1  + $len_chr2  + $len_chr3  + $len_chr4  + $len_chr5; 

my $point_2 = $len_chr1  + $len_chr2 ;
my $point_3 = $len_chr1  + $len_chr2  + $len_chr3 ;
my $point_4 = $len_chr1  + $len_chr2  + $len_chr3  + $len_chr4;

for (my $i = 1; $i <= $num_regions; $i++)
{
	my $random_chr = 	int(rand($total_chr ));
	
	if ($random_chr < $len_chr1)
	{
		my $start = $random_chr + 1;
		my $end = $start + $region_len - 1;
		if ($end > $len_chr1) {$end = $len_chr1}
		
		#print "chr1\t$start\t$end\n";
		printf "chr1\t%0.10u\t%0.10u\n",$start,$end;
		
	}
	elsif ( $random_chr < $point_2)
	{
		my $start = $random_chr - $len_chr1  + 1;
		my $end = $start + $region_len - 1;
		if ($end > $len_chr2) {$end = $len_chr2}
		
		printf "chr2\t%0.10u\t%0.10u\n",$start,$end;
	}
	elsif ( $random_chr < $point_3)
	{
		my $start = $random_chr - $point_2 +1 ;
		my $end = $start + $region_len - 1;
		if ($end > $len_chr3) {$end = $len_chr3}
		
		printf "chr3\t%0.10u\t%0.10u\n",$start,$end;

	}
	elsif ( $random_chr < $point_4)
	{
		my $start = $random_chr - $point_3 +1;
		my $end = $start + $region_len - 1;
		if ($end > $len_chr4) {$end = $len_chr4}
		
		printf "chr4\t%0.10u\t%0.10u\n",$start,$end;

	}
	elsif ( $random_chr < $total_chr)
	{
		my $start = $random_chr - $point_4 +1;
		my $end = $start + $region_len - 1;
		if ($end > $len_chr5) {$end = $len_chr5}
		
		printf "chr5\t%0.10u\t%0.10u\n",$start,$end;

	}
	else
	{
		die "$random_chr >= $total_chr:$!";	
	}
	
}
exit;
=head2 
my $start = $random_chr - $point_ +1;
		my $end = $start + $region_len - 1;
		if ($end > $len_chr) {$end = $len_chr}
		
		print "chr		\t$start\t$end\n";

=cut