#!/usr/bin/perl -w

use strict;

my $usage = "$0 <output_file> <region_num> <region_len>";
die $usage unless (@ARGV == 3);

my $output_file = $ARGV[0];
my $region_num = $ARGV[1];
my $len = $ARGV[2];

`./random_generate_segment_Ath_v1.pl $region_num $len > random`;
#bed file

#`./formated_sort_bed.sh random  random_sorted `;
system("sort -k1,1 -k2,2n random > random_sorted");

`rm -f random`;

`./verify_regions_with_Ecker_data_v6.pl random_sorted > gff`;

`rm -f random_sorted`;

`sort -rnk 9 gff > $output_file `;

`rm -f gff`;

system("head -n50 $output_file | tail -n1");

exit;