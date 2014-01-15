#!/usr/bin/perl -w

#step1 output like
# \#GSM701923_H3K4me2_bowtie_v3m1p4_merged.bam     GSM701924_H3K4me3_bowtie_v3m1p4_merged.bam
# >1:132546-132588
#1       130543-130592   2       22      54      0       0       0       17      0       56      14      170
#2       130593-130642   31      11      45      0       0       0       28      0       43      33      39


#step2.1_sum_of_corresopning_bin_divid_first.pl
# in step 2.0: I add number of each bin first and then divided by input number
# this time(2.1) I divide by input first and then calculate the mean to see wheatehr
# there will be difference in fig

use strict;
use File::Spec;

my $debug = 1;

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <input> <output>\n\n";
die $usage unless(@ARGV == 2);

my $input = shift or die;
my $output = shift or die;

die unless (-e $input);
die if( -e $output);

open(IN, $input) or die "cannot open $input: $!";

die if(-e $output);
open(OUT, ">$output") or die "cannot open $output: $!";

my %sum;

while (<IN>) {
	if (/^#/) {
		s/#//;
		print OUT join("\t", ("order", $_));
		next;
	}
	if (/^>/) {
		next;#code
	}
	chomp;
	my @a = split "\t";
	die $_ unless ($a[1] =~ /(\d+)-(\d+)/);
	
	#for my $i(2..$#a){
		#$sum{$a[0]}->[$i-2] += $a[$i];
	for my $i(2..( $#a - 1 )){
		$sum{$a[0]}->[$i-2] += $a[$i] / ( $a[-1] + 1 );
	}
}

foreach my $order (sort {$a<=>$b} keys %sum){
	print OUT join("\t", ($order, @{$sum{$order}})), "\n";
}


close(IN);
close(OUT);

exit;

sub round {
    my($number) = shift;
    #return int($number + .5);
    return int($number + .5 * ($number <=> 0)); # take care of negative numbers too

}

sub simple_chr{
	my ($chr) = @_;
	if( $chr =~ /chr/i){
		$chr =~  s/chr//i;
	}
	if($chr eq "M" ){
		$chr = "Mt";
	}elsif( $chr eq "C"){
		$chr = "Pt";
	}
	return $chr;
}