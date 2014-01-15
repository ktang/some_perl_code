#!/usr/bin/perl -w

use strict;

my $usage = "\n$0 <input> <cutoff> <col_number> STDOUT\n\n";
die $usage unless (@ARGV == 3);

my $input = shift or die;
my $cutoff = shift or die;
my $col_num = shift or die;

my $index = $col_num - 1;

die unless (-e $input);
open(IN, $input) or die;

my $h = <IN>;
print $h;

#my $cutoff = 30;

while (<IN>){
	my @a = split "\t";
	print if ($a[0] eq "chr");
	
	if($a[$index] eq "NA"){
		next;
	}
	if($a[$index] eq "MANNUAL"){
		#print;
		print STDERR;
		next;

	}
	elsif($a[$index] >= $cutoff and $a[$index] > 10){
		print;
	}
	
}
close IN;

exit;