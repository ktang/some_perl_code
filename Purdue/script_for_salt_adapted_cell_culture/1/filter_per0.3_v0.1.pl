#!/usr/bin/perl -w

use strict;

my $usage = "\n$0 <input> STDOUT\n\n";
die $usage unless (@ARGV == 1);

my $input = shift or die;
my $col_num = 5;

die unless (-e $input);
open(IN, $input) or die;

#chr	pos	ref	percentage	dep	seq	codon_change	AA_change	gene	strand	func
#chr1	2044397	G	90	20	AAaA,AaAA,AaaAAAAAaa	GTG=>ATG	V=>M	AT1G06670	+	nuclear DEIH-boxhelicase


my $h = <IN>;
print $h;

my $cutoff = 30;

while (<IN>){
	my @a = split "\t";
	print if ($a[0] eq "chr");
	if($a[4] eq "MANNUAL"){
		print;
		print STDERR;
	}
	elsif($a[4] >= $cutoff){
		print;
	}
	
}
close IN;

exit;