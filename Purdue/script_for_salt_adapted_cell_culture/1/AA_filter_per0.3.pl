#!/usr/bin/perl -w

use strict;
# 0	 1	2	3		4						5	 6		7						8
#chr2    990     N       15      ^'A^'A^'A^'A^'A^'A^'A^'A^'A^'A^'A^'A^'A^'A^'A   444444444444444 17      AAAA^'A^'A^'A^'A^'A^'A^'A^-A^'A^'A^'A^'A^'A     33334444444444444

#my $usage = "\n$0 <input> <col_number_of_mpileup_seq_1based> STDOUT\n\n";
#die $usage unless (@ARGV == 2);
my $usage = "\n$0 <input> STDOUT\n\n";
die $usage unless (@ARGV == 1);

my $input = shift or die;
#my $col_num = shift or die;
my $col_num = 5;

die unless (-e $input);
open(IN, $input) or die;

#chr	pos	ref	percentage	dep	seq	codon_change	AA_change	gene	strand	func
#chr1	2044397	G	90	20	AAaA,AaAA,AaaAAAAAaa	GTG=>ATG	V=>M	AT1G06670	+	nuclear DEIH-boxhelicase


my $h = <IN>;
#my @h_a = split "\t", $h;
#print join("\t", (@h_a[0..2], "percentage", @h_a[3..$#h_a]) );
print $h;
#my %r;

my $cutoff = 30;

while (<IN>){
	my @a = split "\t";
	#print if ($a[0] eq "chr");
	if($a[3] >= $cutoff){
		print;
	}
	
}
close IN;

exit;