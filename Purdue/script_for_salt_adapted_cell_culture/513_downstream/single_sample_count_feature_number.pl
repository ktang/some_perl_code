#!/usr/bin/perl -w

#step3_categorize_mutation_SNP_v0.2.pl
# modified from 
#/Users/tang58/misc/Zhu_Xiaohong/downstream/src/step3_categorize_mutation_SNP_v0.1.pl

# Apr 17, 2013
			
use strict;

my $debug = 1;
#my %records;

#my %refs;

#input
#chr     pos     ref     mut_WT_0        per_WT_0        mut_WT_150      per_WT_150
#~/misc/Zhu_Xiaohong/Mar27_filter_MAPQ20/mpileup/513loci_downstream_20:59:10_N=517$
#less 513loci_categorize_detail.txt 

my $usage = "$0 \n <input_db_file> STDOUT\n\n";
die $usage unless (@ARGV == 1);
my $input = shift or die "input";
die unless (-e $input);

open (IN, $input) or die;

my $h = <IN>;
chomp $h;
my @tmp = split "\t", $h;
#print $h, "\n";
#chr	pos	ref	mut_WT_0	per_WT_0	dep_WT_0	seq_WT_0	qual_WT_0	type	type_code
#chr1	1807621	G	G=>C	5.88	17	,,.c..,....,...,,	HID#IHDHJJHBFD@B?	IG	7



my %num_records;

while (<IN>){
	chomp;
	my @a = split "\t";
	$num_records{$a[-2]}++;
}

close IN;

for my $k (sort keys %num_records){
#	print join("\t", ($k, $num_records{$k}) ), "\n";
	print $k, "\t";
}
print "\n";
for my $k (sort keys %num_records){
#	print join("\t", ($k, $num_records{$k}) ), "\n";
	print $num_records{$k}, "\t";
}
print "\n";

exit;