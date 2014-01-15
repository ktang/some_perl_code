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
#chr	pos	ref	mut_WT_0	per_WT_0	dep_WT_0	seq_WT_0	qual_WT_0	type	type_code	Gene	Strand	GeneType	GeneAnnot	TE	TEFamily	Intergenic	Promoter_1kb
#chr1	27767199	G	G=>INS	85.71	21	,+1a.+1Aa+1a.+1A,+1a.+1A.+1A.+1Aa.+1A.+1A,+1a.+1A.+1A,+1a,+1a,+1a.+2AA,+1a,a	*D'HAG#B'JI3JJA(:J9?!	IG	7	NONE	NONE	NONE	NONE	NONE	NONE	AT1G73840-AT1G73850	AT1G73840;AT1G73850
#0	1		2	3	4	5	6												7	8

my %num_records;

while (<IN>){
	chomp;
	my @a = split "\t";
	$num_records{$a[8]}++;
}

close IN;

print "name", "\t";

for my $k (sort keys %num_records){
#	print join("\t", ($k, $num_records{$k}) ), "\n";
	print $k, "\t";
}
print "\n";

print $input, "\t";
for my $k (sort keys %num_records){
#	print join("\t", ($k, $num_records{$k}) ), "\n";
	print $num_records{$k}, "\t";
}
print "\n";

exit;