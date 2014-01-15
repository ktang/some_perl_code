#!/usr/bin/perl -w

use strict;
use File::Spec;

my $debug = 0;
#my $cutoff = 10;

#      1 chr
#      2 pos
#      3 ref
#      4 dep_WT_0  3
#      5 seq_WT_0
#      6 qual_WT_0
#      7 dep_WT_150
#      8 seq_WT_150
#      9 qual_WT_150
#     10 dep_ddc_0
#     11 seq_ddc_0
#     12 qual_ddc_0
#     13 dep_ddc_150
#     14 seq_ddc_150
#     15 qual_ddc_150
#     16 dep_nrpe1_0
#     17 seq_nrpe1_0
#     18 qual_nrpe1_0
#     19 dep_nrpe1_150
#     20 seq_nrpe1_150
#     21 qual_nrpe1_150
#     22 dep_017-6
#     23 seq_017-6
#     24 qual_017-6
#     25 dep_072-2
#     26 seq_072-2
#     27 qual_072-2
#     28 dep_109-5
#     29 seq_109-5
#     30 qual_109-5
#     31 dep_Nie129-3  30
#     32 seq_Nie129-3
#     33 qual_Nie129-3
#     34 per_WT_0
#     35 per_WT_150
#     36 per_ddc_0
#     37 per_ddc_150
#     38 per_nrpe1_0
#     39 per_nrpe1_150
#     40 per_017-6
#     41 per_072-2
#     42 per_109-5
#     43 per_Nie129-3

my $usage = "\n$0 \n <input> <outdir> <postfix> <cutoff> \n\n";
die $usage unless (@ARGV == 4);

my $input = shift or die;
#my $output = shift or die;
my $outdir = shift or die;

my $postfix = shift or die;

my $cutoff = shift or die;

die unless (-e $input);
#die if (-e $output);
#open (OUT, ">$output") or die;
open(IN, $input) or die;

my $h = <IN>;
chomp $h;
my @h_a = split "\t", $h;
#print join("\t", (@h_a[0..2], "percentage", @h_a[3..$#h_a]) );

#my @fhr;
#my @fhw;
#for my $i(0..$last_index){
#	open($fhr[$i], "<" , $files[$i]) or die "$files[$i]";
#	open($fhw[$i], ">" , $outputs[$i]) or die $outputs[$i], ": $!";
#	print {$fhw[$i]} $head, "\n";
#}


my @output_file_names;
my @fhw;
#0	1
#33	34

for my $i (0..5){
	my $order = $i + 1;
	my $label = $h_a[$i+33];
	$label =~ s/per_//;
	my $file_name =  $order . "_" . $label . "_". $postfix .".txt";
	my $output = File::Spec->catfile($outdir,$file_name);
	die if (-e $output);
	
	push @output_file_names, $file_name;
	print STDERR $file_name, "\n";
}

if ($debug) {
	exit;#code
}


for my $i (0..5){
	my $output = File::Spec->catfile($outdir, $output_file_names[$i]);
	die if (-e $output);
	open($fhw[$i], ">" , $output) or die $output, ": ", $!;
	print {$fhw[$i]} join("\t", @h_a[0..2, 3*$i + 3 , $i+33, 3*$i + 4, 3*$i + 5 ]), "\n";
}

while (<IN>){
	chomp;
	my @a = split "\t";
	#print if ($a[0] eq "chr");
	my ( $chr, $pos, $ref ) = @a[0..2];
	for my $i (0..5){
		my $dep  = $a[3*$i + 3];
		my $seq  = $a[3*$i + 4];
		my $qual = $a[3*$i + 5];
		my $per  = $a[$i+33];
		if ($per ne "NA" and $per >= $cutoff) {
			print  {$fhw[$i]} join("\t", ($chr, $pos, $ref, $dep, $per, $seq, $qual ) ), "\n";
		}
	}
}
close IN;

exit;
