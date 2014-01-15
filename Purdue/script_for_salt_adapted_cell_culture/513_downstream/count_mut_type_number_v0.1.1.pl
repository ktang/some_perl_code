#!/usr/bin/perl -w

# v0.1.1
# combine del and ins


use strict;
my $debug = 0;

my %lables;
my %nums;
#~/misc/Zhu_Xiaohong/downstream/dep5/step1_all_common_and_other_23:57:04_N=622$ less not_all_common_3404loci.txt  | perl ../../src/count_mut_type_number.pl > not_all_common_3404loci_mutation_type_count.txt 

my $usage = "$0 \n <input_db_file> STDOUT\n\n";
die $usage unless (@ARGV == 1);
my $input = shift or die "input";
die unless (-e $input);

my @head = ("type");

open (IN, $input) or die;

my $h = <IN>;
chomp $h;
my @tmp = split "\t", $h;
for  (my $i = 3; $i<=13; $i+=2){
	push @head, $tmp[$i];
}


while (<IN>){
	chomp;
	my @a = split "\t";
#	my $indel = 0;
	my ( $chr, $pos, $ref ) = @a[0..2];
	my $x = 0;
	for  (my $i = 3; $i<=13; $i+=2){
		$x ++;
		my $lab = $a[$i];
		if( $lab =~ /INS/){
			$lables{ "INS" } = 1;
			$nums{$x}->{"INS" } ++;
		}elsif($lab =~ /DEL/){
			$lables{ "DEL" } = 1;
			$nums{$x}->{"DEL" } ++;
		}
		else{
			$lables{ $lab } = 1;
			$nums{$x}->{$lab} ++;
		}
		  
	}
}

close IN;

print join("\t", @head), "\n";

foreach my $lab ( sort keys %lables){
	print $lab, "\t";
	for my $x (1..5){
		my $n = 0;
		if (defined $nums{$x}->{$lab}  ){
			$n = $nums{$x}->{$lab};
		}
		print $n, "\t";
	}
	
	my $x = 6;
	my $n = 0;
	if (defined $nums{$x}->{$lab}  ){
		$n = $nums{$x}->{$lab};
	}
	print $n, "\n";
}

exit;

