#!/usr/bin/perl -w
use strict;
use File::Spec;
my $debug = 0;


my %hash = (	"AT=>GC" => ["A=>G", "T=>C"],
		"GC=>AT" => ["G=>A", "C=>T"],
		"AT=>CG" => ["A=>C", "T=>G"],
		"AT=>TA" => ["A=>T", "T=>A"],
		"GC=>TA" => ["G=>T", "C=>A"],
		"GC=>CG" => ["G=>C", "C=>G"]
);

my @types = ("AT=>GC", "GC=>AT", "AT=>CG", "AT=>TA", "GC=>TA", "GC=>CG" );

if($debug){
	print STDERR "debug = 1\n\n";
}
my $usage = "$0 \n <input> STDOUT\n\n";
die $usage unless(@ARGV == 1);
my $input = shift or die;

die unless (-e $input);

open(IN, $input) or die;#
my $head = <IN>;
print $head;
my %nums; #h{mut_type}->file = 

while (<IN>) {
	chomp;
	my @a = split "\t";
	$nums{$a[0]} = [@a[1..6] ];
}

my %total;

foreach my $type (@types){
	my ($sub1, $sub2) = @{$hash{$type}};
	for my $i(1..6){
		$total{$type}->[$i-1]  =  $nums{$sub1}->[$i-1] + $nums{$sub2}->[$i-1]
	}
}

my @indel = ();
foreach my $type(keys %nums){
	if ($type =~ /INS|DEL/) {
		for my $i(1..6){
			$indel[$i -1 ] +=  $nums{$type}->[$i-1];
		}
	}
	
}

print join("\t", ("INDEL", @indel)),"\n";

foreach my $type (@types){
	print join("\t", ($type, @{$total{$type}})),"\n"
}

close IN;
exit;