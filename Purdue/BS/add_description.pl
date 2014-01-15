#!/usr/bin/perl -w
#filter total percentage = 0;
use strict;

my $indir = ".";

opendir(DIR, $indir) or die "cannot open $indir";

my @inputs = grep /.xls$/, readdir DIR;

print STDERR join ("\n", @inputs) , "\n\n";

my $des_file = "/Users/tang58/DataBase/TAIR_Col0_genome/TAIR9_functional_descriptions";

my %descs;

open (DES, $des_file) or die "cannnot open $des_file";
while(<DES>){
	chomp;
	my @a = split "\t";
	my $gene = $a[0];# substr($a[0],0,9);
	$descs{$gene} = [@a[1..2]];
}
close(DES);

foreach my $input(@inputs){
	if($input =~ /(\S+)\.xls$/){
		my $output = $1."_annotation.xls";
		
		open (IN, $input);
		
		if (-e $output){
			print STDERR "$output exists\n";
			next;
		}
		
		else{
			open (OUT, ">$output");
			my $head = <IN>;
			chomp $head;
			print OUT $head,"\tType\tShort_description\n";
			while (<IN>){
				chomp;
				my @a = split "\t";
				if (defined $descs{$a[3]}){
					print OUT join("\t", ($_, @{$descs{$a[3]}})), "\n";
				}
				else{
					print STDERR "$a[3] has no description\n";
				}
			}
			close(OUT);
		}#else
		
		close(IN);
	}
}