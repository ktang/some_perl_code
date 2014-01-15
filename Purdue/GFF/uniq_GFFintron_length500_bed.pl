#!/usr/bin/perl -w

#use File::Spec;
#my $outFile_ge9 = File::Spec->catfile($outDir, join(".", @parts) . "_ge9." . $ext);

use strict;
use File::Spec;

my $debug = 0;

if($debug){
	print STDERR "debug = 1\n\n";
}

#  /Users/tang58/DataBase/TAIR10/GFF/Kai/TAIR10_GFF3_intron_Kai.gff

my $length_cutoff = 500;

my $usage = "$0 \n <input_GFF> <output_bed>\n\n";
die $usage unless(@ARGV == 2);

my $input = shift or die;
my $output = shift or die;

die unless (-e $input);
die if( -e $output);

my %records;

my @array;

open(IN, $input) or die "cannot open $input: $!";
while(<IN>){
	chomp;
	my @a = split "\t";
	my ($chr, $start, $end) = ($a[0], $a[3], $a[4]);
	my $length = $a[4] - $a[3] + 1;
	next unless ($length >= $length_cutoff);
	next if ( $a[0] eq  "chrc" or $a[0] eq  "chrm");
	my $label = join("_", ($chr, $start, $end)) ;
	next if ( defined $records{ $label }) ;
	
	#$records{ $label } = $_;
	$records{ $label } = 1;
	push @array, $_;
}
close(IN);

die if(-e $output);
open(OUT, ">>$output") or die "cannot open $output: $!";

print OUT join("\t", ( "chr", "start", "end", "length", "strand", "Gene")), "\n";


#chr1    TAIR10  intron  3914    3995    .       +       .       Parent=AT1G01010.1

#foreach my $label (sort keys %records){
for my $i (0..$#array){
	my $line = $array[$i];
	my @a = split "\t", $line;
	my ($chr, $start, $end, $strand) = ( $a[0], @a[3..4], $a[6]);
	my ($temp, $gene) = split "=", $a[-1];
	my $length = $end - $start + 1;
	print OUT join("\t", ( $chr, $start, $end, $length, $strand, $gene)), "\n";
}
close(OUT);


