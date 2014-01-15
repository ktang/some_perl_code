#!/usr/bin/perl -w

#get_depth_for_gff_file_v0.0.pl
# purpose: for transposon movement detection
use strict;
use File::Spec;

my $debug = 1;


my $usage = "\n $0 \n  <in_gff>  <in_bam> <output>\n\n";

die $usage unless (@ARGV == 3);

my $in_gff_file = shift or die;
my $in_bam_file = shift or die;
my $output_file = shift or die;

die unless (-e $in_bam_file);
die unless (-e $in_gff_file);
die if (-e $output_file);

if ($debug){
	print STDERR "\n\n", join("\n", ( $in_gff_file, $in_bam_file,$output_file )), "\n\n";
	exit;
}

my @gff_lines = ();
read_gff($in_gff_file, \@gff_lines);

open(OUT, ">>$output_file") or die;
print OUT join("\t", ("chr", "start", "end", "length", "dep_sum", "dep_avg", "ID", "strand")), "\n";

foreach my $i(0..$#gff_lines){
	my @a = split "\t", $gff_lines[$i];
	my ($chr, $start, $end) = @a[0,3,4];
	
	my $dep_sum = get_dep_sum($in_bam_file, $chr, $start, $end);
	my $length = $end - $start + 1;
	my $dep_avg = sprintf ("%.3f", $dep_sum/$length);
	
	print OUT join("\t", ($chr, $start, $end, $length, $dep_sum,  $dep_avg, $a[-1], $a[6] )),"\n"; 
}
close OUT;

exit;
#chr1    TAIR10  transposable_element_gene       433031  433819  .       -       .       ID=AT1G02228;Note=transposable_element_gene;Name=AT1G02228;Derives_from=AT1TE01405
#chr1    TAIR10  transposable_element_gene       846664  847739  .       +       .       ID=AT1G03420;Note=transposable_element_gene;Name=AT1G03420;Derives_from=AT1TE02770

#read_gff($in_gff_file, \@gff_lines);
sub read_gff{
	my($file, $ref) = @_;
	die unless (-e $file);
	
	open(IN, $file) or die;
	my $i = -1;
	while (<IN>){
		chomp;
		$i++;
		#my @a = split "\t";
		$ref->[$i] = $_;
	}
	close IN;
}


#get_dep_sum ($in_bam_file, $chr, $start, $end);
sub get_dep_sum{
	my ($file, $chr_sub, $s, $e )  = @_;
	die unless (-e $file);
	my $region = "$chr_sub:$s-$e";
	
	open(IN, "samtools depth -r $region $file |");
	my $sum = 0;
	while (<IN>){
		chomp;
		my @a = split "\t";
		$sum += $a[-1];
	}
	close IN;
	return $sum;
}