#!/usr/bin/perl -w

# calculate the difference between mutants and WT
# note only Col0 #2 is used

use strict;

my $usage = "$0 <in1> <out2> <out>";
die $usage unless(@ARGV == 3);
my ($in1,$in2,$outfile) = @ARGV[0..2];


open (C2, "/Users/kaitang/Desktop/Nimblegen/Sep_20_new_analysis/data/gff_ten/33_50533302_C2_nimblegen_method_sorted.txt")
	or die "cannot open C2";

open (IN1,"$in1")
	or die "cannot open $in1:$!";

open (IN2, "$in2")
	or die "cannot open $in2 :$!";

open (OUT,">$outfile")
	or die "cannot open $outfile :$!";

my ($l1,$l2,$c2);
my (@pts1,@pts2,@pts_c2);

while ($l1 = <IN1>,$l2 = <IN2>, $c2 = <C2>)
{
	chomp $l1;
	chomp $l2;
	chomp $c2;
	@pts1 = split "\t", $l1;
	@pts2 = split "\t", $l2;
	@pts_c2 = split "\t", $c2;
	
	if ( ($pts1[0] eq $pts2[0])and ($pts1[0] eq $pts_c2[0])and
	     ($pts1[3] eq $pts_c2[3]) and ($pts1[3] eq $pts2[3]))
	{
		my $mean = 0.5* ($pts1[5] + $pts2[5]);
		my $diff = $mean - $pts_c2[5];
		$pts1[5] = $diff;
		my $out_line = join "\t",@pts1;
		print OUT "$out_line\n";
	}
	else 
	{die "wrong";}	
}

close(IN1);
close(IN2);
close(OUT);
print STDERR "\a";
exit;