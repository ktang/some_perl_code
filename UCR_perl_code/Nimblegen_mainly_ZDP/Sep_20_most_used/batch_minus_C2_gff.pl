#!/usr/bin/perl -w

# calculate the difference between mutants and WT
# note only Col0 #2 is used

use strict;

my $usage = "$0 <in_dir> <out_dir>";
die $usage unless(@ARGV == 2);
my ($indir, $outdir) = @ARGV[0..1];

opendir (INDIR, $indir) or die "Cannot open dir $indir:$!";
my @files = grep {/\.gff$/} readdir INDIR;

open (C2, "/Users/kaitang/Desktop/Nimblegen/Sep_20_new_analysis/data/gff_ten/33_50533302_C2_nimblegen_method_sorted.txt")
	or die "cannot open C2";

my @c2 = <C2>;
close (C2);

  
foreach my $file(@files)
{
	
	my $i = 0;
	if($file =~ /(\S+)\.gff$/)
	{
		 my $pre = $1;
		 my $output = $pre . "_minusC2.gff";
		 
		 open(OUT, ">$outdir/$output")
		   or die "cannot open output file $output:$!";
		   
		 open (IN, "$indir/$file")
		   or die "cannot open input file $file:$!";
		   
		 while (my $line = <IN>)
		 {
		 	chomp $line;
		 	my @pts = split /\t/,$line;
		 	my $thisC2 = $c2[$i];
		 	$i++;
		 	chomp $thisC2;
		 	my @pts_C2 = split "\t" , $thisC2;
		 	if ( ($pts[0] eq $pts_C2[0] )and ($pts[3] eq $pts_C2[3]))
		 	{
		 		my $diff = $pts[5] - $pts_C2[5];
		 		$pts[5] = $diff;
		 		my $out_line = join "\t",@pts;
		 		print OUT "$out_line\n";
		 	}	
		 	else
		 	{die "wrong order"}
		 }	
	
		close(IN);
		close(OUT);
	} 
}

print STDERR "\a";
exit;