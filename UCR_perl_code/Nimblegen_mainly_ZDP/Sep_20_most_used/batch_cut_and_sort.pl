#!/usr/bin/perl -w
# batch convert fastq file to fasta file
use strict;

my $usage = "$0 <indir> <outdir>";
die $usage unless(@ARGV == 2);
my ($indir, $outdir) = @ARGV[0..1];



opendir (INDIR, $indir) or die "Cannot open dir $indir:$!";
my @files = grep {/\.txt$/} readdir INDIR;
foreach my $file(@files){
     if($file =~ /(\S+)\.txt$/){
		 my $pre = $1;
		 my $output = $pre . "_sorted.txt";
		 my $cmd = "cut -f 1-4,7-12 $indir/$file > temp";
		 
		 print $cmd, "\n";
		 `$cmd`;
		 `sort -rnk 4 temp > $outdir/$output`;
		 `rm temp`;
	 }
}
