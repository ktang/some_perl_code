#!/usr/bin/perl -w
# batch convert fastq file to fasta file
use strict;

my $usage = "$0 <fastq dir> <fasta dir>";
die $usage unless(@ARGV >= 2);
my ($indir, $outdir) = @ARGV[0..1];

my $output_qual = 1;

opendir (INDIR, $indir) or die "Cannot open dir $indir:$!";
my @files = grep {/\.fastq$/} readdir INDIR;
foreach my $file(@files){
     if($file =~ /(\S+)\.fastq$/){
		 my $pre = $1;
		 my $output = $pre . ".fasta";
		 my $cmd = "perl fastq2fasta.pl $indir/$file $outdir/$output";
		 if($output_qual){
		     my $qual = $pre . ".qual";
			 $cmd .= " $outdir/$qual";
		 }
		 print $cmd, "\n";
		 $cmd .= " 2>&1";
		 my $msg = `$cmd 2>&1`;
		 print $msg, "\n";
	 }
}
