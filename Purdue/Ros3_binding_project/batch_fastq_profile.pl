#!/usr/bin/perl -w
# batch convert fastq file to fasta file
use strict;

my $usage = "$0 <indir>";
die $usage unless(@ARGV == 1);
my $indir = $ARGV[0];


opendir (INDIR, $indir) or die "Cannot open dir $indir:$!";
my @files = grep {/\.fastq$/} readdir INDIR;
foreach my $file(@files){
	if($file =~ /(\S+)\.fastq$/){
		my $pre = $1;
		my $output = $pre . "_profile.txt";
		my $cmd = "time ./fastq_profile.pl -i $indir/$file -o $indir/$output";
		print $cmd, "\n";
		`$cmd`;
	}
}
