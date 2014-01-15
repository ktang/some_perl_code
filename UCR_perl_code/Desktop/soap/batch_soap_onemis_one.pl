#!/usr/bin/perl -w
# batch map unique seqs to genome
use strict;

my $usage = "$0 <unique seq dir> <soap out dir> <database>";
die $usage unless(@ARGV >= 3);
my ($indir, $outdir, $db) = @ARGV[0..2];

opendir (INDIR, $indir) or die "Cannot open dir $indir:$!";
my @files = grep {/\.fasta$/} readdir INDIR;
foreach my $file(@files){
     if($file =~ /(\S+)\.fasta$/){
		 my $pre = $1;
		 my $output = $pre . "_vs_cdna.soapout";
		 my $cmd = "/mnt/disk2/kai/soap/soap2.19release/soap -a $indir/$file -D $db -v 1 -M 1 -r 1 -o $outdir/$output";
		 
		 print $cmd, "\n";
		 #$cmd .= " 2>&1";
		 #my $msg = `$cmd 2>&1`;
		 #print $msg, "\n";
		 `$cmd`;
	 }
}
