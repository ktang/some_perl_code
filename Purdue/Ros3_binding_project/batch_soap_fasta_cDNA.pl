#!/usr/bin/perl -w
# 

use strict;

my $usage = "$0 <dir>";
die $usage unless(@ARGV == 1 );
my $indir = $ARGV[0];
#my ($indir, $db) = @ARGV[0..1];
my $db = "~/DataBase/TAIR9_cDNA/ath.cDNA.fa.index";

opendir (INDIR, $indir) or die "Cannot open dir $indir:$!";

my @files = grep {/\.fasta$/} readdir INDIR;



foreach my $file(@files){
	if($file =~/(\S+).fasta$/ ){
		my $pre = $1;    
		my $out = $1."_vs_cDNA_n1r2v2M4.soap";
		if (-e "$indir/$out"){die "$out exists!!!"}
		my $cmd = "soap -a $indir/$file -D $db -n 1 -r 2 -v 2 -M 4 -o $indir/$out 2>&1 | tail";
		print $cmd, "\n";
		system("$cmd");
	}
}

