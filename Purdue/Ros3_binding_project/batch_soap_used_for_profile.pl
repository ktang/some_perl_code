#!/usr/bin/perl -w
# 

use strict;

my $usage = "$0 <dir> <database>";
die $usage unless(@ARGV == 2 );

my ($indir, $db) = @ARGV[0..1];


opendir (INDIR, $indir) or die "Cannot open dir $indir:$!";

my @files = grep {/\.fastq$/} readdir INDIR;



foreach my $file(@files){
     
		 my $cmd = "soap -a $indir/$file -D $db -n 1 -r 2 -v 2 -M 4 -o /dev/null 2>&1 | tail";
		 print $cmd, "\n";
		 system("$cmd");
	 
}

