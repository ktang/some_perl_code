#!/usr/bin/perl -w

use strict;

my $usage = "$0 <dir> <file>";
die $usage unless(@ARGV == 2 );
my $indir = $ARGV[0];
my $infile = $ARGV[1];

#opendir (INDIR, $indir) or die "Cannot open dir $indir:$!";

my $pre = "NONE";
if ($infile =~ /(\S+).soap$/){
	$pre = $1;
}

my $output = $pre."_no_ChrCM.soap";
if (-e "$indir/$output" or $pre eq "NONE") {
	die "$output exists or wrong name!!!\n";
}

open (IN,"$indir/$infile") or die "cannot open $infile";
open (OUT,">$indir/$output") or die "cannot opne $output";

while (<IN>){
	chomp;
	my @a = split "\t";
	if($a[7] =~ /ATC/i or $a[7] =~ /ATM/i){
		next;
	}else{
		print OUT join("\t",@a), "\n";		
	}
}
close(IN);
close(OUT);