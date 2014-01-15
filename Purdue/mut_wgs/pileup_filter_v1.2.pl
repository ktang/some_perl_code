#!/usr/bin/perl -w

=head2
Copyright (C) 2010 Kai Tang
version 1.1 (Feb 16,2011)
this is designed for use after get orignal SNP pileup,
which may contains only one read diff from ref.

usage:

v1.2 two cutoffs
 <input_pileup> <cutoff_No> <cutoff_percentage> STDOUT
depth >= cutoff_no and diff/depth >= cutoff_per

v1.1.pl <input_pileupOrPart> <cutoff> stdout

if cutoff < 1, then means 
diff/depth >= cutoff 
will be output;


if cutoff > 1, then means
diff > cutoff 
will be output.


usage should be
samtools pileup |.pl > filterd.pileup
Nov 25,2010
=cut
use strict;

my $usage = "$0 \n<input_pileup> <cutoff_No> <cutoff_percentage> STDOUT \n";

die $usage unless (@ARGV == 3);
my $file = $ARGV[0];

open (IN, $file)
	or die "cannot open $file:$!";

my $cutoff_no = $ARGV[1];
my $cutoff_per = $ARGV[2];

	while (<IN>){
		chomp;
		my @pts = split /\t/;
		my $plus = ($pts[4]=~ tr/././ );
		my $minus = ($pts[4] =~ tr/,/,/);
		my $same = $plus + $minus;
		my $total = $pts[3];
		my $diff = $total - $same;
	#	if ($same >= $total ){
	#		print STDERR "$_\n";
	#		die ("plus = $plus\tminus = $minus\n")
	#	}
		if ($diff/$total >= $cutoff_per  and $total >= $cutoff_no){
			print join("\t", @pts[0..4]), "\n";
		}
	}
	



exit;