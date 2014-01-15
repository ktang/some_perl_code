#!/usr/bin/perl -w

=head2
Copyright (C) 2010 Kai Tang
version 

usage should be
samtools pileup |.pl > filterd.pileup
Nov 25,2010
=cut

my $dep_cutoff = 5;
my $per_cutoff = 0.5;

use strict;
while (<>){
	chomp;
	my @pts = split /\t/;
	my $plus = ($pts[4]=~ tr/././ );
	my $minus = ($pts[4] =~ tr/,/,/);
	my $total = $plus + $minus;
#	if ($total > $pts[3] ){
#		print STDERR "$_\n";
#		die ("plus = $plus\tminus = $minus\n")
#	}
	if ( $pts[3] >= $dep_cutoff ){
		my $per = $total / $pts[3];
		if($per <= $per_cutoff){
			print $_, "\n";
		}
	}
}

exit;