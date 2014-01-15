#!/usr/bin/perl -w

=head2
Copyright (C) 2010 Kai Tang
version 

consensus
#Chr	pos	ref	depth	A	C	G	T
output
#Chr	pos	ref	depth	A	C	G	T	SNP_base


usage should be
samtools pileup |.pl > filterd.pileup
Nov 25,2010
=cut

my $debug = 0;

use strict;
while (<>){
	chomp;
	my @pts = split /\t/;
	my $max = max(@pts[4..7]);
	my $snp = "N";
	if($debug) {print STDERR "max = $max\t$snp\n";}
	
	if($max == $pts[4]){$snp = "A"}
		if($debug) {print STDERR "max = $max\t$snp\n";}
	if($max == $pts[5]){$snp = "C"}
		if($debug) {print STDERR "max = $max\t$snp\n";}
	if($max == $pts[6]){$snp = "G"}
		if($debug) {print STDERR "max = $max\t$snp\n";}
	if($max == $pts[7]){$snp = "T"}
		if($debug) {print STDERR "max = $max\t$snp\n";}
	print join("\t",(@pts,$snp)),"\n";
}


sub max{
	my @a_sorted = sort {$a<=>$b} @_;
	if ($a_sorted[3] == $a_sorted[2]){
		print STDERR "warning:\n",join("\t",@_), "\n";
	}	
	
	return $a_sorted[3];
}
exit;
