#!/usr/bin/perl -w
# get length distribution of small RNAs
use strict;
use Bio::SeqIO;

if(@ARGV >= 1){
	foreach my $file(@ARGV){
        my $seqin = Bio::SeqIO->new(-file=>$file, -format=>'fasta');
		my %reads;
		while(my $seq = $seqin->next_seq){
			$reads{$seq->length}++;
		}
		my $total = 0;
		print "Length distribution in file $file\n";
		foreach my $len(sort {$a<=>$b} keys %reads){
			print $len, "\t", $reads{$len}, "\n";
			$total += $reads{$len};
		}
		print "Total", "\t", $total, "\n";
		$seqin->close;
	}
}
