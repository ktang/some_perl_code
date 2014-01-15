#!/usr/bin/perl -w
# based on soap output, count how many sequences have at least one exact match
# to the genome (redundancy is taken into account)
use Bio::SeqIO;
use strict;

my $usage = "$0 <soap output> <unique seqs file>";
die $usage unless(@ARGV >= 2);
my ($soap, $unique) = @ARGV[0..1];

my $seqin = Bio::SeqIO->new(-file=>$unique, -format=>'fasta');
my %rdn;
while(my $seq = $seqin->next_seq){
	if(defined $rdn{$seq->id}){
		die "sequence id ", $seq->id, " is seen more than once";
	}else{
		if($seq->desc =~ /RDN=(\d+)/){
			$rdn{$seq->id} = $1;
		}else{
			die "Seq ", $seq->id, " has bad desc: ", $seq->desc;
		}
	}
}

open(SF, $soap) or die "Can't open $soap:$!";
my $matched = 0;
my $unique_match = 0;
my %seqsSeen;
while(<SF>){
	my @temp = split /\t/;
	if(!defined $seqsSeen{$temp[0]}){
		$seqsSeen{$temp[0]} = 1;
		$matched += $rdn{$temp[0]};
		$unique_match++;
	}
}
print STDERR "Number of unique seqs matching genome: $unique_match\n";
print STDERR "Number of total seqs matching genome: $matched\n";
	
