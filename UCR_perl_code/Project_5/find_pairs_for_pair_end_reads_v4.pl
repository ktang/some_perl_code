#!/usr/bin/perl -w

=head2
Copyright (C) 2010 Kai Tang
version 4(just a try) :

ues ref as hash value; 
############ 

remove /1 /2 for fastq, and output the same id


for id which are identical, not /1 /2 ; only store the seq. no quality.
the basic idea is read one file in memary, and index
them in a hash, the read another file and handle,if 2 find pair in 1,
output them in two pair files,delete hash in 1.
Nov 24,2010
=cut

use strict;
use Bio::SeqIO;
use Bio::Seq;

my $debug = 0;

my $usage = "$0 <input_pair1> <input_pair2> <output_pair1> <output_pair2> <SE_out>";
die $usage unless(@ARGV == 5);

my ($input_pair1, $input_pair2, $output_pair1, $output_pair2, $SE_out) = @ARGV[0..4];

my $inseq1 = Bio::SeqIO->new(-file => "<$input_pair1",-format => 'fastq');
my $inseq2 = Bio::SeqIO->new(-file => "<$input_pair2",-format => 'fastq');

my $out1 = Bio::SeqIO->new(-file => ">$output_pair1",-format=>'fasta',-width=> 300);
my $out2 = Bio::SeqIO->new(-file => ">$output_pair2",-format=>'fasta',-width => 300);

my $sigle_out = Bio::SeqIO->new(-file => ">$SE_out",-format=>'fasta',-width=>300);

my %hash;

while (my $seq1 = $inseq1->next_seq)
{
	my ($id) = split /\//,$seq1->id;
	$hash {$id} = \$seq1->seq;
}

$inseq1->close;
my ($pair,$single1,$single2) = (0,0,0);
while (my $seq = $inseq2->next_seq)
{
	my ($id) = split /\//, $seq->id;
	my $outseq2 = Bio::Seq->new(-id => $id,-seq=>$seq->seq);
	if (exists $hash{$id})
	{
		$pair++;
	#	my $id1 = $id."/1";
	#	my $seq_out = Bio::Seq->new(-id => $id1,-seq =>$hash{$id});
		my $seq_out = Bio::Seq->new(-id => $id,-seq =>${$hash{$id}});
		$out1->write_seq($seq_out);
		
		
		$out2->write_seq($outseq2);
		delete $hash{$id};
	
	}
	else {	
		$single2++;
		
		$sigle_out->write_seq($outseq2);
	}
}

$inseq2->close;

foreach my $key (keys %hash)
{
	$single1++;
#	my $id = $key."/1";
	my $seq_out = Bio::Seq->new(-id => $key,-seq =>${$hash{$key}});
	$sigle_out->write_seq($seq_out);
}

$sigle_out->close;
$out1->close;
$out2->close;
print "pairs = $pair\nsingle in $input_pair1 is $single1\nsingle in $input_pair2 is $single2\n";
exit;
