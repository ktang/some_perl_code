#!/usr/bin/perl -w

=head2
Copyright (C) 2010 Kai Tang
version 1 : for id which are identical, not /1 /2  
the basic idea is read one file in memary, and index
them in a hash, the read another file and handle,if 2 find pair in 1,
output them in two pair files,delete hash in 1.
Nov 24,2010
=cut

use strict;
use Bio::SeqIO;

my $usage = "$0 <input_pair1> <input_pair2> <output_pair1> <output_pair2> <SE_out>";
die $usage unless(@ARGV == 5);

my ($input_pair1, $input_pair2, $output_pair1, $output_pair2, $SE_out) = @ARGV[0..4];

my $inseq1 = Bio::SeqIO->new(-file => "<$input_pair1",-format => 'fastq');
my $inseq2 = Bio::SeqIO->new(-file => "<$input_pair2",-format => 'fastq');

if ( (-e $output_pair1) or (-e $output_pair2) or (-e $SE_out))
{die "\n\nat least one of the output_file exists!!\n\n"}

my $out1 = Bio::SeqIO->new(-file => ">$output_pair1",-format=>'fasta',-width=> 300);
my $out2 = Bio::SeqIO->new(-file => ">$output_pair2",-format=>'fasta',-width => 300);

my $sigle_out = Bio::SeqIO->new(-file => ">$SE_out",-format=>'fasta',-width=>300);

my %hash;

while (my $seq1 = $inseq1->next_seq)
{
	$hash {$seq1->id} = $seq1;
}

$inseq1->close;
my ($pair,$single1,$single2) = (0,0,0);
while (my $seq = $inseq2->next_seq)
{
	if (exists $hash{$seq->id})
	{
		$pair++;
		$out1->write_seq($hash{$seq->id});
		$out2->write_seq($seq);
		delete $hash{$seq->id};
	}
	else {	
		$single2++;
		$sigle_out->write_seq($seq);
	}
}
$inseq2->close;
foreach my $key (keys %hash)
{
	$single1++;
	$sigle_out->write_seq($hash{$key});
}

$sigle_out->close;

$out1->close;
$out2->close;
print "pairs = $pair\nsingle in $input_pair1 is $single1\nsingle in $input_pair2 is $single2\n";
exit;
