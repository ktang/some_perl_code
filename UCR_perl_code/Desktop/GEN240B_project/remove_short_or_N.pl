#!/usr/bin/perl -w
# remove the intergenic region which 
# is shorter than 200bp or has many Ns

use strict;
use Bio::SeqIO;

my $usage = "$0 <infile> <outfile> ";

       
my $seq_in = Bio::SeqIO->new('-file' => "<$ARGV[0]",
                             '-format' => 'fasta');

my $seq_out = Bio::SeqIO->new('-file' => ">$ARGV[1]",
                              '-format' => 'fasta');

while (my $inseq = $seq_in->next_seq)
{
	# if the sequence length is shorter than 200
	# or contains continuous more than 10 N, remove
	# the sequence from file
	if($inseq->length <=200 || $inseq->seq =~ /N{10,}/)
	{
		;
	}
	else
	{
		$seq_out->write_seq($inseq);
	}
}


$seq_in->close;
$seq_out->close;
print"\a";
exit;