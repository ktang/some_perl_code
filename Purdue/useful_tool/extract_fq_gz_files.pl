#!/usr/bin/perl -w

use strict;
use Bio::SeqIO;


my $usage = "$0 <infile> <outfile> <width>\n";
die $usage unless (@ARGV == 3);

my $infile = shift or die "input";
my $outfile = shift or die "output";
my $width = shift or die "width";

die unless (-e $infile);
die if (-e $outfile);

my $seq_in = Bio::SeqIO->new('-file'    => $infile,
                             '-format'  => 'fasta');
my $seq_out = Bio::SeqIO->new('-file'   => ">$outfile",
                              '-format' => 'fasta',
                              '-width'  => $width);

while (my $inseq = $seq_in->next_seq) {
	$seq_out->write_seq($inseq);
}
$seq_in->close;
$seq_out->close;
exit;