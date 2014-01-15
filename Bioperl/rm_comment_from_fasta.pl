#!/usr/bin/perl -w

#remove comment after id

use strict;
use Bio::SeqIO;
use Bio::Seq;

my $usage = "$0 <input_fasta>";
die $usage unless(@ARGV == 1);





my $seq_in = new Bio::SeqIO(-format => 'fasta',
             		-file => $ARGV[0]);
             		
my $seq_out = Bio::SeqIO->new(-file => ">$ARGV[0].no_comment",
                             -format => 'fasta');

while( my $seq = $seq_in->next_seq() )
{
	  my $out=Bio::Seq->new (-seq=>$seq->seq(),-id => $seq->id);
	  $seq_out->write_seq($out);
}

$seq_in->close;
$seq_out->close;

exit;