#!/usr/bin/perl -w

use strict;
use Bio::SeqIO;

my $usage = "$0 <in1> <in2> <out>";
die $usage unless(@ARGV >= 3);
my ($input1, $input2, $outFile) = @ARGV[0..2];

my $IN1;

open ($IN1, '<', $input1)
	or die "Cannot open $input1: $!";



my $line;
my %genes;
my $name;

while($line = <$IN1>)
  {
  	chomp($line);
  	$genes{$line} = -1;
  }

close($IN1);

my $in = new Bio::SeqIO(-format => 'fasta',
             		-file => $input2);
             		
 my $seq_out = Bio::SeqIO->new('-file' => ">$outFile",
                              '-format' => 'fasta');

while( my $seq = $in->next_seq() )
{
  $name = $seq->id;
  if (exists $genes{$name})
  {
  	$seq_out->write_seq($seq);
  }
}

print "\a";
exit;

